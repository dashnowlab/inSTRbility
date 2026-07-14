#!/usr/bin/env python3
# Prevent OpenBLAS fork-deadlock: must be set BEFORE numpy/scipy import.
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
"""
calc_instability.py  —  Somatic Tandem-Repeat Instability Analysis
===================================================================

Two-regime genome-wide pipeline:

  Regime 1  NB-Geometric compound model (nbgeom_modelling.py)
            Fast, age-free, FFT-exact. Covers >95 % of WGS loci.

  Regime 2  Two-phase CTMC with age-marginalised inference
            (regime2_marginal.py).  Applied only to loci flagged by
            Regime 1 goodness-of-fit as having heavy-tailed distributions.
            Donor age is optional: supply --donor-age for per-year rates,
            or --tissue-type for a biological age prior.

Six-phase workflow
------------------
  Phase 0  Parse TSV reads + VCF founder lengths
  Phase 1  Regime 1 Pass 1  — estimate global step-size params
  Phase 2  Regime 1 Pass 2  — per-locus inference  (parallelised)
  Phase 3  Flag Regime 2 loci from R1 GoF output
  Phase 4  Regime 2 Pass 1  — estimate global CTMC params (T1, T2, p_exp)
  Phase 5  Regime 2 Pass 2  — per-locus CTMC inference + GoF (parallelised)
  Phase 6  Merge R1 and R2 results into one final TSV

Usage
-----
  # Regime 1 only
  python calc_instability.py -i reads.tsv -v genotypes.vcf -o out.tsv -t 8

  # + Regime 2, tissue age prior (no exact age needed)
  python calc_instability.py -i reads.tsv -v genotypes.vcf -o out.tsv -t 8 \\
      --run-regime2 --tissue-type adult_blood

  # + Regime 2 with exact donor age
  python calc_instability.py -i reads.tsv -v genotypes.vcf -o out.tsv -t 8 \\
      --run-regime2 --donor-age 65

  # + Regime 2 with per-sample age file
  python calc_instability.py -i reads.tsv -v genotypes.vcf -o out.tsv -t 8 \\
      --run-regime2 --age-file ages.tsv

  # + Regime 2 with known repeat preset (skips T2/p_exp estimation)
  python calc_instability.py -i reads.tsv -v genotypes.vcf -o out.tsv -t 8 \\
      --run-regime2 --repeat-preset htt_cag --donor-age 65

Input TSV (one row per read, tab-separated):
  chrom  start  end  motif  read_id  haplotype  length_bp  allele  avg_meth  meth_bases

VCF (optional):
  FORMAT/AL field in base pairs; divided internally by motif length.
  If absent, modal observed length is used as the founder allele.
"""

import argparse
import csv
import gzip
import sys
import time
import traceback
from multiprocessing import Pool
from pathlib import Path

import cyvcf2
import numpy as np
from tqdm import tqdm

from nbgeom_modelling import (
    _modal_length,
    auto_fit,
    estimate_regime1_global,
    Regime1GlobalParams,
)
from regime2_ctmc import CTMCGlobalParams, PRESET_MAP
from regime2_ctmc_v2 import estimate_T1_from_founders, estimate_global_params
from regime2_marginal import (
    TISSUE_AGE_PRIORS,
    _marginalised_pmf,
    fit_locus_marginal,
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="NB-Geometric + CTMC somatic TR instability analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    io = p.add_argument_group("I/O")
    io.add_argument("-i", "--input",   required=True,
                    help="Per-read TSV (chrom start end motif read_id hap "
                         "length_bp allele avg_meth meth_bases)")
    io.add_argument("-v", "--vcf",     default=None,
                    help="VCF with FORMAT/AL germline allele lengths (bp). "
                         "Omit to use modal observed length as founder.")
    io.add_argument("-o", "--output",  default=None,
                    help="Output prefix. Default: <input_stem>_instability")
    io.add_argument("-t", "--threads", type=int, default=1,
                    help="Worker processes (default: 1)")

    r1 = p.add_argument_group("Regime 1 options")
    r1.add_argument("--min-reads", type=int, default=10,
                    help="Min reads per haplotype (default: 10)")
    r1.add_argument("--n-ppp",     type=int, default=300,
                    help="PPP simulations for R1 GoF (default: 300; 0=skip)")
    r1.add_argument("--calibrate",  action="store_true",
                    help="Run Phase 2B: recalibrate q from Phase 2 dispersion "
                         "ratios and rerun Phase 2 with corrected params. "
                         "Recommended for first-run on any new sample.")
    r1.add_argument("--calibrate-min-reads", type=int, default=20,
                    help="Min reads per locus for calibration (default: 20)")

    r2 = p.add_argument_group("Regime 2 options")
    r2.add_argument("--run-regime2",    action="store_true",
                    help="Run CTMC on heavy-tailed loci flagged by R1 GoF")
    r2.add_argument("--tissue-type",    default="adult_blood",
                    choices=list(TISSUE_AGE_PRIORS),
                    help="Tissue-type age prior when --donor-age is absent "
                         "(default: adult_blood)")
    r2.add_argument("--donor-age",      type=float, default=None,
                    help="Exact donor age in years (overrides --tissue-type)")
    r2.add_argument("--age-file",       default=None,
                    help="TSV with columns sample_id and age. "
                         "Matched by input filename stem.")
    r2.add_argument("--repeat-preset",  default=None,
                    choices=list(PRESET_MAP),
                    help="Use a pre-defined CTMC preset (e.g. htt_cag) and "
                         "skip Pass 1 global-param estimation.")
    r2.add_argument("--regime2-min-loci", type=int, default=3,
                    help="Min flagged loci to run R2 (default: 3)")
    r2.add_argument("--regime2-threads",  type=int, default=None,
                    help="Workers for R2 inference. Defaults to --threads.")
    r2.add_argument("--r2-n-ppp",         type=int, default=200,
                    help="PPP simulations for R2 GoF (default: 200; 0=skip)")
    r2.add_argument("--r2-n-quad",        type=int, default=15,
                    help="Age quadrature points for R2 (default: 15)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def get_stratum(motif_len: int) -> int:
    if motif_len <= 3:   return 3
    elif motif_len <= 6: return 6
    return 7


def _founder_from_genotypes(genotypes, hap_idx: int, lengths: list) -> float:
    """Return founder length in repeat units from VCF genotypes or modal."""
    try:
        return float(genotypes[hap_idx])
    except (IndexError, TypeError, ValueError):
        li = np.round(np.asarray(lengths)).astype(int)
        v, c = np.unique(li, return_counts=True)
        return float(v[np.argmax(c)])


def winsorized_mean_abs_delta(
    lengths, founder_length: float, threshold: float = 100.0
) -> float:
    """Mean |Δ| Winsorized at threshold. Matches Handsaker et al. (2024)."""
    d = np.abs(np.asarray(lengths, dtype=float) - founder_length)
    return float(np.minimum(d, threshold).mean())


# ---------------------------------------------------------------------------
# VCF / age-file loading
# ---------------------------------------------------------------------------

def load_genotypes_from_vcf(vcf_path: str) -> dict:
    t0, genotypes = time.time(), {}
    for v in cyvcf2.VCF(vcf_path):
        if v.FILTER not in ("PASS", None):
            continue
        motif = v.INFO.get("MOTIFS", "")
        if not motif:
            continue
        al = v.format("AL")
        if al is not None and len(al):
            key = f"{v.CHROM}:{v.start}-{v.INFO.get('END')+1}_{motif}"
            genotypes[key] = al[0] / len(motif)
    print(f"Loaded {len(genotypes)} VCF genotypes in {time.time()-t0:.1f}s",
          file=sys.stderr)
    return genotypes


def load_age_from_file(age_file: str, sample_id: str) -> float | None:
    try:
        with open(age_file) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row.get("sample_id", "").strip() == sample_id:
                    return float(row["age"])
    except Exception as e:
        print(f"Warning: could not read age file: {e}", file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# Null GoF sentinel (when fit fails or GoF is skipped)
# ---------------------------------------------------------------------------

_NULL_GOF: dict = {
    "ks_stat": float("nan"), "ad_stat": float("nan"),
    "ppp_var": float("nan"), "dispersion_ratio": float("nan"),
    "chisq_p": float("nan"), "fit_quality": float("nan"),
    "fit_quality_adj": float("nan"), "regime2_recommended": float("nan"),
}

_NULL_R2_GOF: dict = {
    "r2_ks_stat":        float("nan"),
    "r2_ppp_var":        float("nan"),
    "r2_tail_coverage":  float("nan"),
    "r2_mean_bias":      float("nan"),
    "r2_fit_quality":    float("nan"),
    "r2_model_warning":  "",
}


# ---------------------------------------------------------------------------
# Output column definitions
# ---------------------------------------------------------------------------

_R1_COLS = [
    # Locus identity
    "chrom", "start", "end", "motif", "haplotype", "regime",
    # Allele summary
    "founder_length", "founder_length_bp",
    "avg_deviation", "winsorized_mean_abs_delta", "n_reads",
    # NB-Geometric posteriors (posterior mean only)
    "r_e", "r_c", "p_plus", "p_minus",
    "mu_expansion", "mu_contraction",
    # Instability indices
    "instability_index", "net_bias", "p_expansion",
    "expansion_rate", "contraction_rate",
    "net_mutation_rate", "p_any_mutation",
    # Regime 1 GoF
    "gof_ks_stat", "gof_ad_stat", "gof_ppp_var",
    "gof_ppp_skewness", "gof_observed_var", "gof_fitted_var",
    "gof_dispersion_ratio", "gof_chisq_stat", "gof_chisq_p",
    "gof_chisq_df", "gof_chisq_n_bins",
    "gof_regime2_recommended", "gof_fit_quality", "gof_fit_quality_adj",
    "elapsed_s",
]

_R2_COLS = [
    # Locus identity
    "chrom", "start", "end", "motif", "haplotype", "regime",
    # Allele summary
    "founder_length", "founder_length_bp",
    "avg_deviation", "winsorized_mean_abs_delta", "n_reads",
    # Age / global CTMC params
    "donor_age", "tissue_type", "T1", "T2", "p_exp",
    # CTMC rate posteriors
    "r1", "r1_ci_lo", "r1_ci_hi",
    "r2", "r2_ci_lo", "r2_ci_hi",
    # Derived rates
    "net_rate_phaseA", "net_rate_phaseA_ci_lo", "net_rate_phaseA_ci_hi",
    "net_rate_phaseB", "net_rate_phaseB_ci_lo", "net_rate_phaseB_ci_hi",
    "rate_acceleration", "rate_acceleration_ci_lo", "rate_acceleration_ci_hi",
    # Distribution summary
    "instability_index", "mean_delta", "phase_B_fraction",
    # Regime 2 GoF
    "r2_ks_stat", "r2_ppp_var",
    "r2_tail_coverage", "r2_mean_bias",
    "r2_fit_quality", "r2_model_warning",
    "elapsed_s",
]


# ---------------------------------------------------------------------------
# Regime 2 goodness-of-fit
# ---------------------------------------------------------------------------

def _compute_r2_gof(
    observed_lengths: np.ndarray,
    founder_length:   float,
    r1_est:           float,
    r2_est:           float,
    T1:               float,
    T2:               float,
    p_exp:            float,
    min_L:            int,
    max_L:            int,
    age_mean:         float,
    age_sd:           float,
    n_quad:           int,
    n_ppp:            int,
) -> dict:
    """
    Goodness-of-fit metrics for the Regime 2 (marginalised CTMC) fit.

    Metrics
    -------
    r2_ks_stat : float
        KS statistic between empirical CDF of observed deltas and the
        fitted marginalised PMF CDF.
        Good: < 0.15.  Poor: > 0.30.

    r2_ppp_var : float
        Posterior predictive p-value on variance. Fraction of n_ppp
        simulated datasets (drawn from the fitted PMF) whose variance
        exceeds the observed variance.
        Well-calibrated: 0.10–0.90.  Systematic misfit: < 0.05 or > 0.95.

    r2_tail_coverage : float
        Fraction of observed reads that fall below the 99th percentile of
        the fitted marginalised PMF. Measures whether the tail of the
        distribution is captured.
        Good: > 0.95.  Poor (tail underestimated): < 0.90.

    r2_mean_bias : float
        (fitted_mean_delta − observed_mean_delta) / observed_std_delta.
        Standardised bias of the fitted mean. Good: |bias| < 0.5.

    r2_fit_quality : float in [0, 1]
        Composite score combining all four metrics.
        > 0.60 = good fit.  < 0.40 = poor fit.

    r2_model_warning : str
        Human-readable warning if the fit is suspect.
        Empty string when fit is acceptable.
    """
    if n_ppp == 0:
        return dict(_NULL_R2_GOF)

    obs     = np.asarray(observed_lengths, dtype=float)
    L0      = int(round(founder_length))
    deltas  = obs - L0
    n       = len(deltas)
    obs_var = float(deltas.var())
    obs_mean = float(deltas.mean())
    obs_std  = max(float(deltas.std()), 1e-6)

    # Fitted PMF
    pmf = _marginalised_pmf(
        L0, r1_est, r2_est, T1, T2, p_exp, min_L, max_L,
        age_mean, age_sd, n_quad,
    )
    lengths_grid = np.arange(min_L, max_L + 1, dtype=float)
    deltas_grid  = lengths_grid - L0
    fitted_mean  = float((deltas_grid * pmf).sum())

    # ── KS statistic ────────────────────────────────────────────────────────
    # Compare empirical CDF of observed deltas vs fitted PMF CDF
    # Use only the delta range covered by both
    pmf_cdf  = np.cumsum(pmf)
    obs_sorted = np.sort(deltas)
    ks_vals  = []
    for d in obs_sorted:
        idx = int(d - min_L)
        if 0 <= idx < len(pmf_cdf):
            emp_cdf = float(np.searchsorted(obs_sorted, d, side="right")) / n
            ks_vals.append(abs(emp_cdf - pmf_cdf[idx]))
    ks_stat = float(max(ks_vals)) if ks_vals else float("nan")

    # ── PPP variance ────────────────────────────────────────────────────────
    rng = np.random.default_rng(seed=0)
    pmf_clipped = np.clip(pmf, 0, None)
    pmf_clipped /= pmf_clipped.sum()
    support = np.arange(min_L, max_L + 1, dtype=int)
    sim_vars_exceed = 0
    for _ in range(n_ppp):
        sim_lengths = rng.choice(support, size=n, p=pmf_clipped)
        sim_deltas  = sim_lengths.astype(float) - L0
        if float(sim_deltas.var()) >= obs_var:
            sim_vars_exceed += 1
    ppp_var = float(sim_vars_exceed) / n_ppp

    # ── Tail coverage ────────────────────────────────────────────────────────
    # 99th percentile of fitted PMF
    cdf_99_idx = int(np.searchsorted(pmf_cdf, 0.99))
    if cdf_99_idx < len(lengths_grid):
        p99_length = float(lengths_grid[cdf_99_idx])
    else:
        p99_length = float(lengths_grid[-1])
    tail_coverage = float(np.mean(obs <= p99_length))

    # ── Mean bias ────────────────────────────────────────────────────────────
    mean_bias = (fitted_mean - obs_mean) / obs_std

    # ── Composite fit quality ─────────────────────────────────────────────
    # KS factor: 1.0 if ks<0.10, linear decay to 0.0 at ks=0.40
    ks_factor = max(0.0, 1.0 - max(ks_stat - 0.10, 0.0) / 0.30)

    # PPP factor: 1.0 if ppp in [0.10,0.90], decay outside
    if 0.10 <= ppp_var <= 0.90:
        ppp_factor = 1.0
    elif ppp_var < 0.10:
        ppp_factor = ppp_var / 0.10
    else:
        ppp_factor = max(0.0, 1.0 - (ppp_var - 0.90) / 0.10)

    # Tail factor: 1.0 if ≥0.95, linear decay below
    tail_factor = min(1.0, tail_coverage / 0.95)

    # Bias factor: 1.0 if |bias|<0.5, linear decay to 0.0 at |bias|=2.0
    bias_factor = max(0.0, 1.0 - max(abs(mean_bias) - 0.5, 0.0) / 1.5)

    fit_quality = float(
        0.35 * ks_factor
        + 0.30 * ppp_factor
        + 0.20 * tail_factor
        + 0.15 * bias_factor
    )

    # ── Warnings ─────────────────────────────────────────────────────────────
    warnings = []
    if tail_coverage < 0.85:
        warnings.append("tail_undercovered")
    if ppp_var < 0.05:
        warnings.append("overdispersed")
    if ppp_var > 0.95:
        warnings.append("underdispersed")
    if abs(mean_bias) > 1.0:
        warnings.append("mean_bias_high")
    if ks_stat > 0.30:
        warnings.append("ks_poor")
    # Special warning for loci where no single-L0 CTMC fits well
    obs_max_delta = float(np.max(np.abs(deltas)))
    if tail_coverage < 0.85 and obs_max_delta > 200:
        warnings.append("heterogeneous_history_likely")

    return {
        "r2_ks_stat":       round(ks_stat,      4),
        "r2_ppp_var":       round(ppp_var,       4),
        "r2_tail_coverage": round(tail_coverage, 4),
        "r2_mean_bias":     round(mean_bias,     4),
        "r2_fit_quality":   round(fit_quality,   4),
        "r2_model_warning": ";".join(warnings),
    }


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def parse_input_tsv(
    input_path: str,
    genotypes:  dict | None,
) -> list:
    """
    Stream-parse the per-read TSV into a list of locus dicts.

    Each dict: {chrom, start, end, motif, haplotypes: {hap: [lengths_ru]},
                genotypes: array|None}
    Lengths are converted to repeat units (bp / motif_size).
    """
    ins = (gzip.open(input_path, "rt")
           if input_path.endswith(".gz")
           else open(input_path, "rt"))

    data:  list = []
    info:  dict = {}
    haps:  dict = {}
    prev_key: str | None = None

    pbar = tqdm(unit="rows", unit_scale=True, ncols=80,
                smoothing=0.1, position=1, desc="Reading TSV")

    for line in ins:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        (chrom, start, end, motif, _read_id,
         haplotype, length_bp, _allele,
         _avg_meth, _meth_bases) = parts

        start  = int(start)
        end    = int(end)
        ml     = len(motif)
        span   = end - start
        pbar.update(1)

        # Filter: motif ≤ 6 bp, span 20–500 bp
        if ml > 6 or span < 20 or span > 500:
            continue

        key = f"{chrom}:{start}-{end}_{motif}"
        info.setdefault(key, {"chrom": chrom, "start": start,
                               "end": end,   "motif": motif})

        hap  = int(haplotype)
        unit = float(length_bp) # / ml

        if key != prev_key and prev_key is not None:
            _flush(data, info, prev_key, haps, genotypes)
            haps = {}
            del info[prev_key]

        haps.setdefault(hap, []).append(unit)
        prev_key = key

    if haps and prev_key is not None:
        _flush(data, info, prev_key, haps, genotypes)

    pbar.close()
    ins.close()
    return data


def _flush(data, info, key, haps, genotypes):
    data.append({
        "chrom":      info[key]["chrom"],
        "start":      info[key]["start"],
        "end":        info[key]["end"],
        "motif":      info[key]["motif"],
        "haplotypes": dict(haps),
        "genotypes":  (genotypes.get(key) if genotypes else None),
    })


# ---------------------------------------------------------------------------
# Phase 1 — Regime 1 global parameter estimation
# ---------------------------------------------------------------------------

def _build_pass1_data(catalog: list, min_reads: int = 10) -> list:
    out = []
    for locus in catalog:
        ml = len(locus["motif"])
        for hap_idx, (hap, lengths) in enumerate(locus["haplotypes"].items()):
            arr = np.asarray(lengths, dtype=float)
            if len(arr) < min_reads:
                continue
            founder = _founder_from_genotypes(
                locus["genotypes"], hap_idx, arr
            )
            # MAD outlier removal
            med = float(np.median(arr))
            mad = max(float(np.median(np.abs(arr - med))),
                      float(arr.std()) + 1e-6)
            arr = arr[np.abs(arr - med) <= 5.0 * mad]
            if len(arr) < min_reads:
                continue
            out.append({
                "deltas":         (np.round(arr) - round(founder)).astype(int),
                "founder_length": founder,
                "locus_id":       f"{locus['chrom']}_{locus['start']}",
                "haplotype":      hap,
                "motif_length":   ml,
            })
    return out


def run_r1_pass1(catalog: list, min_reads: int = 10) -> dict:
    """Regime 1 Pass 1: estimate global step-size params, stratified by motif."""
    loci_data = _build_pass1_data(catalog, min_reads)
    print(f"  R1 Pass 1: {len(loci_data)} haplotype-loci", file=sys.stderr)

    strata: dict = {}
    for e in loci_data:
        s = get_stratum(e["motif_length"])
        strata.setdefault(s, []).append(e)

    gp_map: dict = {}
    labels = {3: "1-3 bp", 6: "4-6 bp", 7: "7+ bp"}
    for s, entries in sorted(strata.items()):
        print(f"    Stratum {labels[s]}: {len(entries)} haplotype-loci",
              file=sys.stderr)
        if len(entries) < 4:
            gp_map[s] = Regime1GlobalParams.from_defaults(s)
        else:
            gp_map[s] = estimate_regime1_global(entries, min_reads=min_reads)
            print(f"    {gp_map[s]}", file=sys.stderr)
    return gp_map


# ---------------------------------------------------------------------------
# Phase 2B — Regime 1 q recalibration (optional)
# ---------------------------------------------------------------------------

def run_r1_calibration(
    gp_map:    dict,
    r1_tsv:    str,
    min_reads: int = 20,
) -> dict:
    """
    Phase 2B: Recalibrate the NegBin q parameter per stratum from the
    dispersion ratios in the Phase 2 output TSV.

    The NB-Geometric model has a q parameter controlling how many mutation
    events occur per cell lineage. If q is wrong, the model's fitted variance
    systematically differs from the observed variance:

      dispersion_ratio = observed_var / fitted_var
      median < 1 → model overpredicts variance → increase q
      median > 1 → model underpredicts variance → decrease q

    This method reads gof_dispersion_ratio from the TSV, computes the
    median per stratum (filtering by n_reads >= min_reads and regime=R1),
    and back-calculates the q that would bring the median to 1.0.

    Calibration is per-stratum because motif length strongly affects the
    somatic mutation rate and thus the appropriate q:
      Short motifs (1-4bp): typically need higher q (fewer events per lineage)
      Long motifs  (5-6bp): typically well-calibrated at the initial q

    Parameters
    ----------
    gp_map : dict
        {stratum: Regime1GlobalParams} from run_r1_pass1().
    r1_tsv : str
        Path to Phase 2 output TSV.
    min_reads : int
        Minimum reads per locus to use in calibration (default: 20).
        Low-depth loci have noisy dispersion ratios.

    Returns
    -------
    dict : calibrated {stratum: Regime1GlobalParams}
    """
    import csv

    # Read dispersion ratios and n_reads per stratum from TSV
    stratum_disp: dict[int, list] = {3: [], 6: [], 7: []}
    try:
        with open(r1_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    n   = int(float(row["n_reads"]))
                    dr  = float(row["gof_dispersion_ratio"])
                    ml  = len(row["motif"])
                    reg = row.get("regime", "R1")
                except (ValueError, KeyError):
                    continue
                if not (n >= min_reads and np.isfinite(dr) and dr > 0
                        and reg == "R1"):
                    continue
                s = get_stratum(ml)
                if s in stratum_disp:
                    stratum_disp[s].append(dr)
    except FileNotFoundError:
        print(f"  Calibration: {r1_tsv} not found, skipping.", file=sys.stderr)
        return gp_map

    labels = {3: "1-3bp", 6: "4-6bp", 7: "7+bp"}
    new_gp_map = {}
    import copy

    for stratum, gp in gp_map.items():
        disp_vals = stratum_disp.get(stratum, [])
        if len(disp_vals) < 5:
            print(f"  Stratum {labels[stratum]}: too few loci ({len(disp_vals)}) "
                  f"for calibration — keeping q={gp.q_e:.3f}", file=sys.stderr)
            new_gp_map[stratum] = gp
            continue

        median_dr = float(np.median(disp_vals))

        # Back-calculate new q.
        # NB variance ∝ (1-q)/q.  We want new_var = old_var * median_dr.
        # => new_rate = old_rate * median_dr  where rate = (1-q)/q
        current_rate = (1.0 - gp.q_e) / gp.q_e
        new_rate     = current_rate * median_dr
        new_q        = float(np.clip(1.0 / (1.0 + new_rate), 0.05, 0.99))

        direction = ("overpredicts" if median_dr < 1 else "underpredicts")
        print(f"  Stratum {labels[stratum]} (n={len(disp_vals)}): "
              f"median disp_ratio={median_dr:.3f} — model {direction} variance "
              f"→ q: {gp.q_e:.3f} → {new_q:.3f}", file=sys.stderr)

        new_gp = copy.deepcopy(gp)
        new_gp.q_e = new_q
        new_gp.q_c = new_q
        new_gp_map[stratum] = new_gp

    return new_gp_map


# ---------------------------------------------------------------------------
# Phase 2 — Regime 1 per-locus worker
# ---------------------------------------------------------------------------

def _r1_worker(
    chunk:     list,
    gp_map:    dict,
    fout:      str,
    thread_id: int,
    n_ppp:     int,
) -> None:
    """Process a chunk of loci with R1 and write TSV rows."""
    import os
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    out = open(fout, "wt")
    if thread_id == 0:
        print("\t".join(_R1_COLS), file=out)

    null_gof = {
        "ks_statistic": float("nan"), "ad_statistic": float("nan"),
        "ppp_variance":  float("nan"), "ppp_skewness":  float("nan"),
        "observed_var":  float("nan"), "fitted_var":    float("nan"),
        "dispersion_ratio": float("nan"),
        "chisq_statistic": float("nan"), "chisq_pvalue":  float("nan"),
        "chisq_df": 0, "chisq_n_bins": 0,
        "regime2_recommended": float("nan"),
        "fit_quality": float("nan"), "fit_quality_adjusted": float("nan"),
    }

    for locus in chunk:
        chrom, start, end, motif = (
            locus["chrom"], locus["start"], locus["end"], locus["motif"]
        )
        ml  = len(motif)
        gp  = gp_map[get_stratum(ml)]
        gt  = locus["genotypes"]

        for hap_idx, (hap, lengths) in enumerate(locus["haplotypes"].items()):
            t0      = time.time()
            founder = _founder_from_genotypes(gt, hap_idx, lengths)

            try:
                result, diag, gof = auto_fit(
                    lengths,
                    founder_length  = founder,
                    regime1_params  = gp,
                    haplotype_label = f"hap{hap}",
                    n_ppp           = n_ppp,
                )
            except Exception as e:
                print(f"WARN R1 {chrom}:{start} hap{hap}: {e}",
                      file=sys.stderr)
                continue

            if gof is None:
                gof = dict(null_gof)
                gof["regime2_recommended"] = int(diag.is_heavy_tailed)

            avg   = float(np.mean(np.abs(
                np.asarray(lengths, float) - result.founder_length
            )))
            wmad  = winsorized_mean_abs_delta(lengths, result.founder_length)

            row = [
                chrom, start, end, motif, hap, "R1",
                round(result.founder_length, 4),
                round(result.founder_length * ml, 1),
                round(avg,  4), round(wmad, 4),
                result.n_reads,
                round(result.r_e[0],            4),
                round(result.r_c[0],            4),
                round(result.p_plus[0],         4),
                round(result.p_minus[0],        4),
                round(result.mu_expansion[0],   4),
                round(result.mu_contraction[0], 4),
                round(result.instability_index[0], 4),
                round(result.net_bias[0],          4),
                round(result.p_expansion[0],       4),
                round(result.expansion_rate[0],    4),
                round(result.contraction_rate[0],  4),
                round(result.net_mutation_rate[0], 4),
                round(result.p_any_mutation[0],    4),
                gof["ks_statistic"],    gof["ad_statistic"],
                gof["ppp_variance"],    gof["ppp_skewness"],
                gof["observed_var"],    gof["fitted_var"],
                gof["dispersion_ratio"],
                gof["chisq_statistic"], gof["chisq_pvalue"],
                gof["chisq_df"],        gof["chisq_n_bins"],
                gof["regime2_recommended"],
                gof["fit_quality"],     gof["fit_quality_adjusted"],
                round(time.time() - t0, 3),
            ]
            print(*row, sep="\t", file=out)

    out.close()


def _run_r1(catalog, gp_map, fout, threads, n_ppp, pbar):
    chunks = np.array_split(catalog, threads)
    fouts  = [fout] + [f"{fout}.r1p{i}" for i in range(1, threads)]

    def _ok(_): pbar.update()
    def _err(e):
        pbar.update()

        print(f"\nR1 worker error:", file=sys.stderr)
        traceback.print_exception(type(e), e, e.__traceback__)

    with Pool(processes=threads) as pool:
        for i in range(threads):
            pool.apply_async(
                _r1_worker,
                args=(chunks[i].tolist(), gp_map, fouts[i], i, n_ppp),
                callback=_ok, error_callback=_err,
            )
        pool.close(); pool.join()

    if threads > 1:
        with open(fout, "at") as f:
            for i in range(1, threads):
                p = f"{fout}.r1p{i}"
                try:
                    f.write(open(p).read())
                    Path(p).unlink()
                except FileNotFoundError:
                    pass


# ---------------------------------------------------------------------------
# Phase 3 — Flag Regime 2 loci
# ---------------------------------------------------------------------------

def identify_regime2_loci(r1_tsv: str, catalog: list) -> tuple[set, dict]:
    """
    Read R1 TSV and return (flagged_locus_keys, catalog_index).
    A locus is flagged when gof_regime2_recommended == 1.
    """
    idx = {
        f"{l['chrom']}:{l['start']}-{l['end']}_{l['motif']}": l
        for l in catalog
    }
    flagged: set = set()
    try:
        with open(r1_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    flag = int(float(row.get("gof_regime2_recommended", 0) or 0))
                except (ValueError, TypeError):
                    flag = 0
                if flag == 1:
                    k = (f"{row['chrom']}:{row['start']}"
                         f"-{row['end']}_{row['motif']}")
                    flagged.add(k)
    except Exception as e:
        print(f"Warning: could not read R1 output: {e}", file=sys.stderr)
    return flagged, idx


# ---------------------------------------------------------------------------
# Phase 4 — Regime 2 global parameter estimation
# ---------------------------------------------------------------------------

def run_r2_pass1(
    flagged:       set,
    catalog_idx:   dict,
    donor_age:     float | None,
    tissue_type:   str,
    preset:        str | None,
    min_loci:      int,
) -> CTMCGlobalParams | None:
    """
    Estimate (or load preset) CTMC global params from flagged loci.
    Returns None when estimation is not possible.
    """
    if preset is not None:
        gp = CTMCGlobalParams.from_preset(preset)
        print(f"  R2 Pass 1: using preset '{preset}'", file=sys.stderr)
        return gp

    if len(flagged) < min_loci:
        print(f"  R2 Pass 1: only {len(flagged)} flagged loci "
              f"(need ≥ {min_loci}). Skipping.", file=sys.stderr)
        return None

    # Resolve age for Pass 1B
    if donor_age is not None:
        age_mean = donor_age
    else:
        age_mean = TISSUE_AGE_PRIORS[tissue_type]["mean"]

    loci_data, founders = [], []
    for key in flagged:
        locus = catalog_idx.get(key)
        if locus is None:
            continue
        for hap_idx, (hap, lengths) in enumerate(locus["haplotypes"].items()):
            if len(lengths) < 5:
                continue
            f = _founder_from_genotypes(locus["genotypes"], hap_idx, lengths)
            founders.append(f)
            loci_data.append({
                "lengths":        np.asarray(lengths, dtype=float),
                "founder_length": f,
                "donor_age":      age_mean,
            })

    if len(loci_data) < min_loci:
        print(f"  R2 Pass 1: too few valid haplotypes ({len(loci_data)}).",
              file=sys.stderr)
        return None

    T1 = estimate_T1_from_founders(founders, percentile=5.0)
    print(f"  R2 Pass 1: T1={T1:.1f} (from {len(founders)} founder lengths)",
          file=sys.stderr)

    gp = estimate_global_params(loci_data, T1=T1, verbose=True)
    print(f"  {gp}", file=sys.stderr)
    return gp


# ---------------------------------------------------------------------------
# Phase 5 — Regime 2 per-locus worker
# ---------------------------------------------------------------------------

def _r2_worker(
    chunk:       list,
    ctmc_gp:     CTMCGlobalParams,
    donor_age:   float | None,
    tissue_type: str,
    fout:        str,
    thread_id:   int,
    n_ppp:       int,
    n_quad:      int,
) -> None:
    """Process a chunk of Regime 2 loci and write TSV rows."""
    import os
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    # TISSUE_AGE_PRIORS and _marginalised_pmf imported at module level
    out = open(fout, "wt")
    if thread_id == 0:
        print("\t".join(_R2_COLS), file=out)

    # Resolve age prior once per worker
    if donor_age is not None:
        age_mean  = float(donor_age)
        age_sd    = 0.5
        rep_age   = age_mean
        ttype_out = f"exact:{donor_age:.1f}"
    else:
        ap       = TISSUE_AGE_PRIORS[tissue_type]
        age_mean = ap["mean"]
        age_sd   = ap["sd"]
        rep_age  = age_mean
        ttype_out = tissue_type

    for item in chunk:
        chrom, start, end, motif = (
            item["chrom"], item["start"], item["end"], item["motif"]
        )
        ml      = len(motif)
        hap     = item["hap"]
        lengths = item["lengths"]
        founder = item["founder"]
        t0      = time.time()

        try:
            result = fit_locus_marginal(
                lengths,
                founder_length  = founder,
                global_params   = ctmc_gp,
                donor_age       = donor_age,
                tissue_type     = tissue_type,
                n_quad          = n_quad,
                r1_n            = 30,
                r2_n            = 12,
                haplotype_label = f"hap{hap}",
            )
        except Exception as e:
            print(f"WARN R2 {chrom}:{start} hap{hap}: {e}", file=sys.stderr)
            continue

        # Compute Regime 2 GoF
        lengths_int = np.round(np.asarray(lengths)).astype(int)
        L0_int      = int(round(founder))
        obs_max     = int(lengths_int.max())
        obs_min     = int(lengths_int.min())
        obs_span    = max(obs_max - L0_int, L0_int - obs_min, 1)
        # min_L must not exceed L0 (L0 must be inside the state space)
        min_L = min(max(obs_min - 10, 0), L0_int)

        # Adaptive cap: cover observed data but limit matrix size
        # Cap chosen so n_states = max_L - min_L + 1 stays manageable.
        # obs_max + 30 buffer, but hard cap at min_L + 150 (151 states max).
        # Reads beyond max_L go into the absorbing sink -- same as Handsaker et al.
        max_L = min(
            obs_max + 30,           # cover observed data with small buffer
            min_L + 150,            # hard cap: n_states ≤ 151
            1000,                   # absolute ceiling
        )

        gof2 = _compute_r2_gof(
            observed_lengths = np.asarray(lengths, float),
            founder_length   = founder,
            r1_est           = result.r1[0],
            r2_est           = result.r2[0],
            T1               = ctmc_gp.T1,
            T2               = ctmc_gp.T2,
            p_exp            = ctmc_gp.p_exp,
            min_L            = min_L,
            max_L            = max_L,
            age_mean         = age_mean,
            age_sd           = age_sd,
            n_quad           = n_quad,
            n_ppp            = n_ppp,
        )

        avg  = float(np.mean(np.abs(np.asarray(lengths, float) - founder)))
        wmad = winsorized_mean_abs_delta(lengths, founder)

        def _v(t, i):
            v = t[i]
            return round(v, 6) if (v is not None and v == v) else ""

        row = [
            chrom, start, end, motif, hap, "R2",
            round(founder,          4),
            round(founder * ml,     1),
            round(avg,  4), round(wmad, 4),
            result.n_reads,
            round(rep_age,          1),
            ttype_out,
            round(result.T1,        2),
            round(result.T2,        2),
            round(result.p_exp,     4),
            _v(result.r1,  0), _v(result.r1,  1), _v(result.r1,  2),
            _v(result.r2,  0), _v(result.r2,  1), _v(result.r2,  2),
            _v(result.net_rate_phaseA, 0),
            _v(result.net_rate_phaseA, 1),
            _v(result.net_rate_phaseA, 2),
            _v(result.net_rate_phaseB, 0),
            _v(result.net_rate_phaseB, 1),
            _v(result.net_rate_phaseB, 2),
            _v(result.rate_acceleration, 0),
            _v(result.rate_acceleration, 1),
            _v(result.rate_acceleration, 2),
            round(result.instability_index[0],  4),
            round(result.mean_delta[0],         4),
            round(result.phase_B_fraction[0],   4),
            gof2["r2_ks_stat"],
            gof2["r2_ppp_var"],
            gof2["r2_tail_coverage"],
            gof2["r2_mean_bias"],
            gof2["r2_fit_quality"],
            gof2["r2_model_warning"],
            round(time.time() - t0, 3),
        ]
        print(*row, sep="\t", file=out)

    out.close()


def _run_r2(r2_items, ctmc_gp, donor_age, tissue_type,
            fout, threads, n_ppp, n_quad, pbar):
    chunks = np.array_split(r2_items, threads)
    fouts  = [fout] + [f"{fout}.r2p{i}" for i in range(1, threads)]

    def _ok(_): pbar.update()
    def _err(e):
        pbar.update()
        print(f"R2 worker error: {e}", file=sys.stderr)

    with Pool(processes=threads) as pool:
        for i in range(threads):
            pool.apply_async(
                _r2_worker,
                args=(chunks[i].tolist(), ctmc_gp, donor_age, tissue_type,
                      fouts[i], i, n_ppp, n_quad),
                callback=_ok, error_callback=_err,
            )
        pool.close(); pool.join()

    if threads > 1:
        with open(fout, "at") as f:
            for i in range(1, threads):
                p = f"{fout}.r2p{i}"
                try:
                    f.write(open(p).read())
                    Path(p).unlink()
                except FileNotFoundError:
                    pass


# ---------------------------------------------------------------------------
# Phase 6 — Merge R1 + R2
# ---------------------------------------------------------------------------

def merge_results(r1_path: str, r2_path: str, out_path: str) -> None:
    """
    Replace R1 rows with R2 rows for flagged loci.
    All columns from both regimes are preserved; missing values are empty.
    """
    # Key = (chrom, start, end, motif, haplotype)
    def _key(row): return (row["chrom"], row["start"],
                            row["end"],   row["motif"], row["haplotype"])

    r2_rows = {}
    with open(r2_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            r2_rows[_key(row)] = row

    shared  = [c for c in _R1_COLS if c in _R2_COLS]
    r1_only = [c for c in _R1_COLS if c not in _R2_COLS]
    r2_only = [c for c in _R2_COLS if c not in _R1_COLS]
    all_cols = shared + r1_only + r2_only

    n_r1 = n_r2 = 0
    with open(out_path, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=all_cols,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        with open(r1_path) as f:
            for r1_row in csv.DictReader(f, delimiter="\t"):
                k = _key(r1_row)
                if k in r2_rows:
                    out_row = {c: "" for c in all_cols}
                    out_row.update(r2_rows[k])
                    n_r2 += 1
                else:
                    out_row = {c: "" for c in all_cols}
                    out_row.update(r1_row)
                    n_r1 += 1
                writer.writerow(out_row)

    print(f"  Merged: {n_r1} R1 rows + {n_r2} R2 rows → {out_path}",
          file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = _parse_args()

    # Output prefix
    prefix = args.output or str(
        Path(args.input).parent / Path(args.input).stem
    )
    r1_out    = f"{prefix}_r1.tsv"
    r2_out    = f"{prefix}_r2.tsv"
    final_out = f"{prefix}_instability.tsv"

    r2_threads = args.regime2_threads or args.threads

    # Donor age resolution
    donor_age: float | None = args.donor_age
    if args.run_regime2 and args.age_file and donor_age is None:
        sample_id = Path(args.input).stem
        donor_age = load_age_from_file(args.age_file, sample_id)
        if donor_age is None:
            print(f"  Sample '{sample_id}' not in age file — "
                  f"using tissue prior '{args.tissue_type}'", file=sys.stderr)

    # ── Phase 0 ─────────────────────────────────────────────────────────────
    print("\nPhase 0: loading inputs...", file=sys.stderr)
    genotypes = load_genotypes_from_vcf(args.vcf) if args.vcf else None
    catalog   = parse_input_tsv(args.input, genotypes)
    print(f"  {len(catalog)} loci loaded", file=sys.stderr)

    # ── Phase 1 ─────────────────────────────────────────────────────────────
    print("\nPhase 1: Regime 1 global param estimation...", file=sys.stderr)
    gp_map = run_r1_pass1(catalog, min_reads=args.min_reads)

    # ── Phase 2 ─────────────────────────────────────────────────────────────
    print("\nPhase 2: Regime 1 per-locus inference...", file=sys.stderr)
    r1_pbar = tqdm(
        unit="chunks" if args.threads > 1 else "loci",
        total=args.threads if args.threads > 1 else len(catalog),
        unit_scale=True, ncols=80, smoothing=0.1,
        position=2, desc="R1",
    )
    _run_r1(catalog, gp_map, r1_out, args.threads, args.n_ppp, r1_pbar)
    r1_pbar.close()
    print(f"  R1 results → {r1_out}", file=sys.stderr)

    # ── Phase 2B: q recalibration (optional) ────────────────────────────────
    if args.calibrate:
        print("\nPhase 2B: recalibrating q from Phase 2 dispersion ratios...",
              file=sys.stderr)
        gp_map_cal = run_r1_calibration(
            gp_map, r1_out, min_reads=args.calibrate_min_reads
        )
        # Check whether any stratum actually changed
        changed = any(
            abs(gp_map_cal[s].q_e - gp_map[s].q_e) > 0.005
            for s in gp_map
        )
        if changed:
            print("  Re-running Phase 2 with calibrated params...",
                  file=sys.stderr)
            r2_pbar_cal = tqdm(
                unit="chunks" if args.threads > 1 else "loci",
                total=args.threads if args.threads > 1 else len(catalog),
                unit_scale=True, ncols=80, smoothing=0.1,
                position=2, desc="R1 (calibrated)",
            )
            _run_r1(catalog, gp_map_cal, r1_out,
                    args.threads, args.n_ppp, r2_pbar_cal)
            r2_pbar_cal.close()
            gp_map = gp_map_cal   # use calibrated params for all downstream phases
            print(f"  Calibrated R1 results → {r1_out}", file=sys.stderr)
        else:
            print("  q unchanged (<0.005 difference) — skipping rerun.",
                  file=sys.stderr)

    if not args.run_regime2:
        import shutil
        shutil.copy(r1_out, final_out)
        print(f"\nDone. Final output: {final_out}", file=sys.stderr)
        print("(Tip: add --run-regime2 to also run the CTMC model "
              "on heavy-tailed loci.)", file=sys.stderr)
        sys.exit(0)

    # ── Phase 3 ─────────────────────────────────────────────────────────────
    print("\nPhase 3: identifying Regime 2 loci...", file=sys.stderr)
    flagged, cat_idx = identify_regime2_loci(r1_out, catalog)
    print(f"  {len(flagged)} loci flagged for Regime 2", file=sys.stderr)

    if len(flagged) < args.regime2_min_loci:
        import shutil
        shutil.copy(r1_out, final_out)
        print(f"  Fewer than {args.regime2_min_loci} flagged loci — "
              f"final output is R1 only: {final_out}", file=sys.stderr)
        sys.exit(0)

    # ── Phase 4 ─────────────────────────────────────────────────────────────
    print("\nPhase 4: Regime 2 global param estimation...", file=sys.stderr)
    ctmc_gp = run_r2_pass1(
        flagged, cat_idx, donor_age, args.tissue_type,
        args.repeat_preset, args.regime2_min_loci,
    )
    if ctmc_gp is None:
        import shutil
        shutil.copy(r1_out, final_out)
        print(f"  R2 estimation failed — final output is R1 only: {final_out}",
              file=sys.stderr)
        sys.exit(0)

    # ── Phase 5 ─────────────────────────────────────────────────────────────
    print("\nPhase 5: Regime 2 per-locus inference + GoF...", file=sys.stderr)

    r2_items = []
    for key in flagged:
        locus = cat_idx.get(key)
        if not locus:
            continue
        for hap_idx, (hap, lengths) in enumerate(locus["haplotypes"].items()):
            if len(lengths) < args.min_reads:
                continue
            f = _founder_from_genotypes(locus["genotypes"], hap_idx, lengths)
            r2_items.append({
                "chrom":   locus["chrom"],
                "start":   locus["start"],
                "end":     locus["end"],
                "motif":   locus["motif"],
                "hap":     hap,
                "lengths": np.asarray(lengths, dtype=float),
                "founder": f,
            })

    r2_pbar = tqdm(
        unit="chunks", total=r2_threads,
        unit_scale=True, ncols=80, smoothing=0.1,
        position=2, desc="R2",
    )
    _run_r2(
        r2_items, ctmc_gp, donor_age, args.tissue_type,
        r2_out, r2_threads, args.r2_n_ppp, args.r2_n_quad, r2_pbar,
    )
    r2_pbar.close()
    print(f"  R2 results → {r2_out}", file=sys.stderr)

    # ── Phase 6 ─────────────────────────────────────────────────────────────
    print("\nPhase 6: merging results...", file=sys.stderr)
    merge_results(r1_out, r2_out, final_out)
    print(f"\nDone. Final output: {final_out}", file=sys.stderr)

# 8668970161  - Domestic
# 12679411037 - International
# 1004028383
# 8008291040
