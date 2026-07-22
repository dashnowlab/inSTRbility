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

Six-phase workflow
------------------
  Phase 0  Parse TSV reads + VCF founder lengths
  Phase 1  Regime 1 Pass 1  — estimate global step-size params
  Phase 2  Regime 1 Pass 2  — per-locus inference  (parallelised)

Usage
-----
  # Regime 1 only
  python calc_instability.py -i reads.tsv -v genotypes.vcf -o out.tsv -t 8

Input TSV (one row per read, tab-separated):
  chrom  start  end  motif  read_id  haplotype  length_bp  allele  avg_meth  meth_bases

VCF (optional):
  FORMAT/AL field in base pairs; divided internally by motif length.
  If absent, modal observed length is used as the founder allele.
"""

import argparse
import gzip
import sys
import time
import traceback
import copy
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from tqdm import tqdm

from nbgeom_modelling   import (auto_fit, estimate_regime1_global, Regime1GlobalParams)
from parse_inputs       import parse_input_tsv, load_genotypes_from_vcf
from calibrate          import recalibrate, recalibrate_distinct
from suggest_re_rc_grid import suggest_grid_ceiling


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
        li   = np.round(np.asarray(lengths)).astype(int)
        v, c = np.unique(li, return_counts=True)
        return float(v[np.argmax(c)])


# ---------------------------------------------------------------------------
# Null GoF sentinel (when fit fails or GoF is skipped)
# ---------------------------------------------------------------------------

_NULL_GOF: dict = {
    "ks_stat": float("nan"), "ad_stat": float("nan"),
    "ppp_var": float("nan"), "dispersion_ratio": float("nan"),
    "chisq_p": float("nan"), "fit_quality": float("nan"),
    "fit_quality_adj": float("nan"), "regime2_recommended": float("nan"),
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


# ---------------------------------------------------------------------------
# Phase 1 — Regime 1 global parameter estimation
# ---------------------------------------------------------------------------

def _build_pass1_data(catalog: list, min_reads: int = 10) -> list:
    """
    Build a list of per-haplotype-locus dicts for Phase 1 global parameter estimation.
    Each dict: {deltas, founder_length, locus_id, haplotype, motif_length}

    @param catalog list of locus dicts from parse_input_tsv()
    @param min_reads minimum reads per haplotype to include in Phase 1

    """
    out = []
    for locus in catalog:
        ml = len(locus["motif"])
        for hap_idx, (hap, lengths) in enumerate(locus["haplotypes"].items()):
            arr = np.asarray(lengths, dtype=float)

            # coverage filter: skip haplotype-loci with too few reads
            if len(arr) < min_reads: continue

            founder = _founder_from_genotypes(locus["genotypes"], hap_idx, arr)
            # MAD outlier removal
            median = float(np.median(arr))
            mad    = max(float(np.median(np.abs(arr - median))), float(arr.std()) + 1e-6)
            arr = arr[np.abs(arr - median) <= 5.0 * mad]

            if len(arr) < min_reads: continue

            out.append({
                "deltas":         (np.round(arr) - round(founder)).astype(int),
                "founder_length": founder,
                "locus_id":       f"{locus['chrom']}_{locus['start']}",
                "haplotype":      hap,
                "motif_length":   ml,
            })

    return out


def run_r1_pass1(catalog: list, min_reads: int = 10) -> dict:
    """
    Regime 1 Pass 1: estimate global step-size params, stratified by motif.

    @param catalog list of locus dicts from parse_input_tsv()
    @param min_reads minimum reads per haplotype to include in Phase 1
    @return dict of {stratum: Regime1GlobalParams}
    """

    loci_data = _build_pass1_data(catalog, min_reads)
    print(f"  R1 Pass 1: {len(loci_data)} haplotype-loci", file=sys.stderr)

    strata: dict = {}
    for e in loci_data:
        s = get_stratum(e["motif_length"])
        strata.setdefault(s, []).append(e)

    gp_map: dict = {}
    labels = {3: "1-3 bp", 6: "4-6 bp", 7: "7+ bp"}
    for s, entries in sorted(strata.items()):
        print(f"    Stratum {labels[s]}: {len(entries)} haplotype-loci", file=sys.stderr)
        if len(entries) < 4:
            gp_map[s] = Regime1GlobalParams.from_defaults(s)
        else:
            gp_map[s] = estimate_regime1_global(entries, min_reads=min_reads)
            print(f"    {gp_map[s]}", file=sys.stderr)
    return gp_map


# ---------------------------------------------------------------------------
# Phase 2 — Regime 1 per-locus worker
# ---------------------------------------------------------------------------

def _r1_worker(chunk: list,
               gp_map: dict,
               fout: str,
               thread_id: int,
               n_ppp: int) -> None:
    """
    Process a chunk of loci with R1 and write TSV rows.

    @param chunk list of locus dicts
    @param gp_map dict of {stratum: Regime1GlobalParams}
    @param fout output TSV file
    @param thread_id thread ID
    @param n_ppp number of permutations per parameter
    """

    out = open(fout, "wt")
    if thread_id == 0:
        print("\t".join(_R1_COLS), file=out)

    null_gof = {
        "ks_statistic":     float("nan"),   "ad_statistic":     float("nan"),
        "ppp_variance":     float("nan"),   "ppp_skewness":     float("nan"),
        "observed_var":     float("nan"),   "fitted_var":       float("nan"),
        "dispersion_ratio": float("nan"),   "chisq_statistic":  float("nan"),
        "chisq_pvalue":     float("nan"),   "chisq_df": 0, "chisq_n_bins": 0, 
        "fit_quality": float("nan"),        "fit_quality_adjusted": float("nan"),
    }

    for locus in chunk:
        chrom, start, end, motif = (locus["chrom"], locus["start"], locus["end"], locus["motif"])

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

            avg   = float(np.mean(np.abs(np.asarray(lengths, float) - result.founder_length)))

            row = [
                chrom, start, end, motif, hap, "R1",
                round(result.founder_length, 4),
                round(result.founder_length * ml, 1),
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


def _run_r1(catalog: list,
            gp_map:  dict,
            fout:    str,
            threads: int,
            n_ppp:   int,
            pbar):
    """
    Run Regime 1 per-locus inference in parallel and write TSV output.

    @param catalog list of locus dicts
    @param gp_map dict of {stratum: Regime1GlobalParams}
    @param fout output TSV file
    @param threads number of worker processes
    @param n_ppp number of permutations per parameter
    @param pbar tqdm progress bar
    """

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


if __name__ == "__main__":
    args = _parse_args()

    # Output prefix
    prefix = args.output or str(Path(args.input).parent / Path(args.input).stem)
    output = f"{prefix}_instr.tsv"

    # ── Phase 0 ─────────────────────────────────────────────────────────────
    print("\nPhase 0: loading inputs...", file=sys.stderr)
    genotypes = load_genotypes_from_vcf(args.vcf) if args.vcf else None
    catalog   = parse_input_tsv(args.input, genotypes)
    print(f"  {len(catalog)} loci loaded", file=sys.stderr)

    # ── Phase 1 ─────────────────────────────────────────────────────────────
    print("\nPhase 1: Regime 1 global param estimation...", file=sys.stderr)
    gp_map = run_r1_pass1(catalog, min_reads=args.min_reads)

    # rec = suggest_grid_ceiling(catalog, gp_map, percentile=99, margin=1.5)
    # for label in rec:
    #     print(label)
    #     for stratum in rec[label]:
    #         print(f"  Stratum: {stratum}")
    #         for key in sorted(list(rec[label][stratum].keys())):
    #             print(f"    {key}:\t{rec[label][stratum][key]}")
    # sys.exit(1)

    # ── Phase 2 ─────────────────────────────────────────────────────────────
    print("\nPhase 2: Regime 1 per-locus inference...", file=sys.stderr)
    r1_pbar = tqdm(
        unit="chunks" if args.threads > 1 else "loci",
        total=args.threads if args.threads > 1 else len(catalog),
        unit_scale=True, ncols=80, smoothing=0.1,
        position=2, desc="R1",
    )
    _run_r1(catalog, gp_map, output, args.threads, args.n_ppp, r1_pbar)
    r1_pbar.close()
    print(f"  R1 results → {output}", file=sys.stderr)

    # ── Phase 2B: q recalibration (optional) ────────────────────────────────
    if args.calibrate:
        print("\nPhase 2B: recalibrating q from Phase 2 dispersion ratios...", file=sys.stderr)

        gp_map_cal = copy.deepcopy(gp_map)
        for stratum in gp_map:
            gp_map_cal[stratum] = recalibrate_distinct(
                results_tsv = output,
                loci_data   = catalog,
                base_params = gp_map[stratum],
                min_reads   = args.calibrate_min_reads,
                verbose     = True,
            )

            if gp_map_cal[stratum] is False:
                gp_map_cal = recalibrate(gp_map, output, min_reads=args.calibrate_min_reads, strata=[stratum])

        # Check whether any stratum actually changed
        changed = any(abs(gp_map_cal[s].q_e - gp_map[s].q_e) > 0.005 for s in gp_map)
        if changed:
            print("  Re-running Phase 2 with calibrated params...", file=sys.stderr)
            r2_pbar_cal = tqdm(
                unit="chunks" if args.threads > 1 else "loci",
                total=args.threads if args.threads > 1 else len(catalog),
                unit_scale=True, ncols=80, smoothing=0.1,
                position=2, desc="R1 (calibrated)",
            )
            _run_r1(catalog, gp_map_cal, output, args.threads, args.n_ppp, r2_pbar_cal)
            r2_pbar_cal.close()
            gp_map = gp_map_cal   # use calibrated params for all downstream phases
            print(f"  Calibrated R1 results → {output}", file=sys.stderr)
        else:
            print("  q unchanged (<0.005 difference) — skipping rerun.", file=sys.stderr)
