#!/usr/bin/env python3
"""
calc_instability.py
-------------------
Genome-wide somatic instability analysis using the NB-Geometric model.

Usage
-----
    python calc_instability.py -i input.tsv -v genotypes.vcf -o output.tsv -t 4

Input TSV columns (tab-separated, one row per read):
    chrom  start  end  motif  read_id  haplotype  length  allele  avg_meth  meth_bases

    - length is in base pairs; the script converts to repeat units (length / motif_size)
    - haplotype is 1 or 2 (integer)

VCF (optional):
    Must contain FORMAT field AL (allele length in base pairs).
    If not provided, the modal observed length per haplotype is used as the founder.

Output TSV columns:
    chrom  start  end  motif  haplotype  founder_length  founder_length_bp
    avg_deviation  n_reads  r_e  r_c  p_plus  p_minus  mu_expansion
    mu_contraction  instability_index  net_bias  p_expansion  expansion_rate
    contraction_rate  net_mutation_rate  p_any_mutation
    ks_statistic  ad_statistic  ppp_variance  ppp_skewness  observed_var
    fitted_var  dispersion_ratio  chisq_statistic  chisq_pvalue  chisq_df
    chisq_n_bins  regime2_recommended  fit_quality  fit_quality_adjusted
"""

import time
import argparse
import sys
import gzip
import math
from pathlib import Path
from multiprocessing import Pool

import cyvcf2
import numpy as np
from tqdm import tqdm

from nbgeom_modelling import (
    _modal_length,
    estimate_regime1_global,
    Regime1GlobalParams,
    auto_fit,
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    p = argparse.ArgumentParser(
        description="NB-Geometric somatic instability analysis"
    )
    p.add_argument(
        "-i", "--input", required=True,
        help="Input TSV file with per-read allele lengths for catalog loci"
    )
    p.add_argument(
        "-v", "--vcf", required=False,
        help="VCF file with genotype information (FORMAT/AL field). "
             "If not provided, the modal observed length is used as founder."
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output TSV file path. Defaults to <input_stem>_instability.tsv"
    )
    p.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of parallel worker processes (default: 1)"
    )
    p.add_argument(
        "--n-ppp", type=int, default=300,
        help="PPP simulations per locus for goodness-of-fit (default: 300). "
             "Set to 0 to skip GoF computation entirely (faster)."
    )
    p.add_argument(
        "--min-reads", type=int, default=10,
        help="Minimum reads per haplotype (default: 10)"
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# VCF loading
# ---------------------------------------------------------------------------

def load_genotypes_from_vcf(vcf_path: str) -> dict:
    """
    Load per-locus genotype allele lengths from a VCF file.

    Keys are formatted as "CHROM:start-END_MOTIF".
    Values are numpy arrays of allele lengths in repeat units
    (base pairs divided by motif length from the AL FORMAT field).

    Only PASS variants are included.
    """
    start_time = time.time()
    genotypes   = {}
    reader      = cyvcf2.VCF(vcf_path)

    for variant in reader:
        if variant.FILTER not in ("PASS", None):
            continue
        motif = variant.INFO.get("MOTIF", "")
        if not motif:
            continue
        motif_len = len(motif)
        key       = f"{variant.CHROM}:{variant.start}-{variant.INFO.get('END')}_{motif}"
        # AL is in base pairs; divide by motif length to get repeat units
        al = variant.format("AL")
        if al is not None and len(al) > 0:
            genotypes[key] = al[0] / motif_len

    print(
        f"Loaded genotypes for {len(genotypes)} loci "
        f"in {time.time() - start_time:.2f}s",
        file=sys.stderr,
    )
    return genotypes


# ---------------------------------------------------------------------------
# Helper: average absolute deviation from founder length
# ---------------------------------------------------------------------------

def average_median_delta(alens, founder_length=None):
    """
    Mean absolute deviation of read lengths from the founder/reference length.

    Parameters
    ----------
    alens : array-like
        Per-read allele lengths in repeat units.
    founder_length : float, optional
        Reference length. If None, the modal observed length is used.

    Returns
    -------
    float : mean |length - founder_length|
    """
    lengths = np.asarray(alens, dtype=float)
    # BUG FIX: was `if founder_length is not None` which always overwrote the
    # supplied founder_length with the modal length. Corrected to `is None`.
    if founder_length is None:
        founder_length = float(_modal_length(lengths))
    deltas = np.abs(lengths - founder_length)
    return float(np.mean(deltas))


# ---------------------------------------------------------------------------
# Null GoF dict — used when auto_fit returns gof=None
# ---------------------------------------------------------------------------

_NULL_GOF = {
    "ks_statistic":         float("nan"),
    "ad_statistic":         float("nan"),
    "ppp_variance":         float("nan"),
    "ppp_skewness":         float("nan"),
    "observed_var":         float("nan"),
    "fitted_var":           float("nan"),
    "dispersion_ratio":     float("nan"),
    "chisq_statistic":      float("nan"),
    "chisq_pvalue":         float("nan"),
    "chisq_df":             0,
    "chisq_n_bins":         0,
    "regime2_recommended":  float("nan"),
    "fit_quality":          float("nan"),
    "fit_quality_adjusted": float("nan"),
}

# Column order for the output TSV
_HEADER = [
    "chrom", "start", "end", "motif", "haplotype",
    "founder_length",      # in repeat units
    "founder_length_bp",   # in base pairs (founder_length × motif_size)
    "avg_deviation",       # mean |read_length - founder_length| in repeat units
    "n_reads",
    # NegBin parameters (posterior mean only)
    "r_e", "r_c", "p_plus", "p_minus",
    "mu_expansion", "mu_contraction",
    # Instability indices (posterior mean only)
    "instability_index", "net_bias", "p_expansion",
    # Mutation rates (posterior mean only)
    "expansion_rate", "contraction_rate", "net_mutation_rate", "p_any_mutation",
    # Goodness of fit
    "ks_statistic", "ad_statistic", "ppp_variance", "ppp_skewness",
    "observed_var", "fitted_var", "dispersion_ratio",
    "chisq_statistic", "chisq_pvalue", "chisq_df", "chisq_n_bins",
    "regime2_recommended", "fit_quality", "fit_quality_adjusted",
    # Diagnostic info
    "diag_n_reads", "diag_founder_length", "diag_modal_length",
    "diag_n_outliers", "diag_outlier_lengths",
    "diag_delta_mean", "diag_delta_std", "diag_delta_max", "diag_delta_min",
    "diag_tail_ratio", "diag_skewness", "diag_is_heavy_tailed", "diag_recommended_regime"
]


# ---------------------------------------------------------------------------
# Worker function (runs in a subprocess)
# ---------------------------------------------------------------------------

def _worker_function(
    chunk_data:              list,
    global_params_by_stratum: dict,
    fout:                    str,
    thread_idx:              int,
    n_ppp:                   int,
):
    """
    Process a chunk of loci and write results to fout.

    Parameters
    ----------
    chunk_data : list of locus dicts
    global_params_by_stratum : {stratum_int: Regime1GlobalParams}
    fout : output file path for this worker
    thread_idx : 0-based worker index (thread 0 writes the header)
    n_ppp : PPP simulations per locus (0 = skip GoF)
    """
    out = open(fout, "wt")

    # Only the first worker writes the header
    if thread_idx == 0:
        print("\t".join(_HEADER), file=out)

    for locus in chunk_data:
        chrom     = locus["chrom"]
        start     = locus["start"]
        end       = locus["end"]
        motif     = locus["motif"]
        motif_len = len(motif)
        stratum   = get_stratum(motif_len)
        gp        = global_params_by_stratum[stratum]
        genotypes = locus["genotypes"]   # repeat units (already divided by motif_len)

        for hap_idx, (hap_key, lengths) in enumerate(locus["haplotypes"].items()):

            # ── Founder length ──────────────────────────────────────────────
            try:
                # genotypes is a numpy array of allele lengths in repeat units
                founder = float(genotypes[hap_idx])
            except (IndexError, TypeError, ValueError):
                # Fallback: modal observed length
                lengths_int = np.round(np.array(lengths)).astype(int)
                vals, cts   = np.unique(lengths_int, return_counts=True)
                founder     = float(vals[np.argmax(cts)])

            # ── Fit ─────────────────────────────────────────────────────────
            try:
                if start == 10440:
                    print(lengths, hap_key)
                result, diag, gof = auto_fit(
                    lengths,
                    founder_length  = founder,
                    regime1_params  = gp,
                    haplotype_label = f"hap{hap_key}",
                    n_ppp           = n_ppp,
                )
            except Exception as e:
                print(
                    f"ERROR fitting {chrom}:{start}-{end} hap{hap_key}: {e}",
                    file=sys.stderr,
                )
                continue

            # ── Guard: gof=None when Regime 2 or fit failed ─────────────────
            if gof is None:
                gof = dict(_NULL_GOF)
                # Use the diagnostic heavy-tail flag as regime2 signal
                gof["regime2_recommended"] = int(diag.is_heavy_tailed)

            # ── Average deviation ────────────────────────────────────────────
            avg_dev = average_median_delta(lengths, founder_length=result.founder_length)

            # ── Write row ────────────────────────────────────────────────────
            # result.r_e etc. are (mean, lo95, hi95) tuples — use [0] for mean
            row = [
                chrom, start, end, motif, hap_key,
                round(result.founder_length, 4),
                round(result.founder_length * motif_len, 1),  # founder in bp
                round(avg_dev, 4),
                result.n_reads,
                # NB parameters — posterior mean only
                round(result.r_e[0],           4),
                round(result.r_c[0],           4),
                round(result.p_plus[0],        4),
                round(result.p_minus[0],       4),
                round(result.mu_expansion[0],  4),
                round(result.mu_contraction[0],4),
                # Instability
                round(result.instability_index[0], 4),
                round(result.net_bias[0],          4),
                round(result.p_expansion[0],       4),
                # Mutation rates
                round(result.expansion_rate[0],    4),
                round(result.contraction_rate[0],  4),
                round(result.net_mutation_rate[0], 4),
                round(result.p_any_mutation[0],    4),
                # GoF
                gof["ks_statistic"],
                gof["ad_statistic"],
                gof["ppp_variance"],
                gof["ppp_skewness"],
                gof["observed_var"],
                gof["fitted_var"],
                gof["dispersion_ratio"],
                gof["chisq_statistic"],
                gof["chisq_pvalue"],
                gof["chisq_df"],
                gof["chisq_n_bins"],
                gof["regime2_recommended"],
                gof["fit_quality"],
                gof["fit_quality_adjusted"],
                diag.n_reads,
                diag.founder_length,
                diag.modal_length,
                diag.n_outliers,
                diag.outlier_lengths,
                diag.delta_mean,
                diag.delta_std,
                diag.delta_max,
                diag.delta_min,
                diag.tail_ratio,
                diag.skewness,
                diag.is_heavy_tailed,
                diag.recommended_regime
            ]
            print(*row, sep="\t", file=out)

    out.close()


# ---------------------------------------------------------------------------
# Multiprocessing coordinator
# ---------------------------------------------------------------------------

def _multiprocess(
    data:                    list,
    global_params_by_stratum: dict,
    fout:                    str,
    threads:                 int,
    n_ppp:                   int,
    progress,
):
    """
    Split data across worker processes and merge results.

    Each worker writes to its own file; thread 0 writes to fout directly,
    threads 1..N write to fout.partN and are merged afterwards.
    """
    chunks = np.array_split(data, threads)
    # Thread 0 → main output file; others → part files
    fouts  = [fout] + [f"{fout}.part{i}" for i in range(1, threads)]

    def on_success(_):
        progress.update()

    def on_error(e):
        progress.update()
        print(f"Worker failed: {e}", file=sys.stderr)

    # BUG FIX: was `args.threads` (NameError); corrected to `threads` parameter
    with Pool(processes=threads) as pool:
        for thread_idx in range(threads):
            pool.apply_async(
                _worker_function,
                args=(
                    chunks[thread_idx].tolist(),
                    global_params_by_stratum,
                    fouts[thread_idx],
                    thread_idx,
                    n_ppp,
                ),
                callback       = on_success,
                error_callback = on_error,
            )
        pool.close()
        pool.join()

    # Merge part files into the main output file (append, no re-header)
    if threads > 1:
        with open(fout, "at") as out:
            for i in range(1, threads):
                part_path = f"{fout}.part{i}"
                try:
                    with open(part_path, "rt") as part:
                        for line in part:
                            out.write(line)
                    Path(part_path).unlink()
                except FileNotFoundError:
                    pass   # worker may have produced no output (empty chunk)


# ---------------------------------------------------------------------------
# Pass 1 helpers
# ---------------------------------------------------------------------------

def build_pass1_data(
    catalog: list,
    min_reads: int = 10,
    outlier_mad_factor: float = 5.0,
) -> list:
    """
    Convert locus dicts into the list-of-dicts format expected by
    estimate_regime1_global().

    Each output dict has 'deltas', 'founder_length', and 'motif_length'.
    Outlier reads (> outlier_mad_factor × MAD from median) are removed
    before computing deltas to prevent artefacts from skewing the
    step-size regression.
    """
    loci_data = []

    for locus in catalog:
        motif_len  = len(locus["motif"])
        genotypes  = locus["genotypes"]    # repeat units (already divided)
        haplotypes = locus["haplotypes"]   # {hap_key: [lengths in repeat units]}

        for hap_idx, (hap_key, lengths_raw) in enumerate(haplotypes.items()):
            lengths = np.array(lengths_raw, dtype=float)

            if len(lengths) < min_reads:
                continue

            # Founder length in repeat units
            try:
                founder_units = float(genotypes[hap_idx])
            except (IndexError, TypeError, ValueError):
                lengths_int   = np.round(lengths).astype(int)
                vals, cts     = np.unique(lengths_int, return_counts=True)
                founder_units = float(vals[np.argmax(cts)])

            # Outlier removal: reads > outlier_mad_factor × MAD from median
            median_l = float(np.median(lengths))
            mad      = float(np.median(np.abs(lengths - median_l)))
            if mad < 1e-6:
                mad = float(lengths.std()) + 1e-6
            keep          = np.abs(lengths - median_l) <= outlier_mad_factor * mad
            lengths_clean = lengths[keep]

            if len(lengths_clean) < min_reads:
                continue

            deltas = (np.round(lengths_clean) - round(founder_units)).astype(int)

            loci_data.append({
                "deltas":         deltas,
                "founder_length": founder_units,
                "locus_id":       f"{locus['chrom']}_{locus['start']}",
                "haplotype":      hap_key,
                "motif_length":   motif_len,
            })

    return loci_data


def get_stratum(motif_length: int) -> int:
    """Map motif length in bp to the nearest stratum key (3, 6, or 7)."""
    if motif_length <= 3:   return 3
    elif motif_length <= 6: return 6
    else:                   return 7


def run_pass1(
    catalog:           list,
    min_reads:         int  = 10,
    stratify_by_motif: bool = True,
) -> dict:
    """
    Estimate global NB-Geometric step-size parameters from the catalog.

    Returns a dict mapping stratum key → Regime1GlobalParams.
    When stratify_by_motif=False, returns {0: GlobalParams}.
    """
    loci_data = build_pass1_data(catalog, min_reads=min_reads)
    print(f"Pass 1: {len(loci_data)} haplotype-loci after filtering", file=sys.stderr)

    if not stratify_by_motif:
        gp = estimate_regime1_global(loci_data, min_reads=min_reads)
        print(gp, file=sys.stderr)
        return {0: gp}

    # Stratify into short / medium / long repeat strata
    strata: dict = {}
    for entry in loci_data:
        s = get_stratum(entry["motif_length"])
        strata.setdefault(s, []).append(entry)

    global_params = {}
    labels = {3: "1-3bp (short)", 6: "4-6bp (medium)", 7: "7+bp (long)"}
    for stratum, entries in sorted(strata.items()):
        print(f"\nStratum {labels[stratum]}: {len(entries)} haplotype-loci",
              file=sys.stderr)
        if len(entries) < 4:
            print("  Too few loci — using defaults", file=sys.stderr)
            global_params[stratum] = Regime1GlobalParams.from_defaults(
                repeat_unit_length=stratum
            )
        else:
            gp = estimate_regime1_global(entries, min_reads=min_reads)
            print(gp, file=sys.stderr)
            global_params[stratum] = gp

    return global_params


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = _parse_args()

    # ── Default output path ─────────────────────────────────────────────────
    if args.output is None:
        args.output = str(
            Path(args.input).parent / (Path(args.input).stem + "_instability.tsv")
        )

    # ── Load VCF genotypes ──────────────────────────────────────────────────
    genotypes = load_genotypes_from_vcf(args.vcf) if args.vcf else None

    # ── Progress bars ───────────────────────────────────────────────────────
    inread_pbar = tqdm(
        unit="rows", unit_scale=True, ncols=80, smoothing=0.1,
        position=1, desc="Reading input",
    )

    if args.threads == 1:
        main_pbar = tqdm(
            unit="loci", unit_scale=True, ncols=80, smoothing=0.1,
            position=3, desc="Processing loci",
        )
        if genotypes is not None:
            main_pbar.total = len(genotypes)
    else:
        main_pbar = tqdm(
            unit="chunks", unit_scale=True, ncols=80, smoothing=0.1,
            total=args.threads, position=3, desc="Processing chunks",
        )

    # ── Read input file ─────────────────────────────────────────────────────
    ins_fh = (
        gzip.open(args.input, "rt")
        if args.input.endswith(".gz")
        else open(args.input, "rt")
    )

    data       = []
    prev_key   = None
    haplotypes: dict = {}
    info:       dict = {}

    for line in ins_fh:
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        chrom, start, end, motif, read_id, haplotype, length, allele, avg_meth, meth_bases = parts

        start        = int(start)
        end          = int(end)
        motif_length = len(motif)
        repeat_length = end - start
        inread_pbar.update(1)

        # Filter out very short/long spans and long motifs (Regime 1 only)
        if motif_length > 6 or repeat_length < 20 or repeat_length > 500:
            continue

        locus_key = f"{chrom}:{start}-{end}_{motif}"

        # Store locus metadata on first encounter
        if locus_key not in info:
            info[locus_key] = {
                "chrom": chrom,
                "start": start,
                "end":   end,
                "motif": motif,
            }

        # Accumulate read lengths in repeat units
        haplotype = int(haplotype)
        units     = float(length) / motif_length

        # When the locus key changes, flush the previous locus
        if locus_key != prev_key and prev_key is not None:
            if haplotypes:
                data.append({
                    "chrom":      info[prev_key]["chrom"],
                    "start":      info[prev_key]["start"],
                    "end":        info[prev_key]["end"],
                    "motif":      info[prev_key]["motif"],
                    "haplotypes": haplotypes,
                    "genotypes":  (
                        genotypes[prev_key]
                        if genotypes and prev_key in genotypes
                        else None
                    ),
                })
            haplotypes = {}
            # BUG FIX: delete info entry AFTER using it, not before the final append
            del info[prev_key]
        haplotypes.setdefault(haplotype, []).append(units)
        prev_key = locus_key

    # Flush the final locus
    if haplotypes and prev_key is not None:
        data.append({
            "chrom":      info[prev_key]["chrom"],
            "start":      info[prev_key]["start"],
            "end":        info[prev_key]["end"],
            "motif":      info[prev_key]["motif"],
            "haplotypes": haplotypes,
            "genotypes":  (
                genotypes[prev_key]
                if genotypes and prev_key in genotypes
                else None
            ),
        })

    inread_pbar.close()
    ins_fh.close()
    print(f"\nLoaded {len(data)} loci for analysis", file=sys.stderr)

    # ── Pass 1: estimate global step-size parameters ────────────────────────
    global_params_by_stratum = run_pass1(data, min_reads=args.min_reads)

    # ── Pass 2: per-locus instability inference ─────────────────────────────
    _multiprocess(
        data,
        global_params_by_stratum,
        args.output,
        args.threads,
        args.n_ppp,
        main_pbar,
    )

    main_pbar.close()
    print(f"\nResults written to {args.output}", file=sys.stderr)