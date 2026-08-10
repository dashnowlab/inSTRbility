#!/usr/bin/env python3
# Prevent OpenBLAS fork-deadlock: must be set BEFORE numpy/scipy import.
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")


import argparse
import sys
import time
import traceback
import copy
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from tqdm import tqdm

from src.model.nbgeom_modelling   import (auto_fit, estimate_nbgeom_global, NBGeomGlobalParams)
from src.model.parse_inputs       import parse_input_tsv, load_genotypes_from_vcf
from src.model.calibrate          import recalibrate, recalibrate_distinct
from src.model.suggest_re_rc_grid import suggest_grid_ceiling
from src.model.utils              import get_stratum, STRATUM_LABELS
from src.model.write              import OUT_COLS


def model_parser(subparsers) -> argparse.Namespace:
    p = subparsers.add_parser(
        "model",
        description="NB-Geometric + CTMC somatic TR instability analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
        help="NB-Geometric + CTMC somatic TR instability analysis",
    )

    p.add_argument("-i", "--input", required=True,                  help="Per-read TSV (chrom start end motif read_id hap "
                                                                         "length_bp allele avg_meth meth_bases)")
    p.add_argument("-v", "--vcf",           type=str, default=None, help="VCF with FORMAT/AL germline allele lengths (bp). "
                                                                         "Omit to use modal observed length as founder.")
    p.add_argument("-o", "--output",        type=str, default=None, help="Output prefix. Default: <input_stem>_instability")
    p.add_argument("-t", "--threads",       type=int, default=1,    help="Worker processes (default: 1)")

    p.add_argument("--min-reads",           type=int, default=10,   help="Min reads per haplotype (default: 10)")
    p.add_argument("--n-ppp",               type=int, default=300,  help="PPP simulations for R1 GoF (default: 300; 0=skip)")

    p.add_argument("--calibrate-min-reads", type=int, default=20,   help="Min reads per locus for calibration (default: 20)")

    p.set_defaults(func=run_model)


def founder_from_genotypes(genotypes, hap_idx, lengths):
    """
    Return founder length in repeat units from VCF genotypes or modal.
    
    @param genotypes list of founder lengths (bp) from VCF, or None
    @param hap_idx index of haplotype to use
    @param lengths list of observed lengths (repeat units) for this haplotype
    @return founder length in repeat units (float)
    """
    try:
        return float(genotypes[hap_idx])
    except (IndexError, TypeError, ValueError):
        li   = np.round(np.asarray(lengths)).astype(int)
        v, c = np.unique(li, return_counts=True)
        return float(v[np.argmax(c)])


def build_pass1_data(catalog, min_reads=10):
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

            founder_length = founder_from_genotypes(locus["genotypes"], hap_idx, arr)

            out.append({
                "deltas":         (np.round(arr) - round(founder_length)).astype(int),
                "founder_length": founder_length,
                "locus_id":       f"{locus['chrom']}-{locus['start']}-{locus['end']}_{locus['motif']}",
                "haplotype":      hap,
                "motif_length":   ml,
            })

    return out


def run_pass1(catalog, min_reads = 10):
    """
    Pass 1: estimate global step-size params, stratified by motif.

    @param catalog list of locus dicts from parse_input_tsv()
    @param min_reads minimum reads per haplotype to include in Phase 1
    @return dict of {stratum: NBGeomGlobalParams}
    """

    loci_data = build_pass1_data(catalog, min_reads)
    print(f"  R1 Pass 1: {len(loci_data)} haplotype-loci", file=sys.stderr)

    strata = {}
    for locus in loci_data:
        s = get_stratum(locus["motif_length"])
        strata.setdefault(s, []).append(locus)

    gp_map = {}
    for s, entries in sorted(strata.items()):
        print(f"    Stratum {STRATUM_LABELS.get(s, f'{s}bp')}: {len(entries)} haplotype-loci", file=sys.stderr)
        if len(entries) < 4:
            gp_map[s] = NBGeomGlobalParams.from_defaults(s)
        else:
            gp_map[s] = estimate_nbgeom_global(entries, min_reads=min_reads)
            print(f"    {gp_map[s]}", file=sys.stderr)
    return gp_map


def r1_worker(chunk, gp_map, fout, thread_id, n_ppp) -> None:
    """
    Process a chunk of loci with R1 and write TSV rows.

    @param chunk list of locus dicts
    @param gp_map dict of {stratum: NBGeomGlobalParams}
    @param fout output TSV file
    @param thread_id thread ID
    @param n_ppp number of permutations per parameter
    """

    out = open(fout, "wt")
    if thread_id == 0:
        print("\t".join(OUT_COLS), file=out)

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
            founder = founder_from_genotypes(gt, hap_idx, lengths)

            # try:
            result, gof = auto_fit(ml, lengths, founder, gp, n_ppp=n_ppp)

            # except Exception as e:
            #     print(f"WARN R1 {chrom}:{start} hap{hap}: {e}",
            #           file=sys.stderr)
            #     continue

            if gof is None:
                gof = dict(null_gof)

            avg   = float(np.mean(np.abs(np.asarray(lengths, float) - result.founder_length)))

            row = [
                chrom, start, end, motif, hap,
                round(result.founder_length, 4),
                round(result.founder_length * ml, 1),
                result.n_reads,
                round(result.r_e[0],            4),
                round(result.r_c[0],            4),
                round(result.q_e,   4),
                round(result.q_c, 4),
                round(result.p_e[0],         4),
                round(result.p_c[0],        4),
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
                gof["fit_quality"],     gof["fit_quality_adjusted"],
                round(time.time() - t0, 3),
            ]
            print(*row, sep="\t", file=out)

    out.close()


def run_r1(catalog, gp_map, fout, threads, n_ppp, pbar):
    """
    Run Regime 1 per-locus inference in parallel and write TSV output.

    @param catalog list of locus dicts
    @param gp_map dict of {stratum: NBGeomGlobalParams}
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
                r1_worker,
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


def run_model(args):
    """
    Run inSTRbility modelling workflow: Phase 0 → Phase 1 → Phase 2 → Phase 3 recalibration
    """

    # Output file
    output = args.output
    if output is None:
        prefix = str(Path(args.input).parent / Path(args.input).stem)
        output = f"{prefix}_instrbility.tsv"
    else:
        output = args.output

    # ── Phase 0 ─────────────────────────────────────────────────────────────
    print("\nPhase 0: loading inputs...", file=sys.stderr)
    genotypes = load_genotypes_from_vcf(args.vcf) if args.vcf else None
    catalog   = parse_input_tsv(args.input, genotypes)
    print(f"  {len(catalog)} loci loaded", file=sys.stderr)

    # ── Phase 1 ─────────────────────────────────────────────────────────────
    print("\nPhase 1: Regime 1 global param estimation...", file=sys.stderr)
    gp_map = run_pass1(catalog, min_reads=args.min_reads)

    # rec = suggest_grid_ceiling(catalog, gp_map, percentile=99, margin=1.5)
    # for label in rec:
    #     for stratum in rec[label]:
    #         n = 0
    #         if rec[label][stratum]["suggested_ceiling"] > 25: n = 100
    #         elif rec[label][stratum]["suggested_ceiling"] > 10: n = 100
    #         elif rec[label][stratum]["suggested_ceiling"] > 0: n = 100
    #         print(label, stratum, rec[label][stratum]["suggested_ceiling"], n, sep="\t")
    #         if label == "r_e":
    #             RE_GRID[stratum] = np.exp(np.linspace(np.log(0.01), np.log(rec[label][stratum]["suggested_ceiling"]), n))
    #         else:
    #             RC_GRID[stratum] = np.exp(np.linspace(np.log(0.01), np.log(rec[label][stratum]["suggested_ceiling"]), n))

    # print("Grids for R1 inference:")
    # for stratum in sorted(RE_GRID):
    #     print(f"  Stratum {STRATUM_LABELS.get(stratum, f'{stratum}bp')}: r_e grid = {min(RE_GRID[stratum])} - {max(RE_GRID[stratum])}")
    # for stratum in sorted(RC_GRID):
    #     print(f"  Stratum {STRATUM_LABELS.get(stratum, f'{stratum}bp')}: r_c grid = {min(RC_GRID[stratum])} - {max(RC_GRID[stratum])}")

    # ── Phase 2 ─────────────────────────────────────────────────────────────
    print("\nPhase 2: Regime 1 per-locus inference...", file=sys.stderr)
    r1_pbar = tqdm(
        unit="chunks" if args.threads > 1 else "loci",
        total=args.threads if args.threads > 1 else len(catalog),
        unit_scale=True, ncols=80, smoothing=0.1,
        position=2, desc="R1",
    )
    run_r1(catalog, gp_map, output, args.threads, args.n_ppp, r1_pbar)
    r1_pbar.close()
    print(f"  R1 results → {output}", file=sys.stderr)

    # ── Phase 2B: q recalibration (optional) ────────────────────────────────
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
        run_r1(catalog, gp_map_cal, output, args.threads, args.n_ppp, r2_pbar_cal)
        r2_pbar_cal.close()
        gp_map = gp_map_cal   # use calibrated params for all downstream phases
        print(f"  Calibrated R1 results → {output}", file=sys.stderr)
    else:
        print("  q unchanged (<0.005 difference) — skipping rerun.", file=sys.stderr)
