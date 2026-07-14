#!/usr/bin/env python3
"""
analyze_locus_across_samples.py — Single-locus instability, fit per sample
============================================================================

Given ONE locus and a manifest of multiple samples, this script:

  1. Streams each sample's per-read TSV and pulls out ONLY the reads at the
     target locus (fast — doesn't build a genome-wide catalog).
  2. Gets that sample's founder (germline) length from its VCF, or falls
     back to the modal observed length if no VCF is given.
  3. Fits Regime 1 (NB-Geometric) independently for each sample, using a
     SHARED set of global step-size parameters (Regime1GlobalParams) so
     the per-sample fits are directly comparable — no genome-wide Pass 1
     is re-run per sample.
  4. For any sample flagged `gof_regime2_recommended` (heavy-tailed,
     Regime 1 doesn't fit well) — IF --run-regime2 is given — continues
     to Regime 2 (two-phase CTMC) for that sample, exactly as Phases 4-5
     of calc_instability.py would, adapted to "flagged samples at one
     locus" instead of "flagged loci in one sample":
       - Global CTMC params (T1, T2, p_exp) come from --repeat-preset if
         given (recommended for known disease repeats, e.g. htt_cag), or
         are estimated from the flagged samples themselves (needs
         --regime2-min-samples of them; same two-stage estimator used
         genome-wide, regime2_ctmc_v2.estimate_global_params).
       - Each flagged sample is then fit independently with
         fit_locus_marginal() + Regime 2 GoF, using its own founder
         length and (optionally) its own donor age.
  5. Writes one row per sample (per haplotype) to a TSV. Samples that
     needed Regime 2 have both their R1 diagnostics and their R2 fit
     (r2_-prefixed columns) so you can see why R1 was rejected and what
     R2 found instead.

IMPORTANT — pooling vs. comparing:
  Each sample keeps its OWN founder length and (if given) its own donor
  age throughout. Reads are never pooled across samples; every fit is
  per-sample. Regime 2's global params (T1/T2/p_exp) are the one thing
  shared across flagged samples when estimated from data (this mirrors
  how the genome-wide pipeline shares them across flagged loci) — use
  --repeat-preset instead if you want zero cross-sample dependency.

Manifest format (TSV, no header needed):
    sample_id   reads_tsv            vcf            donor_age
    sample1     /path/sample1.tsv    /path/sample1.vcf   62
    sample2     /path/sample2.tsv.gz /path/sample2.vcf   -
    sample3     /path/sample3.tsv    -                   -

  vcf column: "-" or "" = no VCF, use modal length as founder.
  donor_age column: optional entirely (3-column manifests are fine).
    "-" or "" = no exact age for this sample -> use --tissue-type prior.

Usage
-----
  # Regime 1 only (as before)
  python analyze_locus_across_samples.py --manifest samples.tsv \
      --chrom chr4 --start 3074876 --end 3074933 --motif CAG \
      --stratum disease_cag --output out.tsv

  # + continue to Regime 2 for flagged samples, using a known preset
  python analyze_locus_across_samples.py --manifest samples.tsv \
      --chrom chr4 --start 3074876 --end 3074933 --motif CAG \
      --stratum disease_cag --run-regime2 --repeat-preset htt_cag \
      --tissue-type adult_blood --output out.tsv

  # + continue to Regime 2, estimating T2/p_exp from the flagged samples
  #   themselves (no preset available for this repeat)
  python analyze_locus_across_samples.py --manifest samples.tsv \
      --chrom chr1 --start 1000 --end 1030 --motif CAG \
      --run-regime2 --regime2-min-samples 3 --output out.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys

import numpy as np

from nbgeom_modelling import analyse_haplotype, Regime1GlobalParams
from calc_instability import load_genotypes_from_vcf, _compute_r2_gof
from regime2_ctmc import CTMCGlobalParams, PRESET_MAP
from regime2_ctmc_v2 import estimate_T1_from_founders, estimate_global_params
from regime2_marginal import TISSUE_AGE_PRIORS, fit_locus_marginal


# ---------------------------------------------------------------------------
# Locus-filtered TSV reader (streams the file, never loads it all into RAM)
# ---------------------------------------------------------------------------

def read_locus_haplotypes(reads_tsv: str, chrom: str, start: int, end: int,
                           motif: str) -> dict[int, list[float]]:
    """
    Stream a per-read TSV and pull out {haplotype: [length_bp, ...]} for
    exactly one locus, matching on (chrom, start, end, motif).

    Input TSV columns (tab-separated, matches calc_instability.py):
        chrom start end motif read_id haplotype length_bp allele avg_meth meth_bases
    """
    opener = gzip.open if reads_tsv.endswith(".gz") else open
    haps: dict[int, list[float]] = {}

    with opener(reads_tsv, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            (row_chrom, row_start, row_end, row_motif, _read_id,
             haplotype, length_bp) = parts[:7]

            if (row_chrom != chrom or int(row_start) != start
                    or int(row_end) != end): # or row_motif != motif):
                continue

            haps.setdefault(int(haplotype), []).append(float(length_bp))

    return haps


def founder_from_vcf_or_modal(genotypes: dict | None, locus_key: str,
                               hap_idx: int, lengths: list[float]) -> float:
    """Founder length in repeat units: VCF genotype if available, else modal."""
    if genotypes is not None:
        gt = genotypes.get(locus_key)
        try:
            return float(gt[hap_idx]) if hasattr(gt, "__len__") else float(gt)
        except (IndexError, TypeError, ValueError):
            pass
    arr = np.round(np.asarray(lengths)).astype(int)
    v, c = np.unique(arr, return_counts=True)
    return float(v[np.argmax(c)])


# ---------------------------------------------------------------------------
# Manifest loading
# ---------------------------------------------------------------------------

def load_manifest(path: str) -> list[dict]:
    rows = []
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for r in reader:
            if not r or r[0].startswith("#"):
                continue
            sample_id, reads_tsv = r[0], r[1]
            vcf = r[2] if len(r) > 2 and r[2] not in ("", "-") else None
            age = None
            if len(r) > 3 and r[3] not in ("", "-"):
                age = float(r[3])
            rows.append({"sample_id": sample_id, "reads_tsv": reads_tsv,
                         "vcf": vcf, "donor_age": age})
    return rows


# ---------------------------------------------------------------------------
# Regime 2 continuation for flagged samples
# ---------------------------------------------------------------------------

def _r2_min_L_max_L(lengths_int: np.ndarray, L0: int) -> tuple[int, int]:
    obs_max = int(lengths_int.max())
    obs_min = int(lengths_int.min())
    obs_span = max(obs_max - L0, L0 - obs_min, 1)
    min_L = min(max(obs_min - 10, 0), L0)          # must not exceed L0
    max_L = min(obs_max + max(30, obs_span // 2), min_L + 150, 1000)
    return min_L, max_L


def run_regime2_for_flagged(flagged_rows: list[dict], preset: str | None,
                             tissue_type: str, global_donor_age: float | None,
                             regime2_min_samples: int, r2_n_ppp: int,
                             r2_n_quad: int) -> None:
    """
    Mutates each row in flagged_rows in place, adding r2_-prefixed keys.

    flagged_rows entries must carry '_lengths' (np.ndarray) and
    '_founder' (float) from the R1 stage, plus optionally 'donor_age'.
    """
    if preset is not None:
        gp2 = CTMCGlobalParams.from_preset(preset)
        print(f"\nRegime 2: using preset '{preset}' "
              f"(T1={gp2.T1:.1f} T2={gp2.T2:.1f} p_exp={gp2.p_exp:.3f})",
              file=sys.stderr)
    else:
        if len(flagged_rows) < regime2_min_samples:
            print(f"\nRegime 2: only {len(flagged_rows)} flagged sample(s) "
                  f"(need >= {regime2_min_samples} to estimate T2/p_exp from "
                  f"data, or supply --repeat-preset). Skipping Regime 2.",
                  file=sys.stderr)
            return

        loci_data, founders = [], []
        for row in flagged_rows:
            age = row.get("donor_age") or global_donor_age or \
                  TISSUE_AGE_PRIORS[tissue_type]["mean"]
            loci_data.append({
                "lengths":        row["_lengths"],
                "founder_length": row["_founder"],
                "donor_age":      age,
            })
            founders.append(row["_founder"])

        T1 = estimate_T1_from_founders(founders, percentile=5.0)
        print(f"\nRegime 2: estimating global params from "
              f"{len(flagged_rows)} flagged sample(s), T1={T1:.1f} ...",
              file=sys.stderr)
        gp2 = estimate_global_params(loci_data, T1=T1, verbose=True)

    for row in flagged_rows:
        lengths = row["_lengths"]
        founder = row["_founder"]
        age = row.get("donor_age") or global_donor_age
        lengths_int = np.round(lengths).astype(int)
        L0 = int(round(founder))

        try:
            result = fit_locus_marginal(
                lengths, founder_length=founder, global_params=gp2,
                donor_age=age, tissue_type=tissue_type,
                n_quad=r2_n_quad, r1_n=30, r2_n=12,
                haplotype_label=row["haplotype"],
            )
        except Exception as e:
            print(f"  WARN Regime 2 fit failed for {row['sample_id']}: {e}",
                  file=sys.stderr)
            continue

        min_L, max_L = _r2_min_L_max_L(lengths_int, L0)
        age_mean = age if age is not None else TISSUE_AGE_PRIORS[tissue_type]["mean"]
        age_sd = 0.5 if age is not None else TISSUE_AGE_PRIORS[tissue_type]["sd"]

        gof2 = _compute_r2_gof(
            observed_lengths=lengths, founder_length=founder,
            r1_est=result.r1[0], r2_est=result.r2[0],
            T1=gp2.T1, T2=gp2.T2, p_exp=gp2.p_exp,
            min_L=min_L, max_L=max_L,
            age_mean=age_mean, age_sd=age_sd,
            n_quad=r2_n_quad, n_ppp=r2_n_ppp,
        )

        row["analysis_regime"] = "R2"
        row["r2_T1"], row["r2_T2"], row["r2_p_exp"] = gp2.T1, gp2.T2, gp2.p_exp
        row["r2_r1"], row["r2_r1_ci_lo"], row["r2_r1_ci_hi"] = result.r1
        row["r2_r2"], row["r2_r2_ci_lo"], row["r2_r2_ci_hi"] = result.r2
        row["r2_net_rate_phaseA"] = result.net_rate_phaseA[0]
        row["r2_net_rate_phaseB"] = result.net_rate_phaseB[0]
        row["r2_rate_acceleration"] = result.rate_acceleration[0]
        row["r2_instability_index"] = result.instability_index[0]
        row["r2_mean_delta"] = result.mean_delta[0]
        row["r2_phase_B_fraction"] = result.phase_B_fraction[0]
        row.update(gof2)

        print(f"  [{row['sample_id']}] R2 fit: r1={result.r1[0]:.4f} "
              f"r2={result.r2[0]:.4f} net_rate_phaseA={result.net_rate_phaseA[0]:.3f} "
              f"r2_fit_quality={gof2['r2_fit_quality']:.2f} "
              f"warning='{gof2['r2_model_warning']}'", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True,
                   help="TSV: sample_id  reads_tsv  vcf  [donor_age]")
    p.add_argument("--chrom",  required=True)
    p.add_argument("--start",  required=True, type=int)
    p.add_argument("--end",    required=True, type=int)
    p.add_argument("--motif",  required=True)
    p.add_argument("--stratum", default=None,
                   choices=["str_short", "str_medium", "str_long", "disease_cag",
                            "disease_cgg", "disease_gaa", "disease_ctg"],
                   help="Named Regime 1 global-params preset shared across all samples.")
    p.add_argument("--repeat-unit-length", type=int, default=None,
                   help="Use from_defaults(repeat_unit_length) instead of a stratum. "
                        "Defaults to len(motif).")
    p.add_argument("--min-reads", type=int, default=10,
                   help="Skip a sample/haplotype with fewer reads than this (default 10)")
    p.add_argument("--n-ppp", type=int, default=300,
                   help="Regime 1 PPP goodness-of-fit simulations (default 300; 0=skip)")
    p.add_argument("--output", required=True)

    r2 = p.add_argument_group("Regime 2 (for samples flagged gof_regime2_recommended)")
    r2.add_argument("--run-regime2", action="store_true",
                     help="Continue to the two-phase CTMC for any flagged sample.")
    r2.add_argument("--repeat-preset", default=None, choices=list(PRESET_MAP),
                     help="Known CTMC preset (e.g. htt_cag). Skips global-param "
                          "estimation entirely -- recommended when available.")
    r2.add_argument("--tissue-type", default="adult_blood",
                     choices=list(TISSUE_AGE_PRIORS),
                     help="Age prior for samples without an exact donor_age "
                          "(default: adult_blood)")
    r2.add_argument("--donor-age", type=float, default=None,
                     help="Exact donor age applied to all samples that don't "
                          "have their own age in the manifest.")
    r2.add_argument("--regime2-min-samples", type=int, default=3,
                     help="Min flagged samples needed to estimate T2/p_exp from "
                          "data when no --repeat-preset is given (default 3).")
    r2.add_argument("--r2-n-ppp", type=int, default=200,
                     help="Regime 2 PPP simulations (default 200; 0=skip)")
    r2.add_argument("--r2-n-quad", type=int, default=15,
                     help="Age quadrature points for Regime 2 (default 15)")
    args = p.parse_args()

    locus_key = f"{args.chrom}:{args.start}-{args.end}_{args.motif}"

    if args.stratum:
        gp = Regime1GlobalParams.from_stratum(args.stratum)
        print(f"Using shared Regime 1 global params: stratum='{args.stratum}'",
              file=sys.stderr)
    else:
        rul = args.repeat_unit_length or len(args.motif)
        gp = Regime1GlobalParams.from_defaults(repeat_unit_length=rul)
        print(f"Using shared Regime 1 global params: "
              f"from_defaults(repeat_unit_length={rul})", file=sys.stderr)

    manifest = load_manifest(args.manifest)
    print(f"{len(manifest)} samples in manifest", file=sys.stderr)

    rows_out = []
    for entry in manifest:
        sid, reads_tsv, vcf = entry["sample_id"], entry["reads_tsv"], entry["vcf"]
        print(f"\n[{sid}] reading {reads_tsv} ...", file=sys.stderr)

        haps = read_locus_haplotypes(reads_tsv, args.chrom, args.start,
                                      args.end, args.motif)
        print('\n\n\n',haps, file=sys.stderr)
        if not haps:
            print(f"  no reads found at locus in {sid} — skipping", file=sys.stderr)
            continue

        genotypes = None
        if vcf:
            genotypes = load_genotypes_from_vcf(vcf)

        for hap_idx, (hap, lengths) in enumerate(haps.items()):
            if len(lengths) < args.min_reads:
                print(f"  hap {hap}: only {len(lengths)} reads (<{args.min_reads}) — skipping",
                      file=sys.stderr)
                continue

            founder = founder_from_vcf_or_modal(genotypes, locus_key, hap, lengths)
            lengths_arr = np.asarray(lengths, dtype=float)

            result = analyse_haplotype(
                lengths=lengths_arr, founder_length=founder, global_params=gp,
                haplotype_label=f"{sid}_hap{hap}", n_ppp=args.n_ppp,
            )
            result["sample_id"] = sid
            result["analysis_regime"] = "R1"
            result["_lengths"] = lengths_arr        # kept for possible Regime 2 use
            result["_founder"] = founder
            result["donor_age"] = entry.get("donor_age")
            rows_out.append(result)
            print(f"  hap {hap}: n_reads={result['n_reads']} founder={founder:.1f} "
                  f"instability_index={result['instability_index']:.3f} "
                  f"net_bias={result['net_bias']:.3f} "
                  f"fit_quality={result['fit_quality']:.2f} "
                  f"gof_regime2_recommended={result['gof_regime2_recommended']}",
                  file=sys.stderr)

    if not rows_out:
        print("\nNo samples produced results — check locus coordinates and --min-reads.",
              file=sys.stderr)
        sys.exit(1)

    # ── Regime 2 continuation ────────────────────────────────────────────
    if args.run_regime2:
        flagged_rows = [r for r in rows_out if r.get("gof_regime2_recommended")]
        if not flagged_rows:
            print("\nRegime 2: no samples flagged gof_regime2_recommended — "
                  "Regime 1 fit was adequate for all samples.", file=sys.stderr)
        else:
            print(f"\n{len(flagged_rows)}/{len(rows_out)} sample-haplotypes "
                  f"flagged for Regime 2: "
                  f"{[r['sample_id'] for r in flagged_rows]}", file=sys.stderr)
            run_regime2_for_flagged(
                flagged_rows, args.repeat_preset, args.tissue_type,
                args.donor_age, args.regime2_min_samples,
                args.r2_n_ppp, args.r2_n_quad,
            )
    elif any(r.get("gof_regime2_recommended") for r in rows_out):
        n_flagged = sum(1 for r in rows_out if r.get("gof_regime2_recommended"))
        print(f"\nNote: {n_flagged} sample(s) flagged gof_regime2_recommended "
              f"(Regime 1 fit is poor for them). Re-run with --run-regime2 "
              f"(and --repeat-preset if this is a known repeat) to fit those "
              f"samples with the CTMC model.", file=sys.stderr)

    # Sort by instability_index, most unstable first (Regime 1 scale;
    # r2_instability_index is on the same Var(Delta) scale for R2 rows)
    rows_out.sort(key=lambda r: r["instability_index"], reverse=True)

    priority_cols = [
        "sample_id", "haplotype", "analysis_regime", "n_reads", "founder_length",
        "instability_index", "instability_index_ci_lo", "instability_index_ci_hi",
        "net_bias", "net_bias_ci_lo", "net_bias_ci_hi",
        "p_expansion", "winsorized_mean_abs_delta",
        "fit_quality", "gof_regime2_recommended",
        "r_e", "r_c", "q_e", "q_c", "p_plus", "p_minus",
        "r2_T1", "r2_T2", "r2_p_exp",
        "r2_r1", "r2_r1_ci_lo", "r2_r1_ci_hi",
        "r2_r2", "r2_r2_ci_lo", "r2_r2_ci_hi",
        "r2_net_rate_phaseA", "r2_net_rate_phaseB", "r2_rate_acceleration",
        "r2_instability_index", "r2_mean_delta", "r2_phase_B_fraction",
        "r2_ks_stat", "r2_ppp_var", "r2_tail_coverage", "r2_mean_bias",
        "r2_fit_quality", "r2_model_warning",
    ]
    all_keys = set()
    for r in rows_out:
        all_keys.update(r.keys())
    all_keys.discard("_lengths"); all_keys.discard("_founder")
    ordered_cols = priority_cols + [c for c in sorted(all_keys) if c not in priority_cols]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ordered_cols, delimiter="\t",
                                 extrasaction="ignore")
        writer.writeheader()
        for row in rows_out:
            writer.writerow(row)

    print(f"\n{len(rows_out)} sample-haplotype rows -> {args.output}", file=sys.stderr)
    print("\nSummary (sorted by instability_index, descending):", file=sys.stderr)
    for r in rows_out:
        regime = r.get("analysis_regime", "R1")
        flag = "  <- gof_regime2_recommended" if r.get("gof_regime2_recommended") else ""
        extra = ""
        if regime == "R2":
            extra = (f"  [R2: net_rate_phaseA={r['r2_net_rate_phaseA']:.3f} "
                      f"r2_fit_quality={r['r2_fit_quality']:.2f}]")
        print(f"  {r['sample_id']:<15} [{regime}] "
              f"instability_index={r['instability_index']:8.3f}  "
              f"net_bias={r['net_bias']:7.3f}  "
              f"winsor_mean|Δ|={r['winsorized_mean_abs_delta']:6.2f}  "
              f"fit_quality={r['fit_quality']:.2f}{flag}{extra}", file=sys.stderr)


if __name__ == "__main__":
    main()