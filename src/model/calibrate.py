import numpy as np
import csv
import sys
import time

from src.model.utils import get_stratum, stratum_labels

"""
calibrate_qe_qc_split.py — Independent q_e / q_c calibration by direction
============================================================================

Regime1GlobalParams.calibrate_from_results() computes ONE dispersion ratio
from the combined Delta = E - C and applies it to both q_e and q_c. This
module splits that into two independent calibrations -- one from the
expansion-only (positive Delta) reads, one from the contraction-only
(negative Delta) reads -- using the exact same back-calculation formula,
just applied twice instead of once.

Why this needs the raw per-locus deltas, not just the results TSV
-------------------------------------------------------------------
calibrate_from_results() only needs `gof_dispersion_ratio`, a single
per-locus scalar already computed from the COMBINED Delta distribution.
There is no column in the output TSV that separates "how much of that
mismatch came from the expansion side" from "how much came from the
contraction side" -- that information was never retained past the
combined variance. So this function needs a second input: the raw
per-locus reads (the same `loci_data` structure fed to
estimate_regime1_global), so it can compute observed variance separately
for the positive and negative subsets of each locus's Delta values.

The fitted r_e, r_c, p_plus, p_minus (used to compute each direction's
FITTED variance for comparison) still come from the results TSV, since
those are the actual per-locus point estimates from Pass 2 -- no need to
refit anything.

Usage
-----
    from calibrate_qe_qc_split import calibrate_from_results_split

    gp_split = calibrate_from_results_split(
        results_tsv="initial_results.tsv",
        loci_data=loci_data,          # same structure as estimate_regime1_global
        base_params=gp_initial,
    )
    print(gp_split.q_e, gp_split.q_c)   # now genuinely different values
"""

def _read_results_tsv(results_tsv: str) -> dict:
    """Same pandas-or-csv-fallback pattern as calibrate_from_results, for consistency."""
    try:
        import pandas as pd
        df = pd.read_csv(results_tsv, sep="\t")
        return {col: df[col].values for col in df.columns}
    except ImportError:
        import csv
        rows = []
        with open(results_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                rows.append(row)
        if not rows:
            return {}
        out = {k: [] for k in rows[0]}
        for row in rows:
            for k, v in row.items():
                try:
                    out[k].append(float(v))
                except (ValueError, TypeError):
                    out[k].append(v)
        return {k: np.array(v) for k, v in out.items()}


def recalibrate_distinct(
    results_tsv: str,
    loci_data: list[dict],
    base_params,                       # Regime1GlobalParams
    min_reads: int = 20,
    min_reads_per_direction: int = 8,
    target_percentile: float = 50,
    verbose: bool = True,
):
    """
    Calibrate q_e and q_c independently, using each locus's expansion-only
    and contraction-only reads separately.

    Parameters
    ----------
    results_tsv : str
        Path to the TSV from write_tsv()/analyse_genome_wide(). Must contain
        'haplotype_label', 'n_reads', 'r_e', 'r_c', 'p_plus', 'p_minus'.
    loci_data : list of dicts, each with:
        'haplotype_label' : str   -- must match a row in results_tsv
        'deltas'          : array -- raw per-read Delta_i = observed - founder
        (same structure as estimate_regime1_global's input, plus the label)
    base_params : Regime1GlobalParams
        Starting parameters. p_plus/p_minus/r are untouched; only q_e, q_c
        are adjusted, independently.
    min_reads : int
        Minimum total reads at a locus to consider it at all (default 20).
    min_reads_per_direction : int
        Minimum reads on ONE side (expansion or contraction) to include
        that side's variance in that direction's calibration (default 8).
        This is deliberately independent per direction -- a locus can
        contribute to the q_e calibration without qualifying for q_c, and
        vice versa, since one direction is very commonly much rarer than
        the other.
    target_percentile : float
        Percentile of the ratio distribution to use (default 50 = median).

    Returns
    -------
    A copy of base_params with q_e and q_c independently recalibrated.
    """
    results = _read_results_tsv(results_tsv)
    required = ["haplotype", "n_reads", "r_e", "r_c", "p_plus", "p_minus"]
    missing = [c for c in required if c not in results]
    if missing:
        raise ValueError(f"results_tsv is missing required columns: {missing}")

    lookup = {}
    for i, label in enumerate(results["haplotype"]):
        lookup[label] = {
            "n_reads": float(results["n_reads"][i]),
            "r_e": float(results["r_e"][i]),
            "r_c": float(results["r_c"][i]),
            "p_plus": float(results["p_plus"][i]),
            "p_minus": float(results["p_minus"][i]),
        }

    qe0, qc0 = base_params.q_e, base_params.q_c
    ratios_e, ratios_c = [], []
    n_loci_e = n_loci_c = 0

    for locus in loci_data:
        label = locus.get("haplotype")
        row = lookup.get(label)
        if row is None or row["n_reads"] < min_reads:
            continue

        deltas = np.asarray(locus["deltas"], dtype=float)
        pos = deltas[deltas > 0]
        neg = -deltas[deltas < 0]

        r_e, p_plus = row["r_e"], row["p_plus"]
        r_c, p_minus = row["r_c"], row["p_minus"]

        # Expansion side
        if len(pos) >= min_reads_per_direction:
            obs_var_e = float(pos.var())
            mu_Ne  = r_e * (1 - qe0) / qe0
            var_Ne = r_e * (1 - qe0) / qe0**2
            fitted_var_e = mu_Ne * (1 - p_plus) / p_plus**2 + var_Ne / p_plus**2
            if fitted_var_e > 0 and np.isfinite(obs_var_e):
                ratios_e.append(obs_var_e / fitted_var_e)
                n_loci_e += 1

        # Contraction side
        if len(neg) >= min_reads_per_direction:
            obs_var_c = float(neg.var())
            mu_Nc  = r_c * (1 - qc0) / qc0
            var_Nc = r_c * (1 - qc0) / qc0**2
            fitted_var_c = mu_Nc * (1 - p_minus) / p_minus**2 + var_Nc / p_minus**2
            if fitted_var_c > 0 and np.isfinite(obs_var_c):
                ratios_c.append(obs_var_c / fitted_var_c)
                n_loci_c += 1

    if n_loci_e < 10:
        # raise UserWarning(
        #     f"Only {n_loci_e} loci have >= {min_reads_per_direction} expansion "
        #     f"reads. Need at least 10 to calibrate q_e -- lower "
        #     f"min_reads_per_direction, or fall back to the combined "
        #     f"calibrate_from_results() for this direction."
        # )
        return False
    
    if n_loci_c < 10:
        # raise UserWarning(
        #     f"Only {n_loci_c} loci have >= {min_reads_per_direction} contraction "
        #     f"reads. Need at least 10 to calibrate q_c -- lower "
        #     f"min_reads_per_direction, or fall back to the combined "
        #     f"calibrate_from_results() for this direction."
        # )
        return False

    def _back_calc(current_q, ratios):
        target_ratio = float(np.percentile(ratios, target_percentile))
        current_rate = (1.0 - current_q) / current_q
        new_rate = current_rate * target_ratio
        new_q = float(np.clip(1.0 / (1.0 + new_rate), 0.05, 0.99))
        return new_q, target_ratio

    new_qe, ratio_e_med = _back_calc(qe0, ratios_e)
    new_qc, ratio_c_med = _back_calc(qc0, ratios_c)

    if verbose:
        print(f"Split calibration:")
        print(f"  Expansion:   {n_loci_e} loci, median ratio={ratio_e_med:.4f}, "
              f"q_e: {qe0:.3f} -> {new_qe:.3f}")
        print(f"  Contraction: {n_loci_c} loci, median ratio={ratio_c_med:.4f}, "
              f"q_c: {qc0:.3f} -> {new_qc:.3f}")
        if abs(new_qe - new_qc) < 0.02:
            print("  Note: q_e and q_c converged to nearly the same value -- "
                  "the symmetric assumption may have been fine for this data.")

    new_params = copy.deepcopy(base_params)
    new_params.q_e = new_qe
    new_params.q_c = new_qc
    return new_params


def recalibrate(
    gp_map:    dict,
    r1_tsv:    str,
    min_reads: int = 20,
    strata:    list | None = None
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
                if strata is not None and s not in strata:
                    continue
                if s in stratum_disp:
                    stratum_disp[s].append(dr)
    except FileNotFoundError:
        print(f"  Calibration: {r1_tsv} not found, skipping.", file=sys.stderr)
        return gp_map


    new_gp_map = {}
    import copy

    for stratum, gp in gp_map.items():
        if strata is not None and stratum not in strata:
            new_gp_map[stratum] = gp
            continue
        disp_vals = stratum_disp.get(stratum, [])
        if len(disp_vals) < 5:
            print(f"  Stratum {stratum_labels.get(stratum, f'{stratum}bp')}: too few loci ({len(disp_vals)}) "
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
        print(f"  Stratum {stratum_labels.get(stratum, f'{stratum}bp')} (n={len(disp_vals)}): "
              f"median disp_ratio={median_dr:.3f} — model {direction} variance "
              f"→ q: {gp.q_e:.3f} → {new_q:.3f}", file=sys.stderr)

        new_gp = copy.deepcopy(gp)
        new_gp.q_e = new_q
        new_gp.q_c = new_q
        new_gp_map[stratum] = new_gp

    return new_gp_map
