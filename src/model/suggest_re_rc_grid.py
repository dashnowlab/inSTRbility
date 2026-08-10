#!/usr/bin/env python3
"""
suggest_re_rc_grid.py — Data-driven r_e / r_c grid ceilings
================================================================

DEFAULT_RE_GRID/DEFAULT_RC_GRID's ceilings (8.0, 6.0) are empirical
defaults tuned for typical STR behaviour, not derived constants. This
tool gives you a fast, closed-form way to check whether they're
appropriate for YOUR data, using the exact relationship:

    Var(E) = r_e * K(p_plus, q_e)      where K is a known constant
                                        once p_plus/q_e are fixed (Pass 1)

This means a rough r_e (and r_c) can be estimated per locus directly from
observed variance -- no grid search needed -- just to see what range your
data actually spans, before committing to an expensive fit.

Two things this gives you
--------------------------
1. suggest_grid_ceiling(): a data-driven ceiling recommendation, from the
   distribution of these closed-form estimates across your loci.
2. check_grid_saturation(): given an ALREADY-FITTED results TSV, flags any
   locus whose fitted r_e/r_c landed suspiciously close to the current
   grid's edge -- the silent-saturation failure mode discussed earlier
   (fit_locus_r1 has no built-in boundary check despite the docstring
   implying one exists).

Usage
-----
    from suggest_re_rc_grid import suggest_grid_ceiling, check_grid_saturation

    # Before fitting: how wide should my grid be?
    rec = suggest_grid_ceiling(loci_data, global_params)
    print(rec)

    # After fitting: did any locus saturate the CURRENT grid?
    flagged = check_grid_saturation("results.tsv", re_grid=..., rc_grid=...)
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from src.model.utils import get_stratum

def _K(p: float, q: float) -> float:
    """Closed form: Var(sum of geometric steps over NB(r,q) events) = r * K(p,q)."""
    return (1 - q) / (p**2 * q**2) * ((1 - p) * q + 1)


def plot_histogram(data, name='default.png'):
    fig, ax = plt.subplots(figsize=(6, 4))
    data = np.asarray(data, dtype=float)
    
    # Histogram of observed data as a probability mass at each integer Delta
    bin_edges = np.arange(data.min() - 1, data.max() + 1, 5)
    ax.hist(data, bins=bin_edges, color="#a9d18e", edgecolor="#68944a", label=f"Observed reads (n={len(data)})")
    # ax.plot(delta_grid, pmf, color="#0b8789", lw=1.8, linestyle="-", label="Fitted NB-Geometric PMF")
    ax.set_xlabel(r"Re values")
    ax.set_ylabel("Frequency")
    ax.legend(fontsize=9)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(0,200)
    # ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'../test_data/{name}', dpi=300)


def suggest_grid_ceiling(
    loci_data: list[dict],
    global_params,                  # Regime1GlobalParams
    percentile: float = 99.0,
    margin: float = 1.5,
    min_reads_per_direction: int = 6,
    winsorize_at: float | None = 100.0,
) -> dict:
    """
    Estimate an appropriate r_e/r_c grid ceiling from your own data, using
    the closed-form Var(E) = r_e * K(p_plus, q_e) relationship (no fitting
    required).
 
    Parameters
    ----------
    loci_data : list of dicts, each with:
        'lengths'         : array-like -- per-read allele lengths
        'founder_length'  : float
    global_params : Regime1GlobalParams
        Supplies p_plus(L0)/p_minus(L0) and q_e/q_c.
    percentile : float
        Which percentile of the per-locus r estimates to size the grid
        around (default 99th -- covers all but the most extreme outlier
        loci; raise to 99.9 if you have thousands of loci and want to
        cover essentially everything).
    margin : float
        Safety multiplier applied on top of that percentile (default 1.5x)
        so the ceiling isn't sitting exactly at the edge of your observed
        range either.
    min_reads_per_direction : int
        Minimum reads in a direction to trust that locus's variance
        estimate (default 8).
    winsorize_at : float or None
        Cap |Delta| at this value (in repeat units, same philosophy as
        winsorized_mean_abs_delta) BEFORE computing per-locus variance.
        This is a method-of-moments estimator with no robustness of its
        own -- a handful of outlier reads (artifacts, or genuine but
        extreme somatic events) can inflate a locus's raw variance, and
        therefore its r-hat, by orders of magnitude with no biological
        meaning at the Regime-1 level. Set to None to disable (matches
        the original, non-robust behaviour).
 
    Returns
    -------
    dict with suggested ceilings, percentile diagnostics (including a
    max/p99 "tail concentration" ratio -- a large ratio means a handful
    of loci are driving the estimate, not a broad population shift), and
    n loci used.
    """

    r_e_estimates, r_c_estimates = {}, {}

    L0 = None

    for locus in loci_data:
        ml = len(locus['motif'])
        stratum = get_stratum(ml)
        gp = global_params[stratum]
        stratum = get_stratum(ml)
        if stratum not in r_e_estimates: r_e_estimates[stratum] = []
        if stratum not in r_c_estimates: r_c_estimates[stratum] = []
        for h,hap in enumerate(sorted(list(locus["haplotypes"].keys()))):
            lengths = np.asarray(locus["haplotypes"][hap], dtype=float)
            if locus["genotypes"] is None or locus["genotypes"][h] is None:
                values, counts = np.unique(lengths, return_counts=True)
                L0 = values[np.argmax(counts)]
            else: L0 = locus["genotypes"][h]
            deltas = np.round(lengths).astype(int) - int(round(L0))
 
            pos = deltas[deltas > 0]
            neg = -deltas[deltas < 0]
            if winsorize_at is not None:
                pos = np.minimum(pos, winsorize_at)
                neg = np.minimum(neg, winsorize_at)
    
            p_plus  = gp.p_plus(L0)
            p_minus = gp.p_minus(L0)
            qe, qc  = gp.q_e, gp.q_c
    
            if len(pos) >= min_reads_per_direction:
                var_pos = float(pos.var())
                r_e_hat = var_pos / _K(p_plus, qe)
                if r_e_hat > 100:
                    f = round(locus["genotypes"][h], 3)
                    l = [round(x, 3) for x in locus["haplotypes"][hap] if round(x, 3) != f]
                if np.isfinite(r_e_hat) and r_e_hat > 0:
                    r_e_estimates[stratum].append(r_e_hat)
    
            if len(neg) >= min_reads_per_direction:
                var_neg = float(neg.var())
                r_c_hat = var_neg / _K(p_minus, qc)
                if np.isfinite(r_c_hat) and r_c_hat > 0:
                    r_c_estimates[stratum].append(r_c_hat)


        def _summarize(estimates, current_default, label):
            if len(estimates) < 10:
                return {
                    "label": label, "n_loci": len(estimates),
                    "suggested_ceiling": None,
                    "note": f"Only {len(estimates)} loci with enough reads -- "
                            f"too few to recommend a ceiling confidently.",
                }
            arr = np.array(estimates)
            pctile_val = float(np.percentile(arr, percentile))
            suggested = pctile_val * margin
            current_pctile_rank = float((arr < current_default).mean() * 100)
            p50 = float(np.percentile(arr, 50))
            p90 = float(np.percentile(arr, 90))
            max_val = float(arr.max())
            tail_ratio = max_val / pctile_val if pctile_val > 0 else float("inf")
            n_above_10x_p99 = int((arr > 10 * pctile_val).sum())
            interpretation = (
                "OUTLIER-DRIVEN: a small number of loci sit far beyond the rest -- "
                "investigate/exclude them (check gof_regime2_recommended) rather "
                "than sizing the grid to cover them."
                if tail_ratio > 20 else
                "Broadly elevated -- the shift looks like it applies across many "
                "loci, not just a few extreme ones; likely a genuine grid-sizing issue."
            )
            return {
                "label": label,
                "n_loci": len(estimates),
                "p50": p50, "p90": p90,
                f"p{percentile:g}_estimate": pctile_val,
                "max_estimate_in_data": max_val,
                "tail_concentration_ratio (max/p{:g})".format(percentile): tail_ratio,
                "n_loci_over_10x_the_percentile": n_above_10x_p99,
                "suggested_ceiling": suggested,
                "current_default": current_default,
                "current_default_covers_pctile": current_pctile_rank,
                "interpretation": interpretation,
            }

    r_e_summary, r_c_summary = {}, {}
    strata = set(list(r_e_estimates.keys()) + list(r_c_estimates.keys()))
    for stratum in strata:
        if stratum in r_e_estimates: r_e_summary[stratum] = _summarize(r_e_estimates[stratum], 8.0, "r_e")
        if stratum in r_c_estimates: r_c_summary[stratum] = _summarize(r_c_estimates[stratum], 6.0, "r_c")

    return {"r_e": r_e_summary, "r_c": r_c_summary}


def check_grid_saturation(
    results_tsv: str,
    re_grid: np.ndarray,
    rc_grid: np.ndarray,
    edge_tolerance: float = 0.05,
) -> list[dict]:
    """
    Flag loci whose FITTED r_e or r_c landed within `edge_tolerance` (as a
    fraction of the grid's range, in log space) of the grid's upper edge --
    the silent-saturation failure mode: the posterior mean got pulled
    toward the ceiling because the true value is unconstrained above it,
    and nothing in fit_locus_r1 warns you this happened.

    Parameters
    ----------
    results_tsv : str
        Path to a results TSV (must have 'haplotype_label', 'r_e', 'r_c').
    re_grid, rc_grid : np.ndarray
        The EXACT grids used when the fit was run (e.g. DEFAULT_RE_GRID,
        or whatever you passed to fit_locus_r1).
    edge_tolerance : float
        Fraction of the log-range counted as "near the edge" (default 0.05
        = top 5% of the grid's log-span).

    Returns
    -------
    List of dicts for flagged loci: {haplotype_label, r_e, r_c, flagged_on}
    """
    import csv

    log_re_max, log_re_min = np.log(re_grid.max()), np.log(re_grid.min())
    log_rc_max, log_rc_min = np.log(rc_grid.max()), np.log(rc_grid.min())
    re_edge = re_grid.max() / np.exp(edge_tolerance * (log_re_max - log_re_min))
    rc_edge = rc_grid.max() / np.exp(edge_tolerance * (log_rc_max - log_rc_min))

    flagged = []
    with open(results_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                r_e, r_c = float(row["r_e"]), float(row["r_c"])
            except (KeyError, ValueError):
                continue
            hit_e = r_e >= re_edge
            hit_c = r_c >= rc_edge
            if hit_e or hit_c:
                flagged.append({
                    "haplotype_label": row.get("haplotype_label", "?"),
                    "r_e": r_e, "r_c": r_c,
                    "flagged_on": ("r_e" if hit_e else "") + ("+r_c" if hit_c else ""),
                })
    return flagged


if __name__ == "__main__":
    # Self-test with synthetic data spanning a range of true r values,
    # including some deliberately ABOVE the default ceiling of 8.0
    rng = np.random.default_rng(3)

    class GP:
        q_e, q_c = 0.5, 0.5
        def p_plus(self, L0):  return 0.4
        def p_minus(self, L0): return 0.5

    gp = GP()
    true_rs = np.concatenate([
        rng.uniform(0.05, 2.0, 40),
        rng.uniform(2.0, 15.0, 15),   # some genuinely need r_e > 8
    ])
    loci_data = []
    for i, r_true in enumerate(true_rs):
        n = rng.integers(30, 80)
        Ne = rng.negative_binomial(r_true, 0.5, size=n)
        E = np.array([rng.geometric(0.4, size=max(ne, 0)).sum() if ne > 0 else 0 for ne in Ne])
        Nc = rng.negative_binomial(0.8, 0.5, size=n)
        C = np.array([rng.geometric(0.5, size=max(nc, 0)).sum() if nc > 0 else 0 for nc in Nc])
        lengths = 44 + E - C
        loci_data.append({"lengths": lengths, "founder_length": 44.0})

    rec = suggest_grid_ceiling(loci_data, gp, percentile=95, margin=1.5)
    print("r_e:", rec["r_e"])
    print("r_c:", rec["r_c"])
