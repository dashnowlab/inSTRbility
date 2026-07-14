#!/usr/bin/env python3
"""
plot_regime2_fit.py — Overlay the fitted Regime 2 (two-phase CTMC) distribution
================================================================================

Same idea as plot_regime1_fit.py, but for the CTMC model. Unlike Regime 1's
fitted PMF (buried inline inside _compute_gof), Regime 2 already exposes
reusable PMF functions:

    regime2_ctmc._ctmc_pmf(L0, age, r1, r2, T1, T2, p_exp, min_L, max_L)
        -- PMF at a single, exactly-known age (matrix exponential of Q).

    regime2_marginal._marginalised_pmf(L0, r1, r2, T1, T2, p_exp, min_L,
                                        max_L, age_mean, age_sd, n_quad)
        -- PMF averaged over a Gamma(age_mean, age_sd) age prior, via
           quadrature. This is what fit_locus_marginal() actually fits
           against, so it's the correct curve to overlay when donor age
           is uncertain (the normal case for tissue-type priors).

This module wraps whichever is appropriate and adds T1/T2 threshold
markers, since those are the mechanistically interesting reference points
in Regime 2 (Regime 1 has no analogous length-dependent threshold).

Usage
-----
    from plot_regime2_fit import plot_regime2_fit
    from regime2_ctmc import CTMCGlobalParams

    gp2 = CTMCGlobalParams.from_preset("htt_cag")

    # Let it fit for you (age uncertain -> tissue-type prior):
    plot_regime2_fit(lengths, founder_length, global_params=gp2,
                      tissue_type="adult_blood", out_path="locus_r2_fit.png")

    # With an exactly-known donor age:
    plot_regime2_fit(lengths, founder_length, global_params=gp2,
                      donor_age=63, out_path="locus_r2_fit.png")

    # Reuse a result you already have (CTMCResult, or a dict with
    # r2_-prefixed keys like analyze_locus_across_samples.py produces):
    plot_regime2_fit(lengths, founder_length, global_params=gp2,
                      result=my_result, out_path="locus_r2_fit.png")
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from regime2_ctmc import _ctmc_pmf
from regime2_marginal import _marginalised_pmf, TISSUE_AGE_PRIORS, _age_prior_params


# ---------------------------------------------------------------------------
# State-space bounds -- must match fit_locus_marginal's own adaptive logic,
# or the PMF you plot won't be over the same grid the fit actually used.
# ---------------------------------------------------------------------------

def _state_space(lengths_int: np.ndarray, L0: int, max_n_states: int = 80) -> tuple[int, int]:
    obs_max  = int(lengths_int.max())
    obs_min  = int(lengths_int.min())
    obs_span = max(obs_max - L0, L0 - obs_min, 1)
    min_L    = min(max(obs_min - 10, 0), L0)
    max_L    = min(min_L + max_n_states - 1,
                   obs_max + max(30, obs_span // 2),
                   1000)
    return min_L, max_L


def _get(result, *names):
    """
    Pull a value off `result` trying each name in turn, working whether
    result is a CTMCResult dataclass (attrs, possibly (mean, lo, hi) tuples),
    or a flat dict (possibly with r2_-prefixed keys, as produced by
    analyze_locus_across_samples.py).
    """
    for name in names:
        val = None
        if hasattr(result, name):
            val = getattr(result, name)
        elif isinstance(result, dict) and name in result:
            val = result[name]
        if val is not None:
            return val[0] if isinstance(val, (tuple, list, np.ndarray)) else val
    raise KeyError(f"None of {names} found on result ({type(result)})")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_regime2_fit(
    lengths, founder_length: float,
    global_params=None,
    result=None,
    donor_age: Optional[float] = None,
    tissue_type: str = "adult_blood",
    n_quad: int = 15,
    max_n_states: int = 80,
    out_path: Optional[str] = None,
    title: Optional[str] = None,
    show_cdf: bool = False,
    mark_thresholds: bool = True,
):
    """
    Overlay the fitted two-phase CTMC distribution on a histogram of
    observed Delta values.

    Parameters
    ----------
    lengths : array-like
        Per-read allele lengths (repeat units).
    founder_length : float
        Germline / founder length L0.
    global_params : CTMCGlobalParams
        Supplies T1, T2, p_exp. Required (either to fit fresh, or -- even
        when `result` is given -- because T1/T2/p_exp may not be attached
        to `result` depending on how it was built).
    result : CTMCResult or dict, optional
        A previous fit (from fit_locus_ctmc / fit_locus_marginal, or a row
        from analyze_locus_across_samples.py's output with r2_ columns).
        If not given, fits fresh via fit_locus_marginal().
    donor_age : float, optional
        Exact donor age. If omitted, uses the `tissue_type` prior (the
        same behaviour as fit_locus_marginal).
    tissue_type : str
        Age prior when donor_age is not given.
    n_quad : int
        Age-quadrature points for the marginalised PMF (default 15;
        matches the resolution used elsewhere, not fit_locus_marginal's
        faster default of 3 -- for a plot you want the smooth version).
    max_n_states : int
        Must match whatever was used at fit time (default 80, matching
        fit_locus_marginal's default) so the PMF grid lines up with the
        state space the fit actually used.
    out_path, title, show_cdf : see plot_regime1_fit.py.
    mark_thresholds : bool
        If True, draws vertical dashed lines at T1-L0 and T2-L0 (the
        length-dependent phase thresholds, in Delta space).

    Returns
    -------
    (fig, ax) or (fig, (ax1, ax2)) if show_cdf=True
    """
    if global_params is None:
        raise ValueError("global_params is always required (supplies T1/T2/p_exp).")

    lengths = np.asarray(lengths, dtype=float)
    lengths_int = np.round(lengths).astype(int)
    L0 = int(round(founder_length))
    deltas = lengths_int - L0

    if result is None:
        from regime2_marginal import fit_locus_marginal
        result = fit_locus_marginal(
            lengths, founder_length=founder_length, global_params=global_params,
            donor_age=donor_age, tissue_type=tissue_type, n_quad=3,
        )

    T1    = _get(result, "T1", "r2_T1") if hasattr(result, "T1") or (
             isinstance(result, dict) and ("T1" in result or "r2_T1" in result)
             ) else global_params.T1
    T2    = _get(result, "T2", "r2_T2") if hasattr(result, "T2") or (
             isinstance(result, dict) and ("T2" in result or "r2_T2" in result)
             ) else global_params.T2
    p_exp = _get(result, "p_exp", "r2_p_exp") if hasattr(result, "p_exp") or (
             isinstance(result, dict) and ("p_exp" in result or "r2_p_exp" in result)
             ) else global_params.p_exp
    r1_est = _get(result, "r1", "r2_r1")
    r2_est = _get(result, "r2", "r2_r2")

    min_L, max_L = _state_space(lengths_int, L0, max_n_states)

    if donor_age is not None:
        age_mean, age_sd = float(donor_age), 0.5
    else:
        age_mean, age_sd = _age_prior_params(tissue_type)

    pmf = _marginalised_pmf(L0, r1_est, r2_est, T1, T2, p_exp,
                             min_L, max_L, age_mean, age_sd, n_quad=n_quad)
    lengths_grid = np.arange(min_L, max_L + 1)
    delta_grid = lengths_grid - L0

    fit_q = None
    for name in ("r2_fit_quality", "gof_fit_quality", "fit_quality"):
        if isinstance(result, dict) and name in result:
            fit_q = result[name]
            break
    warning = None
    for name in ("r2_model_warning", "gof_model_warning", "model_warning"):
        if isinstance(result, dict) and name in result:
            warning = result[name]
            break

    subtitle_bits = [f"n={len(lengths)}", f"L0={L0}",
                      f"r1={r1_est:.4f}", f"r2={r2_est:.4f}"]
    if fit_q is not None:
        subtitle_bits.append(f"fit_quality={fit_q:.2f}")
    if warning:
        subtitle_bits.append(f"[{warning}]")
    default_title = "Regime 2 (two-phase CTMC) fit  (" + "  ".join(subtitle_bits) + ")"

    if show_cdf:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    else:
        fig, ax1 = plt.subplots(figsize=(7, 4.2))

    bin_edges = np.arange(delta_grid.min() - 0.5, delta_grid.max() + 1.5, 1)
    ax1.hist(deltas, bins=bin_edges, density=True, alpha=0.5,
              color="steelblue", label=f"Observed reads (n={len(deltas)})")
    ax1.plot(delta_grid, pmf, color="darkorange", lw=1.8,
              label="Fitted CTMC PMF (age-marginalised)")
    if mark_thresholds:
        if min_L <= T1 <= max_L:
            ax1.axvline(T1 - L0, color="gray", ls="--", lw=1, alpha=0.7, label=f"T1={T1:.0f}")
        if min_L <= T2 <= max_L:
            ax1.axvline(T2 - L0, color="black", ls="--", lw=1, alpha=0.7, label=f"T2={T2:.0f}")
    ax1.set_xlabel(r"$\Delta$ (repeat units from founder)")
    ax1.set_ylabel("Probability")
    ax1.legend(fontsize=8)
    ax1.set_title(title or default_title, fontsize=9.5)

    if show_cdf:
        obs_sorted = np.sort(deltas)
        emp_cdf = np.arange(1, len(obs_sorted) + 1) / len(obs_sorted)
        fitted_cdf = np.cumsum(pmf)
        ax2.step(obs_sorted, emp_cdf, where="post", color="steelblue",
                  label="Empirical CDF")
        ax2.plot(delta_grid, fitted_cdf, color="darkorange", lw=1.8,
                  label="Fitted CDF")
        ks_val = None
        for name in ("r2_ks_stat", "gof_ks_statistic"):
            if isinstance(result, dict) and name in result:
                ks_val = result[name]; break
        ax2.set_title(f"KS statistic = {ks_val:.3f}" if ks_val is not None
                      else "Age-marginalised CDF comparison", fontsize=9.5)
        ax2.set_xlabel(r"$\Delta$ (repeat units from founder)")
        ax2.set_ylabel("Cumulative probability")
        ax2.legend(fontsize=8)

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
        print(f"Saved: {out_path}")

    return (fig, (ax1, ax2)) if show_cdf else (fig, ax1)


if __name__ == "__main__":
    from regime2_ctmc import CTMCGlobalParams

    rng = np.random.default_rng(7)
    L0 = 44
    # Bimodal-ish: bulk near founder + a subset that crossed T2 and expanded fast
    bulk = np.round(rng.normal(L0 + 3, 2.5, 45)).astype(int)
    tail = np.round(rng.normal(L0 + 60, 15, 12)).astype(int)
    lengths = np.concatenate([bulk, tail])

    gp2 = CTMCGlobalParams.from_preset("htt_cag")
    fig, axes = plot_regime2_fit(lengths, L0, global_params=gp2,
                                  donor_age=60, out_path="regime2_fit_selftest.png",
                                  show_cdf=True)
    print("Self-test complete.")
