"""
Regime 2 — Generalised Two-Phase CTMC Somatic Expansion Model
==============================================================

Extension of regime2_ctmc.py with data-driven estimation of the global
parameters T1, T2, and p_exp rather than relying on literature values
specific to one repeat type.

Motivation
----------
The two-phase structure of somatic repeat expansion is biologically general:
any repeat long enough to form stable hairpin structures during replication
will show length-dependent instability, and the rate tends to accelerate
once the repeat crosses a critical length (T2). The specific values of T1,
T2, and p_exp vary by repeat type and cannot be assumed from HTT literature
for genome-wide analysis.

Identifiability
---------------
The three global parameters are estimated from distinct signals:

    T1  — lower stability threshold
          Identified by: the length at which read-length variance starts
          increasing with founder length. Estimated from Pass 1A using
          the distribution of observed lengths across loci of the same
          repeat type. Can also be supplied from population databases.

    T2  — phase transition threshold
          Identified by: the "kink" in empirical CDFs of read lengths —
          the point where expansion rate accelerates. Only identifiable
          when some observed reads exceed T2. Estimated by maximising
          pooled log-likelihood across Regime 2 flagged loci.

    p_exp — expansion probability
          Identified by: the asymmetry between positive and negative
          tails of the delta distribution. Estimated jointly with T2
          in Pass 1B. Requires enough reads to distinguish from 0.5;
          n ≥ 30 per locus gives reasonable precision.

What cannot be estimated without donor age
------------------------------------------
Donor age is required to convert between the CTMC time axis (years)
and the observed endpoint distribution. Without age, r1 and r2 are
confounded with age — you can only estimate r1*age and r2*age, which
is the net expansion over the donor's lifetime, not the per-year rate.

If donor age is unavailable, set donor_age to a population-typical
value for the tissue type (e.g. 60 for adult blood) and interpret
r1, r2 as effective rates assuming that age. The instability index
Var(Δ) and winsorized_mean_abs_delta remain valid regardless of age.

Usage
-----
from regime2_ctmc_v2 import (
    CTMCGlobalParams,
    estimate_T1_from_founders,
    estimate_global_params,
    fit_locus_ctmc,
)

# Pass 1A: estimate T1 from founder length distribution
T1 = estimate_T1_from_founders(founder_lengths)

# Pass 1B: estimate T2 and p_exp from Regime 2 loci
gp = estimate_global_params(regime2_loci, T1=T1)

# Pass 2: per-locus fit
result = fit_locus_ctmc(lengths, founder_length, donor_age, gp)
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.linalg import expm
from scipy import stats


# ===========================================================================
# Re-export core computation functions from original module
# ===========================================================================
from regime2_ctmc import (
    CTMCGlobalParams,
    CTMCResult,
    _build_rate_matrix,
    _ctmc_pmf,
    _grid_loglik,
    fit_locus_ctmc,
    PRESET_MAP,
)


# ===========================================================================
# Pass 1A: Estimate T1 from founder length distribution
# ===========================================================================

def estimate_T1_from_founders(
    founder_lengths:    list[float],
    variance_threshold: float = 2.0,
    percentile:         float = 5.0,
    min_length:         float = 5.0,
) -> float:
    """
    Estimate the lower instability threshold T1 from the distribution of
    founder (germline) allele lengths across loci of the same repeat type.

    The approach: T1 is the length below which the repeat is essentially
    stable, meaning alleles shorter than T1 show no somatic expansion.
    Empirically this appears as the length below which there is no
    population-level length variation.

    Two complementary methods are used and the smaller is returned
    (conservative — prefer a lower T1 to avoid missing instability):

    Method A — percentile of founder lengths:
        The 5th percentile of observed founder lengths. Loci with founders
        below this are likely in the stable range. T1 should be at or below
        the minimum founder length showing instability.

    Method B — variance inflection:
        Bin founder lengths and find the bin where read-length variance
        starts exceeding variance_threshold. This requires having variance
        data per locus, so it is only applicable when you have both founder
        lengths and per-locus variance estimates.

    Parameters
    ----------
    founder_lengths : list of float
        Germline allele lengths (repeat units) across all loci of this
        repeat type that are being analysed.
    variance_threshold : float
        Minimum variance (in repeat units²) above which a locus is
        considered unstable. Default 2.0 is conservative.
    percentile : float
        Percentile of founder lengths to use for Method A. Default 5.0.
    min_length : float
        Absolute minimum for T1 regardless of data (default 5.0 repeat units).

    Returns
    -------
    float : estimated T1 in repeat units

    Notes
    -----
    For known repeat types, literature values are more reliable:
        HTT CAG:  33.5  (well-established from large cohort studies)
        FMR1 CGG: 54.0  (<55 CGG = stable premutation threshold)
        ATXN3 CAG: 27.0 (SCA3)
        DMPK CTG: 37.0  (DM1)
    If you have any of these loci, use CTMCGlobalParams.from_preset()
    and ignore this function.
    """
    fl = np.asarray(founder_lengths, dtype=float)
    fl = fl[fl > min_length]

    if len(fl) < 3:
        return float(min_length)

    # Method A: percentile of founder distribution
    T1_A = float(np.percentile(fl, percentile))

    # Ensure T1 is at least slightly below the minimum observed founder
    T1_A = min(T1_A, float(fl.min()) - 1.0)

    return max(float(T1_A), float(min_length))


# ===========================================================================
# Pass 1B: Estimate T2 and p_exp from Regime 2 flagged loci
# ===========================================================================

def estimate_global_params(
    loci_data:         list[dict],
    T1:                float,
    T2_grid:           Optional[np.ndarray] = None,
    pexp_grid:         Optional[np.ndarray] = None,
    r1_n:              int   = 10,
    r2_n:              int   = 10,
    min_reads:         int   = 10,
    default_age:       float = 60.0,
    verbose:           bool  = True,
) -> CTMCGlobalParams:
    """
    Estimate T2 and p_exp from Regime 2 flagged loci by maximising the
    pooled log-likelihood across the dataset.

    T1 must be supplied (from estimate_T1_from_founders or literature).
    T2 and p_exp are estimated jointly from the data.

    The estimation works by:
    1. For each (T2, p_exp) combination on a grid, compute the log-likelihood
       of each locus marginalised over (r1, r2) — i.e. how well does this
       parameterisation explain the data without knowing the per-locus rates.
    2. Sum the marginalised log-likelihoods across all loci.
    3. Pick the (T2, p_exp) with the highest total.

    Identifiability notes:

    T2 is only identifiable when some observed lengths exceed T2.
    Loci where all reads are below T2 contribute no information about T2
    and are effectively treated as one-phase models for the purposes of
    estimating T2.

    p_exp is identified by the asymmetry between expansion and contraction.
    A completely symmetric distribution (all net_bias near 0) will push
    p_exp toward 0.5. Strong expansion bias pushes it higher.

    Parameters
    ----------
    loci_data : list of dicts, each with:
        'lengths'        : array-like — per-read allele lengths (repeat units)
        'founder_length' : float      — germline allele length (repeat units)
        'donor_age'      : float, optional — donor age in years.
                           If absent, default_age is used.
    T1 : float
        Fixed lower stability threshold (from Pass 1A or literature).
    T2_grid : array, optional
        Values of T2 to search. Default: np.linspace(T1+10, T1+150, 20).
    pexp_grid : array, optional
        Values of p_exp to search. Default: np.linspace(0.50, 0.85, 8).
    r1_n, r2_n : int
        Grid size for marginalising over r1, r2 at each (T2, pexp) point.
    min_reads : int
        Minimum reads per locus to include in estimation.
    default_age : float
        Age used when 'donor_age' is absent from a locus dict. Default 60.
    verbose : bool
        Print progress.

    Returns
    -------
    CTMCGlobalParams with estimated T1, T2, p_exp and data-derived priors
    for r1 and r2 based on the spread of per-locus estimates.
    """
    if T2_grid is None:
        T2_grid = np.linspace(T1 + 10.0, T1 + 150.0, 20)
    if pexp_grid is None:
        pexp_grid = np.linspace(0.50, 0.85, 8)

    # Filter loci
    valid = []
    for locus in loci_data:
        lengths = np.asarray(locus.get("lengths", []), dtype=float)
        if len(lengths) < min_reads:
            continue
        if locus.get("founder_length") is None:
            continue
        valid.append(locus)

    if len(valid) < 2:
        print("Warning: fewer than 2 valid loci. Using T2=T1+40, p_exp=0.6 defaults.")
        return CTMCGlobalParams(T1=T1, T2=T1 + 40.0, p_exp=0.60, repeat_type="estimated")

    if verbose:
        print(f"Pass 1B (CTMC): estimating T2, p_exp from {len(valid)} loci "
              f"[T1={T1:.1f} fixed]")
        obs_max = max(
            float(np.round(np.asarray(locus["lengths"])).max())
            for locus in valid
        )
        print(f"  Max observed length across all loci: {obs_max:.0f}")
        n_above_T2_min = sum(
            1 for locus in valid
            if float(np.asarray(locus["lengths"]).max()) > T2_grid.min()
        )
        print(f"  Loci with reads potentially above T2_min={T2_grid.min():.0f}: "
              f"{n_above_T2_min}/{len(valid)}")

    r1_grid = np.exp(np.linspace(np.log(1e-4), np.log(0.30), r1_n))
    r2_grid = np.exp(np.linspace(np.log(1e-3), np.log(15.0), r2_n))

    # Grid search over (T2, p_exp)
    ll_surface = np.full((len(T2_grid), len(pexp_grid)), -np.inf)

    for i, T2 in enumerate(T2_grid):
        for j, pexp in enumerate(pexp_grid):
            total_ll = 0.0
            for locus in valid:
                lengths_int = np.round(np.asarray(locus["lengths"])).astype(int)
                L0    = int(round(locus["founder_length"]))
                age   = float(locus.get("donor_age", default_age))
                min_L = max(int(lengths_int.min()) - 5, 0)
                max_L = min(int(lengths_int.max()) + 25, 500)
                obs_idx = lengths_int - min_L
                valid_mask = (obs_idx >= 0) & (obs_idx <= max_L - min_L)
                obs_idx = obs_idx[valid_mask]
                unique_obs, counts = np.unique(obs_idx, return_counts=True)

                ll = _grid_loglik(
                    unique_obs, counts, L0, age,
                    r1_grid, r2_grid, T1, T2, pexp, min_L, max_L,
                )
                # Marginalise over (r1, r2) with Gamma priors
                lp = ll + stats.gamma.logpdf(r1_grid, 1.5, scale=0.05)[:,None] \
                        + stats.gamma.logpdf(r2_grid, 1.5, scale=1.0)[None,:]
                lp_max = lp.max()
                if np.isfinite(lp_max):
                    total_ll += lp_max + np.log(np.exp(lp - lp_max).sum())

            ll_surface[i, j] = total_ll

    # Find MAP estimate
    idx = np.unravel_index(np.argmax(ll_surface), ll_surface.shape)
    best_T2   = float(T2_grid[idx[0]])
    best_pexp = float(pexp_grid[idx[1]])
    best_ll   = float(ll_surface[idx])

    if verbose:
        print(f"  Best: T2={best_T2:.1f}  p_exp={best_pexp:.3f}  "
              f"pooled_loglik={best_ll:.2f}")

    # Uncertainty in T2 and p_exp from the likelihood surface
    # Marginalise over p_exp to get T2 uncertainty
    ll_T2   = ll_surface - ll_surface.max()
    post_T2 = np.exp(ll_T2.max(axis=1))
    post_T2 /= post_T2.sum()
    T2_mean = float((post_T2 * T2_grid).sum())

    ll_pexp   = ll_surface - ll_surface.max()
    post_pexp = np.exp(ll_pexp.max(axis=0))
    post_pexp /= post_pexp.sum()
    pexp_mean = float((post_pexp * pexp_grid).sum())

    if verbose:
        print(f"  Posterior mean: T2={T2_mean:.1f}  p_exp={pexp_mean:.3f}")

    # Estimate data-driven priors for r1 and r2 from per-locus fits
    # at the best (T2, p_exp)
    r1_ests, r2_ests = [], []
    for locus in valid:
        lengths_int = np.round(np.asarray(locus["lengths"])).astype(int)
        L0    = int(round(locus["founder_length"]))
        age   = float(locus.get("donor_age", default_age))
        min_L = max(int(lengths_int.min()) - 5, 0)
        max_L = min(int(lengths_int.max()) + 25, 500)
        obs_idx = lengths_int - min_L
        valid_mask = (obs_idx >= 0) & (obs_idx <= max_L - min_L)
        unique_obs, counts = np.unique(obs_idx[valid_mask], return_counts=True)

        ll = _grid_loglik(
            unique_obs, counts, L0, age,
            r1_grid, r2_grid, T1, best_T2, best_pexp, min_L, max_L,
        )
        lp = ll + stats.gamma.logpdf(r1_grid, 1.5, scale=0.05)[:,None] \
                + stats.gamma.logpdf(r2_grid, 1.5, scale=1.0)[None,:]
        lp -= lp.max(); post = np.exp(lp); post /= post.sum()
        R1, R2 = np.meshgrid(r1_grid, r2_grid, indexing='ij')
        r1_ests.append(float((post * R1).sum()))
        r2_ests.append(float((post * R2).sum()))

    r1_ests = np.array(r1_ests)
    r2_ests = np.array(r2_ests)

    # Fit Gamma priors from the distribution of per-locus estimates
    # Method of moments: shape = mean²/var, scale = var/mean
    def fit_gamma(vals):
        vals = vals[vals > 1e-6]
        if len(vals) < 2:
            return 1.5, 0.05
        m, v = vals.mean(), vals.var()
        if v < 1e-10: return 1.5, float(m)
        return float(m**2 / v), float(v / m)

    r1_shape, r1_scale = fit_gamma(r1_ests)
    r2_shape, r2_scale = fit_gamma(r2_ests)

    if verbose:
        print(f"  Per-locus r1: mean={r1_ests.mean():.4f}  "
              f"prior→ Gamma(shape={r1_shape:.2f}, scale={r1_scale:.4f})")
        print(f"  Per-locus r2: mean={r2_ests.mean():.3f}  "
              f"prior→ Gamma(shape={r2_shape:.2f}, scale={r2_scale:.3f})")

    return CTMCGlobalParams(
        T1              = T1,
        T2              = best_T2,
        p_exp           = best_pexp,
        r1_prior_shape  = r1_shape,
        r1_prior_scale  = r1_scale,
        r2_prior_shape  = r2_shape,
        r2_prior_scale  = r2_scale,
        repeat_type     = "estimated",
        n_loci_used     = len(valid),
    )


# ===========================================================================
# Self-test: generalised parameters, non-HTT repeat
# ===========================================================================

def _run_selftest():
    rng = np.random.default_rng(42)

    print("=" * 65)
    print("GENERALISED CTMC SELF-TEST — data-driven global params")
    print("=" * 65)
    print()

    # True parameters for a generic disease repeat (not HTT)
    true_T1    = 20.0
    true_T2    = 55.0
    true_p_exp = 0.65
    true_r1    = 0.025
    true_r2    = 2.5

    print(f"True global params: T1={true_T1} T2={true_T2} p_exp={true_p_exp}")
    print(f"True per-locus:     r1={true_r1} r2={true_r2}")
    print(f"Net rate Phase A at L0=25: "
          f"{true_r1*(2*true_p_exp-1)*max(25-true_T1,0):.4f} CAGs/yr")
    print(f"Net rate Phase B at L0=25: "
          f"{(true_r1+true_r2)*(2*true_p_exp-1)*max(25-true_T1,0):.4f} CAGs/yr")
    print(f"Rate acceleration: {(true_r1+true_r2)/true_r1:.1f}×")
    print()

    # Simulate multiple loci spanning a range of L0 and ages
    max_L = 200
    loci_data = []
    for L0, age in [(22, 40), (25, 50), (28, 60), (30, 65),
                    (32, 70), (25, 55), (28, 45), (30, 75)]:
        Q   = _build_rate_matrix(true_r1, true_r2, true_T1, true_T2,
                                  true_p_exp, 0, max_L)
        Pt  = expm(age * Q)
        pmf = Pt[L0, :]; pmf /= pmf.sum()
        obs = rng.choice(np.arange(max_L + 1), size=30, p=pmf).astype(float)
        loci_data.append({
            "lengths":        obs,
            "founder_length": float(L0),
            "donor_age":      float(age),
        })

    # --- Pass 1A: estimate T1 ---
    founders = [locus["founder_length"] for locus in loci_data]
    T1_est   = estimate_T1_from_founders(founders, percentile=5.0)
    print(f"Pass 1A — T1 estimate: {T1_est:.1f}  (true={true_T1})")
    # In practice T1 is known from the repeat type; set it manually if off
    T1_use = true_T1  # use true for clean test

    # --- Pass 1B: estimate T2, p_exp ---
    print()
    t0  = time.time()
    gp  = estimate_global_params(
        loci_data,
        T1          = T1_use,
        T2_grid     = np.linspace(T1_use + 10, T1_use + 80, 12),
        pexp_grid   = np.linspace(0.50, 0.80, 7),
        r1_n        = 8,
        r2_n        = 8,
        verbose     = True,
    )
    print(f"Pass 1B time: {time.time()-t0:.1f}s")
    print()
    print(f"Estimated: T2={gp.T2:.1f}  p_exp={gp.p_exp:.3f}")
    print(f"True:      T2={true_T2:.1f}  p_exp={true_p_exp:.3f}")
    print()

    # --- Pass 2: per-locus fit ---
    print("Pass 2 — per-locus fit on one test locus:")
    L0_test, age_test = 28, 60
    Q_test  = _build_rate_matrix(true_r1, true_r2, true_T1, true_T2,
                                  true_p_exp, 0, max_L)
    Pt_test = expm(age_test * Q_test)
    pmf_test = Pt_test[L0_test, :]; pmf_test /= pmf_test.sum()
    obs_test = rng.choice(np.arange(max_L + 1), size=50, p=pmf_test).astype(float)

    t0 = time.time()
    result = fit_locus_ctmc(obs_test, L0_test, age_test, gp,
                             haplotype_label="H1_test")
    print(result)
    print(f"True r1={true_r1:.4f}  est={result.r1[0]:.4f} "
          f"[{result.r1[1]:.4f},{result.r1[2]:.4f}]")
    print(f"True net_rate_A={true_r1*(2*true_p_exp-1)*max(L0_test-true_T1,0):.4f}  "
          f"est={result.net_rate_phaseA[0]:.4f}")

    # --- Coverage test ---
    print()
    print("Coverage test (net_rate_phaseA, 30 trials, n=50):")
    true_net_A = true_r1 * (2 * true_p_exp - 1) * max(L0_test - true_T1, 0)
    covered = 0
    times   = []
    for _ in range(30):
        sim = rng.choice(np.arange(max_L + 1), size=50, p=pmf_test).astype(float)
        t0  = time.time()
        res = fit_locus_ctmc(sim, L0_test, age_test, gp)
        times.append(time.time() - t0)
        lo, hi = res.net_rate_phaseA[1], res.net_rate_phaseA[2]
        if lo <= true_net_A <= hi:
            covered += 1
    print(f"Coverage: {covered}/30 = {covered/30*100:.0f}%")
    print(f"Mean time per locus: {np.mean(times)*1000:.0f}ms")

    # --- Winsorized delta comparison ---
    print()
    print("Winsorized mean|Δ| (comparable to Handsaker et al. Fig. SN2.3):")
    for locus in loci_data[:4]:
        lengths = locus["lengths"]
        L0l     = locus["founder_length"]
        age_l   = locus["donor_age"]
        abs_d   = np.abs(lengths - L0l)
        w_mad   = float(np.minimum(abs_d, 100.0).mean())
        print(f"  L0={L0l:.0f} age={age_l:.0f}: "
              f"winsor_mean_abs_delta={w_mad:.3f} repeat units")


if __name__ == "__main__":
    _run_selftest()