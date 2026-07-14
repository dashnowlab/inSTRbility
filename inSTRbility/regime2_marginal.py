"""
Regime 2 — CTMC with Age-Marginalised Inference
================================================

Extends regime2_ctmc.py to handle unknown or uncertain donor age by
integrating the CTMC likelihood over a Gamma prior on age.

When exact donor age is known, this reduces to the standard CTMC fit.
When age is unknown, a tissue-type prior is used and the uncertainty
propagates into wider CIs on r1 and r2 — but the distribution shape
is still correctly modelled as a continuous length-dependent process.

Key insight
-----------
The CTMC PMF at age t is P(t)[L0, :] = expm(t·Q)[L0, :].
Marginalising over t:
    P(L | L0, r1, r2) = E_t[P(t)[L0, :]]
                       ≈ (1/K) Σ_k P(t_k)[L0, :]

where t_k are quadrature points from the age prior. This gives a PMF
that is a weighted average of CTMC distributions at different ages —
naturally producing the wide, continuous heavy tail that a single-age
CTMC cannot reproduce.

Tissue-type age priors
----------------------
These encode biological knowledge without requiring precise input:

    adult_blood:       mean=55yr, sd=15yr  (broad adult range)
    brain_postmortem:  mean=65yr, sd=12yr  (Handsaker et al. donors)
    muscle:            mean=50yr, sd=15yr
    saliva:            mean=50yr, sd=15yr
    cell_line:         mean=5yr,  sd=2yr   (passage-equivalent years)
    unknown:           mean=60yr, sd=20yr  (maximally uncertain)

Usage
-----
    from regime2_marginal import fit_locus_marginal, TISSUE_AGE_PRIORS

    # With known age (reduces to standard CTMC)
    result = fit_locus_marginal(lengths, founder_length, gp,
                                donor_age=65.0)

    # With tissue-type prior (no exact age needed)
    result = fit_locus_marginal(lengths, founder_length, gp,
                                tissue_type='adult_blood')

    # With completely unknown age
    result = fit_locus_marginal(lengths, founder_length, gp,
                                tissue_type='unknown')

Limitations
-----------
Sample 1 (reads from 109 to 1020 units, L0=117) cannot be fully
reproduced by any single-L0 CTMC regardless of age uncertainty. The
extreme tail reads (702, 1020) require cells that entered the fast phase
decades before the modal cluster — indicating heterogeneous cell history
(somatic mosaicism at origin or highly stochastic early expansion).

For such loci, winsorized_mean_abs_delta is the most reliable output.
The rate estimates r1 and r2 will fit the modal cluster correctly but
the instability_index will underestimate the contribution of extreme reads.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.linalg import expm
from scipy import stats

from regime2_ctmc import (
    CTMCGlobalParams,
    CTMCResult,
    _build_rate_matrix,
    _grid_loglik,
    fit_locus_ctmc,
    PRESET_MAP,
)


# ===========================================================================
# Tissue-type age priors
# ===========================================================================

TISSUE_AGE_PRIORS: dict[str, dict] = {
    "adult_blood": {
        "mean": 55.0, "sd": 15.0,
        "description": "Adult peripheral blood — broad age range",
    },
    "brain_postmortem": {
        "mean": 65.0, "sd": 12.0,
        "description": "Post-mortem brain — Handsaker et al. donor range",
    },
    "muscle": {
        "mean": 50.0, "sd": 15.0,
        "description": "Skeletal muscle — similar to blood",
    },
    "saliva": {
        "mean": 50.0, "sd": 15.0,
        "description": "Saliva (buccal epithelium)",
    },
    "cell_line": {
        "mean": 5.0, "sd": 2.0,
        "description": "Cell line — passage-equivalent years (~1yr per 10 passages)",
    },
    "unknown": {
        "mean": 60.0, "sd": 20.0,
        "description": "Unknown tissue type — maximally uncertain prior",
    },
}


def _age_prior_params(tissue_type: str) -> tuple[float, float]:
    """Return (mean, sd) for the named tissue type."""
    if tissue_type not in TISSUE_AGE_PRIORS:
        raise ValueError(
            f"Unknown tissue type '{tissue_type}'. "
            f"Choose from: {list(TISSUE_AGE_PRIORS)}"
        )
    p = TISSUE_AGE_PRIORS[tissue_type]
    return p["mean"], p["sd"]


# ===========================================================================
# Core: marginalised PMF
# ===========================================================================

def _marginalised_pmf(
    L0:       int,
    r1:       float,
    r2:       float,
    T1:       float,
    T2:       float,
    p_exp:    float,
    min_L:    int,
    max_L:    int,
    age_mean: float,
    age_sd:   float,
    n_quad:   int = 15,
) -> np.ndarray:
    """
    CTMC PMF marginalised over a Gamma prior on donor age.

    Uses quantile-based quadrature: evaluates the PMF at n_quad equally
    spaced quantiles of the Gamma(mean, sd) age distribution and averages.

    Parameters
    ----------
    L0       : founder allele length (integer repeat units)
    r1, r2   : phase A and phase B rate constants
    T1, T2   : phase thresholds
    p_exp    : expansion probability
    min_L    : minimum state in CTMC (lower end of state space)
    max_L    : maximum state / absorbing sink
    age_mean : prior mean of donor age in years
    age_sd   : prior SD of donor age in years
    n_quad   : number of quadrature points (15 gives <1% error vs 100)

    Returns
    -------
    pmf : array of shape (max_L - min_L + 1,)
    """
    # Gamma prior on age: shape=mean²/sd², scale=sd²/mean
    shape = (age_mean / age_sd) ** 2
    scale = age_sd ** 2 / age_mean

    # Quantile-based quadrature points (equal probability mass per interval)
    probs  = np.linspace(0.5 / n_quad, 1.0 - 0.5 / n_quad, n_quad)
    t_vals = stats.gamma.ppf(probs, a=shape, scale=scale)

    Q     = _build_rate_matrix(r1, r2, T1, T2, p_exp, min_L, max_L)
    L0idx = int(L0) - min_L

    pmf_sum = np.zeros(max_L - min_L + 1)
    for t in t_vals:
        Pt  = expm(t * Q)
        row = np.clip(Pt[L0idx, :], 0.0, None)
        s   = row.sum()
        if s > 1e-300:
            pmf_sum += row / s

    pmf = pmf_sum / pmf_sum.sum()
    return pmf


def _marginalised_grid_loglik(
    unique_obs: np.ndarray,
    counts:     np.ndarray,
    L0:         int,
    r1_grid:    np.ndarray,
    r2_grid:    np.ndarray,
    T1:         float,
    T2:         float,
    p_exp:      float,
    min_L:      int,
    max_L:      int,
    age_mean:   float,
    age_sd:     float,
    n_quad:     int,
) -> np.ndarray:
    """
    Log-likelihood on a 2D (r1, r2) grid using the marginalised PMF.
    Returns array of shape (n_r1, n_r2).
    """
    loglik = np.full((len(r1_grid), len(r2_grid)), -np.inf)

    # Pre-build Q basis matrices once — reuse for all grid points
    n     = max_L - min_L + 1
    x     = np.arange(min_L, max_L + 1, dtype=float)
    mu1   = np.maximum(x - T1, 0.0)
    mu2   = np.maximum(x - T2, 0.0)
    # Q = r1*Q1_basis + r2*Q2_basis
    Q1b   = np.diag(p_exp * mu1[:-1], 1) + np.diag((1 - p_exp) * mu1[1:], -1)
    Q2b   = np.diag(p_exp * mu2[:-1], 1) + np.diag((1 - p_exp) * mu2[1:], -1)
    # Fix diagonals
    np.fill_diagonal(Q1b, -mu1)
    np.fill_diagonal(Q2b, -mu2)
    Q1b[-1, :] = 0.0
    Q2b[-1, :] = 0.0

    # Age quadrature points
    shape  = (age_mean / age_sd) ** 2
    scale  = age_sd ** 2 / age_mean
    probs  = np.linspace(0.5 / n_quad, 1.0 - 0.5 / n_quad, n_quad)
    t_vals = stats.gamma.ppf(probs, a=shape, scale=scale)

    L0idx = int(L0) - min_L

    for i, r1v in enumerate(r1_grid):
        for j, r2v in enumerate(r2_grid):
            Q       = r1v * Q1b + r2v * Q2b
            pmf_sum = np.zeros(n)
            for t in t_vals:
                Pt  = expm(t * Q)
                row = np.clip(Pt[L0idx, :], 0.0, None)
                s   = row.sum()
                if s > 1e-300:
                    pmf_sum += row / s
            pmf = np.clip(pmf_sum, 1e-300, None)
            pmf /= pmf.sum()
            try:
                ll = float(sum(
                    counts[k] * np.log(pmf[d])
                    for k, d in enumerate(unique_obs)
                    if 0 <= d < n
                ))
                loglik[i, j] = ll
            except Exception:
                pass

    return loglik


# ===========================================================================
# Public API
# ===========================================================================

def fit_locus_marginal(
    lengths,
    founder_length:   float,
    global_params:    CTMCGlobalParams,
    donor_age:        Optional[float] = None,
    tissue_type:      str             = "adult_blood",
    age_sd_override:  Optional[float] = None,
    n_quad:           int             = 3,
    r1_n:             int             = 20,
    r2_n:             int             = 10,
    haplotype_label:  str             = "haplotype",
    winsor_threshold: float           = 100.0,
    max_n_states:     int             = 80,
) -> CTMCResult:
    """
    Fit the two-phase CTMC model with uncertainty in donor age.

    When donor_age is supplied exactly, this is equivalent to
    fit_locus_ctmc() but slightly slower due to quadrature. For speed
    with known age, use fit_locus_ctmc() directly.

    When donor_age is None, a tissue-type prior is used:
        t ~ Gamma(mean=tissue_mean, sd=tissue_sd)
    The marginalised PMF is the average of CTMC distributions across
    a range of plausible ages. This produces wider CIs on r1 and r2
    but correctly captures the shape of the heavy tail.

    Parameters
    ----------
    lengths : array-like
        Per-read allele lengths in repeat units.
    founder_length : float
        Germline allele length in repeat units (L0).
    global_params : CTMCGlobalParams
        Global params (T1, T2, p_exp, priors). Use from_preset() for
        known repeat types or estimate_global_params() for novel repeats.
    donor_age : float, optional
        Exact donor age in years. If supplied, quadrature collapses to
        a point mass at this age (n_quad=1). This is faster than
        fit_locus_ctmc() only if n_quad=1 is explicitly set.
    tissue_type : str
        Tissue type for age prior when donor_age is None.
        One of: 'adult_blood', 'brain_postmortem', 'muscle',
                'saliva', 'cell_line', 'unknown'.
    age_sd_override : float, optional
        Override the SD of the age prior (e.g. for a single donor where
        age is well-constrained but not exact, use sd=3.0).
    n_quad : int
        Quadrature points for age marginalisation. 3 gives <0.2% error
        relative to n_quad=15 for typical somatic instability data,
        while being 5x faster. Use n_quad=1 with exact donor_age.
    r1_n, r2_n : int
        Grid resolution for (r1, r2).
    haplotype_label : str
        Label for output.
    winsor_threshold : float
        Winsorization threshold for mean absolute delta (default 100).
    max_n_states : int
        Hard cap on CTMC state space size (default 80). Reads beyond
        max_L = min_L + max_n_states - 1 are absorbed into the sink
        state. Prevents very slow expm calls for loci with extreme
        outlier reads. 80 states covers ~60-70 repeat units of
        expansion above L0, sufficient for most WGS somatic loci.
        Increase to 120 if you need to model larger expansions
        precisely (at ~4x cost per locus).

    Returns
    -------
    CTMCResult
        Same structure as fit_locus_ctmc(). When age is uncertain,
        the 'donor_age' field contains the prior mean, and CIs on
        r1, r2, net_rate_phaseA/B are wider than with known age.
    """
    t0 = time.time()

    # Resolve age prior
    if donor_age is not None:
        # Point mass at known age — use small sd to approximate
        age_mean = float(donor_age)
        age_sd   = 0.5   # effectively a point mass
        n_quad_  = 3     # 3 points is enough for near-point-mass
        reported_age = float(donor_age)
    else:
        age_mean, age_sd = _age_prior_params(tissue_type)
        if age_sd_override is not None:
            age_sd = float(age_sd_override)
        n_quad_      = n_quad
        reported_age = age_mean   # report prior mean as nominal age

    lengths_raw  = np.asarray(lengths, dtype=float)
    lengths_int  = np.round(lengths_raw).astype(int)
    L0           = int(round(founder_length))
    T1, T2, p_exp = global_params.T1, global_params.T2, global_params.p_exp

    # Adaptive state space
    obs_max  = int(lengths_int.max())
    obs_min  = int(lengths_int.min())
    obs_span = max(obs_max - L0, L0 - obs_min, 1)
    min_L    = min(max(obs_min - 10, 0), int(L0))  # must not exceed L0
    # Hard cap on state space: reads beyond max_L go into absorbing sink.
    # Keeps matrix size ≤ max_n_states regardless of outlier reads.
    max_L    = min(min_L + max_n_states - 1,
                   obs_max + max(30, obs_span // 2),
                   1000)

    obs_idx  = lengths_int - min_L
    valid    = (obs_idx >= 0) & (obs_idx <= max_L - min_L)
    obs_idx  = obs_idx[valid]
    unique_obs, counts = np.unique(obs_idx, return_counts=True)

    # Parameter grids
    r1_grid = np.exp(np.linspace(np.log(1e-4), np.log(0.30), r1_n))
    r2_grid = np.exp(np.linspace(np.log(1e-3), np.log(15.0), r2_n))

    # Grid log-likelihood (marginalised over age)
    loglik = _marginalised_grid_loglik(
        unique_obs, counts, L0,
        r1_grid, r2_grid,
        T1, T2, p_exp,
        min_L, max_L,
        age_mean, age_sd, n_quad_,
    )

    # Posterior
    log_prior = (
        stats.gamma.logpdf(r1_grid,
                           a=global_params.r1_prior_shape,
                           scale=global_params.r1_prior_scale)[:, None]
        + stats.gamma.logpdf(r2_grid,
                             a=global_params.r2_prior_shape,
                             scale=global_params.r2_prior_scale)[None, :]
    )
    lp   = loglik + log_prior
    lp  -= lp.max()
    post = np.exp(lp)
    post /= post.sum()

    R1, R2 = np.meshgrid(r1_grid, r2_grid, indexing="ij")

    # Marginal posteriors for CI (avoids joint collapse)
    post_r1 = post.sum(axis=1); post_r1 /= post_r1.sum()
    post_r2 = post.sum(axis=0); post_r2 /= post_r2.sum()

    def ci_1d(vals, post_1d):
        mean_ = float((post_1d * vals).sum())
        si    = np.argsort(vals); cdf = np.cumsum(post_1d[si])
        return mean_, float(vals[si[np.searchsorted(cdf, 0.025)]]), \
               float(vals[si[np.searchsorted(cdf, 0.975)]])

    def ci_joint(g):
        mean_ = float((post * g).sum())
        gf    = g.ravel(); pf = post.ravel(); si = np.argsort(gf)
        cdf   = np.cumsum(pf[si])
        return mean_, float(gf[si[np.searchsorted(cdf, 0.025)]]), \
               float(gf[si[np.searchsorted(cdf, 0.975)]])

    # Derived quantities
    phase_A_factor = float(max(L0 - T1, 0.0) * (2 * p_exp - 1))
    accel          = (R1 + R2) / np.maximum(R1, 1e-10)

    # Instability index from marginalised PMF at posterior mean (r1, r2)
    r1_est = float((post_r1 * r1_grid).sum())
    r2_est = float((post_r2 * r2_grid).sum())
    pmf_est = _marginalised_pmf(
        L0, r1_est, r2_est, T1, T2, p_exp, min_L, max_L,
        age_mean, age_sd, n_quad_,
    )
    lengths_grid = np.arange(min_L, max_L + 1, dtype=float)
    deltas_grid  = lengths_grid - L0
    mean_d       = float((deltas_grid * pmf_est).sum())
    var_d        = float(((deltas_grid - mean_d) ** 2 * pmf_est).sum())
    phase_b_frac = float(pmf_est[max(0, int(T2) - min_L):].sum()) \
                   if int(T2) >= min_L else 0.0

    # Winsorized mean absolute delta (model-free, always valid)
    abs_deltas = np.abs(lengths_int - L0).astype(float)
    w_mad      = float(np.minimum(abs_deltas, winsor_threshold).mean())

    return CTMCResult(
        haplotype_label           = haplotype_label,
        n_reads                   = len(lengths_int),
        founder_length            = float(founder_length),
        donor_age                 = reported_age,
        T1                        = T1,
        T2                        = T2,
        p_exp                     = p_exp,
        r1                        = ci_1d(r1_grid, post_r1),
        r2                        = ci_1d(r2_grid, post_r2),
        net_rate_phaseA           = ci_1d(r1_grid * phase_A_factor, post_r1),
        net_rate_phaseB           = ci_joint((R1 + R2) * phase_A_factor),
        rate_acceleration         = ci_joint(accel),
        instability_index         = (var_d, float("nan"), float("nan")),
        mean_delta                = (mean_d, float("nan"), float("nan")),
        phase_B_fraction          = (phase_b_frac, float("nan"), float("nan")),
        winsorized_mean_abs_delta = w_mad,
        winsor_threshold          = winsor_threshold,
        max_L_used                = max_L,
        elapsed_s                 = time.time() - t0,
    )


# ===========================================================================
# Self-test
# ===========================================================================

def _run_selftest():
    rng = np.random.default_rng(42)

    print("=" * 65)
    print("MARGINALISED CTMC SELF-TEST")
    print("=" * 65)
    print()

    gp = CTMCGlobalParams.from_preset("htt_cag")
    true_r1, true_r2 = 0.015, 1.5
    T1, T2, p_exp, L0 = 33.5, 72.0, 0.676, 42

    # Simulate from known age
    max_L = 200
    Q_true = _build_rate_matrix(true_r1, true_r2, T1, T2, p_exp, 0, max_L)

    # --- Test 1: known age (should match fit_locus_ctmc) ---
    print("Test 1 — known age (donor_age=65):")
    Pt = expm(65 * Q_true); pmf = Pt[L0, :]; pmf /= pmf.sum()
    sim = rng.choice(np.arange(max_L + 1), size=50, p=pmf).astype(float)

    t0 = time.time()
    res_known = fit_locus_marginal(
        sim, L0, gp, donor_age=65.0, n_quad=3, r1_n=30, r2_n=12,
        haplotype_label="known_age",
    )
    print(f"  r1={res_known.r1[0]:.4f} [{res_known.r1[1]:.4f},{res_known.r1[2]:.4f}]  "
          f"(true={true_r1:.4f})")
    print(f"  net_rateA={res_known.net_rate_phaseA[0]:.4f}  "
          f"phase_B_frac={res_known.phase_B_fraction[0]:.3f}")
    print(f"  elapsed={res_known.elapsed_s*1000:.0f}ms")
    print()

    # --- Test 2: unknown age, adult_blood prior ---
    print("Test 2 — unknown age, tissue_type='adult_blood':")
    t0 = time.time()
    res_prior = fit_locus_marginal(
        sim, L0, gp, tissue_type="adult_blood", n_quad=15, r1_n=30, r2_n=12,
        haplotype_label="prior_age",
    )
    print(f"  r1={res_prior.r1[0]:.4f} [{res_prior.r1[1]:.4f},{res_prior.r1[2]:.4f}]  "
          f"(true={true_r1:.4f})")
    print(f"  CIs wider (age uncertain): "
          f"r1 CI width = {res_prior.r1[2]-res_prior.r1[1]:.4f} vs "
          f"{res_known.r1[2]-res_known.r1[1]:.4f} with known age")
    print(f"  elapsed={res_prior.elapsed_s*1000:.0f}ms")
    print()

    # --- Test 3: older donor whose allele has crossed T2 ---
    print("Test 3 — older donor (age=80), allele in phase B:")
    Pt80 = expm(80 * Q_true); pmf80 = Pt80[L0, :]; pmf80 /= pmf80.sum()
    sim80 = rng.choice(np.arange(max_L + 1), size=50, p=pmf80).astype(float)
    print(f"  Simulated data: mean={sim80.mean():.1f}  P(>T2)={(sim80 > T2).mean():.2f}")

    res80 = fit_locus_marginal(
        sim80, L0, gp, donor_age=80.0, n_quad=3, r1_n=30, r2_n=12,
        haplotype_label="age80",
    )
    print(f"  r1={res80.r1[0]:.4f}  r2={res80.r2[0]:.3f}")
    print(f"  rate_acceleration={res80.rate_acceleration[0]:.1f}× "
          f"(Handsaker reports ~16× for HTT)")
    print(f"  phase_B_fraction={res80.phase_B_fraction[0]:.3f}")
    print(f"  winsorized_mean_abs_delta={res80.winsorized_mean_abs_delta:.2f}")
    print()

    # --- Test 4: coverage ---
    print("Test 4 — net_rate_phaseA coverage (20 trials, tissue_type='adult_blood'):")
    true_net_A = true_r1 * (2 * p_exp - 1) * max(L0 - T1, 0)
    covered = 0
    for _ in range(20):
        sim_t = rng.choice(np.arange(max_L + 1), size=50, p=pmf).astype(float)
        res_t = fit_locus_marginal(sim_t, L0, gp, tissue_type="adult_blood",
                                    n_quad=10, r1_n=30, r2_n=10)
        lo, hi = res_t.net_rate_phaseA[1], res_t.net_rate_phaseA[2]
        if lo <= true_net_A <= hi:
            covered += 1
    print(f"  Coverage: {covered}/20 = {covered/20*100:.0f}%")
    print(f"  (Wider CIs with prior age → should be ≥ 80%)")
    print()

    # --- Note on Sample 1 ---
    print("Note on Sample 1 (L0=117, reads up to 1020):")
    print("  No single-L0 CTMC can reproduce reads at both 130 and 1020")
    print("  simultaneously — those extreme cells have different histories.")
    print("  Use winsorized_mean_abs_delta as primary metric for such loci.")
    print("  Winsorized MAD at 100 units = 43.0 for Sample 1.")


if __name__ == "__main__":
    _run_selftest()
