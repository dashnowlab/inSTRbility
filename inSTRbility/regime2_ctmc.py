"""
Regime 2 — Two-Phase CTMC Somatic Expansion Model
===================================================

Replaces the previous MDN-based Regime 2 with a CTMC (continuous-time
Markov chain) implementation following Handsaker et al. (2024) Cell.

Biological model
----------------
The HTT CAG repeat (and other disease repeats) undergoes length-dependent
somatic expansion. The mutation rate accelerates sharply at a second
threshold (T2 ≈ 72-80 CAGs), creating two phases:

    Phase A (T1 < L < T2): slow expansion
        μ(L) = r1 * (L − T1)
        net rate ≈ 3.5% CAGs/year at typical HTT alleles

    Phase B (L ≥ T2): fast expansion
        μ(L) = r1*(L−T1) + r2*(L−T2)
        net rate ≈ 57.6% CAGs/year — 16× faster than Phase A

The CTMC evolves each cell's repeat length over the donor's lifetime:
    Q[L, L+1] = p_exp * μ(L)     — expansion
    Q[L, L-1] = (1−p_exp) * μ(L) — contraction
    P(t) = expm(t · Q)

The PMF at age t is the L0-th row of P(t): P(t)[L0, :].

Global parameters (per repeat type, not per sample)
----------------------------------------------------
    T1     — lower instability threshold (fixed by biology)
    T2     — phase transition threshold (estimated globally)
    p_exp  — expansion probability (remarkably stable across donors)

Per-locus parameters (estimated per sample per locus)
-----------------------------------------------------
    r1     — phase A rate constant
    r2     — phase B rate constant

Key outputs
-----------
    net_rate_phaseA  — r1*(2*p_exp−1)*(L0−T1) CAGs/year in slow phase
    net_rate_phaseB  — (r1+r2)*(2*p_exp−1)*(L0−T1) CAGs/year in fast phase
    phase_B_fraction — fraction of cells predicted to have crossed T2 at donor age
    instability_index — Var(Δ) at donor age under fitted model

Comparison with previous Regime 2
-----------------------------------
Old: MDN-based amortised inference over 6 sigmoid parameters (a0,a1,b0,b1,c0,c1)
     Time parameterised as discrete NegBin events
     No donor age input
     No Winsorized instability index

New: Exact CTMC inference (matrix exponential)
     Time parameterised in years (donor age as explicit input)
     Two-phase rate function matching Handsaker et al.
     p_exp fixed globally (biologically grounded)
     Winsorized mean absolute delta reported alongside Var(Δ)

Reference
---------
Handsaker RE et al. (2024). Long somatic DNA-repeat expansion drives
neurodegeneration in Huntington's disease. Cell.
Section 4: Repeat expansion dynamics.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.linalg import expm
from scipy import stats


# ===========================================================================
# Global parameter presets — biologically grounded from Handsaker et al.
# ===========================================================================

# HTT CAG repeat (HD) — directly from paper section 4.12
HTT_CAG_PARAMS = dict(
    T1    = 33.5,    # lower instability threshold (repeats < T1 are stable)
    T2    = 72.0,    # phase A→B transition threshold (72.2 ± 2.53 in paper)
    p_exp = 0.676,   # expansion probability (0.676 ± 0.011 across donors)
)

# FMR1 CGG repeat — T1 from known biology, T2/p_exp less constrained
# Premutation (55-200): use with caution, limited published CTMC data
FMR1_CGG_PARAMS = dict(
    T1    = 54.0,    # <55 CGG repeats are generally stable
    T2    = 90.0,    # approximate phase transition (less well-characterised)
    p_exp = 0.70,    # expansion-biased; approximate
)

# Generic / unknown repeat type — conservative defaults
GENERIC_PARAMS = dict(
    T1    = 35.0,
    T2    = 80.0,
    p_exp = 0.65,
)

PRESET_MAP = {
    "htt_cag":   HTT_CAG_PARAMS,
    "fmr1_cgg":  FMR1_CGG_PARAMS,
    "generic":   GENERIC_PARAMS,
}


# ===========================================================================
# Global parameter dataclass
# ===========================================================================

@dataclass
class CTMCGlobalParams:
    """
    Global parameters for the two-phase CTMC somatic expansion model.

    These parameters describe the biology of a specific repeat type and are
    shared across all samples and loci of that type. They should be estimated
    once from a reference dataset (ideally multiple donors spanning a range of
    ages and allele lengths) and then fixed for per-locus inference.

    Parameters
    ----------
    T1 : float
        Lower instability threshold. Repeats shorter than T1 are essentially
        stable (mutation rate = 0). For HTT CAG: 33.5 (Handsaker et al.).

    T2 : float
        Phase transition threshold. Above T2, the mutation rate accelerates
        sharply (Phase B). For HTT CAG: ~72 CAGs (Handsaker et al.).

    p_exp : float
        Expansion probability at each mutation event (constant across lengths
        and donors for a given repeat type). For HTT CAG: 0.676.
        Determines the net expansion bias: net_rate = μ(L) × (2·p_exp − 1).

    r1_prior_shape, r1_prior_scale : float
        Gamma prior parameters for r1 (phase A rate constant).
        Default: shape=1.5, scale=0.02 (weakly informative, centred ~0.03).

    r2_prior_shape, r2_prior_scale : float
        Gamma prior parameters for r2 (phase B rate constant).
        Default: shape=1.5, scale=1.0 (wider prior for faster rate).

    Notes
    -----
    The net expansion rate in each phase at repeat length L is:
        Phase A: r1 × (2·p_exp − 1) × max(L − T1, 0)  [CAGs/year]
        Phase B: (r1+r2) × (2·p_exp − 1) × max(L − T1, 0)  [CAGs/year]

    The ratio (r1+r2)/r1 is the phase B / phase A rate acceleration.
    Handsaker et al. find this ratio ≈ 16 for HTT CAG across donors.
    """
    T1:               float = 33.5
    T2:               float = 72.0
    p_exp:            float = 0.676
    r1_prior_shape:   float = 1.5
    r1_prior_scale:   float = 0.02
    r2_prior_shape:   float = 1.5
    r2_prior_scale:   float = 1.0
    repeat_type:      str   = "generic"
    n_loci_used:      int   = 0

    @classmethod
    def from_preset(cls, preset: str) -> "CTMCGlobalParams":
        """
        Construct from a named preset.

        Parameters
        ----------
        preset : str
            One of: 'htt_cag', 'fmr1_cgg', 'generic'.

        Examples
        --------
        gp = CTMCGlobalParams.from_preset('htt_cag')
        """
        if preset not in PRESET_MAP:
            raise ValueError(
                f"Unknown preset '{preset}'. Choose from: {list(PRESET_MAP)}"
            )
        params = PRESET_MAP[preset]
        return cls(**params, repeat_type=preset)

    def net_rate_phaseA(self, L: float) -> float:
        """Net expansion rate in phase A at repeat length L (CAGs/year)."""
        return float(self.r1_prior_scale * (2 * self.p_exp - 1)
                     * max(L - self.T1, 0))

    def __repr__(self) -> str:
        return (
            f"CTMCGlobalParams(\n"
            f"  repeat_type = '{self.repeat_type}'\n"
            f"  T1={self.T1:.1f}  T2={self.T2:.1f}  p_exp={self.p_exp:.3f}\n"
            f"  Phase A net rate at L=40: "
            f"{self.r1_prior_scale*(2*self.p_exp-1)*max(40-self.T1,0):.4f} CAGs/year "
            f"(at prior mean r1={self.r1_prior_scale:.3f})\n"
            f"  n_loci_used = {self.n_loci_used}\n)"
        )


# ===========================================================================
# Result dataclass
# ===========================================================================

@dataclass
class CTMCResult:
    """
    Posterior summary for the two-phase CTMC model.

    All tuple fields are (posterior_mean, ci_lo_95, ci_hi_95).
    """
    haplotype_label:   str
    n_reads:           int
    founder_length:    float
    donor_age:         float    # years — required input

    # Global params used
    T1:     float
    T2:     float
    p_exp:  float

    # Per-locus posteriors
    r1:              tuple[float, float, float]  # phase A rate constant
    r2:              tuple[float, float, float]  # phase B rate constant

    # Derived expansion rates (CAGs/year)
    net_rate_phaseA: tuple[float, float, float]
    net_rate_phaseB: tuple[float, float, float]
    rate_acceleration: tuple[float, float, float]  # (r1+r2)/r1

    # Distribution summary at donor age
    instability_index:   tuple[float, float, float]  # Var(Δ) under fitted model
    mean_delta:          tuple[float, float, float]  # Mean(Δ) under fitted model
    phase_B_fraction:    tuple[float, float, float]  # P(length > T2) at donor age
    winsorized_mean_abs_delta: float   # mean|Δ| Winsorized at winsor_threshold

    # Metadata
    winsor_threshold:  float
    max_L_used:        int
    elapsed_s:         float

    def __repr__(self) -> str:
        def f(n, v):
            return f"  {n:<28s}= {v[0]:.4f}  95%CI [{v[1]:.4f}, {v[2]:.4f}]"
        lines = [
            f"CTMCResult(haplotype='{self.haplotype_label}', "
            f"n={self.n_reads}, L0={self.founder_length:.1f}, age={self.donor_age:.0f}yr)",
            f"  T1={self.T1:.1f}  T2={self.T2:.1f}  p_exp={self.p_exp:.3f}",
            "",
            "  --- Rate parameters ---",
            f("r1 (phase A rate)",          self.r1),
            f("r2 (phase B rate)",          self.r2),
            f("rate_acceleration (B/A)",    self.rate_acceleration),
            "",
            "  --- Net expansion rates (CAGs/year at L0) ---",
            f("net_rate_phaseA",            self.net_rate_phaseA),
            f"    [r1·(2p−1)·(L0−T1)]",
            f("net_rate_phaseB",            self.net_rate_phaseB),
            f"    [(r1+r2)·(2p−1)·(L0−T1)]",
            "",
            "  --- Distribution at donor age ---",
            f("instability_index",          self.instability_index),
            f"    [Var(Δ) under fitted PMF]",
            f("mean_delta",                 self.mean_delta),
            f"    [Mean(Δ) under fitted PMF]",
            f("phase_B_fraction",           self.phase_B_fraction),
            f"    [P(length > T2={self.T2:.0f}) at age {self.donor_age:.0f}]",
            f"  winsorized_mean_abs_delta = "
            f"{self.winsorized_mean_abs_delta:.4f}"
            f"  [mean|Δ| Winsorized at {self.winsor_threshold:.0f} units]",
            f"\n  elapsed = {self.elapsed_s*1000:.1f} ms",
        ]
        return "\n".join(lines)

    def to_dict(self) -> dict:
        d = {
            "regime":          "2_ctmc",
            "haplotype":       self.haplotype_label,
            "n_reads":         self.n_reads,
            "founder_length":  self.founder_length,
            "donor_age":       self.donor_age,
            "T1":              self.T1,
            "T2":              self.T2,
            "p_exp":           self.p_exp,
            "winsor_threshold": self.winsor_threshold,
            "winsorized_mean_abs_delta": self.winsorized_mean_abs_delta,
            "max_L_used":      self.max_L_used,
        }
        for field_name in [
            "r1", "r2", "rate_acceleration",
            "net_rate_phaseA", "net_rate_phaseB",
            "instability_index", "mean_delta", "phase_B_fraction",
        ]:
            v = getattr(self, field_name)
            d[field_name]            = v[0]
            d[field_name + "_ci_lo"] = v[1]
            d[field_name + "_ci_hi"] = v[2]
        d["elapsed_s"] = self.elapsed_s
        return d


# ===========================================================================
# Core CTMC computation
# ===========================================================================

def _build_rate_matrix(
    r1: float,
    r2: float,
    T1: float,
    T2: float,
    p_exp: float,
    min_L: int,
    max_L: int,
) -> np.ndarray:
    """
    Build the CTMC rate matrix Q for the TwoPhaseLinear model.

    The state space is [min_L, ..., max_L]. max_L is an absorbing sink.

    Q[x, x+1] = p_exp   × μ(x)    expansion
    Q[x, x-1] = (1-p_exp) × μ(x)  contraction
    Q[x, x]   = −μ(x)              diagonal

    where μ(x) = r1·max(x−T1, 0) + r2·max(x−T2, 0).
    """
    n  = max_L - min_L + 1
    x  = np.arange(min_L, max_L + 1, dtype=float)
    mu = r1 * np.maximum(x - T1, 0.0) + r2 * np.maximum(x - T2, 0.0)

    # Upper diagonal: expansion Q[i, i+1] = p_exp * mu[i]
    # Lower diagonal: contraction Q[i, i-1] = (1-p_exp) * mu[i]
    Q = np.diag(p_exp * mu[:-1], 1) + np.diag((1 - p_exp) * mu[1:], -1)
    np.fill_diagonal(Q, -mu)  # row sums = 0 (before sink adjustment)

    # max_L is absorbing: all outgoing rates removed
    Q[-1, :] = 0.0
    return Q


def _ctmc_pmf(
    L0: int,
    age: float,
    r1: float,
    r2: float,
    T1: float,
    T2: float,
    p_exp: float,
    min_L: int,
    max_L: int,
) -> np.ndarray:
    """
    Compute the PMF of repeat lengths at donor age for a single (r1, r2) pair.

    Returns array of shape (max_L - min_L + 1,) giving P(length = k) for
    k in [min_L, ..., max_L].
    """
    Q  = _build_rate_matrix(r1, r2, T1, T2, p_exp, min_L, max_L)
    Pt = expm(age * Q)
    L0_idx = int(L0) - min_L
    pmf    = np.clip(Pt[L0_idx, :], 0.0, None)
    pmf   /= pmf.sum()
    return pmf


def _grid_loglik(
    observed_idx: np.ndarray,
    counts:        np.ndarray,
    L0:            int,
    age:           float,
    r1_grid:       np.ndarray,
    r2_grid:       np.ndarray,
    T1:            float,
    T2:            float,
    p_exp:         float,
    min_L:         int,
    max_L:         int,
) -> np.ndarray:
    """
    Log-likelihood on a 2D grid of (r1, r2) values.

    Returns array of shape (n_r1, n_r2).
    """
    n_r1, n_r2 = len(r1_grid), len(r2_grid)
    loglik     = np.full((n_r1, n_r2), -np.inf)
    L0_idx     = int(L0) - min_L

    for i, r1v in enumerate(r1_grid):
        for j, r2v in enumerate(r2_grid):
            try:
                Q  = _build_rate_matrix(r1v, r2v, T1, T2, p_exp, min_L, max_L)
                Pt = expm(age * Q)
                pmf = np.clip(Pt[L0_idx, :], 1e-300, None)
                pmf /= pmf.sum()
                # Log-likelihood: sum of count * log(pmf) at observed lengths
                ll = float(sum(
                    counts[k] * np.log(pmf[d])
                    for k, d in enumerate(observed_idx)
                    if d < len(pmf)
                ))
                loglik[i, j] = ll
            except Exception:
                pass  # expm can fail for extreme parameter values

    return loglik


# ===========================================================================
# Public API: fit one locus
# ===========================================================================

def fit_locus_ctmc(
    lengths,
    founder_length:    float,
    donor_age:         float,
    global_params:     CTMCGlobalParams,
    haplotype_label:   str   = "haplotype",
    r1_n:              int   = 12,
    r2_n:              int   = 12,
    winsor_threshold:  float = 100.0,
    seed:              int   = 0,
) -> CTMCResult:
    """
    Fit the two-phase CTMC somatic expansion model for one haplotype.

    This is the Regime 2 replacement, following Handsaker et al. (2024).
    It requires donor age as an explicit input — the CTMC evolves repeat
    lengths in continuous time (years), making age a fundamental parameter
    rather than an optional annotation.

    Parameters
    ----------
    lengths : array-like
        Per-read allele lengths in repeat units (already divided by motif size).
        Float values are rounded to the nearest integer.
    founder_length : float
        Germline/inherited repeat length in repeat units (L0).
    donor_age : float
        Donor age in years at the time of tissue sampling. For post-mortem
        brain samples this is age at death. This is required — the CTMC
        cannot fit without it.
    global_params : CTMCGlobalParams
        Global parameters (T1, T2, p_exp, priors) for this repeat type.
        Use CTMCGlobalParams.from_preset('htt_cag') for HTT CAG repeats.
    haplotype_label : str
        Label for output.
    r1_n, r2_n : int
        Grid resolution for r1 and r2. Default 12×12 = 144 grid points.
        Increase for higher precision (costs proportionally more time).
    winsor_threshold : float
        Winsorization threshold for mean absolute delta (default 100 repeat
        units, matching Handsaker et al. Figure SN2.3). Set to np.inf to
        disable Winsorization.
    seed : int
        RNG seed (not used in current implementation, reserved for future
        Monte Carlo uncertainty propagation).

    Returns
    -------
    CTMCResult
        All quantities with 95% posterior CIs from the 2D grid posterior.

    Notes
    -----
    The key outputs to report are:

    net_rate_phaseA : expansion rate in the slow phase (CAGs/year at L0).
        This is the primary per-donor instability measure for most HD donors
        where the majority of cells have not yet crossed T2.

    net_rate_phaseB : expansion rate in the fast phase (CAGs/year at L0).
        Highly uncertain unless many cells have been observed past T2.

    rate_acceleration : ratio of phase B to phase A rate (~16× for HTT).
        A well-constrained value here indicates cells spanning both phases.

    phase_B_fraction : fraction of cells predicted past T2 at donor age.
        Near 0 = mostly phase A (early disease / young donor).
        Near 1 = mostly phase B (advanced disease / older donor).

    winsorized_mean_abs_delta : directly comparable to somatic instability
        index values reported by Handsaker et al. (Fig. SN2.3).
    """
    t0 = time.time()

    lengths     = np.asarray(lengths, dtype=float)
    lengths_int = np.round(lengths).astype(int)
    L0          = int(round(founder_length))

    T1, T2, p_exp = global_params.T1, global_params.T2, global_params.p_exp

    # ── Adaptive state space ─────────────────────────────────────────────────
    # Set max_L to just cover observed data + a buffer proportional to the
    # observed spread. This keeps matrix size small for typical loci.
    obs_max  = int(lengths_int.max())
    obs_min  = int(lengths_int.min())
    obs_span = max(obs_max - L0, L0 - obs_min, 1)
    max_L    = min(obs_max + max(30, obs_span // 2), 1000)
    min_L    = max(obs_min - 10, 0)

    # Observed lengths as indices into [min_L, ..., max_L]
    obs_idx = lengths_int - min_L
    valid   = (obs_idx >= 0) & (obs_idx <= max_L - min_L)
    obs_idx = obs_idx[valid]
    unique_obs, counts = np.unique(obs_idx, return_counts=True)

    # ── Parameter grids ──────────────────────────────────────────────────────
    # r1: phase A rate constant — typically small (0.001 to 0.1)
    # r2: phase B rate constant — can be much larger (0.1 to 10)
    # Log-spaced grids give fine resolution at small values where most
    # loci live, while still covering high-rate loci.
    # r1 uses more grid points because the posterior is sharp (well-identified
    # even when no cells have crossed T2). Coarse r1 grids collapse CIs.
    r1_grid = np.exp(np.linspace(np.log(1e-4), np.log(0.30), max(r1_n, 30)))
    r2_grid = np.exp(np.linspace(np.log(1e-3), np.log(15.0), r2_n))

    # ── Grid log-likelihood ──────────────────────────────────────────────────
    loglik = _grid_loglik(
        unique_obs, counts, L0, donor_age,
        r1_grid, r2_grid, T1, T2, p_exp, min_L, max_L,
    )

    # ── Posterior ────────────────────────────────────────────────────────────
    log_prior = (
        stats.gamma.logpdf(
            r1_grid,
            a=global_params.r1_prior_shape,
            scale=global_params.r1_prior_scale,
        )[:, None]
        + stats.gamma.logpdf(
            r2_grid,
            a=global_params.r2_prior_shape,
            scale=global_params.r2_prior_scale,
        )[None, :]
    )
    lp   = loglik + log_prior
    lp  -= lp.max()
    post = np.exp(lp)
    post /= post.sum()

    R1, R2 = np.meshgrid(r1_grid, r2_grid, indexing="ij")

    # Marginal posteriors — sum over the other axis
    # r1 posterior is well-identified (sharp) even without phase B cells.
    # r2 posterior may be flat when no cells have crossed T2.
    # CIs from the joint (r1,r2) grid collapse when one axis is flat
    # because all mass concentrates on one grid point.
    # Using marginal posteriors gives correct CIs for each parameter.
    post_r1 = post.sum(axis=1)   # (n_r1,)  marginal over r2
    post_r2 = post.sum(axis=0)   # (n_r2,)  marginal over r1
    post_r1 /= post_r1.sum()
    post_r2 /= post_r2.sum()

    def ci_1d(vals: np.ndarray, post_1d: np.ndarray) -> tuple[float, float, float]:
        """CI from a 1D marginal posterior."""
        mean_ = float((post_1d * vals).sum())
        si    = np.argsort(vals)
        cdf   = np.cumsum(post_1d[si])
        lo95  = float(vals[si[np.searchsorted(cdf, 0.025)]])
        hi95  = float(vals[si[np.searchsorted(cdf, 0.975)]])
        return mean_, lo95, hi95

    def ci_joint(g: np.ndarray) -> tuple[float, float, float]:
        """CI from the joint 2D posterior (for derived quantities)."""
        mean_ = float((post * g).sum())
        gf    = g.ravel()
        pf    = post.ravel()
        si    = np.argsort(gf)
        cdf   = np.cumsum(pf[si])
        lo95  = float(gf[si[np.searchsorted(cdf, 0.025)]])
        hi95  = float(gf[si[np.searchsorted(cdf, 0.975)]])
        return mean_, lo95, hi95

    # ── Derived quantities ───────────────────────────────────────────────────
    # Net expansion rates at L0 (CAGs/year)
    phase_A_factor = float(max(L0 - T1, 0.0) * (2 * p_exp - 1))
    net_rateA_1d = r1_grid * phase_A_factor           # (n_r1,)
    net_rateB_1d = r2_grid * phase_A_factor            # (n_r2,) — only r2 component
    accel        = (R1 + R2) / np.maximum(R1, 1e-10)  # (n_r1, n_r2)

    # Instability index and mean delta from fitted PMF at posterior mean
    r1_est = float((post * R1).sum())
    r2_est = float((post * R2).sum())
    pmf_est = _ctmc_pmf(L0, donor_age, r1_est, r2_est, T1, T2, p_exp, min_L, max_L)
    lengths_grid = np.arange(min_L, max_L + 1, dtype=float)
    deltas_grid  = lengths_grid - L0
    mean_d  = float((deltas_grid * pmf_est).sum())
    var_d   = float(((deltas_grid - mean_d) ** 2 * pmf_est).sum())
    phase_b_frac = float(pmf_est[int(T2) - min_L:].sum()) if int(T2) >= min_L else 0.0

    # Winsorized mean absolute delta (matches Handsaker et al. definition)
    # Computed directly from observed reads, not from the fitted model
    abs_deltas    = np.abs(lengths_int - L0).astype(float)
    abs_deltas_w  = np.minimum(abs_deltas, winsor_threshold)
    winsor_mad    = float(abs_deltas_w.mean())

    return CTMCResult(
        haplotype_label          = haplotype_label,
        n_reads                  = len(lengths_int),
        founder_length           = float(founder_length),
        donor_age                = float(donor_age),
        T1                       = T1,
        T2                       = T2,
        p_exp                    = p_exp,
        r1                       = ci_1d(r1_grid, post_r1),
        r2                       = ci_1d(r2_grid, post_r2),
        net_rate_phaseA          = ci_1d(net_rateA_1d, post_r1),
        net_rate_phaseB          = ci_joint((R1 + R2) * phase_A_factor),
        rate_acceleration        = ci_joint(accel),
        instability_index        = (var_d, float("nan"), float("nan")),
        mean_delta               = (mean_d, float("nan"), float("nan")),
        phase_B_fraction         = (phase_b_frac, float("nan"), float("nan")),
        winsorized_mean_abs_delta = winsor_mad,
        winsor_threshold         = winsor_threshold,
        max_L_used               = max_L,
        elapsed_s                = time.time() - t0,
    )


# ===========================================================================
# Pass 1: estimate T2 and p_exp from multiple donors/loci
# ===========================================================================

def estimate_ctmc_global(
    loci_data: list[dict],
    T1:        float = 33.5,
    T2_grid:   Optional[np.ndarray] = None,
    pexp_grid: Optional[np.ndarray] = None,
    r1_n:      int   = 10,
    r2_n:      int   = 10,
) -> CTMCGlobalParams:
    """
    Estimate T2 and p_exp from multiple donors/loci by maximising the
    pooled log-likelihood across the dataset.

    Uses a grid search over (T2, p_exp) — these are typically stable
    across donors for a given repeat type, so pooling maximises power.

    Parameters
    ----------
    loci_data : list of dicts, each with:
        'lengths'        : array of observed lengths (repeat units)
        'founder_length' : float (repeat units)
        'donor_age'      : float (years)
    T1 : float
        Fixed lower threshold (not estimated — determined by biology).
    T2_grid : array, optional
        Values of T2 to search. Default: np.linspace(50, 120, 15).
    pexp_grid : array, optional
        Values of p_exp to search. Default: np.linspace(0.55, 0.85, 7).
    r1_n, r2_n : int
        Grid resolution for per-locus rate fitting at each (T2, pexp) point.

    Returns
    -------
    CTMCGlobalParams with estimated T2, p_exp and priors derived from
    the distribution of per-locus r1, r2 estimates.
    """
    if T2_grid is None:
        T2_grid   = np.linspace(50.0, 120.0, 15)
    if pexp_grid is None:
        pexp_grid = np.linspace(0.55, 0.85, 7)

    # Filter to loci with sufficient reads and donor age
    valid = []
    for locus in loci_data:
        if (len(locus.get("lengths", [])) >= 5 and
                locus.get("donor_age") is not None and
                locus.get("founder_length") is not None):
            valid.append(locus)

    if len(valid) < 2:
        print("Warning: fewer than 2 valid loci for Pass 1. Using defaults.")
        return CTMCGlobalParams(T1=T1, repeat_type="estimated")

    print(f"Pass 1 (CTMC): estimating T2, p_exp from {len(valid)} loci...")

    best_ll  = -np.inf
    best_T2  = float(T2_grid[len(T2_grid)//2])
    best_pexp = 0.676

    r1_grid = np.exp(np.linspace(np.log(1e-4), np.log(0.30), r1_n))
    r2_grid = np.exp(np.linspace(np.log(1e-3), np.log(15.0), r2_n))

    for T2 in T2_grid:
        for pexp in pexp_grid:
            total_ll = 0.0
            for locus in valid:
                lengths_int = np.round(np.asarray(locus["lengths"])).astype(int)
                L0          = int(round(locus["founder_length"]))
                age         = float(locus["donor_age"])
                obs_max  = int(lengths_int.max())
                obs_min  = int(lengths_int.min())
                obs_span = max(obs_max - L0, L0 - obs_min, 1)
                max_L = min(obs_max + max(30, obs_span//2), 500)
                min_L = max(obs_min - 10, 0)
                obs_idx = lengths_int - min_L
                valid_mask = (obs_idx >= 0) & (obs_idx <= max_L - min_L)
                obs_idx = obs_idx[valid_mask]
                unique_obs, counts = np.unique(obs_idx, return_counts=True)

                ll = _grid_loglik(
                    unique_obs, counts, L0, age,
                    r1_grid, r2_grid, T1, T2, pexp, min_L, max_L,
                )
                # Marginalise over r1, r2 to get locus likelihood
                ll_max = ll.max()
                if np.isfinite(ll_max):
                    total_ll += ll_max + np.log(np.exp(ll - ll_max).sum())

            if total_ll > best_ll:
                best_ll   = total_ll
                best_T2   = float(T2)
                best_pexp = float(pexp)

    print(f"  Best: T2={best_T2:.1f}  p_exp={best_pexp:.3f}  loglik={best_ll:.2f}")

    return CTMCGlobalParams(
        T1           = T1,
        T2           = best_T2,
        p_exp        = best_pexp,
        repeat_type  = "estimated",
        n_loci_used  = len(valid),
    )


# ===========================================================================
# Self-test
# ===========================================================================

def _run_selftest():
    import time
    rng = np.random.default_rng(42)

    print("=" * 65)
    print("CTMC REGIME 2 SELF-TEST — TwoPhaseLinear model")
    print("=" * 65)

    # Parameters from Handsaker et al.
    true_r1, true_r2 = 0.012, 1.5
    T1, T2, p_exp    = 33.5, 72.0, 0.676
    L0, age          = 42, 65

    gp = CTMCGlobalParams.from_preset("htt_cag")
    print(gp)
    print()

    # Simulate data from the true model
    max_L = 200
    Q_true = _build_rate_matrix(true_r1, true_r2, T1, T2, p_exp, 0, max_L)
    Pt     = expm(age * Q_true)
    pmf_true = Pt[L0, :]; pmf_true /= pmf_true.sum()

    lengths_sim = rng.choice(np.arange(max_L + 1), size=50, p=pmf_true).astype(float)

    print(f"Simulated reads (n=50, L0={L0}, age={age}):")
    print(f"  mean={lengths_sim.mean():.1f}  max={lengths_sim.max():.0f}"
          f"  P(>T2)={(lengths_sim > T2).mean():.2f}")
    print()

    # Fit
    t0 = time.time()
    result = fit_locus_ctmc(lengths_sim, L0, age, gp, haplotype_label="H1")
    print(result)

    print()
    print(f"True r1={true_r1:.4f}  estimated={result.r1[0]:.4f}")
    print(f"True r2={true_r2:.3f}   estimated={result.r2[0]:.3f}")
    print(f"True net rate phase A: "
          f"{true_r1*(2*p_exp-1)*max(L0-T1,0):.4f} CAGs/yr  "
          f"estimated: {result.net_rate_phaseA[0]:.4f}")
    print(f"Winsorized mean|Δ|: {result.winsorized_mean_abs_delta:.3f} "
          f"(comparable to Handsaker et al. Fig. SN2.3)")

    # Coverage test
    print()
    print("--- Net rate phase A: CI coverage (30 trials, n=50) ---")
    true_net_A = true_r1 * (2*p_exp-1) * max(L0-T1, 0)
    covered = 0
    times   = []
    for _ in range(30):
        sim = rng.choice(np.arange(max_L+1), size=50, p=pmf_true).astype(float)
        res = fit_locus_ctmc(sim, L0, age, gp)
        times.append(res.elapsed_s)
        lo, hi = res.net_rate_phaseA[1], res.net_rate_phaseA[2]
        if lo <= true_net_A <= hi:
            covered += 1
    print(f"Coverage: {covered}/30 = {covered/30*100:.0f}%")
    print(f"Mean time: {np.mean(times)*1000:.0f}ms")

    # Test with older donor where phase B cells are visible
    print()
    print("--- Older donor (age=80) where Phase B is populated ---")
    age2 = 80
    Pt2      = expm(age2 * Q_true)
    pmf2     = Pt2[L0, :]; pmf2 /= pmf2.sum()
    lengths2 = rng.choice(np.arange(max_L+1), size=50, p=pmf2).astype(float)
    print(f"Simulated reads (age={age2}): mean={lengths2.mean():.0f}"
          f"  max={lengths2.max():.0f}  P(>T2)={(lengths2>T2).mean():.2f}")
    result2 = fit_locus_ctmc(lengths2, L0, age2, gp, haplotype_label="H1_old")
    print(f"r1: {result2.r1[0]:.4f} [{result2.r1[1]:.4f},{result2.r1[2]:.4f}]  "
          f"(true={true_r1:.4f})")
    print(f"r2: {result2.r2[0]:.3f} [{result2.r2[1]:.3f},{result2.r2[2]:.3f}]   "
          f"(true={true_r2:.3f})")
    print(f"Rate acceleration: {result2.rate_acceleration[0]:.1f}×  "
          f"(paper reports ~16× for HTT)")
    print(f"Phase B fraction: {result2.phase_B_fraction[0]:.3f}")


if __name__ == "__main__":
    _run_selftest()