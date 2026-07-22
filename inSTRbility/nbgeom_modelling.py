"""
Somatic Instability Index — NB-Geometric Length-Dependent Model
===============================================================

Unified model covering two biological regimes of tandem repeat instability,
both using Negative-Binomial mutation counts and Geometric step sizes
(biologically motivated: slippage probability decays exponentially with step
size, giving a geometric step-size distribution; NB count captures
cell-to-cell heterogeneity in mutation rate).

Both regimes are run in parallel for comparison when compare_models() is called.

─────────────────────────────────────────────────────────────────────────────
REGIME 1 — Small-to-moderate somatic excursions (|Δ| ≪ L₀)
─────────────────────────────────────────────────────────────────────────────
Valid for most WGS loci where typical somatic length changes are small
relative to the founder allele length. Step-size distribution is evaluated
at L₀ and held constant (linearisation approximation).

Model:
  E = total expansion ~ NegBin(r·α,  p₊)   [NB(r_e, p₊) of Geom(p₊) steps]
  C = total contraction ~ NegBin(r·(1−α), p₋)
  Δ = E − C  ~  Generalised Skellam(r_e, p₊, r_c, p₋)

Parameters:
  Global (Pass 1):  p₊(L₀), p₋(L₀) — step-size params as functions of L₀
                    q_e, q_c         — NegBin probability for count distribution
  Per-locus Pass 2: r — total event dispersion
                    α — expansion fraction ∈ (0,1)

PMF: computed via FFT of the characteristic function.

Instability index (closed form):
  Var(Δ) = μ_{N_e}·(1−p₊)/p₊² + σ²_{N_e}/p₊²
          + μ_{N_c}·(1−p₋)/p₋² + σ²_{N_c}/p₋²

Python API
──────────
  from nb_geom_instability import (
      fit_locus, compare_models,
      Regime1GlobalParams, stimate_regime1_global,
  )

  # Fit both regimes and compare
  r1 = compare_models(lengths, founder_length=80, n_steps=5, regime1_params=gp1)
  print(r1)

Dependencies: numpy, scipy, torch (for Regime 2 MDN)
              pymc (LONG_READ_WGS_LARGE only, not implemented in this version)
"""

from __future__ import annotations

import json
import time
import warnings
import logging
from dataclasses import dataclass, asdict, field
from typing import Optional

import math
import numpy as np
from scipy import stats
from scipy.stats import skew as _skew, kurtosis as _kurt


logging.getLogger("pymc").setLevel(logging.ERROR)


# ===========================================================================
# Constants and grid defaults
# ===========================================================================

# Regime 1 per-locus grid: (r_e, r_c) with narrow (p₊, p₋) grid around global
# Log-spaced r grids give fine resolution at small values (most genome-wide
# loci have few mutations per lineage, r ~ 0.01-0.5) while still covering
# high-instability loci (r up to 8). Linear spacing would waste most grid
# points in the high-r region where few loci actually live.
# 25 × 25 × 7 × 7 = 30,625 grid points, ~150ms per locus
DEFAULT_RE_GRID  = np.exp(np.linspace(np.log(0.01), np.log(8.0), 25))  # log-spaced
DEFAULT_RC_GRID  = np.exp(np.linspace(np.log(0.01), np.log(6.0), 25))  # log-spaced
DEFAULT_PP_N     = 7     # number of p₊ grid points around global estimate
DEFAULT_PM_N     = 7     # number of p₋ grid points around global estimate
DEFAULT_P_LOG_SD = 0.2   # prior SD in log-odds space for p₊ and p₋ grids


# ===========================================================================
# Global parameter dataclasses
# ===========================================================================

@dataclass
class Regime1GlobalParams:
    """
    Global parameters for Regime 1 (NB-Geometric compound).

    Step-size parameters estimated from pooled high-depth reads across loci.
    These are locus-biology parameters (how repeat slippage mechanics work)
    not sample-level parameters, so they are shared across samples.

    p₊(L₀): geometric parameter for expansion steps.
            Smaller p₊ → larger mean expansion step (mean = 1/p₊).
            Can vary with L₀ via log-odds linear model:
            logit(p₊(L₀)) = lp0 + lp1·L₀

    p₋(L₀): geometric parameter for contraction steps.
            Similarly parameterised.

    q_e, q_c: NegBin probability for the expansion and contraction
              count distributions. Estimated from the shape of the
              count distribution across cells.
    """
    # Log-odds linear model for p₊: logit(p₊) = lp0 + lp1·L₀
    lp0:  float = 0.0    # logit(p₊) intercept
    lp1:  float = 0.0    # logit(p₊) slope (negative → bigger steps at longer alleles)
    # Log-odds linear model for p₋
    lm0:  float = 0.0
    lm1:  float = 0.0
    # NegBin probability for count distributions
    q_e:  float = 0.5    # expansion count NB probability
    q_c:  float = 0.5    # contraction count NB probability
    # Diagnostics
    n_loci_used: int = 0

    def p_plus(self, L0: float) -> float:
        """Expansion step-size geometric parameter at founder length L0."""
        return float(_sigmoid_scalar(self.lp0 + self.lp1 * L0))

    def p_minus(self, L0: float) -> float:
        """Contraction step-size geometric parameter at founder length L0."""
        return float(_sigmoid_scalar(self.lm0 + self.lm1 * L0))

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)

    @classmethod
    def from_json(cls, s: str) -> "Regime1GlobalParams":
        return cls(**json.loads(s))

    @classmethod
    def from_defaults(
        cls,
        repeat_unit_length: int = 3,
        q_e: float = 0.5,
        q_c: float = 0.5,
    ) -> "Regime1GlobalParams":
        """
        Biologically motivated defaults WITHOUT running Pass 1.

        Use when you want to run per-locus inference before estimating global
        params, or as a quick sanity-check starting point. The narrow p₊/p₋
        floating grid in Pass 2 corrects for any bias here, so mild
        initialisation errors do not badly affect per-locus posteriors.

        repeat_unit_length : int
            Length of the repeat unit in base pairs.
              1-3 bp (STR):           p₊ ≈ 0.60, p₋ ≈ 0.70  (mean steps ~1.7 / ~1.4)
              4-6 bp (medium STR):    p₊ ≈ 0.45, p₋ ≈ 0.55  (mean steps ~2.2 / ~1.8)
              7+ bp (long/disease):   p₊ ≈ 0.30, p₋ ≈ 0.45  (mean steps ~3.3 / ~2.2)
        q_e, q_c : float
            NegBin probability for count distributions (default 0.5).
        """
        def logit(p): return math.log(p / (1.0 - p))
        if repeat_unit_length <= 3:
            pp, pm = 0.60, 0.70
        elif repeat_unit_length <= 6:
            pp, pm = 0.45, 0.55
        else:
            pp, pm = 0.30, 0.45
        return cls(lp0=logit(pp), lp1=0.0, lm0=logit(pm), lm1=0.0,
                   q_e=q_e, q_c=q_c, n_loci_used=0)

    @classmethod
    def from_stratum(
        cls,
        stratum: str,
        q_e: float = 0.5,
        q_c: float = 0.5,
    ) -> "Regime1GlobalParams":
        """
        Named-stratum defaults for genome-wide catalogs.

        stratum options
        ---------------
        'str_short'   : di/trinucleotide STRs (AC, AT, CAG etc.)
        'str_medium'  : tetra/pentanucleotide STRs
        'str_long'    : hexanucleotide and longer STRs
        'disease_cag' : CAG repeats (HD, SCA)
        'disease_cgg' : CGG repeats (FMR1 premutation; use Regime 2 for full mutation)
        'disease_gaa' : GAA repeats (FRDA)
        'disease_ctg' : CTG repeats (DM1)
        """
        STRATA = {
            "str_short":   (0.60, 0.70),
            "str_medium":  (0.45, 0.55),
            "str_long":    (0.30, 0.45),
            "disease_cag": (0.35, 0.50),
            "disease_cgg": (0.30, 0.45),
            "disease_gaa": (0.25, 0.40),
            "disease_ctg": (0.35, 0.50),
        }
        if stratum not in STRATA:
            raise ValueError(
                f"Unknown stratum '{stratum}'. Choose from: {list(STRATA)}"
            )
        def logit(p): return math.log(p / (1.0 - p))
        pp, pm = STRATA[stratum]
        return cls(lp0=logit(pp), lp1=0.0, lm0=logit(pm), lm1=0.0,
                   q_e=q_e, q_c=q_c, n_loci_used=0)

    @classmethod
    def calibrate_from_results(
        cls,
        results_tsv: str,
        base_params: "Regime1GlobalParams",
        min_reads: int = 20,
        target_percentile: float = 50,
    ) -> "Regime1GlobalParams":
        """
        Adjust q_e and q_c from observed dispersion ratios in a results TSV.

        When the model is systematically overdispersed (median dispersion_ratio
        < 1, median PPP_var > 0.9), the NegBin count distribution is predicting
        too many mutation events per lineage. This method back-calculates the q
        value that aligns the model's expected variance with the typical observed
        variance in your data.

        The adjustment is derived from:
            dispersion_ratio = observed_var / fitted_var
                             ≈ (actual mutation rate) / (assumed mutation rate)

        If dispersion_ratio = 0.19 (your case), the model is predicting ~5×
        too many events. Increasing q from 0.5 to ~0.83 reduces the expected
        count by that factor, realigning the model with your data.

        Parameters
        ----------
        results_tsv : str
            Path to TSV from write_tsv() or analyse_genome_wide().
        base_params : Regime1GlobalParams
            Starting parameters. p₊ and p₋ are kept; only q is adjusted.
        min_reads : int
            Only use loci with at least this many reads (default 20).
            Low-depth loci have unstable dispersion ratio estimates.
        target_percentile : float
            Percentile of dispersion ratio to target (default 50 = median).
            Using the median is robust to the extreme outliers at the tails
            (very stable loci and very unstable loci).

        Returns
        -------
        Regime1GlobalParams with adjusted q_e and q_c.

        Usage
        -----
        # Run once on your initial results, then rerun analysis with new params
        gp_initial = Regime1GlobalParams.from_defaults(repeat_unit_length=3)
        results    = analyse_genome_wide(catalog, gp_initial, ...)
        write_tsv(results, 'initial_results.tsv')

        gp_calibrated = Regime1GlobalParams.calibrate_from_results(
            'initial_results.tsv', gp_initial
        )
        print(gp_calibrated)   # check new q value

        results2 = analyse_genome_wide(catalog, gp_calibrated, ...)
        write_tsv(results2, 'calibrated_results.tsv')
        """
        import copy
        try:
            import pandas as pd
            df = pd.read_csv(results_tsv, sep="\t")
        except ImportError:
            import csv
            rows = []
            with open(results_tsv) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    rows.append(row)
            df_data = {k: [] for k in rows[0]} if rows else {}
            for row in rows:
                for k, v in row.items():
                    try: df_data[k].append(float(v))
                    except (ValueError, TypeError): df_data[k].append(v)
            class _DF:
                def __init__(self, d): self._d = d
                def __getitem__(self, k): return np.array(self._d[k])
                def __len__(self): return len(next(iter(self._d.values())))
            df = _DF(df_data)

        try:
            n_reads = np.array(df["n_reads"], dtype=float)
            disp    = np.array(df["gof_dispersion_ratio"], dtype=float)
        except (KeyError, TypeError):
            raise ValueError(
                "results_tsv must contain 'n_reads' and 'gof_dispersion_ratio' columns. "
                "Ensure write_tsv() was used to generate the file."
            )

        # Filter to reliable loci
        valid = (n_reads >= min_reads) & np.isfinite(disp) & (disp > 0)
        if valid.sum() < 10:
            raise ValueError(
                f"Only {valid.sum()} loci pass filters (n_reads >= {min_reads}, "
                f"valid dispersion_ratio). Need at least 10 for calibration."
            )

        target_ratio = float(np.percentile(disp[valid], target_percentile))
        current_q    = base_params.q_e

        # Back-calculate q from dispersion ratio.
        # At fixed r and p, Var(Delta) ∝ (1-q)/q  (the NB overdispersion).
        # If observed_var = target_ratio × fitted_var, we need to scale
        # (1-q)/q by target_ratio to match:
        #   new_rate = current_rate × target_ratio
        #   new_q    = 1 / (1 + new_rate)
        current_rate = (1.0 - current_q) / current_q
        new_rate     = current_rate * target_ratio
        new_q        = float(np.clip(1.0 / (1.0 + new_rate), 0.05, 0.99))

        n_overdispersed   = int((disp[valid] > 1.5).sum())
        n_underdispersed  = int((disp[valid] < 0.67).sum())
        n_good            = int(valid.sum()) - n_overdispersed - n_underdispersed

        print(f"Calibration summary ({int(valid.sum())} loci with n >= {min_reads}):")
        print(f"  Dispersion ratio at {target_percentile}th percentile: {target_ratio:.4f}")
        print(f"  {'Model overdispersed' if target_ratio < 1 else 'Model underdispersed'} "
              f"({'fitted_var > observed_var' if target_ratio < 1 else 'fitted_var < observed_var'})")
        print(f"  Loci breakdown: {n_good} well-fit, "
              f"{n_overdispersed} underdispersed, {n_underdispersed} overdispersed")
        print(f"  q adjusted: {current_q:.3f} → {new_q:.3f}")
        print(f"  Effect: mean mutation count per lineage "
              f"{'reduced' if new_q > current_q else 'increased'} "
              f"by factor {target_ratio:.2f}")
        print(f"  Expected PPP_var after calibration: median should shift toward 0.5")

        new_params     = copy.deepcopy(base_params)
        new_params.q_e = new_q
        new_params.q_c = new_q
        return new_params

    def __repr__(self) -> str:
        return (
            f"Regime1GlobalParams(\n"
            f"  p₊(L0) = sigmoid({self.lp0:.4f} + {self.lp1:.5f}·L0)\n"
            f"  p₋(L0) = sigmoid({self.lm0:.4f} + {self.lm1:.5f}·L0)\n"
            f"  q_e={self.q_e:.3f}, q_c={self.q_c:.3f}\n"
            f"  n_loci={self.n_loci_used}\n)"
        )


# ===========================================================================
# Result dataclasses
# ===========================================================================

@dataclass
class Regime1Result:
    """
    Posterior summary for Regime 1 (NB-Geometric compound) inference.
    All quantities as (posterior_mean, ci_lo_95, ci_hi_95).
    """
    haplotype_label: str
    n_reads:         int
    founder_length:  float
    p_plus_prior:    float   # global p₊ at this L₀
    p_minus_prior:   float   # global p₋ at this L₀
    q_e:             float
    q_c:             float

    # Per-locus posteriors
    r_e:               tuple[float, float, float]  # expansion NB dispersion
    r_c:               tuple[float, float, float]  # contraction NB dispersion
    p_plus:            tuple[float, float, float]  # expansion step-size param
    p_minus:           tuple[float, float, float]  # contraction step-size param
    mu_expansion:      tuple[float, float, float]  # mean total expansion = r_e*(1-p₊)/p₊
    mu_contraction:    tuple[float, float, float]  # mean total contraction
    instability_index: tuple[float, float, float]  # Var(Δ) — closed form
    net_bias:          tuple[float, float, float]  # Mean(Δ) = μ_E - μ_C
    p_expansion:       tuple[float, float, float]  # P(net_bias > 0)

    # ── Mutation rates ────────────────────────────────────────────────────────
    # All rates are per cell lineage (accumulated over all divisions from zygote
    # to sampled tissue). Divide by G (mitotic generations) for per-division rates.
    expansion_rate:          tuple[float, float, float]  # mean repeat units gained
    contraction_rate:        tuple[float, float, float]  # mean repeat units lost
    net_mutation_rate:       tuple[float, float, float]  # net directional drift (signed)
    p_any_mutation:          tuple[float, float, float]  # P(≥1 slippage event occurred)
    # Per-division rates (None if mitotic_generations not supplied to fit_locus_r1)
    expansion_rate_per_div:   Optional[tuple[float, float, float]]
    contraction_rate_per_div: Optional[tuple[float, float, float]]
    net_rate_per_div:         Optional[tuple[float, float, float]]
    mitotic_generations:      Optional[float]   # G supplied by caller

    elapsed_s: float

    def __repr__(self) -> str:
        def f(n, v): return f"  {n:<26s}= {v[0]:.4f}  95%CI [{v[1]:.4f}, {v[2]:.4f}]"
        G_str = (f"{self.mitotic_generations:.0f}" if self.mitotic_generations
                 else "not supplied")
        lines = [
            f"Regime1Result(haplotype='{self.haplotype_label}', n={self.n_reads}, L0={self.founder_length:.1f})",
            f"  p₊ prior={self.p_plus_prior:.4f}  p₋ prior={self.p_minus_prior:.4f}"
            f"  q_e={self.q_e:.3f}  q_c={self.q_c:.3f}",
            "",
            "  --- NegBin expansion/contraction parameters ---",
            f("r_e (exp dispersion)",     self.r_e),
            f("r_c (con dispersion)",     self.r_c),
            f("p₊ (exp step size)",       self.p_plus),
            f("p₋ (con step size)",       self.p_minus),
            f("μ_expansion",              self.mu_expansion),
            f("μ_contraction",            self.mu_contraction),
            "",
            "  --- Instability indices ---",
            f("instability_index",        self.instability_index),
            f"    [Var(Δ) = Var(E) + Var(C)]",
            f("net_bias",                 self.net_bias),
            f"    [Mean(Δ) = μ_E − μ_C]",
            f("p_expansion",              self.p_expansion),
            "",
            f"  --- Mutation rates (per lineage, G={G_str} divisions) ---",
            f("expansion_rate",           self.expansion_rate),
            f"    [repeat units gained per lineage]",
            f("contraction_rate",         self.contraction_rate),
            f"    [repeat units lost per lineage]",
            f("net_mutation_rate",        self.net_mutation_rate),
            f"    [net drift per lineage, signed]",
            f("p_any_mutation",           self.p_any_mutation),
            f"    [P(≥1 slippage event in lineage)]",
        ]
        if self.expansion_rate_per_div is not None:
            lines += [
                "",
                f"  --- Per-division rates (G={G_str}) ---",
                f("exp_rate/division",        self.expansion_rate_per_div),
                f("con_rate/division",        self.contraction_rate_per_div),
                f("net_rate/division",        self.net_rate_per_div),
            ]
        lines.append(f"\n  elapsed = {self.elapsed_s*1000:.1f} ms")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        d = {"regime": "1", "haplotype": self.haplotype_label,
             "n_reads": self.n_reads, "founder_length": self.founder_length,
             "p_plus_prior": self.p_plus_prior, "p_minus_prior": self.p_minus_prior,
             "mitotic_generations": self.mitotic_generations}
        per_lineage = ["r_e","r_c","p_plus","p_minus","mu_expansion","mu_contraction",
                       "instability_index","net_bias","p_expansion",
                       "expansion_rate","contraction_rate","net_mutation_rate",
                       "p_any_mutation"]
        for field_name in per_lineage:
            v = getattr(self, field_name)
            d[field_name]          = v[0]
            d[field_name+"_ci_lo"] = v[1]
            d[field_name+"_ci_hi"] = v[2]
        per_div = ["expansion_rate_per_div","contraction_rate_per_div","net_rate_per_div"]
        for field_name in per_div:
            v = getattr(self, field_name)
            if v is not None:
                d[field_name]          = v[0]
                d[field_name+"_ci_lo"] = v[1]
                d[field_name+"_ci_hi"] = v[2]
        d["elapsed_s"] = self.elapsed_s
        return d


# ===========================================================================
# Shared utilities
# ===========================================================================

def _sigmoid_scalar(x: float) -> float:
    """Numerically stable sigmoid for scalar values."""
    return float(1 / (1 + np.exp(-np.clip(x, -30, 30))))


def _sigmoid_vec(x: np.ndarray) -> np.ndarray:
    """Numerically stable sigmoid for arrays."""
    return 1 / (1 + np.exp(-np.clip(x, -30, 30)))


def _summary_stats(deltas: np.ndarray, L0: float) -> np.ndarray:
    """
    Compute summary statistics of the Δ distribution for MDN inference.
    Returns 5D vector: [mean, log(1+var), skewness, clipped_kurtosis, L₀/100]
    L₀/100 conditions the MDN on the founder length, allowing one trained
    network to serve all allele lengths.
    """
    d = np.asarray(deltas, dtype=float)
    v = float(d.var())
    return np.array([
        float(d.mean()),
        float(np.log1p(v)),
        float(_skew(d))     if v > 0.01 else 0.0,
        float(np.clip(_kurt(d), -5, 5)) if v > 0.01 else 0.0,
        float(L0) / 100.0,
    ], dtype=np.float32)


# ===========================================================================
# REGIME 1 — NB-Geometric compound via FFT
# ===========================================================================

def _nb_geom_grid_loglik_r1(deltas: np.ndarray, re_grid: np.ndarray, rc_grid: np.ndarray,
                            pp: float, pm: float, qe: float, qc: float, N: Optional[int] = None) -> np.ndarray:
    """
    Vectorised log-likelihood for the NB-Geometric compound (Regime 1)
    over a 2D grid of (r_e, r_c) with p₊, p₋, q_e, q_c fixed.

    Characteristic function of Δ = E − C:
        φ_Δ(t) = [q_e / (1 − (1−q_e)·M_{S+}(e^{−it}))]^{r_e}
                · [q_c / (1 − (1−q_c)·M_{S−}(e^{+it}))]^{r_c}

    where M_S(z) = p·z/(1−(1−p)·z) is the MGF of Geometric(p).

    Sign convention: expansion uses e^{−it} (positive steps in the IFFT),
    contraction uses e^{+it} (negative steps).

    Returns loglik of shape (n_re, n_rc).
    """
    deltas = np.asarray(deltas, dtype=int)
    unique_d, counts = np.unique(deltas, return_counts=True)

    # FFT length: must span the support of Δ without aliasing
    d_abs = max(int(np.abs(deltas).max()), 1)
    if N is None:
        N = min(int(2 ** np.ceil(np.log2(d_abs * 6 + 16))), 512)

    t   = 2 * np.pi * np.arange(N) / N
    emi = np.exp(-1j * t)   # e^{−it}: expansion direction
    ei  = np.exp( 1j * t)   # e^{+it}: contraction direction

    # FFT indices for observed deltas (negative deltas wrap around)
    obs_idx = unique_d % N

    # MGF of Geometric(p) evaluated at e^{−it} and e^{+it}
    M_sp = pp * emi / (1 - (1 - pp) * emi)   # (N,) expansion steps
    M_sm = pm * ei  / (1 - (1 - pm) * ei)    # (N,) contraction steps

    # log PGF of NB count distribution composed with step MGF
    # log G_{NB(r,q)}(M_S) = r · log(q / (1 − (1−q)·M_S))
    log_pgf_e = np.log(np.clip(qe / (1 - (1 - qe) * M_sp), 1e-300, None))  # (N,)
    log_pgf_c = np.log(np.clip(qc / (1 - (1 - qc) * M_sm), 1e-300, None))  # (N,)

    # CF for each r_e and r_c value: outer product over grid × frequency
    cf_exp = np.exp(np.outer(re_grid, log_pgf_e))   # (n_re, N)
    cf_con = np.exp(np.outer(rc_grid, log_pgf_c))   # (n_rc, N)

    # Combined CF: (n_re, n_rc, N)
    combined = cf_exp[:, None, :] * cf_con[None, :, :]

    # IFFT to get PMF for all (r_e, r_c) simultaneously
    pmf_all = np.clip(np.real(np.fft.ifft(combined, axis=2)), 1e-300, None)
    pmf_all /= pmf_all.sum(axis=2, keepdims=True)

    # Extract PMF at observed delta values and accumulate log-likelihood
    pmf_obs = pmf_all[:, :, obs_idx]   # (n_re, n_rc, n_unique)
    return (np.log(pmf_obs) * counts[None, None, :]).sum(axis=2)   # (n_re, n_rc)


def _regime1_grid_posterior(
    deltas: np.ndarray,
    re_grid: np.ndarray,
    rc_grid: np.ndarray,
    pp_grid: np.ndarray,
    pm_grid: np.ndarray,
    pp_prior_mean: float,
    pm_prior_mean: float,
    qe: float,
    qc: float,
    pp_log_sd: float = DEFAULT_P_LOG_SD,
    pm_log_sd: float = DEFAULT_P_LOG_SD,
) -> np.ndarray:
    """
    4D grid posterior over (r_e, r_c, p₊, p₋) for Regime 1.

    p₊ and p₋ float on a narrow grid around the global estimate,
    allowing per-locus data to correct for any bias in the global estimate.

    Priors:
        r_e, r_c ~ Gamma(1.5, scale=1.0)   — moderate event rates
        logit(p₊) ~ Normal(logit(prior), σ²) — log-odds Normal prior
        logit(p₋) ~ Normal(logit(prior), σ²)

    Returns post of shape (n_re, n_rc, n_pp, n_pm), sums to 1.
    """
    n_re, n_rc, n_pp, n_pm = len(re_grid), len(rc_grid), len(pp_grid), len(pm_grid)
    loglik = np.full((n_re, n_rc, n_pp, n_pm), -np.inf)

    # Loop over p₊ values (outer), vectorise over (r_e, r_c, p₋) (inner)
    for ip, pp in enumerate(pp_grid):
        for im, pm in enumerate(pm_grid):
            loglik[:, :, ip, im] = _nb_geom_grid_loglik_r1(
                deltas, re_grid, rc_grid, pp, pm, qe, qc
            )

    # Log-Normal prior on logit(p₊) and logit(p₋) — keeps p values in (0,1)
    def logit(p): return np.log(p / (1 - p + 1e-10) + 1e-10)

    log_prior = (
        stats.gamma.logpdf(re_grid, a=1.5, scale=1.0)[:, None, None, None]
        + stats.gamma.logpdf(rc_grid, a=1.5, scale=1.0)[None, :, None, None]
        + stats.norm.logpdf(logit(pp_grid), logit(pp_prior_mean), pp_log_sd
                            )[None, None, :, None]
        + stats.norm.logpdf(logit(pm_grid), logit(pm_prior_mean), pm_log_sd
                            )[None, None, None, :]
    )

    log_post = loglik + log_prior
    log_post -= log_post.max()
    post = np.exp(log_post)
    post /= post.sum()
    return post


def _regime1_posterior_summaries(
    post: np.ndarray,
    re_grid: np.ndarray,
    rc_grid: np.ndarray,
    pp_grid: np.ndarray,
    pm_grid: np.ndarray,
    qe: float,
    qc: float,
) -> dict[str, tuple[float, float, float]]:
    """
    Compute posterior mean and 95% CI for all Regime 1 quantities.

    Instability index (closed form):
        Var(Δ) = Var(E) + Var(C)
        Var(E) = μ_{N_e}·(1−p₊)/p₊² + σ²_{N_e}/p₊²
        where μ_{N_e} = r_e·(1−q_e)/q_e  and  σ²_{N_e} = r_e·(1−q_e)/q_e²
    """
    RE, RC, PP, PM = np.meshgrid(re_grid, rc_grid, pp_grid, pm_grid, indexing="ij")

    # NegBin moments for expansion and contraction count distributions
    mu_Ne  = RE * (1 - qe) / qe         # mean N_e (expansion events)
    var_Ne = RE * (1 - qe) / qe ** 2    # variance N_e
    mu_Nc  = RC * (1 - qc) / qc         # mean N_c (contraction events)
    var_Nc = RC * (1 - qc) / qc ** 2    # variance N_c

    # Mean total expansion E and contraction C (each is a sum of N Geometric steps)
    # E[E] = E[N_e] · E[S₊] = μ_{N_e} · (1−p₊)/p₊  (Geometric mean = (1-p)/p for support 0,1,2..)
    # Actually Geometric with support {1,2,...}: E[S] = 1/p
    # Here we use support {1,2,...} so E[S₊] = 1/p₊
    mu_E = mu_Ne / PP       # mean total expansion
    mu_C = mu_Nc / PM       # mean total contraction

    # Var(E) = E[N_e]·Var(S₊) + Var(N_e)·E[S₊]²
    # Var(S₊) = (1−p₊)/p₊²  (Geometric variance, support {1,2,...})
    var_E = mu_Ne * (1 - PP) / PP ** 2 + var_Ne / PP ** 2
    var_C = mu_Nc * (1 - PM) / PM ** 2 + var_Nc / PM ** 2

    # P(at least one slippage event occurred in a lineage)
    # = 1 - P(N_e=0) * P(N_c=0)
    # For NegBin(r,q): P(N=0) = (q / (2-q))^r  [PMF at k=0]
    # This is the probability a lineage has accumulated ZERO mutations —
    # complement is the probability of at least one observed slippage event.
    p_zero_e = (qe / (2.0 - qe + 1e-10)) ** RE
    p_zero_c = (qc / (2.0 - qc + 1e-10)) ** RC
    p_any_mut = 1.0 - p_zero_e * p_zero_c

    derived = {
        "r_e":               RE,
        "r_c":               RC,
        "p_plus":            PP,
        "p_minus":           PM,
        "mu_expansion":      mu_E,
        "mu_contraction":    mu_C,
        "instability_index": var_E + var_C,            # Var(Δ) = Var(E) + Var(C)
        "net_bias":          mu_E - mu_C,              # Mean(Δ)
        "p_expansion":       (mu_E > mu_C).astype(float),
        # Mutation rates (per lineage = accumulated over all cell divisions)
        "expansion_rate":    mu_E,                     # = μ_{N_e}/p₊  repeat units gained
        "contraction_rate":  mu_C,                     # = μ_{N_c}/p₋  repeat units lost
        "net_mutation_rate": mu_E - mu_C,              # signed net drift per lineage
        "p_any_mutation":    p_any_mut,                # P(≥1 slippage event)
    }

    post_flat = post.ravel()
    results = {}
    for name, grid in derived.items():
        mean_ = float((post * grid).sum())
        gf    = grid.ravel()
        si    = np.argsort(gf)
        cdf   = np.cumsum(post_flat[si])
        lo95  = float(gf[si[np.searchsorted(cdf, 0.025)]])
        hi95  = float(gf[si[np.searchsorted(cdf, 0.975)]])
        results[name] = (mean_, lo95, hi95)
    return results


# ===========================================================================
# REGIME 1 — Pass 1: Global parameter estimation
# ===========================================================================

def estimate_regime1_global(
    loci_data: list[dict],
    min_reads: int = 10,
    n_pp_grid: int = 20,
    n_pm_grid: int = 20,
) -> Regime1GlobalParams:
    """
    Estimate global Regime 1 parameters by pooling across loci.

    For each locus, the empirical mean and variance of Δ constrain p₊, p₋,
    and the ratio r_e/r_c. Pooling across many loci with different L₀ values
    allows regression of logit(p) on L₀ to estimate the length-dependence slopes.

    Parameters
    ----------
    loci_data : list of dicts with keys:
        'deltas'         : int array  — Δ_i = observed - founder
        'founder_length' : float      — L₀

    Method
    ------
    For each locus: estimate p₊, p₋ from the method of moments using the
    asymmetry between positive and negative tails of the Δ distribution.
    Then regress logit(p₊) and logit(p₋) on L₀ across loci.
    """
    L0s      = []
    pp_ests  = []   # per-locus p₊ estimates
    pm_ests  = []   # per-locus p₋ estimates
    weights  = []

    for locus in loci_data:
        deltas = np.asarray(locus["deltas"], dtype=float)
        L0     = float(locus["founder_length"])
        n      = len(deltas)
        if n < min_reads:
            continue

        # Separate positive (expansion) and negative (contraction) deltas
        pos = deltas[deltas > 0]
        neg = -deltas[deltas < 0]   # flip to positive

        if len(pos) < 3 or len(neg) < 3:
            continue

        # Method of moments: for Geometric(p), mean = 1/p
        # Use reciprocal of mean step size as estimate of p
        pp_est = np.clip(1.0 / max(pos.mean(), 0.5), 0.05, 0.95)
        pm_est = np.clip(1.0 / max(neg.mean(), 0.5), 0.05, 0.95)

        L0s.append(L0)
        pp_ests.append(pp_est)
        pm_ests.append(pm_est)
        weights.append(np.sqrt(n))

    if len(L0s) < 4:
        # Not enough loci — return default (no length dependence)
        return Regime1GlobalParams(n_loci_used=len(L0s))

    L0s     = np.array(L0s)
    pp_ests = np.array(pp_ests)
    pm_ests = np.array(pm_ests)
    W       = np.array(weights)
    X       = np.column_stack([np.ones(len(L0s)), L0s])

    def logit(p): return np.log(p / (1 - p + 1e-10) + 1e-10)

    def wls(y):
        Xw = X * W[:, None]; yw = y * W
        return np.linalg.lstsq(Xw, yw, rcond=None)[0]

    coef_p = wls(logit(pp_ests))
    coef_m = wls(logit(pm_ests))

    return Regime1GlobalParams(
        lp0=float(coef_p[0]), lp1=float(coef_p[1]),
        lm0=float(coef_m[0]), lm1=float(coef_m[1]),
        q_e=0.5, q_c=0.5,   # default — update from richer data if available
        n_loci_used=len(L0s),
    )


# ===========================================================================
# Public API
# ===========================================================================

def fit_locus_r1(
    lengths,
    founder_length: float,
    global_params: Regime1GlobalParams,
    haplotype_label: str = "haplotype",
    re_grid: Optional[np.ndarray] = None,
    rc_grid: Optional[np.ndarray] = None,
    n_pp: int = DEFAULT_PP_N,
    n_pm: int = DEFAULT_PM_N,
    mitotic_generations: Optional[float] = None,
    seed: int = 0,
) -> Regime1Result:
    """
    Fit Regime 1 (NB-Geometric compound) for one phased haplotype.

    Parameters
    ----------
    lengths : array-like of int or float
        Per-read allele lengths for this haplotype.
    founder_length : float
        Germline / founder allele length L₀.
    global_params : Regime1GlobalParams
        p₊(L₀), p₋(L₀), q_e, q_c from Pass 1.
    haplotype_label : str
        Label for output.
    re_grid, rc_grid : optional custom grids for r_e, r_c.
    n_pp, n_pm : number of p₊, p₋ grid points around global estimate.
    mitotic_generations : float, optional
        Number of mitotic cell divisions from zygote to sampled tissue (G).
        If supplied, per-division mutation rates are computed in addition to
        per-lineage rates. Typical values:
            Blood (adult):         50-70
            Colon epithelium:      1000-2000
            Neurons (adult):       ~0  (post-mitotic)
            Fibroblasts (adult):   ~40-60
        If None, only per-lineage rates are reported.
    """
    t0 = time.time()

    lengths     = np.asarray(lengths)
    lengths_int = np.round(lengths).astype(int)
    mode_int    = int(round(founder_length))
    deltas      = (lengths_int - mode_int).astype(int)

    re_grid = re_grid if re_grid is not None else DEFAULT_RE_GRID
    rc_grid = rc_grid if rc_grid is not None else DEFAULT_RC_GRID

    # Global p₊ and p₋ at this L₀
    pp_prior = global_params.p_plus(founder_length)
    pm_prior = global_params.p_minus(founder_length)
    qe       = global_params.q_e
    qc       = global_params.q_c

    # Build narrow grids around global estimates in logit space
    def logit(p): return np.log(p / (1 - p + 1e-10) + 1e-10)
    logit_pp_centre = logit(pp_prior)
    logit_pm_centre = logit(pm_prior)
    pp_grid = _sigmoid_vec(np.linspace(
        logit_pp_centre - 2 * DEFAULT_P_LOG_SD,
        logit_pp_centre + 2 * DEFAULT_P_LOG_SD, n_pp))
    pm_grid = _sigmoid_vec(np.linspace(
        logit_pm_centre - 2 * DEFAULT_P_LOG_SD,
        logit_pm_centre + 2 * DEFAULT_P_LOG_SD, n_pm))

    post = _regime1_grid_posterior(
        deltas, re_grid, rc_grid, pp_grid, pm_grid,
        pp_prior_mean=pp_prior, pm_prior_mean=pm_prior,
        qe=qe, qc=qc,
    )
    s = _regime1_posterior_summaries(post, re_grid, rc_grid, pp_grid, pm_grid, qe, qc)

    # Compute per-division rates if G supplied
    def _scale(tup, G):
        if G is None or G <= 0:
            return None
        return (tup[0]/G, tup[1]/G, tup[2]/G)

    return Regime1Result(
        haplotype_label          = haplotype_label,
        n_reads                  = len(lengths_int),
        founder_length           = founder_length,
        p_plus_prior             = pp_prior,
        p_minus_prior            = pm_prior,
        q_e                      = qe,
        q_c                      = qc,
        r_e                      = s["r_e"],
        r_c                      = s["r_c"],
        p_plus                   = s["p_plus"],
        p_minus                  = s["p_minus"],
        mu_expansion             = s["mu_expansion"],
        mu_contraction           = s["mu_contraction"],
        instability_index        = s["instability_index"],
        net_bias                 = s["net_bias"],
        p_expansion              = s["p_expansion"],
        expansion_rate           = s["expansion_rate"],
        contraction_rate         = s["contraction_rate"],
        net_mutation_rate        = s["net_mutation_rate"],
        p_any_mutation           = s["p_any_mutation"],
        expansion_rate_per_div   = _scale(s["expansion_rate"],   mitotic_generations),
        contraction_rate_per_div = _scale(s["contraction_rate"], mitotic_generations),
        net_rate_per_div         = _scale(s["net_mutation_rate"], mitotic_generations),
        mitotic_generations      = mitotic_generations,
        elapsed_s                = time.time() - t0,
    )


# ===========================================================================
# AUTO-DETECTION AND REGIME ROUTING
# ===========================================================================

@dataclass
class DistributionDiagnostics:
    """
    Pre-fit diagnostics computed from raw reads before model selection.
    Used to route loci to the correct regime and flag data quality issues.
    """
    n_reads:          int
    founder_length:   float
    modal_length:     float
    n_outliers:       int         # reads > 3 MAD from median
    outlier_lengths:  list        # the actual outlier values
    delta_mean:       float
    delta_std:        float
    delta_max:        float
    delta_min:        float
    tail_ratio:       float       # max(|delta|) / std(delta): > 5 → Regime 2
    skewness:         float
    is_heavy_tailed:  bool        # tail_ratio > 5
    recommended_regime: int       # 1 or 2
    warnings:         list        # list of warning strings

    def __repr__(self) -> str:
        lines = [
            f"DistributionDiagnostics(n={self.n_reads}, L0={self.founder_length:.1f})",
            f"  modal_length      = {self.modal_length:.1f}",
            f"  delta range       = [{self.delta_min:.0f}, {self.delta_max:.0f}]",
            f"  delta mean/std    = {self.delta_mean:.1f} / {self.delta_std:.1f}",
            f"  tail_ratio        = {self.tail_ratio:.2f}  "
            f"(max|Δ|/std, >5 → Regime 2)",
            f"  skewness          = {self.skewness:.2f}",
            f"  n_outliers        = {self.n_outliers}",
            f"  is_heavy_tailed   = {self.is_heavy_tailed}",
            f"  recommended_regime= {self.recommended_regime}",
        ]
        if self.outlier_lengths:
            lines.append(f"  outlier_lengths   = {self.outlier_lengths}")
        for w in self.warnings:
            lines.append(f"  WARNING: {w}")
        return "\n".join(lines)


def diagnose_distribution(
    lengths,
    founder_length: float,
    tail_ratio_threshold: float = 5.0,
    outlier_mad_factor:   float = 5.0,
) -> DistributionDiagnostics:
    """
    Compute pre-fit diagnostics on a haplotype's read-length distribution
    and recommend the appropriate regime.

    The key diagnostic is the **tail ratio**: max(|Δ|) / std(Δ).
    - Ratio ≤ 5: distribution consistent with a single NB-Geometric process.
                 Regime 1 is appropriate.
    - Ratio > 5: distribution has extreme outliers or a heavy tail that
                 a single NB-Geometric process cannot explain — the maximum
                 observed length change is more than 5 standard deviations
                 from the centre, which is essentially impossible under a
                 light-tailed NB-Geometric model. Regime 2 is appropriate.

    Additionally flags:
    - Individual outlier reads (>5 MAD from median) that may represent
      reads from a different allele leaking through phasing, or mapping
      artefacts — these should be reviewed before fitting.
    - Very low read counts (n < 5) where neither regime gives reliable results.
    - Collapsed posteriors (detected in fit_locus_r1 via boundary check).

    Parameters
    ----------
    lengths : array-like
        Per-read allele lengths for one phased haplotype.
    founder_length : float
        Germline/founder allele length L₀.
    tail_ratio_threshold : float
        Ratio above which Regime 2 is recommended (default 5.0).
    outlier_mad_factor : float
        Reads more than this many MADs from the median are flagged
        as potential outliers (default 5.0, conservative).

    Returns
    -------
    DistributionDiagnostics
    """
    from scipy.stats import skew as _skew_fn
    lengths     = np.asarray(lengths)
    lengths_int = np.round(lengths).astype(int)
    n           = len(lengths_int)
    warnings_   = []

    # Modal length
    vals, cts = np.unique(lengths_int, return_counts=True)
    modal_length = float(vals[np.argmax(cts)])

    # Deltas
    deltas = (lengths - founder_length).astype(float)

    # Basic stats
    d_mean = float(deltas.mean())
    d_std  = float(deltas.std()) if n > 1 else 0.0
    d_max  = float(deltas.max())
    d_min  = float(deltas.min())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        d_skew = float(_skew_fn(deltas)) if n > 2 else 0.0

    # Tail ratio: how many std deviations is the largest |delta|?
    tail_ratio = float(max(abs(d_max), abs(d_min)) / max(d_std, 1e-6))

    # Outlier detection via MAD
    median_d = float(np.median(deltas))
    mad      = float(np.median(np.abs(deltas - median_d)))
    if mad < 1e-6:
        mad = d_std  # fallback for very tight distributions
    outlier_mask   = np.abs(deltas - median_d) > outlier_mad_factor * mad
    outlier_lengths = sorted(lengths_int[outlier_mask].tolist())
    n_outliers      = int(outlier_mask.sum())

    # Warnings
    if n < 5:
        warnings_.append(f"Very low read depth (n={n}). "
                         "Results unreliable for either regime.")
    if n_outliers > 0:
        warnings_.append(
            f"{n_outliers} read(s) flagged as outliers "
            f"(>{outlier_mad_factor}×MAD from median): {outlier_lengths}. "
            "Consider whether these are from a second allele or mapping artefact. "
            "Remove with filter_outliers=True in fit_locus_r1 if confirmed artefacts."
        )
    if tail_ratio > tail_ratio_threshold:
        warnings_.append(
            f"Heavy tail detected (tail_ratio={tail_ratio:.1f} > {tail_ratio_threshold}). "
            "A single NB-Geometric process cannot explain this distribution. "
            "Regime 2 (length-dependent Markov chain) is required."
        )
    if d_std < 1.0 and n > 5:
        warnings_.append(
            "Near-zero variance in delta distribution. "
            "Locus may be completely stable — instability index will be near zero."
        )

    is_heavy_tailed     = tail_ratio > tail_ratio_threshold
    recommended_regime  = 2 if is_heavy_tailed else 1

    return DistributionDiagnostics(
        n_reads           = n,
        founder_length    = founder_length,
        modal_length      = modal_length,
        n_outliers        = n_outliers,
        outlier_lengths   = outlier_lengths,
        delta_mean        = d_mean,
        delta_std         = d_std,
        delta_max         = d_max,
        delta_min         = d_min,
        tail_ratio        = tail_ratio,
        skewness          = d_skew,
        is_heavy_tailed   = is_heavy_tailed,
        recommended_regime= recommended_regime,
        warnings          = warnings_,
    )


def filter_outliers(
    lengths,
    founder_length: float,
    mad_factor: float = 5.0,
) -> tuple[np.ndarray, list]:
    """
    Remove outlier reads before fitting Regime 1.

    Outliers are reads more than mad_factor × MAD from the median delta.
    Returns (filtered_lengths, removed_lengths).

    Use this when diagnose_distribution() flags individual outlier reads
    that are suspected phasing errors or mapping artefacts — NOT when the
    entire tail is flagged as heavy (that requires Regime 2, not filtering).

    Parameters
    ----------
    lengths : array-like
        Per-read allele lengths.
    founder_length : float
        Germline allele length L₀.
    mad_factor : float
        Outlier threshold in MAD units (default 5.0).

    Returns
    -------
    (filtered_lengths, removed_lengths)
    """
    lengths  = np.asarray(lengths)
    deltas   = lengths - founder_length
    median_d = float(np.median(deltas))
    mad      = float(np.median(np.abs(deltas - median_d)))
    if mad < 1e-6:
        mad = float(deltas.std())
    keep = np.abs(deltas - median_d) <= mad_factor * mad
    return lengths[keep], sorted(lengths[~keep].tolist())


def auto_fit(
    lengths,
    founder_length: float,
    regime1_params: Regime1GlobalParams,
    haplotype_label: str = "haplotype",
    tail_ratio_threshold: float = 5.0,
    outlier_mad_factor: float = 5.0,
    filter_outliers_r1: bool = False,
    mitotic_generations: Optional[float] = None,
    seed: int = 0,
    **kwargs,
) -> tuple[object, "DistributionDiagnostics", Optional[dict]]:
    """
    Automatically diagnose the read-length distribution and route to the
    appropriate regime.

    This is the recommended entry point for genome-wide use. It:
    1. Runs diagnose_distribution() to assess the data
    2. If heavy-tailed (tail_ratio > threshold)
    3. Otherwise routes to fit_locus_r1(), optionally filtering outliers first
    4. Returns (result, diagnostics) so the routing decision is always visible

    Parameters
    ----------
    lengths : array-like
        Per-read allele lengths for one phased haplotype.
    founder_length : float
        Germline allele length L₀.
    regime1_params : Regime1GlobalParams
        Required. Regime 1 global parameters.
    haplotype_label : str
        Label for output.
    tail_ratio_threshold : float
        Threshold for Regime 2 routing (default 5.0).
    outlier_mad_factor : float
        MAD factor for outlier flagging (default 5.0).
    filter_outliers_r1 : bool
        If True and running Regime 1, remove outlier reads before fitting.
        Use only when outliers are confirmed artefacts, not genuine biology.
    mitotic_generations : float, optional
        G for per-division mutation rates (Regime 1 only).
    **kwargs
        Additional arguments forwarded to fit_locus_r1.

    Returns
    -------
    (result, diagnostics, gof) :
        result      : Regime1Result
        diagnostics : DistributionDiagnostics
        gof         : dict or None
            Goodness-of-fit metrics for Regime 1 fits (same keys as
            _compute_gof: ks_statistic, ad_statistic, ppp_variance,
            ppp_skewness, observed_var, fitted_var, dispersion_ratio,
            chisq_statistic, chisq_pvalue, chisq_df, chisq_n_bins, fit_quality,
            fit_quality_adjusted).

    Note: n_ppp (number of PPP simulations, default 300) can be passed
    via **kwargs to control GoF computation speed.
    """
    diag = diagnose_distribution(
        lengths, founder_length,
        tail_ratio_threshold=tail_ratio_threshold,
        outlier_mad_factor=outlier_mad_factor,
    )

    lengths_fit = np.asarray(lengths)
    gof         = None   # populated only for Regime 1 fits

    # Extract n_ppp before forwarding kwargs to fit_locus_r1
    # so it does not cause an unexpected keyword argument error
    n_ppp = kwargs.pop("n_ppp", 300)

    # Optionally filter outliers before Regime 1
    if filter_outliers_r1 and diag.n_outliers > 0:
        lengths_fit, removed = filter_outliers(
            lengths_fit, founder_length, mad_factor=outlier_mad_factor
        )
        diag.warnings.append(
            f"Removed {len(removed)} outlier reads before Regime 1 fit: {removed}"
        )

    result = fit_locus_r1(
        lengths_fit, founder_length, regime1_params,
        haplotype_label=haplotype_label,
        mitotic_generations=mitotic_generations,
        seed=seed, **kwargs
    )

    # Goodness-of-fit for Regime 1 fits
    # n_ppp=0 means skip GoF entirely (useful for speed when not needed)
    if n_ppp == 0:
        gof = None
    else:
        try:
            lengths_int = np.round(lengths_fit).astype(int)
            deltas      = (lengths_int - int(round(result.founder_length))).astype(int)
            gof = _compute_gof(
                deltas,
                r_e = result.r_e[0],
                r_c = result.r_c[0],
                pp  = result.p_plus[0],
                pm  = result.p_minus[0],
                qe  = regime1_params.q_e,
                qc  = regime1_params.q_c,
                n_ppp = n_ppp,
                seed  = seed + 1,
            )
        except Exception as e:
            diag.warnings.append(f"GoF computation failed: {e}")

    return result, diag, gof


# ===========================================================================
# GENOME-WIDE REGIME 1 ANALYSIS
# ===========================================================================

@dataclass
class LocusResult:
    """
    Complete Regime 1 results for one locus (both haplotypes).
    Includes per-haplotype instability estimates and goodness-of-fit metrics.
    """
    locus_id:  str
    h1:        Optional[Regime1Result]
    h2:        Optional[Regime1Result]
    h1_gof:    Optional[dict]   # goodness-of-fit metrics for H1
    h2_gof:    Optional[dict]   # goodness-of-fit metrics for H2
    h1_diag:   Optional[DistributionDiagnostics]
    h2_diag:   Optional[DistributionDiagnostics]
    h1_error:  Optional[str]    # error message if H1 fit failed
    h2_error:  Optional[str]    # error message if H2 fit failed
    elapsed_s: float

    def to_rows(self) -> list[dict]:
        """
        Flatten to a list of dicts (one per haplotype) for TSV output.
        Returns up to 2 rows per locus.
        """
        rows = []
        for hap, result, gof, diag, err in [
            ("H1", self.h1, self.h1_gof, self.h1_diag, self.h1_error),
            ("H2", self.h2, self.h2_gof, self.h2_diag, self.h2_error),
        ]:
            row = {"locus_id": self.locus_id, "haplotype": hap}
            if err:
                row["error"] = err
                rows.append(row)
                continue
            if result is None:
                continue

            # Core result fields
            row.update(result.to_dict())

            # Goodness-of-fit
            if gof:
                for k, v in gof.items():
                    row[f"gof_{k}"] = round(v, 6) if isinstance(v, float) else v

            # Diagnostics
            if diag:
                row["diag_n_reads"]        = diag.n_reads
                row["diag_tail_ratio"]     = round(diag.tail_ratio, 3)
                row["diag_skewness"]       = round(diag.skewness, 3)
                row["diag_n_outliers"]     = diag.n_outliers
                row["diag_is_heavy_tailed"]= int(diag.is_heavy_tailed)
                row["diag_recommended_regime"] = diag.recommended_regime

            rows.append(row)
        return rows


def _compute_gof(
    deltas: np.ndarray,
    r_e:    float,
    r_c:    float,
    pp:     float,
    pm:     float,
    qe:     float,
    qc:     float,
    n_ppp:  int = 300,
    seed:   int = 0,
) -> dict:
    """
    Compute goodness-of-fit metrics for a Regime 1 NB-Geometric fit.

    Metrics
    -------
    ks_statistic : float [0, 1]
        Kolmogorov-Smirnov statistic: max|F_empirical - F_fitted| over all
        observed delta values. Tests whether the bulk of the distribution
        is well described. Lower is better.

    ad_statistic : float [0, ∞)
        Anderson-Darling statistic: a weighted integral of (F_emp - F_fit)²,
        with higher weight in the tails. More sensitive than KS to tail misfit.
        Values below ~5 indicate good fit; values above 100 indicate severe
        tail misfit (as in Sample 1). Lower is better.

    ppp_variance : float [0, 1]
        Posterior predictive p-value using variance as the test statistic.
        Fraction of simulated datasets from the fitted model whose variance
        exceeds the observed variance. Values near 0.5 indicate good fit.
        Values near 0 indicate the model is systematically underdispersed —
        the primary signal that Regime 2 is needed.

    ppp_skewness : float [0, 1]
        Same but using skewness as the test statistic. Values near 0 indicate
        the model underestimates the right-skewness of the data.

    dispersion_ratio : float
        observed_var / fitted_var. Values near 1.0 = good fit.
        Values >> 1.0 = model underdispersed (Regime 2 needed).
        Values << 1.0 = model overdispersed (unlikely with NB-Geometric).

    chisq_statistic : float
        Pearson Chi-Square statistic computed over quantile-based bins of
        the fitted distribution. Bins are defined by equal-probability
        quantiles of the fitted PMF, so each bin has the same expected
        count under the model. Adjacent bins are merged until all expected
        counts are ≥ 5 (the standard Chi-Square validity requirement).
        Lower is better. Unlike AD and KS, Chi-Square has a calibrated
        null distribution (chi-squared with known df), making its p-value
        formally interpretable.

    chisq_pvalue : float [0, 1]
        P-value under the chi-squared null distribution with
        df = n_bins_used - 2 - 1 (subtracting 2 for the per-locus
        estimated parameters r_e and r_c, and 1 for the constraint that
        counts sum to n). Values below 0.05 indicate statistically
        significant misfit. The most formally rigorous of the GoF metrics
        for reporting in publications.

    chisq_df : int
        Degrees of freedom used for the chi-squared p-value.

    chisq_n_bins : int
        Number of bins after merging. Fewer bins = less powerful test.
        If < 3, chi-square results are unreliable (reported as NaN).
    """
    rng    = np.random.default_rng(seed)
    deltas = np.asarray(deltas, dtype=int)
    n      = len(deltas)

    # Degenerate case: all reads at the same length (zero variance).
    # The model fits trivially but GoF metrics are undefined.
    if np.all(deltas == deltas[0]):
        return {
            "ks_statistic":        0.0,
            "ad_statistic":        0.0,
            "ppp_variance":        0.5,
            "ppp_skewness":        0.5,
            "observed_var":        0.0,
            "fitted_var":          0.0,
            "dispersion_ratio":    1.0,
            "chisq_statistic":     float("nan"),
            "chisq_pvalue":        float("nan"),
            "chisq_df":            0,
            "chisq_n_bins":        0,
            "fit_quality":         float("nan"),
            "fit_quality_adjusted":float("nan"),
        }

    # Fitted PMF at posterior mean parameters.
    # Extend the range well beyond the observed deltas so the FFT grid
    # captures enough probability mass for reliable sampling in the PPP step.
    # Minimum span of 40 ensures ks_vals is never empty even for tight
    # distributions where lo and hi would otherwise be equal or very close.
    obs_span = max(int(deltas.max()) - int(deltas.min()), 1)
    lo = int(deltas.min()) - max(5, obs_span)
    hi = int(deltas.max()) + max(5, obs_span)
    N  = min(int(2**np.ceil(np.log2((hi - lo) * 3 + 32))), 8192)
    N  = max(N, 64)

    t  = 2 * np.pi * np.arange(N) / N
    emi, ei = np.exp(-1j * t), np.exp(1j * t)
    M_sp = pp * emi / (1 - (1 - pp) * emi)
    M_sm = pm * ei  / (1 - (1 - pm) * ei)
    log_pgf_e = np.log(np.clip(qe / (1 - (1 - qe) * M_sp), 1e-300, None))
    log_pgf_c = np.log(np.clip(qc / (1 - (1 - qc) * M_sm), 1e-300, None))
    cf      = np.exp(r_e * log_pgf_e) * np.exp(r_c * log_pgf_c)
    pmf_raw = np.real(np.fft.ifft(cf))
    dg      = np.arange(N); dg[dg >= N // 2] -= N
    order   = np.argsort(dg); dg = dg[order]; pmf_raw = pmf_raw[order]
    mask    = (dg >= lo) & (dg <= hi)
    ks_vals = dg[mask]
    pmf     = np.clip(pmf_raw[mask], 1e-300, None)

    # Guard: if the mask produced no values fall back to the full sorted grid.
    if len(ks_vals) == 0:
        ks_vals = dg
        pmf     = np.clip(pmf_raw, 1e-300, None)

    pmf /= pmf.sum()
    cdf  = np.cumsum(pmf)

    # Fitted mean and variance
    fitted_mean = float(np.sum(ks_vals * pmf))
    fitted_var  = float(np.sum(ks_vals ** 2 * pmf) - fitted_mean ** 2)

    # KS statistic
    unique_d = np.sort(np.unique(deltas))
    ecdf     = np.array([(deltas <= d).mean() for d in unique_d])
    cdf_dict = dict(zip(ks_vals.tolist(), cdf.tolist()))
    fcdf     = np.array([cdf_dict.get(int(d), 1.0) for d in unique_d])
    ks_stat  = float(np.max(np.abs(ecdf - fcdf)))

    # Anderson-Darling statistic (tail-weighted)
    eps     = 1e-10
    F       = np.clip(fcdf, eps, 1 - eps)
    ad_stat = float(n * np.mean((ecdf - F) ** 2 / (F * (1 - F))))

    # Posterior predictive p-value
    obs_var  = float(np.var(deltas))
    obs_skew = float(stats.skew(deltas.astype(float))) if n > 2 else 0.0
    sim_vars  = np.zeros(n_ppp)
    sim_skews = np.zeros(n_ppp)
    for i in range(n_ppp):
        sim_d        = rng.choice(ks_vals, size=n, p=pmf)
        sim_vars[i]  = float(np.var(sim_d))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            sim_skews[i] = float(stats.skew(sim_d.astype(float))) if n > 2 else 0.0

    if n_ppp == 0:
        ppp_var  = 0.5   # neutral value when PPP not computed
        ppp_skew = 0.5
    else:
        ppp_var  = float((sim_vars  >= obs_var).mean())
        ppp_skew = float((sim_skews >= obs_skew).mean())

    disp_ratio = obs_var / max(fitted_var, 1e-6)

    # ── Chi-Square GoF with quantile-based binning ────────────────────────
    # Bin edges are quantiles of the FITTED distribution so each bin has
    # equal expected probability under the model. Bins with expected count
    # < 5 are merged with their smallest neighbour until the requirement
    # is met, then chi-square is computed with df = n_bins - 2 - 1.
    chisq_stat = float('nan')
    chisq_pval = float('nan')
    chisq_df   = 0
    chisq_bins = 0

    if len(ks_vals) >= 3:
        # Adaptive bin count: target ~5 reads per bin, min 5 bins max 10
        # With 5 bins and 2 estimated per-locus params, df = 5-2-1 = 2 (minimum valid)
        n_bins_target = min(10, max(5, n // 5))
        fitted_cdf    = np.cumsum(pmf)

        # Quantile bin edges from fitted CDF
        quantile_probs = np.linspace(0, 1, n_bins_target + 1)
        bin_edges = []
        for q in quantile_probs:
            idx = int(np.searchsorted(fitted_cdf, q))
            idx = min(idx, len(ks_vals) - 1)
            bin_edges.append(int(ks_vals[idx]))
        bin_edges = sorted(set(bin_edges))

        if len(bin_edges) >= 3:
            inner = bin_edges[1:-1]          # interior cut points
            n_b   = len(inner) + 1           # number of bins

            # Observed counts per bin
            obs_bins = np.digitize(deltas, inner)   # 0..n_b-1
            obs_c = np.array(
                [(obs_bins == b).sum() for b in range(n_b)], dtype=float
            )

            # Expected counts per bin from fitted PMF
            exp_c = np.zeros(n_b)
            for b in range(n_b):
                if b == 0:
                    mask_b = ks_vals <= inner[0]
                elif b == n_b - 1:
                    mask_b = ks_vals > inner[-1]
                else:
                    mask_b = (ks_vals > inner[b-1]) & (ks_vals <= inner[b])
                exp_c[b] = float(pmf[mask_b].sum()) * n

            # Merge bins with expected count < 5 into smallest neighbour
            # until all bins satisfy the requirement or only 2 bins remain
            while len(exp_c) > 2 and exp_c.min() < 5.0:
                idx_min = int(np.argmin(exp_c))
                if idx_min == 0:
                    merge = 1
                elif idx_min == len(exp_c) - 1:
                    merge = len(exp_c) - 2
                else:
                    merge = (idx_min - 1
                             if exp_c[idx_min-1] <= exp_c[idx_min+1]
                             else idx_min + 1)
                lo_m = min(idx_min, merge)
                hi_m = max(idx_min, merge)
                exp_c = np.concatenate([
                    exp_c[:lo_m],
                    [exp_c[lo_m] + exp_c[hi_m]],
                    exp_c[hi_m+1:]
                ])
                obs_c = np.concatenate([
                    obs_c[:lo_m],
                    [obs_c[lo_m] + obs_c[hi_m]],
                    obs_c[hi_m+1:]
                ])

            chisq_bins = len(exp_c)
            if chisq_bins >= 2:
                valid_b   = exp_c > 0
                chisq_stat = float(
                    np.sum((obs_c[valid_b] - exp_c[valid_b])**2
                           / exp_c[valid_b])
                )
                # df = n_bins - n_estimated_params - 1
                # n_estimated_params = 2 (r_e and r_c; p+/p- come from global)
                # df = n_bins - n_estimated_params - 1
                # Flag as unreliable (NaN) when fewer than 3 bins remain
                # because df <= 0 gives a degenerate test
                raw_df = chisq_bins - 2 - 1
                # Require df >= 2: with df=1 the test has very low power
                # and is prone to false positives from single bin imbalances
                if raw_df < 2:
                    chisq_stat = float('nan')
                    chisq_pval = float('nan')
                    chisq_df   = 0
                else:
                    chisq_df   = raw_df
                    chisq_pval = float(
                        1 - stats.chi2.cdf(chisq_stat, chisq_df)
                    )

    # When the model is OVERdispersed (disp_ratio < 1), chi-square misfit
    # means the global step-size priors are too loose for this locus -- a
    # calibration issue.
    chisq_flag = (
        not np.isnan(chisq_pval) and
        chisq_df >= 2               and
        chisq_bins >= 3             and
        chisq_pval < 0.05           and
        disp_ratio > 1.0               # only flag underdispersed model
    )

    # ── Asymmetric fit quality score ─────────────────────────────────────
    # Standard formula 2*min(PPP,1-PPP) penalises mild overdispersion
    # (PPP_var 0.7-0.9) which is common and expected for stable loci in
    # genome-wide catalogs — the prior prevents r from going to zero so
    # the model is always slightly wider than a perfectly stable locus.
    # This asymmetric formula treats mild overdispersion as acceptable
    # while still penalising severe overdispersion and all underdispersion.
    #
    # Underdispersion (PPP_var 0 → 0.5): linear 0 → 1.0
    # Mild overdispersion (PPP_var 0.5 → 0.90): stays at 1.0
    # Severe overdispersion (PPP_var 0.90 → 1.0): drops 1.0 → 0.5
    if ppp_var <= 0.5:
        fq_base = 2.0 * ppp_var                              # 0 at PPP=0, 1 at PPP=0.5
    elif ppp_var <= 0.90:
        fq_base = 1.0                                        # acceptable overdispersion
    else:
        fq_base = 1.0 - 0.5 * (ppp_var - 0.90) / 0.10      # gentle penalty 1.0→0.5

    # AD penalty: tail misfit degrades score
    # AD < 5: no penalty.  AD 5-100: linear to 0.5×.  AD > 100: severe.
    if ad_stat <= 5.0:
        ad_factor = 1.0
    elif ad_stat <= 100.0:
        ad_factor = 1.0 - 0.5 * (ad_stat - 5.0) / 95.0
    else:
        ad_factor = max(0.5 * float(np.exp(-0.01 * (ad_stat - 100.0))), 0.0)

    fit_quality = float(np.clip(fq_base * ad_factor, 0.0, 1.0))

    # Chi-square adjusted score: only penalise when underdispersed and
    # chi-square is reliable (df >= 2, bins >= 3).
    # Overdispersed chi-square misfit is a calibration issue.
    if (not np.isnan(chisq_pval) and
            chisq_df >= 2 and chisq_bins >= 3 and
            chisq_pval < 0.05 and ppp_var < 0.5):
        fq_adj = fit_quality * min(1.0, chisq_pval / 0.05)
    else:
        fq_adj = fit_quality

    return {
        "ks_statistic":         round(ks_stat,   6),
        "ad_statistic":         round(ad_stat,   4),
        "ppp_variance":         round(ppp_var,   4),
        "ppp_skewness":         round(ppp_skew,  4),
        "observed_var":         round(obs_var,   4),
        "fitted_var":           round(fitted_var, 4),
        "dispersion_ratio":     round(disp_ratio, 4),
        "chisq_statistic":      round(chisq_stat, 4) if not np.isnan(chisq_stat) else float('nan'),
        "chisq_pvalue":         round(chisq_pval, 6) if not np.isnan(chisq_pval) else float('nan'),
        "chisq_df":             chisq_df,
        "chisq_n_bins":         chisq_bins,
        "fit_quality":          round(fit_quality, 4),
        "fit_quality_adjusted": round(fq_adj,       4),
    }


def analyse_genome_wide(
    data: dict,
    global_params: Regime1GlobalParams,
    founder_lengths: Optional[dict] = None,
    mitotic_generations: Optional[float] = None,
    n_ppp: int = 300,
    tail_ratio_threshold: float = 5.0,
    min_reads: int = 5,
    seed: int = 0,
    verbose: bool = True,
) -> list[LocusResult]:
    """
    Run Regime 1 instability analysis across a genome-wide catalog of loci
    for a single sample.

    Parameters
    ----------
    data : dict
        Input data in the format:
            {
                "locus_id": {
                    "H1": [length1, length2, ...],
                    "H2": [length1, length2, ...],
                },
                ...
            }
        Haplotype keys can be any strings (e.g. "H1"/"H2", "mat"/"pat",
        "hap1"/"hap2"). Float lengths are rounded to nearest integer.

    global_params : Regime1GlobalParams
        Global step-size parameters from Pass 1 (or from_defaults() /
        from_stratum() for a quick start).

    founder_lengths : dict, optional
        Germline allele lengths per locus/haplotype:
            {"locus_id": {"H1": 40.0, "H2": 55.0}, ...}
        If not supplied, the modal observed length is used per haplotype.
        Supplying known germline lengths (e.g. from blood or parental data)
        gives more accurate instability estimates.

    mitotic_generations : float, optional
        G — number of mitotic cell divisions from zygote to sampled tissue.
        If supplied, per-division mutation rates are computed in addition to
        per-lineage rates. Typical values:
            Blood (adult):    50-70
            Colon epithelium: 1000-2000
            Neurons (adult):  ~0 (post-mitotic)
        If None, only per-lineage rates are reported.

    n_ppp : int
        Simulations per locus for the posterior predictive p-value.
        Default 300: fast (~10ms overhead per locus) and sufficient to
        distinguish ppp_variance=0 (complete failure) from ppp≈0.5 (good fit).
        Increase to 1000+ for publication-quality GoF estimates.

    tail_ratio_threshold : float
        max|Δ|/std(Δ) above which a haplotype is flagged as heavy-tailed.
        Default 5.0.

    min_reads : int
        Minimum reads per haplotype to attempt fitting. Haplotypes with
        fewer reads get an error entry rather than a failed fit.

    seed : int
        Global RNG seed. Per-locus seeds are derived from this to ensure
        reproducibility without correlation between loci.

    verbose : bool
        Print progress (locus count).

    Returns
    -------
    list[LocusResult]
        One LocusResult per locus. Call .to_rows() on each to get TSV-ready
        dicts, or use write_tsv() to write all results directly.
    """
    results = []
    n_loci  = len(data)

    for i, (locus_id, haplotypes) in enumerate(data.items()):
        t0     = time.time()
        locus_seed = seed + i * 1000   # unique per-locus seed

        h1_result = h2_result = None
        h1_gof    = h2_gof    = None
        h1_diag   = h2_diag   = None
        h1_error  = h2_error  = None

        hap_keys = list(haplotypes.keys())   # e.g. ["H1", "H2"]

        for j, hap_key in enumerate(hap_keys[:2]):   # at most 2 haplotypes
            lengths_raw = np.asarray(haplotypes[hap_key])
            hap_seed    = locus_seed + j * 100

            # Look up founder length if provided
            fl = None
            if founder_lengths and locus_id in founder_lengths:
                fl = founder_lengths[locus_id].get(hap_key)

            # Pre-fit diagnostics
            try:
                diag = diagnose_distribution(
                    lengths_raw,
                    founder_length=fl if fl is not None else float(
                        np.median(lengths_raw)   # use median as placeholder
                    ),
                    tail_ratio_threshold=tail_ratio_threshold,
                )
            except Exception as e:
                diag = None

            # Check minimum reads
            if len(lengths_raw) < min_reads:
                err = (f"Insufficient reads: n={len(lengths_raw)} < "
                       f"min_reads={min_reads}")
                if j == 0: h1_error = err; h1_diag = diag
                else:       h2_error = err; h2_diag = diag
                continue

            # Fit Regime 1
            try:
                result = fit_locus_r1(
                    lengths_raw,
                    founder_length=fl,   # None → uses modal length
                    global_params=global_params,
                    haplotype_label=hap_key,
                    mitotic_generations=mitotic_generations,
                    seed=hap_seed,
                )

                # Update diag with the actual founder length used
                if diag is not None:
                    diag = diagnose_distribution(
                        lengths_raw,
                        founder_length=result.founder_length,
                        tail_ratio_threshold=tail_ratio_threshold,
                    )

                # Goodness of fit at posterior mean parameters
                lengths_int = np.round(lengths_raw).astype(int)
                deltas      = (lengths_int - int(round(result.founder_length))).astype(int)
                gof = _compute_gof(
                    deltas,
                    r_e=result.r_e[0],
                    r_c=result.r_c[0],
                    pp=result.p_plus[0],
                    pm=result.p_minus[0],
                    qe=global_params.q_e,
                    qc=global_params.q_c,
                    n_ppp=n_ppp,
                    seed=hap_seed + 50,
                )

                if j == 0: h1_result=result; h1_gof=gof; h1_diag=diag
                else:       h2_result=result; h2_gof=gof; h2_diag=diag

            except Exception as e:
                err = str(e)
                if j == 0: h1_error=err; h1_diag=diag
                else:       h2_error=err; h2_diag=diag

        results.append(LocusResult(
            locus_id  = locus_id,
            h1        = h1_result,
            h2        = h2_result,
            h1_gof    = h1_gof,
            h2_gof    = h2_gof,
            h1_diag   = h1_diag,
            h2_diag   = h2_diag,
            h1_error  = h1_error,
            h2_error  = h2_error,
            elapsed_s = time.time() - t0,
        ))

        if verbose and (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{n_loci} loci...")

    if verbose:
        print(f"Done. {n_loci} loci processed.")

    return results


def write_tsv(
    results: list[LocusResult],
    output_path: str,
    include_per_div: bool = False,
) -> None:
    """
    Write genome-wide results to a TSV file.

    Parameters
    ----------
    results : list[LocusResult]
        Output of analyse_genome_wide().
    output_path : str
        Path to write the TSV file.
    include_per_div : bool
        If True, include per-division mutation rate columns in output.
        Only meaningful when mitotic_generations was supplied.

    Output columns (in order)
    -------------------------
    locus_id, haplotype, n_reads, founder_length,
    -- NegBin parameters --
    r_e, r_e_ci_lo, r_e_ci_hi,
    r_c, r_c_ci_lo, r_c_ci_hi,
    p_plus, p_plus_ci_lo, p_plus_ci_hi,
    p_minus, p_minus_ci_lo, p_minus_ci_hi,
    -- Instability --
    instability_index, instability_index_ci_lo, instability_index_ci_hi,
    net_bias, net_bias_ci_lo, net_bias_ci_hi,
    p_expansion, p_expansion_ci_lo, p_expansion_ci_hi,
    -- Mutation rates --
    expansion_rate, expansion_rate_ci_lo, expansion_rate_ci_hi,
    contraction_rate, contraction_rate_ci_lo, contraction_rate_ci_hi,
    net_mutation_rate, net_mutation_rate_ci_lo, net_mutation_rate_ci_hi,
    p_any_mutation, p_any_mutation_ci_lo, p_any_mutation_ci_hi,
    -- Goodness of fit --
    gof_ks_statistic, gof_ad_statistic,
    gof_ppp_variance, gof_ppp_skewness,
    gof_dispersion_ratio,
    -- Diagnostics --
    diag_tail_ratio, diag_skewness, diag_n_outliers, diag_is_heavy_tailed,
    -- Metadata --
    elapsed_s, error
    """
    import csv

    all_rows = []
    for lr in results:
        for row in lr.to_rows():
            # Remove per-division columns if not requested
            if not include_per_div:
                for col in ["expansion_rate_per_div", "expansion_rate_per_div_ci_lo",
                            "expansion_rate_per_div_ci_hi", "contraction_rate_per_div",
                            "contraction_rate_per_div_ci_lo","contraction_rate_per_div_ci_hi",
                            "net_rate_per_div","net_rate_per_div_ci_lo",
                            "net_rate_per_div_ci_hi","mitotic_generations"]:
                    row.pop(col, None)
            all_rows.append(row)

    if not all_rows:
        return

    # Canonical column order: put key columns first, GoF last
    priority = [
        "locus_id", "haplotype", "n_reads", "founder_length",
        "r_e", "r_e_ci_lo", "r_e_ci_hi",
        "r_c", "r_c_ci_lo", "r_c_ci_hi",
        "p_plus", "p_plus_ci_lo", "p_plus_ci_hi",
        "p_minus", "p_minus_ci_lo", "p_minus_ci_hi",
        "instability_index", "instability_index_ci_lo", "instability_index_ci_hi",
        "net_bias", "net_bias_ci_lo", "net_bias_ci_hi",
        "p_expansion", "p_expansion_ci_lo", "p_expansion_ci_hi",
        "expansion_rate", "expansion_rate_ci_lo", "expansion_rate_ci_hi",
        "contraction_rate","contraction_rate_ci_lo","contraction_rate_ci_hi",
        "net_mutation_rate","net_mutation_rate_ci_lo","net_mutation_rate_ci_hi",
        "p_any_mutation","p_any_mutation_ci_lo","p_any_mutation_ci_hi",
        "gof_ks_statistic","gof_ad_statistic",
        "gof_ppp_variance","gof_ppp_skewness",
        "gof_dispersion_ratio",
        "gof_chisq_statistic","gof_chisq_pvalue",
        "gof_chisq_df","gof_chisq_n_bins",
        "diag_tail_ratio","diag_skewness","diag_n_outliers","diag_is_heavy_tailed",
        "diag_recommended_regime",
        "elapsed_s","error",
    ]
    all_keys = set()
    for row in all_rows:
        all_keys.update(row.keys())
    ordered = [c for c in priority if c in all_keys]
    ordered += sorted(all_keys - set(ordered))   # any remaining columns

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ordered, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow({k: row.get(k, "") for k in ordered})

    print(f"Written {len(all_rows)} rows to {output_path}")


# ===========================================================================
# SINGLE HAPLOTYPE CONVENIENCE FUNCTION
# ===========================================================================

def analyse_haplotype(
    lengths,
    founder_length: float,
    global_params: Optional[Regime1GlobalParams] = None,
    repeat_unit_length: int = 3,
    stratum: Optional[str] = None,
    haplotype_label: str = "haplotype",
    mitotic_generations: Optional[float] = None,
    n_ppp: int = 300,
    tail_ratio_threshold: float = 5.0,
    seed: int = 0,
) -> dict:
    """
    Analyse one phased haplotype at one locus using Regime 1.

    The simplest entry point — takes a list of allele lengths and the
    founder length, returns a flat dict of all results including
    goodness-of-fit metrics.

    Parameters
    ----------
    lengths : list or array-like of float or int
        Per-read allele lengths for this haplotype, e.g. [40.3, 41.0, 41.7].
        Float values are rounded to the nearest integer repeat unit.
    founder_length : float
        Germline / founder allele length (L₀). This is the reference point
        from which all length changes (Δ) are measured. Use the known
        germline length from blood or parental data if available; otherwise
        the modal observed length is a reasonable approximation.
    global_params : Regime1GlobalParams, optional
        Pre-estimated global step-size parameters from Pass 1. If not
        supplied, defaults are constructed automatically from
        repeat_unit_length or stratum.
    repeat_unit_length : int
        Length of the repeat motif in base pairs (used only when
        global_params is None and stratum is None). Controls default
        step-size priors:
            1-3 bp → p₊ ≈ 0.60, p₋ ≈ 0.70  (small steps, STRs)
            4-6 bp → p₊ ≈ 0.45, p₋ ≈ 0.55
            7+ bp  → p₊ ≈ 0.30, p₋ ≈ 0.45  (larger steps, disease loci)
    stratum : str, optional
        Named repeat stratum (overrides repeat_unit_length when supplied).
        Options: 'str_short', 'str_medium', 'str_long',
                 'disease_cag', 'disease_cgg', 'disease_gaa', 'disease_ctg'.
    haplotype_label : str
        Label written into the output dict (default "haplotype").
    mitotic_generations : float, optional
        G — number of mitotic cell divisions from zygote to sampled tissue.
        When supplied, per-division mutation rates are included in the output.
        Typical values: blood (adult) 50-70, colon epithelium 1000-2000.
    n_ppp : int
        Simulations for the posterior predictive p-value (default 300).
        Set to 0 to skip goodness-of-fit computation entirely (faster).
    tail_ratio_threshold : float
        max|Δ|/std(Δ) above which a haplotype is flagged as heavy-tailed (default 5.0).
    seed : int
        RNG seed for reproducibility.

    Returns
    -------
    dict with the following keys:

    Identity
        haplotype, n_reads, founder_length

    NegBin parameters (posterior mean + 95% CI)
        r_e, r_e_ci_lo, r_e_ci_hi
        r_c, r_c_ci_lo, r_c_ci_hi
        p_plus, p_plus_ci_lo, p_plus_ci_hi
        p_minus, p_minus_ci_lo, p_minus_ci_hi
        mu_N_expansion  (mean expansion events per lineage)
        mu_N_contraction

    Instability
        instability_index, instability_index_ci_lo, instability_index_ci_hi
        net_bias, net_bias_ci_lo, net_bias_ci_hi
        p_expansion, p_expansion_ci_lo, p_expansion_ci_hi

    Mutation rates (per lineage)
        expansion_rate, expansion_rate_ci_lo, expansion_rate_ci_hi
        contraction_rate, contraction_rate_ci_lo, contraction_rate_ci_hi
        net_mutation_rate, net_mutation_rate_ci_lo, net_mutation_rate_ci_hi
        p_any_mutation, p_any_mutation_ci_lo, p_any_mutation_ci_hi

    Per-division rates (only when mitotic_generations is supplied)
        expansion_rate_per_div, expansion_rate_per_div_ci_lo/hi
        contraction_rate_per_div, contraction_rate_per_div_ci_lo/hi
        net_rate_per_div, net_rate_per_div_ci_lo/hi
        mitotic_generations

    Goodness of fit (only when n_ppp > 0)
        gof_ks_statistic      — KS distance between empirical and fitted CDF
        gof_ad_statistic      — Anderson-Darling statistic (tail-sensitive)
        gof_ppp_variance      — P(simulated var ≥ observed var); near 0 = bad fit
        gof_ppp_skewness      — P(simulated skew ≥ observed skew)
        gof_dispersion_ratio  — observed_var / fitted_var; near 1 = good fit

    Diagnostics
        diag_tail_ratio         — max|Δ|/std(Δ); > 5 suggests
        diag_skewness           — skewness of the Δ distribution
        diag_n_outliers         — reads flagged as potential outliers
        diag_is_heavy_tailed    — 1 if tail_ratio > threshold
        diag_warnings           — list of warning strings

    Timing
        elapsed_s

    Examples
    --------
    Basic usage:
        result = analyse_haplotype([40, 41, 41, 42, 40, 43], founder_length=40)
        print(result['instability_index'])

    With known germline length and tissue type:
        result = analyse_haplotype(
            lengths=[55, 56, 57, 58, 56, 55, 60, 54],
            founder_length=55.0,
            stratum='disease_cag',
            mitotic_generations=60.0,   # adult blood
        )

    Checking fit quality before trusting results:
        elif result['gof_ppp_variance'] > 0.1:
            print('Acceptable fit')
        else:
            print('Good fit')
    """
    t0 = time.time()

    # ── Build global params if not supplied ───────────────────────────────
    if global_params is None:
        if stratum is not None:
            global_params = Regime1GlobalParams.from_stratum(stratum)
        else:
            global_params = Regime1GlobalParams.from_defaults(repeat_unit_length)

    # ── Pre-fit diagnostics ───────────────────────────────────────────────
    diag = diagnose_distribution(
        lengths,
        founder_length=founder_length,
        tail_ratio_threshold=tail_ratio_threshold,
    )

    # ── Fit Regime 1 ──────────────────────────────────────────────────────
    result = fit_locus_r1(
        lengths,
        founder_length=founder_length,
        global_params=global_params,
        haplotype_label=haplotype_label,
        mitotic_generations=mitotic_generations,
        seed=seed,
    )

    # ── Goodness of fit ───────────────────────────────────────────────────
    gof = {}
    if n_ppp > 0:
        lengths_int = np.round(np.asarray(lengths)).astype(int)
        deltas      = (lengths_int - int(round(founder_length))).astype(int)
        gof = _compute_gof(
            deltas,
            r_e = result.r_e[0],
            r_c = result.r_c[0],
            pp  = result.p_plus[0],
            pm  = result.p_minus[0],
            qe  = global_params.q_e,
            qc  = global_params.q_c,
            n_ppp = n_ppp,
            seed  = seed + 1,
        )

    # ── Assemble flat output dict ─────────────────────────────────────────
    def _unpack(tup, name):
        """Unpack a (mean, lo, hi) tuple into three named keys."""
        return {
            name:             tup[0],
            name + "_ci_lo":  tup[1],
            name + "_ci_hi":  tup[2],
        }

    out = {
        "haplotype":      haplotype_label,
        "n_reads":        result.n_reads,
        "founder_length": result.founder_length,
    }

    # NegBin parameters
    out.update(_unpack(result.r_e,           "r_e"))
    out.update(_unpack(result.r_c,           "r_c"))
    out["q_e"] = global_params.q_e
    out["q_c"] = global_params.q_c
    out.update(_unpack(result.p_plus,        "p_plus"))
    out.update(_unpack(result.p_minus,       "p_minus"))
    out.update(_unpack(result.mu_expansion,  "mu_N_expansion"))
    out.update(_unpack(result.mu_contraction,"mu_N_contraction"))

    # Instability
    out.update(_unpack(result.instability_index, "instability_index"))
    out.update(_unpack(result.net_bias,          "net_bias"))
    out.update(_unpack(result.p_expansion,       "p_expansion"))

    # Mutation rates (per lineage)
    out.update(_unpack(result.expansion_rate,    "expansion_rate"))
    out.update(_unpack(result.contraction_rate,  "contraction_rate"))
    out.update(_unpack(result.net_mutation_rate, "net_mutation_rate"))
    out.update(_unpack(result.p_any_mutation,    "p_any_mutation"))

    # Per-division rates (only when G was supplied)
    if mitotic_generations is not None:
        out["mitotic_generations"] = mitotic_generations
        if result.expansion_rate_per_div is not None:
            out.update(_unpack(result.expansion_rate_per_div,   "expansion_rate_per_div"))
            out.update(_unpack(result.contraction_rate_per_div, "contraction_rate_per_div"))
            out.update(_unpack(result.net_rate_per_div,         "net_rate_per_div"))

    # Goodness of fit
    for k, v in gof.items():
        out[f"gof_{k}"] = v

    # Composite fit quality scores
    ppp_v = gof.get("ppp_variance", 0.5) if gof else 0.5
    chisq_p = gof.get("chisq_pvalue", float("nan")) if gof else float("nan")
    chisq_df_val = gof.get("chisq_df", 0) if gof else 0

    # Use the same asymmetric formula as _compute_gof:
    # mild overdispersion (PPP 0.5-0.9) is acceptable for stable loci
    if gof:
        out["fit_quality"]          = gof.get("fit_quality", 0.0)
        out["fit_quality_adjusted"] = gof.get("fit_quality_adjusted", 0.0)
    else:
        out["fit_quality"]          = 0.0
        out["fit_quality_adjusted"] = 0.0

    # Winsorized mean absolute delta (matches Handsaker et al. definition)
    # Directly from observed reads, not from the model
    abs_deltas = np.abs(
        np.round(np.asarray(lengths)).astype(int) - int(round(founder_length))
    ).astype(float)
    winsor_threshold = 100.0   # 100 repeat units, matching Handsaker et al. Fig. SN2.3
    out["winsorized_mean_abs_delta"] = round(
        float(np.minimum(abs_deltas, winsor_threshold).mean()), 4
    )
    out["winsor_threshold"] = winsor_threshold

    # Diagnostics
    out["diag_tail_ratio"]      = round(diag.tail_ratio, 3)
    out["diag_skewness"]        = round(diag.skewness, 3)
    out["diag_n_outliers"]      = diag.n_outliers
    out["diag_is_heavy_tailed"] = int(diag.is_heavy_tailed)
    out["diag_warnings"]        = diag.warnings   # list of strings

    out["elapsed_s"] = round(time.time() - t0, 3)

    return out