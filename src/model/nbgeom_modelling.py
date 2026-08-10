from __future__ import annotations

import time
import warnings
import logging

import numpy as np
from   scipy import stats
from   scipy.stats import skew as _skew, kurtosis as _kurt

from src.model.gof import compute_gof
from src.model.structs import NBGeomGlobalParams, NBGeomResult
from src.model.utils import get_stratum, STRATUM_LABELS

logging.getLogger("pymc").setLevel(logging.ERROR)


# Default grid for (r_e, r_c) — log-spaced from 0.01 to 12.0, 50 points each for
# for likelihood evaluation.
DEFAULT_RE_GRID  = np.exp(np.linspace(np.log(0.01), np.log(12.0), 50))  # log-spaced
DEFAULT_RC_GRID  = np.exp(np.linspace(np.log(0.01), np.log(12.0), 50))  # log-spaced
DEFAULT_PE_N     = 10     # number of p_e grid points around global estimate
DEFAULT_PC_N     = 10     # number of p_c grid points around global estimate
DEFAULT_P_LOG_SD = 0.2   # prior SD in log-odds space for p_e and p_c grids


def sigmoid_vec(x):
    """
    Numerically stable sigmoid for arrays.

    @param x: input array-like
    @return: array of sigmoid(x) values, clipped to avoid overflow
    """
    return 1 / (1 + np.exp(-np.clip(x, -30, 30)))


def summary_stats(deltas, L0):
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


def estimate_nbgeom_global(loci_data, min_reads=10):
    """
    Estimate global NBGeom parameters by pooling across loci.

    For each locus, the empirical mean and variance of Δ constrain p₊, p₋,
    and the ratio r_e/r_c. Pooling across many loci with different L₀ values
    allows regression of logit(p) on L₀ to estimate the length-dependence slopes.

    @param loci_data: list of dicts, each with keys:
        'deltas'         : int array  — Δ_i = observed - founder
        'founder_length' : float      — L₀
    @param min_reads: minimum number of reads per locus to include in regression
    @param n_pe_grid: number of p_e grid points around global estimate (for per-locus correction)
    @param n_pc_grid: number of p_c grid points around global estimate (for per-locus correction)
    @return: NBGeomGlobalParams with estimated lp0, lp1, lm0, lm1, q_e, q_c, n_loci_used

    For each locus: estimate p₊, p₋ from the method of moments using the
    asymmetry between positive and negative tails of the Δ distribution.
    Then regress logit(p₊) and logit(p₋) on L₀ across loci.
    """
    L0s      = []
    pe_ests  = []   # per-locus p_e estimates
    pc_ests  = []   # per-locus p_c estimates
    weights  = []

    for locus in loci_data:
        deltas = np.asarray(locus["deltas"], dtype=float)
        L0     = float(locus["founder_length"])
        n      = len(deltas)
        if n < min_reads: continue

        # Separate positive (expansion) and negative (contraction) deltas
        pos =  deltas[deltas > 0]
        neg = -deltas[deltas < 0]   # flip to positive

        # If either side has fewer than 3 data points, skip this locus
        if len(pos) < 3 or len(neg) < 3:
            continue

        # Method of moments: for Geometric(p), mean = 1/p
        # Use reciprocal of mean step size as estimate of p
        pe_est = np.clip(1.0 / max(pos.mean(), 0.5), 0.05, 0.95)
        pc_est = np.clip(1.0 / max(neg.mean(), 0.5), 0.05, 0.95)

        L0s.append(L0)
        pe_ests.append(pe_est)
        pc_ests.append(pc_est)
        weights.append(np.sqrt(n))

    if len(L0s) < 4:
        # Not enough loci — return default (no length dependence)
        return NBGeomGlobalParams(n_loci_used=len(L0s))

    L0s     = np.array(L0s)
    pe_ests = np.array(pe_ests)
    pc_ests = np.array(pc_ests)
    W       = np.array(weights)
    X       = np.column_stack([np.ones(len(L0s)), L0s])

    # logit function converting p in (0,1) to log-odds in (-inf, +inf)
    def logit(p): return np.log(p / (1 - p + 1e-10) + 1e-10)

    def wls(y):
        Xw = X * W[:, None]; yw = y * W
        return np.linalg.lstsq(Xw, yw, rcond=None)[0]

    coef_p = wls(logit(pe_ests))
    coef_m = wls(logit(pc_ests))

    return NBGeomGlobalParams(
        lp0=float(coef_p[0]), lp1=float(coef_p[1]),
        lm0=float(coef_m[0]), lm1=float(coef_m[1]),
        q_e=0.35, q_c=0.35,   # default — update from richer data if available
        n_loci_used=len(L0s),
    )


def nbgeom_grid_loglik(deltas, re_grid, rc_grid, pe, pc, qe, qc, N=None):
    """
    Vectorised log-likelihood for the NB-Geometric compound
    over a 2D grid of (r_e, r_c) with p₊, p₋, q_e, q_c fixed.

    @param deltas: array of observed Δ values (lengths - founder)
    @param re_grid: array of r_e values (expansion NegBin size)
    @param rc_grid: array of r_c values (contraction NegBin size)
    @param pe: p_e (Geometric success probability for expansion steps)
    @param pc: p_c (Geometric success probability for contraction steps)
    @param qe: q_e (NegBin success probability for expansion counts)
    @param qc: q_c (NegBin success probability for contraction counts)
    @param N: FFT length (must be large enough to avoid aliasing)

    @return: 2D array of log-likelihoods with shape (len(re_grid), len(rc_grid))
    """
    deltas = np.asarray(deltas, dtype=int)
    unique_d, counts = np.unique(deltas, return_counts=True)

    # FFT length: must span the support of Δ without aliasing
    d_abs = max(int(np.abs(deltas).max()), 1)
    if N is None:
        N = min(int(2 ** np.ceil(np.log2(d_abs * 6 + 16))), 4096)

    # the characteristic function is evaluated at N equally spaced points on the unit circle
    t   = 2 * np.pi * np.arange(N) / N

    ei_e = np.exp(-1j * t)   # e^{−it}: expansion direction
    ei_c = np.exp( 1j * t)   # e^{+it}: contraction direction

    # FFT indices for observed deltas (negative deltas wrap around)
    obs_idx = unique_d % N

    # Characteristic function evaluated at e^{−it} and e^{+it}
    M_sp = pe * ei_e / (1 - (1 - pe) * ei_e)   # (N,) expansion steps
    M_sm = pc * ei_c / (1 - (1 - pc) * ei_c)    # (N,) contraction steps

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


def nbgeom_grid_posterior(deltas, re_grid, rc_grid, pe_grid, pc_grid, pe_prior_mean, pc_prior_mean,
                          qe, qc, pe_log_sd = DEFAULT_P_LOG_SD, pc_log_sd = DEFAULT_P_LOG_SD):
    """
    4D grid posterior over (r_e, r_c, p_e, p_c) for NBGeom model.

    p_e and p_c float on a narrow grid around the global estimate, allowing per-locus
    data to correct for any bias in the global estimate.

    Priors:
        r_e, r_c  ~ Gamma(1.5, scale=1.0)    — moderate event rates
        logit(p_e) ~ Normal(logit(prior), σ²) — log-odds Normal prior
        logit(p_c) ~ Normal(logit(prior), σ²)

    @param deltas: array of observed Δ values (lengths - founder)
    @param re_grid: array of r_e values (expansion NegBin size)
    @param rc_grid: array of r_c values (contraction NegBin size)
    @param pe_grid: array of p_e values (Geometric success probability for expansion steps)
    @param pc_grid: array of p_c values (Geometric success probability for contraction steps)
    @param pe_prior_mean: global estimate of p_e (from Pass 1)
    @param pc_prior_mean: global estimate of p_c (from Pass 1)
    @param qe: q_e (NegBin success probability for expansion counts)
    @param qc: q_c (NegBin success probability for contraction counts)
    @param pe_log_sd: prior SD in log-odds space for p_e grid
    @param pc_log_sd: prior SD in log-odds space for p_c grid
    @return: 4D array of posterior probabilities with shape (len(re_grid), len(rc_grid), len(pe_grid), len(pc_grid))
    """
    n_re, n_rc, n_pe, n_pc = len(re_grid), len(rc_grid), len(pe_grid), len(pc_grid)
    loglik = np.full((n_re, n_rc, n_pe, n_pc), -np.inf)

    # Loop over p₊ values (outer), vectorise over (r_e, r_c, p₋) (inner)
    for ip, pe in enumerate(pe_grid):
        for im, pc in enumerate(pc_grid):
            loglik[:, :, ip, im] = nbgeom_grid_loglik(deltas, re_grid, rc_grid, pe, pc, qe, qc)

    # Log-Normal prior on logit(p₊) and logit(p₋) — keeps p values in (0,1)
    def logit(p): return np.log(p / (1 - p + 1e-10) + 1e-10)

    log_prior = (
        stats.gamma.logpdf(re_grid, a=1.5, scale=1.0)[:, None, None, None]
        + stats.gamma.logpdf(rc_grid, a=1.5, scale=1.0)[None, :, None, None]
        + stats.norm.logpdf(logit(pe_grid), logit(pe_prior_mean), pe_log_sd
                            )[None, None, :, None]
        + stats.norm.logpdf(logit(pc_grid), logit(pc_prior_mean), pc_log_sd
                            )[None, None, None, :]
    )

    log_post = loglik + log_prior
    log_post -= log_post.max()
    post = np.exp(log_post)
    post /= post.sum()
    return post


def nbgeom_posterior_summaries(post, re_grid, rc_grid, pe_grid, pc_grid, qe, qc):
    """
    Compute posterior mean and 95% CI for all NBGeom model quantities.

    @param post: 4D posterior array over (r_e, r_c, p₊, p₋)
    @param re_grid: array of r_e values (expansion NegBin size)
    @param rc_grid: array of r_c values (contraction NegBin size)
    @param pe_grid: array of p_e values (Geometric success probability for expansion steps)
    @param pc_grid: array of p_c values (Geometric success probability for contraction steps)
    @return: dict of {name: (mean, lo95, hi95)} for each quantity:
    
    
    Instability index (closed form):
        Var(Δ) = Var(E) + Var(C)
        Var(E) = μ_{N_e}·(1−p₊)/p₊² + σ²_{N_e}/p₊²
        where μ_{N_e} = r_e·(1−q_e)/q_e  and  σ²_{N_e} = r_e·(1−q_e)/q_e²
    """
    RE, RC, PE, PC = np.meshgrid(re_grid, rc_grid, pe_grid, pc_grid, indexing="ij")

    # NegBin moments for expansion and contraction count distributions
    mu_Ne  = RE * (1 - qe) / qe         # mean N_e (expansion events)
    var_Ne = RE * (1 - qe) / qe ** 2    # variance N_e
    mu_Nc  = RC * (1 - qc) / qc         # mean N_c (contraction events)
    var_Nc = RC * (1 - qc) / qc ** 2    # variance N_c

    # Mean total expansion E and contraction C (each is a sum of N Geometric steps)
    # E[E] = E[N_e] · E[S₊] = μ_{N_e} · (1−p₊)/p₊  (Geometric mean = (1-p)/p for support 0,1,2..)
    # Actually Geometric with support {1,2,...}: E[S] = 1/p
    # Here we use support {1,2,...} so E[S₊] = 1/p₊
    mu_E = mu_Ne / PE       # mean total expansion
    mu_C = mu_Nc / PC       # mean total contraction

    # Var(E) = E[N_e]·Var(S₊) + Var(N_e)·E[S₊]²
    # Var(S₊) = (1−p₊)/p₊²  (Geometric variance, support {1,2,...})
    var_E = mu_Ne * (1 - PE) / PE ** 2 + var_Ne / PE ** 2
    var_C = mu_Nc * (1 - PC) / PC ** 2 + var_Nc / PC ** 2

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
        "p_e":               PE,
        "p_c":               PC,
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
# Public API
# ===========================================================================

def fit_locus_r1(
    motif_length, lengths,
    founder_length: float,
    global_params: NBGeomGlobalParams,
    haplotype_label: str = "haplotype",
    n_pe: int = DEFAULT_PE_N,
    n_pc: int = DEFAULT_PC_N,
    mitotic_generations: Optional[float] = None,
    seed: int = 0,
) -> NBGeomResult:
    """
    Fit NB-Geometric compound model for one phased haplotype.

    @param lengths: array-like of allele lengths for this haplotype
    @param founder_length: germline / founder allele length L₀
    @param global_params: NBGeomGlobalParams with p₊(L₀), p₋(L₀), q_e, q_c from Pass 1
    @param haplotype_label: label for output
    @param n_pe, n_pc: number of p₊, p₋ grid points around global estimate
    @param mitotic_generations: number of mitotic cell divisions from zygote to sampled tissue (G).
        If supplied, per-division mutation rates are computed in addition to per-lineage rates.
        Typical values:
            Blood (adult):         50-70
            Colon epithelium:      1000-2000
            Neurons (adult):       ~0  (post-mitotic)
            Fibroblasts (adult):   ~40-60
        If None, only per-lineage rates are reported.
    @param seed: random seed for reproducibility (not used in current implementation)
    @return: NBGeomResult with posterior summaries for this haplotype
    """

    t0 = time.time()

    lengths     = np.asarray(lengths)
    lengths_int = np.round(lengths).astype(int)
    mode_int    = int(round(founder_length))
    deltas      = (lengths_int - mode_int).astype(int)

    re_grid = DEFAULT_RE_GRID
    rc_grid = DEFAULT_RC_GRID

    # Global p₊ and p₋ at this L₀
    pe_prior = global_params.p_e(founder_length)
    pc_prior = global_params.p_c(founder_length)
    qe       = global_params.q_e
    qc       = global_params.q_c

    # Build narrow grids around global estimates in logit space
    def logit(p): return np.log(p / (1 - p + 1e-10) + 1e-10)
    logit_pe_centre = logit(pe_prior)
    logit_pc_centre = logit(pc_prior)
    pe_grid = sigmoid_vec(np.linspace(
        logit_pe_centre - 2 * DEFAULT_P_LOG_SD,
        logit_pe_centre + 2 * DEFAULT_P_LOG_SD, n_pe))
    pc_grid = sigmoid_vec(np.linspace(
        logit_pc_centre - 2 * DEFAULT_P_LOG_SD,
        logit_pc_centre + 2 * DEFAULT_P_LOG_SD, n_pc))

    post = nbgeom_grid_posterior(deltas, re_grid, rc_grid, pe_grid, pc_grid,
                                 pe_prior_mean=pe_prior, pc_prior_mean=pc_prior,
                                 qe=qe, qc=qc)
    s = nbgeom_posterior_summaries(post, re_grid, rc_grid, pe_grid, pc_grid, qe, qc)

    # Compute per-division rates if G supplied
    def _scale(tup, G):
        if G is None or G <= 0:
            return None
        return (tup[0]/G, tup[1]/G, tup[2]/G)

    return NBGeomResult(
        haplotype_label          = haplotype_label,
        n_reads                  = len(lengths_int),
        founder_length           = founder_length,
        p_e_prior                = pe_prior,
        p_c_prior                = pc_prior,
        q_e                      = qe,
        q_c                      = qc,
        r_e                      = s["r_e"],
        r_c                      = s["r_c"],
        p_e                      = s["p_e"],
        p_c                      = s["p_c"],
        mu_expansion             = s["mu_expansion"],
        mu_contraction           = s["mu_contraction"],
        instability_index        = s["instability_index"],
        net_bias                 = s["net_bias"],
        p_expansion              = s["p_expansion"],
        expansion_rate           = s["expansion_rate"],
        contraction_rate         = s["contraction_rate"],
        net_mutation_rate        = s["net_mutation_rate"],
        p_any_mutation           = s["p_any_mutation"],
        expansion_rate_per_div   = _scale(s["expansion_rate"],    mitotic_generations),
        contraction_rate_per_div = _scale(s["contraction_rate"],  mitotic_generations),
        net_rate_per_div         = _scale(s["net_mutation_rate"], mitotic_generations),
        mitotic_generations      = mitotic_generations,
        elapsed_s                = time.time() - t0,
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


def auto_fit(motif_length, lengths, founder_length, nbgeom_params, tail_ratio_threshold = 5.0,
             outlier_mad_factor = 5.0, filter_outliers = False, mitotic_generations = None, seed = 0,
             **kwargs):
    """
    Automatically diagnose the read-length distribution and route to the
    appropriate regime.

    This is the recommended entry point for genome-wide use. It:
    1. Runs diagnose_distribution() to assess the data
    2. If heavy-tailed (tail_ratio > threshold)
    3. Otherwise routes to fit_locus_r1(), optionally filtering outliers first
    4. Returns (result, diagnostics) so the routing decision is always visible

    @param lengths: array-like of int or float, per-read allele lengths
    @param founder_length: float, germline allele length L₀
    @param nbgeom_params: NBGeomGlobalParams, required global parameters
    @param haplotype_label: str, label for output
    @param outlier_mad_factor: float, MAD factor for outlier flagging (default 5.0)
    @param tail_ratio_threshold: float, threshold for determining heavy-tailed distributions (default 5.0)
    @param filter_outliers: bool, if True and running Regime 1, remove outlier reads before fitting
    @param mitotic_generations: float, optional, G for per-division mutation rates (Regime 1 only)
    @param seed: int, random seed for reproducibility
    @param **kwargs: additional arguments forwarded to fit_locus_r1

    @returns: tuple of (result, diagnostics, gof)
    - result: NBGeomResult
    - diagnostics: DistributionDiagnostics object with pre-fit diagnostics
    - gof: dict of goodness-of-fit metrics for NBGeom model

    Returns
    -------
    (result, diagnostics, gof) :
        result      : NBGeomResult
        diagnostics : DistributionDiagnostics
        gof         : dict or None
            Goodness-of-fit metrics for NBGeom model (same keys as
            _compute_gof: ks_statistic, ad_statistic, ppp_variance,
            ppp_skewness, observed_var, fitted_var, dispersion_ratio,
            chisq_statistic, chisq_pvalue, chisq_df, chisq_n_bins, fit_quality,
            fit_quality_adjusted).

    Note: n_ppp (number of PPP simulations, default 300) can be passed
    via **kwargs to control GoF computation speed.
    """

    lengths = np.asarray(lengths)
    gof     = None   # populated only for Regime 1 fits

    # Extract n_ppp before forwarding kwargs to fit_locus_r1
    # so it does not cause an unexpected keyword argument error
    n_ppp = kwargs.pop("n_ppp", 300)

    # Optionally filter outliers before Regime 1
    if filter_outliers:
        lengths, removed = filter_outliers(lengths, founder_length, mad_factor=outlier_mad_factor)

    result = fit_locus_r1(motif_length, lengths, founder_length, nbgeom_params,
                          mitotic_generations=mitotic_generations, seed=seed, **kwargs)

    # Goodness-of-fit
    if n_ppp == 0:
        # n_ppp=0 means skip GoF entirely (useful for speed when not needed)
        gof = None
    else:
        try:
            lengths_int = np.round(lengths).astype(int)
            deltas      = (lengths_int - int(round(result.founder_length))).astype(int)
            gof = compute_gof(
                deltas,
                r_e = result.r_e[0],
                r_c = result.r_c[0],
                pe  = result.p_e[0],
                pc  = result.p_c[0],
                qe  = nbgeom_params.q_e,
                qc  = nbgeom_params.q_c,
                n_ppp = n_ppp,
                seed  = seed + 1,
            )
        except Exception as e:
            logging.warning(f"GoF computation failed: {e}")
            gof = None
            # Ignore GoF computation failures

    return result, gof
