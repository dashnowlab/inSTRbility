import numpy as np
from   scipy import stats
import warnings


def compute_gof(
    deltas: np.ndarray,
    r_e:    float,
    r_c:    float,
    pe:     float,
    pc:     float,
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
    M_sp = pe * emi / (1 - (1 - pe) * emi)
    M_sm = pc * ei  / (1 - (1 - pc) * ei)
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
