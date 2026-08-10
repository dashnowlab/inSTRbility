from src.model.structs import LocusResult


OUT_COLS = [
    # Locus identity
    "chrom", "start", "end", "motif", "haplotype",
    # Allele summary
    "founder_length", "founder_length_bp", "n_reads",
    # NB-Geometric posteriors (posterior mean only)
    "r_e", "r_c", "q_e", "q_c", "p_e", "p_c",
    "mu_expansion", "mu_contraction",
    # Instability indices
    "instability_index", "net_bias", "p_expansion",
    "expansion_rate", "contraction_rate",
    "net_mutation_rate", "p_any_mutation",
    # Regime 1 GoF
    "gof_ks_stat", "gof_ad_stat", "gof_ppp_var",
    "gof_ppp_skewness", "gof_observed_var", "gof_fitted_var",
    "gof_dispersion_ratio", "gof_chisq_stat", "gof_chisq_p",
    "gof_chisq_df", "gof_chisq_n_bins", "gof_fit_quality", "gof_fit_quality_adj",
    "elapsed_s",
]


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
