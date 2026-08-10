import numpy as np
import math

from dataclasses import dataclass
from typing import Optional

# ===========================================================================
# Global parameter dataclasses
# ===========================================================================


def sigmoid_scalar(x):
    """
    Numerically stable sigmoid for scalar values.
    
    @param x: input float value
    @return: sigmoid(x) = 1 / (1 + exp(-x)), clipped to avoid overflow
    """
    return float(1 / (1 + np.exp(-np.clip(x, -30, 30))))


@dataclass
class NBGeomGlobalParams:
    """
    Global parameters for Regime 1 (NB-Geometric compound).

    Step-size parameters estimated from pooled high-depth reads across loci.
    These are locus-biology parameters (how repeat slippage mechanics work)
    not sample-level parameters, so they are shared across samples.

    p_e(L₀): geometric parameter for expansion steps.
             Smaller p_e → larger mean expansion step (mean = 1/p_e).
             Can vary with L₀ via log-odds linear model:
             logit(p_e(L₀)) = lp0 + lp1·L₀

    p_c(L₀): geometric parameter for contraction steps.
             Similarly parameterised.

    q_e, q_c: NegBin probability for the expansion and contraction
              count distributions. Estimated from the shape of the
              count distribution across cells.
    """
    # Log-odds linear model for p₊: logit(p₊) = lp0 + lp1·L₀
    lp0:  float = 0.0    # logit(p₊) intercept
    lp1:  float = 0.0    # logit(p₊) slope (negative → bigger steps at longer alleles)
    # Log-odds linear model for p₋: logit(p₋) = lm0 + lm1·L₀
    lm0:  float = 0.0
    lm1:  float = 0.0
    # NegBin probability for count distributions
    q_e:  float = 0.5    # expansion count NB probability
    q_c:  float = 0.5    # contraction count NB probability
    # Diagnostics
    n_loci_used: int = 0


    def p_e(self, L0):
        """Expansion step-size geometric parameter at founder length L0."""
        return float(sigmoid_scalar(self.lp0 + self.lp1 * L0))


    def p_c(self, L0):
        """Contraction step-size geometric parameter at founder length L0."""
        return float(sigmoid_scalar(self.lm0 + self.lm1 * L0))


    @classmethod
    def from_defaults(cls, motif_length=3, q_e=0.5, q_c=0.5):
        """
        Biologically motivated defaults WITHOUT running Pass 1.

        Use when you want to run per-locus inference before estimating global
        params, or as a quick sanity-check starting point. The narrow p₊/p₋
        floating grid in Pass 2 corrects for any bias here, so mild
        initialisation errors do not badly affect per-locus posteriors.

        motif_length : int
            Length of the repeat unit in base pairs.
              1-3 bp (STR):           p₊ ≈ 0.60, p₋ ≈ 0.70  (mean steps ~1.7 / ~1.4)
              4-6 bp (medium STR):    p₊ ≈ 0.45, p₋ ≈ 0.55  (mean steps ~2.2 / ~1.8)
              7+ bp (long/disease):   p₊ ≈ 0.30, p₋ ≈ 0.45  (mean steps ~3.3 / ~2.2)
        q_e, q_c : float
            NegBin probability for count distributions (default 0.5).
        """
        def logit(p): return math.log(p / (1.0 - p))

        if motif_length <= 3:
            pe, pc = 0.60, 0.70

        elif motif_length <= 6:
            pe, pc = 0.45, 0.55

        else:
            pe, pc = 0.30, 0.45

        return cls(lp0=logit(pe), lp1=0.0, lm0=logit(pc), lm1=0.0, q_e=q_e, q_c=q_c, n_loci_used=0)


    @classmethod
    def from_stratum(cls, stratum, q_e=0.5, q_c=0.5):
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
        pe, pc = STRATA[stratum]
        return cls(lp0=logit(pe), lp1=0.0, lm0=logit(pc), lm1=0.0,
                   q_e=q_e, q_c=q_c, n_loci_used=0)


    @classmethod
    def calibrate_from_results(cls, results_tsv, base_params, min_reads=20, target_percentile=50):
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

        @param results_tsv      Path to TSV file written by write_tsv().
        @param base_params      NBGeomGlobalParams used to generate the results.
        @param min_reads        Minimum reads per locus to include in calibration.
        @param target_percentile Percentile of dispersion_ratio to target (default 50).

        @return NBGeomGlobalParams with adjusted q_e and q_c.
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
            f"NBGeomGlobalParams(\n"
            f"  p₊(L0) = sigmoid({self.lp0:.4f} + {self.lp1:.5f}·L0)\n"
            f"  p₋(L0) = sigmoid({self.lm0:.4f} + {self.lm1:.5f}·L0)\n"
            f"  q_e={self.q_e:.3f}, q_c={self.q_c:.3f}\n"
            f"  n_loci={self.n_loci_used}\n)"
        )



@dataclass
class NBGeomResult:
    """
    Posterior summary for NB-Geometric compound) inference.
    All quantities as (posterior_mean, ci_lo_95, ci_hi_95).
    """
    haplotype_label: str
    n_reads:         int
    founder_length:  float
    p_e_prior:    float   # global p₊ at this L₀
    p_c_prior:   float   # global p₋ at this L₀
    q_e:             float
    q_c:             float

    # Per-locus posteriors
    r_e:               tuple[float, float, float]  # expansion NB dispersion
    r_c:               tuple[float, float, float]  # contraction NB dispersion
    p_e:            tuple[float, float, float]  # expansion step-size param
    p_c:           tuple[float, float, float]  # contraction step-size param
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
            f"NBGeomResult(haplotype='{self.haplotype_label}', n={self.n_reads}, L0={self.founder_length:.1f})",
            f"  p₊ prior={self.p_e_prior:.4f}  p₋ prior={self.p_c_prior:.4f}"
            f"  q_e={self.q_e:.3f}  q_c={self.q_c:.3f}",
            "",
            "  --- NegBin expansion/contraction parameters ---",
            f("r_e (exp dispersion)",     self.r_e),
            f("r_c (con dispersion)",     self.r_c),
            f("p₊ (exp step size)",       self.p_e),
            f("p₋ (con step size)",       self.p_c),
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
             "p_e_prior": self.p_e_prior, "p_c_prior": self.p_c_prior,
             "mitotic_generations": self.mitotic_generations}
        per_lineage = ["r_e","r_c","p_e","p_c","mu_expansion","mu_contraction",
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



@dataclass
class LocusResult:
    """
    Complete Regime 1 results for one locus (both haplotypes).
    Includes per-haplotype instability estimates and goodness-of-fit metrics.
    """
    locus_id:  str
    h1:        Optional[NBGeomResult]
    h2:        Optional[NBGeomResult]
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
