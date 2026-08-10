from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.optimize import minimize


def _expit(x: np.ndarray | float) -> np.ndarray | float:
    x = np.asarray(x)
    return 1.0 / (1.0 + np.exp(-np.clip(x, -60, 60)))


def _logit(p: np.ndarray | float) -> np.ndarray | float:
    p = np.asarray(p)
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return np.log(p / (1.0 - p))


def _next_pow2(n: int) -> int:
    n = max(1, int(n))
    return 1 if n == 1 else 2 ** int(np.ceil(np.log2(n)))


def simulate_compound_nb_geom(
    n: int,
    r_e: float,
    r_c: float,
    q_e: float,
    q_c: float,
    p_plus: float,
    p_minus: float,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """
    Simulate Delta = total_expansion - total_contraction.

    NB parameterization:
        P(N=k) = C(k+r-1, k) * (1-q)^k * q^r,  k = 0,1,2,...
    Geometric parameterization:
        S ~ Geom(p) on {1,2,3,...}, mean = 1/p
    """
    if rng is None:
        rng = np.random.default_rng()

    out = np.empty(n, dtype=int)
    for i in range(n):
        ne = rng.negative_binomial(r_e, q_e)
        nc = rng.negative_binomial(r_c, q_c)

        e = rng.geometric(p_plus, size=ne).sum() if ne > 0 else 0
        c = rng.geometric(p_minus, size=nc).sum() if nc > 0 else 0
        out[i] = e - c

    return out


def _delta_pmf_fft(
    r_e: float,
    r_c: float,
    q_e: float,
    q_c: float,
    p_plus: float,
    p_minus: float,
    n_fft: int,
) -> np.ndarray:
    """
    Return PMF of Delta on the FFT grid:
        pmf[k] corresponds to Delta = k if 0 <= k < N//2,
        and Delta = k-N if N//2 <= k < N.
    """
    N = int(n_fft)
    t = 2.0 * np.pi * np.arange(N) / N

    z_minus = np.exp(-1j * t)  # expansion side
    z_plus = np.exp(1j * t)    # contraction side

    # Geometric PGF, support {1,2,...}
    m_plus = p_plus * z_minus / (1.0 - (1.0 - p_plus) * z_minus)
    m_minus = p_minus * z_plus / (1.0 - (1.0 - p_minus) * z_plus)

    # NB PGF composed with geometric PGF
    g_e = q_e / (1.0 - (1.0 - q_e) * m_plus)
    g_c = q_c / (1.0 - (1.0 - q_c) * m_minus)

    cf = np.exp(r_e * np.log(np.clip(g_e, 1e-300, None))) * np.exp(
        r_c * np.log(np.clip(g_c, 1e-300, None))
    )

    pmf = np.real(np.fft.ifft(cf))
    pmf = np.clip(pmf, 0.0, None)

    s = pmf.sum()
    if not np.isfinite(s) or s <= 0:
        raise FloatingPointError("FFT PMF became invalid.")
    pmf /= s
    return pmf


def exact_delta_pmf(
    r_e: float,
    r_c: float,
    q_e: float,
    q_c: float,
    p_plus: float,
    p_minus: float,
    lo: int,
    hi: int,
    pad: int = 8,
    n_fft: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Exact PMF on integer grid [lo, hi] using the FFT representation.
    """
    span = max(abs(lo), abs(hi)) + pad
    if n_fft is None:
        # generous support to reduce circular aliasing
        n_fft = _next_pow2(max(256, 8 * span + 128))
        n_fft = min(n_fft, 65536)

    pmf = _delta_pmf_fft(r_e, r_c, q_e, q_c, p_plus, p_minus, n_fft=n_fft)

    grid = np.arange(n_fft)
    grid[grid >= n_fft // 2] -= n_fft

    order = np.argsort(grid)
    grid = grid[order]
    pmf = pmf[order]

    mask = (grid >= lo) & (grid <= hi)
    grid = grid[mask]
    pmf = pmf[mask]
    pmf = pmf / pmf.sum()

    return grid, pmf


@dataclass
class CompoundNBGeomFit:
    r_e: float
    r_c: float
    q_e: float
    q_c: float
    p_plus: float
    p_minus: float
    nll: float
    success: bool
    message: str
    n_fft: int


def fit_compound_nb_geom(
    deltas: np.ndarray,
    *,
    init_r_e: float = 1.0,
    init_r_c: float = 1.0,
    init_q_e: float = 0.5,
    init_q_c: float = 0.5,
    init_p_plus: float = 0.5,
    init_p_minus: float = 0.5,
    pad: int = 8,
    n_fft: Optional[int] = None,
    n_starts: int = 5,
    seed: int = 1,
) -> CompoundNBGeomFit:
    """
    Fit the full compound NB-Geometric model by maximum likelihood.

    Parameters are all free:
        r_e, r_c > 0
        q_e, q_c in (0,1)
        p_plus, p_minus in (0,1)

    Notes:
      - This uses exact FFT-based likelihood.
      - For best stability, deltas should be integer-valued.
      - If your support is wide, increase pad or n_fft.
    """
    x = np.asarray(np.rint(deltas), dtype=int)
    if x.size == 0:
        raise ValueError("deltas is empty")

    lo = int(x.min()) - pad
    hi = int(x.max()) + pad

    span = max(abs(lo), abs(hi))
    if n_fft is None:
        n_fft = _next_pow2(max(256, 8 * span + 128))
        n_fft = min(n_fft, 65536)

    # Precompute unique counts for speed
    vals, counts = np.unique(x, return_counts=True)
    idx = vals % n_fft

    def nll_from_theta(theta: np.ndarray) -> float:
        log_r_e, log_r_c, logit_q_e, logit_q_c, logit_p_plus, logit_p_minus = theta

        r_e = float(np.exp(log_r_e))
        r_c = float(np.exp(log_r_c))
        q_e = float(_expit(logit_q_e))
        q_c = float(_expit(logit_q_c))
        p_plus = float(_expit(logit_p_plus))
        p_minus = float(_expit(logit_p_minus))

        try:
            pmf = _delta_pmf_fft(r_e, r_c, q_e, q_c, p_plus, p_minus, n_fft=n_fft)
        except FloatingPointError:
            return np.inf

        prob = np.clip(pmf[idx], 1e-300, 1.0)
        return float(-(counts * np.log(prob)).sum())

    base = np.array(
        [
            np.log(max(init_r_e, 1e-8)),
            np.log(max(init_r_c, 1e-8)),
            _logit(init_q_e),
            _logit(init_q_c),
            _logit(init_p_plus),
            _logit(init_p_minus),
        ],
        dtype=float,
    )

    rng = np.random.default_rng(seed)
    starts = [base]
    for _ in range(max(0, n_starts - 1)):
        jitter = rng.normal(
            loc=0.0,
            scale=np.array([0.5, 0.5, 0.7, 0.7, 0.5, 0.5]),
        )
        starts.append(base + jitter)

    best = None
    for start in starts:
        res = minimize(
            nll_from_theta,
            start,
            method="L-BFGS-B",
            options={"maxiter": 500, "ftol": 1e-9},
        )
        if best is None or res.fun < best.fun:
            best = res

    theta = best.x
    return CompoundNBGeomFit(
        r_e=float(np.exp(theta[0])),
        r_c=float(np.exp(theta[1])),
        q_e=float(_expit(theta[2])),
        q_c=float(_expit(theta[3])),
        p_plus=float(_expit(theta[4])),
        p_minus=float(_expit(theta[5])),
        nll=float(best.fun),
        success=bool(best.success),
        message=str(best.message),
        n_fft=int(n_fft),
    )


# Example:
if __name__ == "__main__":
    rng = np.random.default_rng(0)

    true_params = dict(
        r_e=1.8,
        r_c=0.9,
        q_e=0.55,
        q_c=0.65,
        p_plus=0.35,
        p_minus=0.45,
    )

    deltas = simulate_compound_nb_geom(30, rng=rng, **true_params)

    print(len(deltas))

    fit = fit_compound_nb_geom(
        deltas,
        init_r_e=1.0,
        init_r_c=1.0,
        init_q_e=0.55,
        init_q_c=0.65,
        init_p_plus=0.35,
        init_p_minus=0.45,
        n_starts=6,
    )

    print(fit)