
"""
basket_cds.py

Correlated Nth-to-Default Basket CDS valuation and (optional) risk-capital loss analytics.

Core model:
- Each obligor i defaults with intensity (hazard) lambda_i(t), piecewise constant between quoted CDS maturities.
- Dependence across obligors is introduced via a copula applied to default times:
    * Gaussian copula
    * Student-t copula

This implementation is intentionally "production-shaped":
- deterministic discounting (flat or curve),
- robust CDS bootstrap via root finding,
- vectorized Monte Carlo,
- clean separation between calibration, simulation, and pricing.

Author: refactor from an R prototype into Python
License: MIT (recommended if you publish this)
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq
from scipy.stats import norm, t as student_t


Array = NDArray[np.float64]


def discount_factors_flat(r: float, times: Array) -> Array:
    """Flat continuously-compounded discount factors."""
    return np.exp(-r * times)


@dataclass(frozen=True)
class HazardCurve:
    """
    Piecewise-constant hazard curve.

    segments_end: increasing times (years) defining the segment end points.
    hazard: hazard rates (1/year) for each segment, same length as segments_end.
    """
    segments_end: Array
    hazard: Array

    def survival(self, times: Array) -> Array:
        """Survival S(t) = exp(-∫0^t λ(s) ds) under piecewise-constant λ."""
        times = np.asarray(times, dtype=float)
        H = cumulative_hazard(times, self.segments_end, self.hazard)
        return np.exp(-H)

    def inverse_default_time(self, u: Array) -> Array:
        """
        Map U~Unif(0,1) to default time τ by τ = H^{-1}(-ln U).
        """
        u = np.asarray(u, dtype=float)
        x = -np.log(np.clip(u, 1e-16, 1.0))
        return invert_cumulative_hazard(x, self.segments_end, self.hazard)


def cumulative_hazard(times: Array, seg_end: Array, hazard: Array) -> Array:
    """Compute H(t)=∫0^t λ(s) ds for piecewise-constant hazard."""
    times = np.asarray(times, dtype=float)
    seg_end = np.asarray(seg_end, dtype=float)
    hazard = np.asarray(hazard, dtype=float)

    # prepend segment start points
    seg_start = np.concatenate([[0.0], seg_end[:-1]])
    # cumulative hazard at segment ends
    seg_H_end = np.cumsum(hazard * (seg_end - seg_start))
    # for each time t, find segment index
    idx = np.searchsorted(seg_end, times, side="right")  # 0..nseg
    H = np.zeros_like(times, dtype=float)

    # t in segment k: H = H_end(k-1) + hazard[k]*(t - seg_start[k])
    in_range = idx > 0
    k = idx[in_range] - 1
    H_prev = np.where(k > 0, seg_H_end[k - 1], 0.0)
    H[in_range] = H_prev + hazard[k] * (times[in_range] - seg_start[k])

    # t beyond last segment: extend last hazard
    beyond = idx == len(seg_end)
    if np.any(beyond):
        H_last = seg_H_end[-1]
        H[beyond] = H_last + hazard[-1] * (times[beyond] - seg_end[-1])

    return H


def invert_cumulative_hazard(x: Array, seg_end: Array, hazard: Array) -> Array:
    """Invert H(t)=x for piecewise-constant hazard. Returns t."""
    x = np.asarray(x, dtype=float)
    seg_end = np.asarray(seg_end, dtype=float)
    hazard = np.asarray(hazard, dtype=float)

    seg_start = np.concatenate([[0.0], seg_end[:-1]])
    seg_dur = seg_end - seg_start
    seg_H = hazard * seg_dur
    seg_H_end = np.cumsum(seg_H)

    # find first segment where cumulative hazard exceeds x
    idx = np.searchsorted(seg_H_end, x, side="right")  # 0..nseg
    t_out = np.zeros_like(x, dtype=float)

    # within known segments
    in_seg = idx < len(seg_end)
    k = idx[in_seg]
    H_prev = np.where(k > 0, seg_H_end[k - 1], 0.0)
    # avoid divide-by-zero hazard
    lam = np.maximum(hazard[k], 1e-16)
    t_out[in_seg] = seg_start[k] + (x[in_seg] - H_prev) / lam

    # beyond last segment: extend last hazard
    beyond = ~in_seg
    if np.any(beyond):
        lam_last = max(float(hazard[-1]), 1e-16)
        t_out[beyond] = seg_end[-1] + (x[beyond] - seg_H_end[-1]) / lam_last

    return t_out


def cds_par_spread_from_survival(
    S: Array,
    pay_times: Array,
    df: Array,
    recovery: float,
    accrual_on_default: bool = True,
) -> float:
    """
    Given survival at pay_times, compute par spread (decimal) under standard discrete approximation.

    Protection PV: (1-R) * Σ DF(t_i) * (S(t_{i-1}) - S(t_i))
    Premium PV per unit spread:
        Σ DF(t_i) * Δ * S(t_i) + (optional) Σ DF(t_i) * 0.5*Δ*(S(t_{i-1}) - S(t_i))
    """
    dt = np.diff(np.concatenate([[0.0], pay_times]))
    S_prev = np.concatenate([[1.0], S[:-1]])

    prot = (1.0 - recovery) * np.sum(df * (S_prev - S))
    prem = np.sum(df * dt * S)

    if accrual_on_default:
        prem += np.sum(df * 0.5 * dt * (S_prev - S))

    return prot / prem


def bootstrap_hazard_from_cds(
    maturities: Array,
    spreads_bps: Array,
    r: float,
    recovery: float,
    pay_freq: int = 4,
) -> HazardCurve:
    """
    Bootstrap a piecewise-constant hazard curve that matches CDS par spreads at each maturity.

    Inputs:
    - maturities: e.g., [1,2,3,4,5]
    - spreads_bps: par spreads in bps at each maturity
    - r: flat rate (cc)
    - recovery: constant recovery
    - pay_freq: payments per year (4 for quarterly)
    """
    maturities = np.asarray(maturities, dtype=float)
    spreads = np.asarray(spreads_bps, dtype=float) * 1e-4  # bps -> decimal
    assert maturities.ndim == 1 and spreads.ndim == 1 and len(maturities) == len(spreads)

    seg_end = maturities.copy()
    hazard = np.zeros_like(seg_end)

    # master payment grid up to max maturity
    dt = 1.0 / pay_freq
    pay_times_all = np.arange(dt, maturities[-1] + 1e-12, dt)

    for k, (Tk, sk) in enumerate(zip(maturities, spreads)):
        # payment times up to Tk
        pay_times = pay_times_all[pay_times_all <= Tk + 1e-12]
        df = discount_factors_flat(r, pay_times)

        # survival from previous segments, and trial hazard for current segment
        def par_spread_given_lambda(lam: float) -> float:
            hazard_trial = hazard.copy()
            hazard_trial[k] = lam
            curve = HazardCurve(seg_end=seg_end[: k + 1], hazard=hazard_trial[: k + 1])
            S = curve.survival(pay_times)
            return cds_par_spread_from_survival(S, pay_times, df, recovery)

        # solve par_spread_given_lambda(lam) - sk = 0
        # hazard typically in [1e-6, 5] for IG->HY; widen if needed
        f = lambda lam: par_spread_given_lambda(lam) - sk
        lo, hi = 1e-8, 10.0
        # ensure bracket; expand if required
        f_lo, f_hi = f(lo), f(hi)
        while f_lo * f_hi > 0:
            hi *= 2.0
            if hi > 200:
                raise RuntimeError("Failed to bracket hazard; check inputs.")
            f_hi = f(hi)

        hazard[k] = brentq(f, lo, hi, maxiter=200)

    return HazardCurve(segments_end=seg_end, hazard=hazard)


def _chol_psd(corr: Array, jitter: float = 1e-12) -> Array:
    """Cholesky with diagonal jitter for near-PSD matrices."""
    corr = np.asarray(corr, dtype=float)
    if corr.shape[0] != corr.shape[1]:
        raise ValueError("corr must be square")
    # symmetry
    corr = 0.5 * (corr + corr.T)
    # add jitter until Cholesky succeeds
    for j in range(10):
        try:
            return np.linalg.cholesky(corr + (jitter * (10**j)) * np.eye(corr.shape[0]))
        except np.linalg.LinAlgError:
            continue
    # final fallback: eigenvalue cleaning
    w, v = np.linalg.eigh(corr)
    w = np.clip(w, 1e-10, None)
    corr_clean = (v * w) @ v.T
    corr_clean = corr_clean / np.sqrt(np.outer(np.diag(corr_clean), np.diag(corr_clean)))
    return np.linalg.cholesky(corr_clean)


def sample_copula_uniforms(
    n_sims: int,
    corr: Array,
    copula: str = "gaussian",
    t_df: int = 6,
    seed: Optional[int] = 1234,
) -> Array:
    """
    Draw correlated uniforms U in [0,1]^d from a copula.

    copula:
        - "gaussian"
        - "t"
    """
    rng = np.random.default_rng(seed)
    d = corr.shape[0]
    L = _chol_psd(corr)

    Z = rng.standard_normal((n_sims, d)) @ L.T

    if copula.lower() in {"gaussian", "normal"}:
        return norm.cdf(Z)

    if copula.lower() in {"t", "student", "student-t"}:
        # multivariate t via scale mixture: Z / sqrt(V/df), V ~ chi2(df)
        V = rng.chisquare(df=t_df, size=n_sims) / float(t_df)
        T = Z / np.sqrt(V)[:, None]
        return student_t.cdf(T, df=t_df)

    raise ValueError(f"Unknown copula={copula}")


def simulate_default_times(
    curves: List[HazardCurve],
    U: Array,
) -> Array:
    """Convert copula uniforms to default times, shape (n_sims, d)."""
    n_sims, d = U.shape
    if d != len(curves):
        raise ValueError("U dimension must match number of curves")
    taus = np.zeros_like(U, dtype=float)
    for i, curve in enumerate(curves):
        taus[:, i] = curve.inverse_default_time(U[:, i])
    return taus


@dataclass(frozen=True)
class BasketCDSContract:
    maturity: float
    k: int  # kth-to-default
    notional: float
    recovery: float
    pay_freq: int = 4


def price_nth_to_default_mc(
    contract: BasketCDSContract,
    curves: List[HazardCurve],
    corr: Array,
    r: float,
    copula: str = "gaussian",
    t_df: int = 6,
    n_sims: int = 50_000,
    seed: int = 1234,
) -> Dict[str, float]:
    """
    Monte Carlo pricing of kth-to-default basket CDS.

    Returns:
        par_spread_bps, pv_protection, pv_rpv01
    """
    d = len(curves)
    if contract.k < 1 or contract.k > d:
        raise ValueError("k must be between 1 and number of names")

    dt = 1.0 / contract.pay_freq
    pay_times = np.arange(dt, contract.maturity + 1e-12, dt)
    df_pay = discount_factors_flat(r, pay_times)

    U = sample_copula_uniforms(n_sims, corr, copula=copula, t_df=t_df, seed=seed)
    taus = simulate_default_times(curves, U)
    taus_sorted = np.sort(taus, axis=1)
    tau_k = taus_sorted[:, contract.k - 1]

    # Protection leg: pay LGD on kth default if it occurs before maturity
    lgd = 1.0 - contract.recovery
    name_notional = contract.notional / d
    df_tau = np.exp(-r * np.minimum(tau_k, contract.maturity))
    prot_pay = lgd * name_notional * (tau_k <= contract.maturity) * df_tau
    pv_prot = float(np.mean(prot_pay))

    # Premium PV01 per 1bp: payments until min(tau_k, T) + accrued on default
    # For each path, indicator of whether payment date occurs before default
    alive = pay_times[None, :] <= np.minimum(tau_k[:, None], contract.maturity)
    # scheduled premium accruals
    rpv01_sched = np.sum(df_pay[None, :] * dt * alive, axis=1) * contract.notional
    # accrued premium on default (approx): from last payment to default
    # find last payment time strictly before default
    last_idx = np.sum(pay_times[None, :] < np.minimum(tau_k[:, None], contract.maturity), axis=1) - 1
    last_idx = np.clip(last_idx, -1, len(pay_times) - 1)
    last_pay = np.where(last_idx >= 0, pay_times[last_idx], 0.0)
    accr = (tau_k <= contract.maturity) * (tau_k - last_pay)
    rpv01_accr = contract.notional * np.exp(-r * tau_k) * accr
    pv_rpv01 = float(np.mean(rpv01_sched + rpv01_accr))

    par_spread = pv_prot / pv_rpv01  # decimal
    return {
        "par_spread_bps": 1e4 * par_spread,
        "pv_protection": pv_prot,
        "pv_rpv01": pv_rpv01,
    }


def loss_distribution_at_horizon(
    contract: BasketCDSContract,
    curves: List[HazardCurve],
    corr: Array,
    r: float,
    horizon: float = 1.0,
    copula: str = "gaussian",
    t_df: int = 6,
    n_sims: int = 200_000,
    seed: int = 1234,
) -> Dict[str, float]:
    """
    A simple risk-capital oriented loss distribution for the protection seller:
    - Loss is discounted protection payment if kth default occurs within horizon.
    - Premium inflows are ignored (conservative for seller capital) unless you incorporate MtM.

    This is not a full IM/CSA MtM simulation; it is a clean starting point for EC-style capital.
    """
    d = len(curves)
    U = sample_copula_uniforms(n_sims, corr, copula=copula, t_df=t_df, seed=seed)
    taus = simulate_default_times(curves, U)
    tau_k = np.sort(taus, axis=1)[:, contract.k - 1]

    lgd = 1.0 - contract.recovery
    name_notional = contract.notional / d
    loss = lgd * name_notional * (tau_k <= horizon) * np.exp(-r * tau_k)

    loss = loss.astype(float)
    el = float(np.mean(loss))
    var99 = float(np.quantile(loss, 0.99))
    var995 = float(np.quantile(loss, 0.995))

    def es(alpha: float) -> float:
        q = np.quantile(loss, alpha)
        tail = loss[loss >= q]
        return float(np.mean(tail)) if len(tail) else float(q)

    return {
        "expected_loss": el,
        "var_99": var99,
        "var_995": var995,
        "es_99": es(0.99),
        "es_995": es(0.995),
    }
