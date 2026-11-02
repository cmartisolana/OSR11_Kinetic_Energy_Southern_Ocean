"""Trend statistics for the reviewer-2 revision (plan item A1).

One framework for every reported trend:
  * deseasonalise monthly series (remove the monthly climatology);
  * Theil-Sen slope [Sen 1968];
  * 95 % CI and two-sided p-value from an AR-sieve Monte Carlo of the detrended
    residuals [Buhlmann 1997; cf. the AR(6) Monte Carlo of Zhang et al. 2021]:
    an AR(p) model (p <= 12 for monthly, <= 3 for annual data, chosen by AIC)
    is fitted to the residuals, B surrogate series are simulated, and the
    Theil-Sen slope of each gives the null distribution. CI = slope minus the
    2.5/97.5 % quantiles of the null slopes.
    (A circular moving-block bootstrap was tested first and rejected: in
    synthetic AR(1) tests, phi = 0.7-0.9, n = 300, it gave 13-18 % false
    positives at alpha = 0.05; the sieve gave ~5-7 %.)
  * secondary: Hamed & Rao (1998) modified Mann-Kendall (pymannkendall),
    reported for comparison only (also liberal for strongly persistent series);
  * maps: Benjamini-Hochberg false discovery rate (Wilks 2016, alpha_FDR = 0.10).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

try:
    import pymannkendall as _pmk
except Exception:  # pragma: no cover
    _pmk = None


# ---------------------------------------------------------------- basics
def theil_sen(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Median of pairwise slopes; intercept = median(y - b x)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = x.size
    if n < 3:
        return np.nan, np.nan
    i, j = np.triu_indices(n, 1)
    dx = x[j] - x[i]
    ok = dx != 0
    b = np.median((y[j] - y[i])[ok] / dx[ok])
    a = np.median(y - b * x)
    return float(b), float(a)


def deseasonalise(y: np.ndarray, month: np.ndarray) -> np.ndarray:
    """Remove the mean of each calendar month (anomalies keep the overall mean)."""
    y = np.asarray(y, float).copy()
    month = np.asarray(month)
    overall = np.nanmean(y)
    for mo in np.unique(month):
        k = month == mo
        y[k] = y[k] - np.nanmean(y[k]) + overall
    return y


def efolding_time(r: np.ndarray, max_lag: int | None = None) -> float:
    """First lag at which the sample ACF of r drops below 1/e (in samples)."""
    r = np.asarray(r, float)
    r = r[np.isfinite(r)] - np.nanmean(r)
    n = r.size
    max_lag = max_lag or n // 3
    den = np.sum(r * r)
    if den == 0:
        return 1.0
    for k in range(1, max_lag):
        if np.sum(r[:-k] * r[k:]) / den < np.exp(-1):
            return float(k)
    return float(max_lag)


def block_length(r: np.ndarray, monthly: bool = True) -> int:
    tau = efolding_time(r)
    n = int(np.isfinite(r).sum())
    lo = 12 if monthly else 2
    return int(np.clip(np.ceil(2 * tau), lo, max(lo, n // 4)))


def _circular_block_indices(n: int, L: int, rng: np.random.Generator) -> np.ndarray:
    nb = int(np.ceil(n / L))
    starts = rng.integers(0, n, nb)
    idx = (starts[:, None] + np.arange(L)[None, :]) % n
    return idx.ravel()[:n]


@dataclass
class Trend:
    slope: float
    ci_lo: float
    ci_hi: float
    p_boot: float
    p_hamed_rao: float
    n: int
    ar_order: int
    intercept: float

    def as_dict(self) -> dict:
        return asdict(self)


def _fit_ar(res: np.ndarray, maxp: int):
    """AIC-selected AR(p) (OLS, no constant) fitted to mean-removed residuals."""
    from statsmodels.tsa.ar_model import AutoReg
    r = res - res.mean()
    best = None
    for p in range(1, maxp + 1):
        if r.size <= 3 * p:
            break
        m = AutoReg(r, lags=p, trend="n").fit()
        if best is None or m.aic < best[0]:
            best = (m.aic, p, m)
    if best is None:
        return np.zeros(0), float(np.std(r))
    _, p, m = best
    return np.asarray(m.params, float), float(np.std(m.resid, ddof=p))


def _simulate_ar(phis: np.ndarray, sig: float, n: int, B: int, rng, burn: int = 200) -> np.ndarray:
    p = phis.size
    z = np.zeros((B, n + burn))
    e = rng.normal(0.0, sig, (B, n + burn))
    for t in range(n + burn):
        acc = e[:, t]
        for k in range(1, p + 1):
            if t - k >= 0:
                acc = acc + phis[k - 1] * z[:, t - k]
        z[:, t] = acc
    return z[:, burn:]


def _ts_many(x: np.ndarray, Y: np.ndarray, chunk: int = 100) -> np.ndarray:
    i, j = np.triu_indices(x.size, 1)
    dx = x[j] - x[i]
    ok = dx != 0
    i, j, dx = i[ok], j[ok], dx[ok]
    out = np.empty(Y.shape[0])
    for a in range(0, Y.shape[0], chunk):
        Yc = Y[a:a + chunk]
        out[a:a + chunk] = np.median((Yc[:, j] - Yc[:, i]) / dx, axis=1)
    return out


def trend(x: np.ndarray, y: np.ndarray, *, month: np.ndarray | None = None,
          monthly: bool = True, n_boot: int = 1000, seed: int = 0) -> Trend:
    """Theil-Sen slope with AR-sieve Monte Carlo CI and p-value.

    If `month` is given the series is deseasonalised first. NaNs are dropped
    (the gaps are closed; acceptable for the few-percent gaps in these data).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if month is not None:
        y = deseasonalise(y, month)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(x.size)
    b, a = theil_sen(x, y)
    if not np.isfinite(b) or n < 8:
        return Trend(np.nan, np.nan, np.nan, np.nan, np.nan, n, 0, np.nan)
    res = y - (a + b * x)
    phis, sig = _fit_ar(res, 12 if monthly else 3)
    rng = np.random.default_rng(seed)
    sims = _simulate_ar(phis, sig, n, n_boot, rng)
    null = _ts_many(x, sims)
    q_lo, q_hi = np.percentile(null, [2.5, 97.5])
    p = (1 + np.sum(np.abs(null) >= abs(b))) / (1 + null.size)
    p_hr = np.nan
    if _pmk is not None and n >= 8:
        try:
            p_hr = float(_pmk.hamed_rao_modification_test(y).p)
        except Exception:
            p_hr = np.nan
    return Trend(float(b), float(b - q_hi), float(b - q_lo), float(p), p_hr, n,
                 int(phis.size), float(a))


# ---------------------------------------------------------------- maps
def bh_fdr(p: np.ndarray, alpha: float = 0.10) -> np.ndarray:
    """Benjamini-Hochberg: boolean mask of discoveries (NaNs -> False)."""
    p = np.asarray(p, float)
    flat = p.ravel()
    ok = np.isfinite(flat)
    pv = flat[ok]
    m = pv.size
    out = np.zeros(flat.shape, bool)
    if m == 0:
        return out.reshape(p.shape)
    order = np.argsort(pv)
    thresh = alpha * np.arange(1, m + 1) / m
    passed = pv[order] <= thresh
    if passed.any():
        kmax = np.max(np.nonzero(passed)[0])
        sel = np.zeros(m, bool)
        sel[order[: kmax + 1]] = True
        tmp = out[ok]
        tmp[:] = sel
        out[ok] = tmp
    return out.reshape(p.shape)


def _null_sd(x, y, monthly, n_boot, rng):
    b, a = theil_sen(x, y)
    res = y - (a + b * x)
    phis, sig = _fit_ar(res, 12 if monthly else 3)
    sims = _simulate_ar(phis, sig, x.size, n_boot, rng)
    return b, float(np.std(_ts_many(x, sims))), phis, sig


def trend_calibrated(x: np.ndarray, y: np.ndarray, *, month: np.ndarray | None = None,
                     monthly: bool = True, n_outer: int = 200, n_boot: int = 300,
                     seed: int = 0) -> dict:
    """Sieve trend plus a double-bootstrap calibrated p-value.

    Short persistent records make any re-estimated noise model too optimistic
    (in tests the plain sieve gave 7-13 % false positives at nominal 5 %).
    The statistic z = slope / sd(null slopes from the re-estimated AR model) is
    computed for the data and for n_outer null series simulated from the AR
    model fitted to the data (each with its own re-estimated AR model);
    p_cal = fraction of |z_null| >= |z_obs|; CI_cal = slope +/- q95(|z_null|) * sd.
    """
    import os
    if os.environ.get("RSTATS_QUICK"):          # for pipeline smoke tests only
        n_outer, n_boot = 5, 50
    T = trend(x, y, month=month, monthly=monthly, n_boot=n_boot, seed=seed)
    xx = np.asarray(x, float)
    yy = np.asarray(y, float)
    if month is not None:
        yy = deseasonalise(yy, month)
    m = np.isfinite(xx) & np.isfinite(yy)
    xx, yy = xx[m], yy[m]
    rng = np.random.default_rng(seed + 1)
    b, sd, phis, sig = _null_sd(xx, yy, monthly, n_boot, rng)
    if not np.isfinite(sd) or sd == 0:
        out = T.as_dict(); out.update(p_cal=np.nan, z=np.nan, z_crit_cal=np.nan, ci_cal_lo=np.nan, ci_cal_hi=np.nan)
        return out
    z_obs = b / sd
    sims = _simulate_ar(phis, sig, xx.size, n_outer, rng)
    zn = np.empty(n_outer)
    for i, s_ in enumerate(sims):
        bi, sdi, _, _ = _null_sd(xx, s_, monthly, n_boot, rng)
        zn[i] = bi / sdi if sdi > 0 else np.nan
    zn = zn[np.isfinite(zn)]
    p_cal = (1 + np.sum(np.abs(zn) >= abs(z_obs))) / (1 + zn.size)
    zc = float(np.quantile(np.abs(zn), 0.95))
    out = T.as_dict()
    out.update(p_cal=float(p_cal), z=float(z_obs), z_crit_cal=zc,
               ci_cal_lo=float(b - zc * sd), ci_cal_hi=float(b + zc * sd))
    return out


# ---------------------------------------------------------------- v2 test (2026-10-13)
# Bias-corrected AR sieve. The AR model fitted to the detrended residuals of a short persistent
# record underestimates persistence (small-sample bias of AR estimates, increased by removing the
# fitted trend), so surrogates simulated from it are too smooth-free and the test rejects too often
# (stats_size_check.py: 15-19 % false positives for the plain sieve, 6-11 % for trend_calibrated).
# Fix: bootstrap bias correction of the AR coefficients (Kilian 1998, Rev. Econ. Stat. 80, 218-230),
# with the trend removal included in the bias simulation and Kilian's stationarity adjustment.
# trend()/trend_calibrated() above are unchanged, so earlier outputs remain reproducible.
def _ar_design(R: np.ndarray, p: int):
    """R: (B, n) demeaned series -> X (B, n-p, p), Y (B, n-p) for OLS AR(p) without constant."""
    n = R.shape[1]
    X = np.stack([R[:, p - k:n - k] for k in range(1, p + 1)], axis=2)
    return X, R[:, p:]


def _ols_ar_many(R: np.ndarray, p: int) -> np.ndarray:
    X, Y = _ar_design(R, p)
    XtX = np.einsum("bti,btj->bij", X, X)
    XtY = np.einsum("bti,bt->bi", X, Y)
    return np.linalg.solve(XtX, XtY[..., None])[..., 0]


def _is_stationary(phis: np.ndarray) -> bool:
    p = phis.size
    if p == 0:
        return True
    C = np.zeros((p, p)); C[0] = phis
    if p > 1:
        C[1:, :-1] = np.eye(p - 1)
    return bool(np.max(np.abs(np.linalg.eigvals(C))) < 1.0)


def _bias_correct_ar(phis: np.ndarray, sig: float, x: np.ndarray, B: int, rng) -> np.ndarray:
    """Kilian (1998) bias correction; the bias simulation repeats the trend removal (OLS on x)."""
    p = phis.size
    if p == 0 or not _is_stationary(phis):
        return phis
    sims = _simulate_ar(phis, sig, x.size, B, rng)
    X = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(X, sims.T, rcond=None)
    res = sims - (X @ coef).T
    res = res - res.mean(1, keepdims=True)
    bias = _ols_ar_many(res, p).mean(0) - phis
    for delta in np.arange(1.0, -0.001, -0.01):
        phi_bc = phis - delta * bias
        if _is_stationary(phi_bc):
            return phi_bc
    return phis


def _null_bc(x, y, monthly, n_boot, n_bias, rng):
    b, a = theil_sen(x, y)
    res = y - (a + b * x)
    phis, sig = _fit_ar(res, 12 if monthly else 3)
    phi_bc = _bias_correct_ar(phis, sig, x, n_bias, rng)
    null = _ts_many(x, _simulate_ar(phi_bc, sig, x.size, n_boot, rng))
    return b, a, null, phis, phi_bc, sig


def trend_bc(x: np.ndarray, y: np.ndarray, *, month: np.ndarray | None = None, monthly: bool = True,
             n_boot: int = 2000, n_bias: int = 200, calibrate: bool = False, n_outer: int = 100,
             seed: int = 0) -> dict:
    """Theil-Sen slope; CI and two-sided p from a bias-corrected AR-sieve null.

    calibrate=True adds the double-bootstrap calibration of trend_calibrated (z = slope / sd of the
    null slopes; outer nulls simulated from the bias-corrected model, each re-estimated the same way).
    Returns slope, ci_lo, ci_hi, p (and z, p_cal, ci_cal_* if calibrated), ar_order, phi_sum, phi_sum_bc.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    if month is not None:
        y = deseasonalise(y, month)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    rng = np.random.default_rng(seed)
    b, a, null, phis, phi_bc, sig = _null_bc(x, y, monthly, n_boot, n_bias, rng)
    q_lo, q_hi = np.percentile(null, [2.5, 97.5])
    p = (1 + np.sum(np.abs(null) >= abs(b))) / (1 + null.size)
    out = dict(slope=float(b), intercept=float(a), ci_lo=float(b - q_hi), ci_hi=float(b - q_lo), p=float(p),
               n=int(x.size), ar_order=int(phis.size), phi_sum=float(phis.sum()), phi_sum_bc=float(phi_bc.sum()),
               sd_null=float(np.std(null)))
    if calibrate:
        sd = float(np.std(null)); z_obs = b / sd
        nb_in = min(n_boot, 300)
        zn = []
        for s_ in _simulate_ar(phi_bc, sig, x.size, n_outer, rng):
            bi, _, nulli, *_ = _null_bc(x, s_, monthly, nb_in, n_bias, rng)
            sdi = np.std(nulli)
            if sdi > 0:
                zn.append(bi / sdi)
        zn = np.abs(np.array(zn))
        zc = float(np.quantile(zn, 0.95))
        out.update(z=float(z_obs), p_cal=float((1 + np.sum(zn >= abs(z_obs))) / (1 + zn.size)), z_crit_cal=zc,
                   ci_cal_lo=float(b - zc * sd), ci_cal_hi=float(b + zc * sd))
    return out


# ---------------------------------------------------------------- v2b: bias-corrected + calibrated, fast
# Size checks (stats_size_check_bc.py) showed that bias correction alone fixes the persistence but
# not the test size (~9-10 % at nominal 5 % for AR(1) 0.5-0.8, n = 300): the uncertainty of the
# estimated noise model also has to be propagated. trend_bcz therefore calibrates the statistic
# z = Theil-Sen slope / sd, where sd is the exact OLS-slope sd under the bias-corrected AR model
# (Toeplitz autocovariance), against outer null series simulated from the bias-corrected model and
# processed identically (deseasonalised, AR re-selected and bias-corrected). The analytic sd makes
# n_outer = 500-1000 affordable (p floor 0.002-0.001). The calibration absorbs the small TS/OLS
# efficiency difference because data and nulls use the same statistic.
def _fit_ar_fast(r: np.ndarray, maxp: int):
    """AIC-selected AR(p), OLS without constant, common estimation sample (hold-back = maxp)."""
    r = r - r.mean(); n = r.size
    maxp = int(min(maxp, max(1, n // 3 - 1)))
    L = np.column_stack([r[maxp - k:n - k] for k in range(1, maxp + 1)]); Y = r[maxp:]
    best = None
    for p in range(1, maxp + 1):
        X = L[:, :p]
        ph, *_ = np.linalg.lstsq(X, Y, rcond=None)
        s2 = np.mean((Y - X @ ph) ** 2)
        aic = Y.size * np.log(s2) + 2 * p
        if best is None or aic < best[0]:
            best = (aic, p)
    p = best[1]
    X = np.column_stack([r[p - k:n - k] for k in range(1, p + 1)]); Yp = r[p:]
    ph, *_ = np.linalg.lstsq(X, Yp, rcond=None)
    sig = float(np.std(Yp - X @ ph, ddof=p))
    return np.asarray(ph, float), sig


def _ols_slope_sd(phis: np.ndarray, sig: float, x: np.ndarray) -> float:
    from statsmodels.tsa.arima_process import arma_acovf
    from scipy.linalg import toeplitz
    g = arma_acovf(np.r_[1.0, -phis], np.array([1.0]), nobs=x.size, sigma2=sig ** 2)
    xc = x - x.mean(); c = xc / np.sum(xc ** 2)
    return float(np.sqrt(c @ toeplitz(g) @ c))


def _shrink_to_stationary(phis: np.ndarray, rmax: float = 0.99) -> np.ndarray:
    """Scale the AR roots so that the largest companion eigenvalue modulus is rmax (no-op if stationary)."""
    p = phis.size
    if p == 0 or _is_stationary(phis):
        return phis
    C = np.zeros((p, p)); C[0] = phis
    if p > 1:
        C[1:, :-1] = np.eye(p - 1)
    lam = float(np.max(np.abs(np.linalg.eigvals(C))))
    return phis * (rmax / lam) ** np.arange(1, p + 1)


def _bias_correct_ar_iter(phis, sig, x, B, rng, n_iter=6):
    """Iterated bootstrap bias correction: find phi with E[phi_hat | phi] = phi_hat_obs (mean-unbiased).
    Used for annual series (n = 15-25), where the one-step Kilian correction is far too small."""
    target = phis.copy(); cur = _bias_correct_ar(phis, sig, x, B, rng)
    X = np.column_stack([np.ones_like(x), x])
    for _ in range(n_iter):
        cur = _shrink_to_stationary(cur)
        sims = _simulate_ar(cur, sig, x.size, B, rng)
        coef, *_ = np.linalg.lstsq(X, sims.T, rcond=None)
        res = sims - (X @ coef).T; res = res - res.mean(1, keepdims=True)
        cur = cur + (target - _ols_ar_many(res, cur.size).mean(0))
    return _shrink_to_stationary(cur)


import os as _os
# Monthly AR order for trend_bcz: a fixed order (default 12) or "aic". The independent size check on a
# process with the observed (Drake-box EKE) autocorrelation found 9.5 % false positives with AIC:
# AIC picks AR(2) and misses the small high-lag coefficients that carry part of the low-frequency
# persistence. A fixed AR(12) keeps them; the bias correction and the calibration absorb the cost.
MONTHLY_ORDER = _os.environ.get("RSTATS_MONTHLY_ORDER", "both")   # "both" | "aic" | "<int>"


def _fit_ar_order(r: np.ndarray, p: int):
    """OLS AR(p) without constant, full sample."""
    r = r - r.mean(); n = r.size
    X = np.column_stack([r[p - k:n - k] for k in range(1, p + 1)]); Y = r[p:]
    ph, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return np.asarray(ph, float), float(np.std(Y - X @ ph, ddof=p))


def _z_stat(x, y, monthly, n_bias, rng):
    """Statistic of trend_bcz. Monthly: AIC AR(p<=12), one-step Kilian correction. Annual: AR(1),
    iterated bias correction (2026-10-13, after the independent check found annual size 7-12 %).
    Non-stationary models are shrunk to max root modulus 0.99 (identically for data and nulls)."""
    b, a = theil_sen(x, y)
    if monthly:
        if MONTHLY_ORDER == "aic":
            phis, sig = _fit_ar_fast(y - (a + b * x), 12)
        else:
            phis, sig = _fit_ar_order(y - (a + b * x), int(MONTHLY_ORDER))
        phi_bc = _shrink_to_stationary(_bias_correct_ar(phis, sig, x, n_bias, rng))
    else:
        phis, sig = _fit_ar_fast(y - (a + b * x), 1)
        phi_bc = _bias_correct_ar_iter(phis, sig, x, max(n_bias, 200), rng)
    sd = _ols_slope_sd(phi_bc, sig, x)
    return b, a, sd, phis, phi_bc, sig


def trend_bcz(x: np.ndarray, y: np.ndarray, *, month: np.ndarray | None = None, monthly: bool = True,
              n_outer: int = 300, n_bias: int = 100, seed: int = 0) -> dict:
    """Reported test (2026-10-13). Monthly series: the calibrated z test is run with two noise models,
    AIC-selected order and fixed AR(12); p = max of the two, CI = union of the two intervals.
    Measured false-positive rates at nominal 5 % (stats_size_check_bc.py, n = 300):
      AIC order:  AR(1) 0.8 -> 5.3 % (300 series); Drake-EKE-like AR(12) process -> 9.5 % (200)
      AR(12):     AR(1) 0.8 -> 8 % (100);           Drake-EKE-like process -> 3 % (100)
    so the max rule is bounded by the smaller of the two in each case. Annual series: AR(1) with
    iterated bias correction (size 5.0 % at n = 25, 9.5 % at n = 15, AR(1) 0.5, 200 series)."""
    global MONTHLY_ORDER
    if monthly and MONTHLY_ORDER == "both":
        outs = []
        for order in ("aic", "12"):
            MONTHLY_ORDER = order
            try:
                outs.append(_trend_bcz_one(x, y, month=month, monthly=True, n_outer=n_outer, n_bias=n_bias, seed=seed))
            finally:
                MONTHLY_ORDER = "both"
        o = dict(outs[np.argmax([q["p"] for q in outs])])
        o.update(p_aic=outs[0]["p"], p_ar12=outs[1]["p"], z_aic=outs[0]["z"], z_ar12=outs[1]["z"],
                 ci_lo=min(q["ci_lo"] for q in outs), ci_hi=max(q["ci_hi"] for q in outs), p=max(q["p"] for q in outs))
        return o
    return _trend_bcz_one(x, y, month=month, monthly=monthly, n_outer=n_outer, n_bias=n_bias, seed=seed)


def _trend_bcz_one(x: np.ndarray, y: np.ndarray, *, month: np.ndarray | None = None, monthly: bool = True,
                   n_outer: int = 300, n_bias: int = 100, seed: int = 0) -> dict:
    """Theil-Sen slope; p and 95 % CI from the calibrated z statistic (see block comment above)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    mo = None if month is None else np.asarray(month)
    if mo is not None:
        y = deseasonalise(y, mo)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if mo is not None:
        mo = mo[m]
    rng = np.random.default_rng(seed)
    b, a, sd, phis, phi_bc, sig = _z_stat(x, y, monthly, n_bias, rng)
    z = b / sd
    zn = np.empty(n_outer)
    for i, s_ in enumerate(_simulate_ar(phi_bc, sig, x.size, n_outer, rng, burn=500)):
        if mo is not None:
            s_ = deseasonalise(s_, mo)
        bi, _, sdi, *_ = _z_stat(x, s_, monthly, n_bias, rng)
        zn[i] = bi / sdi
    zn = np.abs(zn[np.isfinite(zn)])
    zc = float(np.quantile(zn, 0.95))
    return dict(slope=float(b), intercept=float(a), z=float(z), p=float((1 + np.sum(zn >= abs(z))) / (1 + zn.size)),
                z_crit=zc, ci_lo=float(b - zc * sd), ci_hi=float(b + zc * sd), sd=sd, n=int(x.size),
                ar_order=int(phis.size), phi_sum=float(phis.sum()), phi_sum_bc=float(phi_bc.sum()))
