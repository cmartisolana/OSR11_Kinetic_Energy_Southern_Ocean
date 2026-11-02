"""Trend utilities.

This repository standardizes trend estimation on:
- Theil–Sen slope (scipy.stats.theilslopes) for robust linear trend estimates.
- Modified Mann–Kendall (Yue–Wang) for significance testing when available.

The goal is to avoid mixing multiple slope estimators (e.g., OLS vs Theil–Sen)
across notebooks and scripts.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy import stats


def theil_sen_slope(
    x: np.ndarray,
    y: np.ndarray,
    *,
    alpha: float = 0.05,
) -> Tuple[float, float, float, float]:
    """Compute Theil–Sen slope with confidence interval.

    Parameters
    ----------
    x, y : array-like
        Predictor (typically decimal year) and response.
    alpha : float
        Significance level for the CI; default 0.05 gives 95% CI.

    Returns
    -------
    slope, intercept, slope_lo, slope_hi
        If insufficient finite points, returns NaNs.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if int(m.sum()) < 3:
        return np.nan, np.nan, np.nan, np.nan

    res = stats.theilslopes(y[m], x[m], alpha=1.0 - float(alpha))
    return float(res.slope), float(res.intercept), float(res.low_slope), float(res.high_slope)


def _basic_mann_kendall(y: np.ndarray) -> Tuple[float, float]:
    """Fallback (non-modified) Mann–Kendall test.

    This is used only if `pymannkendall` is not available.
    Returns Kendall's tau and a two-sided p-value using a normal approximation
    with tie correction.
    """
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    n = int(y.size)
    if n < 4:
        return np.nan, np.nan

    s = 0.0
    for i in range(n - 1):
        s += np.sign(y[i + 1 :] - y[i]).sum()

    # Tie correction
    _, counts = np.unique(y, return_counts=True)
    tie_term = float(np.sum(counts * (counts - 1) * (2 * counts + 5)))
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0
    if var_s <= 0:
        return np.nan, np.nan

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0

    p = 2.0 * stats.norm.sf(np.abs(z))
    tau = s / (0.5 * n * (n - 1))
    return float(tau), float(p)


def mann_kendall(y: np.ndarray) -> Tuple[float, float]:
    """Modified Mann–Kendall (Yue–Wang) when available; fallback otherwise.

    Returns
    -------
    tau, p_value
    """
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if int(y.size) < 4:
        return np.nan, np.nan

    try:
        import pymannkendall as mk  # type: ignore

        res = mk.yue_wang_modification_test(y)
        return float(res.Tau), float(res.p)
    except Exception:
        return _basic_mann_kendall(y)
