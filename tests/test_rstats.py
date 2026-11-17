"""Trend test on synthetic series with a known slope (Theil-Sen estimator and the calibrated test of SI Sect. S1.7)."""
import numpy as np
import pytest


def _series(slope_per_year, phi=0.6, n=300, seed=1, amp=1.0):
    rng = np.random.default_rng(seed)
    t = 2001 + (np.arange(n) + 0.5) / 12
    month = np.tile(np.arange(1, 13), n // 12)
    e = np.zeros(n)
    for k in range(1, n):
        e[k] = phi * e[k - 1] + rng.normal(0, amp)
    season = 0.5 * np.sin(2 * np.pi * month / 12)
    return t, month, 10.0 + slope_per_year * (t - t[0]) + season + e


def test_theil_sen_exact_line():
    from osr11 import rstats as R
    x = np.arange(50.0)
    b, a = R.theil_sen(x, 3.0 * x - 7.0)
    assert b == pytest.approx(3.0) and a == pytest.approx(-7.0)


def test_theil_sen_recovers_slope_with_ar1_noise():
    from osr11 import rstats as R
    t, month, y = _series(0.2)
    b, _ = R.theil_sen(t, R.deseasonalise(y, month))
    assert b == pytest.approx(0.2, abs=0.05)


def test_trend_bcz_detects_strong_trend_and_not_null():
    """Calibrated test (fast settings): a clear trend is significant; a pure AR(1) series is not."""
    from osr11 import rstats as R
    t, month, y = _series(0.3, phi=0.5, amp=0.5)
    r = R.trend_bcz(t, y, month=month, monthly=True, n_outer=60, n_bias=30, seed=0)
    assert r["ci_lo"] <= 0.3 <= r["ci_hi"]
    assert r["p"] < 0.05
    t, month, y0 = _series(0.0, phi=0.5, amp=0.5, seed=7)
    r0 = R.trend_bcz(t, y0, month=month, monthly=True, n_outer=60, n_bias=30, seed=0)
    assert r0["ci_lo"] <= 0.0 <= r0["ci_hi"]
    assert r0["p"] > 0.05


def test_bh_fdr_controls_discoveries():
    from osr11 import rstats as R
    p = np.array([0.001, 0.002, 0.5, 0.8, 0.03, 0.9])
    sig = R.bh_fdr(p, alpha=0.10)
    assert sig[0] and sig[1] and not sig[2] and not sig[5]
