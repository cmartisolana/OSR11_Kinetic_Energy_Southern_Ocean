"""x10 (9 Nov 2026): 95 % intervals of the power-law exponents of Fig. 4d (x9: response 2004-2025, wind 2003-2024) by parametric
surrogates, as for the trends (Sect. S1.7): the response is modelled as b*ln tau(y-1) + c*t + AR(1) noise (lag-1 coefficient
bias-corrected, phi + (1 + 3 phi)/n, Kendall 1954), tau is replaced by phase-randomised surrogates of the observed ln tau, the exponent
is re-estimated as in the analysis (trends retained; both series detrended), and the interval is the basic bootstrap interval
b_hat - f*quantiles(b_sim - b_hat), f = t(0.975, n-3)/z(0.975) = 1.07 for the estimated noise level. Coverage is verified by a nested Monte Carlo (outer synthetic records with the true phi).
The x9 pairs bootstrap, block bootstrap and HAC intervals covered only 76-92 % (x9_coverage.csv).
Outputs: x10/x10_intervals.csv, x10/x10_coverage.csv, x10/x10_log.txt

Origin  : x10_sensitivity_intervals.py (9 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : 95 % intervals of the exponents of Fig. 4d (SI S1.6, S5.2).
Inputs  : outputs/f4e/f4e_series.csv (04g)
Run     : python scripts/04_wind_argo/04j_exponent_intervals.py   (no raw data; minutes)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd
from osr11 import paths as PA
OUT = PA.out("x10")
S = pd.read_csv(PA.OUT / "f4e" / "f4e_series.csv")
yrs = np.arange(2001, 2026); ann = lambda c: np.array([S[c][S.year == y].mean() for y in yrs])
W = ann("wind"); X = {k: ann(k) for k in ("KE", "EKE", "slope")}; c4 = yrs >= 2004
sel_y = yrs >= 2004; sel_w = (yrs >= 2003) & (yrs <= 2024); lw = np.log(W / W.mean()); wx = lw[sel_w]; n = wx.size; t = np.arange(n)
def slope_v(x, y):
    xc = x - x.mean(-1, keepdims=True); return (xc * (y - y.mean(-1, keepdims=True))).sum(-1) / (xc * xc).sum(-1)
def detr_v(y):
    tc = t - t.mean(); return y - y.mean(-1, keepdims=True) - np.outer(np.atleast_2d(y)[:, :] @ tc / (tc @ tc), tc).reshape(y.shape)
F0 = np.fft.rfft(lw - lw.mean()); AMP = np.abs(F0)
def surr_v(m, rng):
    ph = rng.uniform(0, 2 * np.pi, (m, F0.size)); ph[:, 0] = 0
    if lw.size % 2 == 0: ph[:, -1] = 0
    return np.fft.irfft(AMP * np.exp(1j * ph), n=lw.size, axis=-1)[:, sel_w] + lw.mean()
def ar1_v(m, phi, sig, rng):
    e = np.empty((m, n)); e[:, 0] = rng.normal(0, sig / np.sqrt(1 - phi ** 2), m)
    for i in range(1, n): e[:, i] = phi * e[:, i - 1] + rng.normal(0, sig, m)
    return e
def fit(ws, ys):
    A = np.c_[np.ones(n), ws, t]; coef, *_ = np.linalg.lstsq(A, ys, rcond=None); e = ys - A @ coef
    phi = np.corrcoef(e[1:], e[:-1])[0, 1]; sig = e.std(ddof=3) * np.sqrt(1 - min(phi, 0.95) ** 2)
    return coef, phi, sig
def interval(ws, ys, rng, m=4000):
    coef, phi, sig = fit(ws, ys); phic = float(np.clip(phi + (1 + 3 * phi) / n, -0.9, 0.95))
    b = float(slope_v(ws, ys)); bd = float(slope_v(detr_v(ws), detr_v(ys)))
    sw = surr_v(m, rng); e = ar1_v(m, phic, sig, rng)
    ysim = b * sw + coef[2] * t + e                       # trend-retained exponent: truth b_hat
    bs = slope_v(sw, ysim); lo, hi = b - TF * np.percentile(bs - b, [97.5, 2.5])
    ysd = bd * sw + coef[2] * t + e                        # detrended exponent: truth bd_hat
    bsd = slope_v(detr_v(sw), detr_v(ysd)); lod, hid = bd - TF * np.percentile(bsd - bd, [97.5, 2.5])
    return b, lo, hi, bd, lod, hid, phi, phic
from scipy import stats
TF = stats.t.ppf(0.975, n - 3) / stats.norm.ppf(0.975)   # small-sample factor for the estimated noise level (1.07)
rng = np.random.default_rng(20261109); rows, log = [], []
P = lambda s: (print(s, flush=True), log.append(s))
for k, v in X.items():
    yy = np.log(v / np.nanmean(v[c4]))[sel_y]
    b, lo, hi, bd, lod, hid, phi, phic = interval(wx, yy, rng)
    rows.append(dict(series=k, b=b, b_lo=lo, b_hi=hi, b_detr=bd, bd_lo=lod, bd_hi=hid, resid_ar1=phi, resid_ar1_corrected=phic))
    P(f"{k}: b = {b:.2f} [{lo:.2f}, {hi:.2f}] | detrended {bd:.2f} [{lod:.2f}, {hid:.2f}] | residual AR(1) {phi:.2f} -> {phic:.2f}")
pd.DataFrame(rows).to_csv(OUT / "x10_intervals.csv", index=False)
# nested coverage: outer synthetic records generated with the observed (uncorrected) fit as truth
cov = []
for R_, k in zip(rows, X):
    yy = np.log(X[k] / np.nanmean(X[k][c4]))[sel_y]; coef, phi, sig = fit(wx, yy)
    NO = 400; sw = surr_v(NO, rng); e = ar1_v(NO, float(np.clip(phi, -0.9, 0.95)), sig, rng)
    hit = hitd = 0
    for j in range(NO):
        ys = R_["b"] * sw[j] + coef[2] * t + e[j]
        _, lo, hi, _, _, _, _, _ = interval(sw[j], ys, rng, m=1000); hit += lo <= R_["b"] <= hi
        ysd = R_["b_detr"] * sw[j] + coef[2] * t + e[j]
        _, _, _, _, lod, hid, _, _ = interval(sw[j], ysd, rng, m=1000); hitd += lod <= R_["b_detr"] <= hid
    cov.append(dict(series=k, coverage=hit / NO, coverage_detrended=hitd / NO, n_outer=NO))
    P(f"  nested coverage {k}: trends retained {hit/NO:.3f}, detrended {hitd/NO:.3f} (n = {NO}; MC s.e. ~0.011)")
pd.DataFrame(cov).to_csv(OUT / "x10_coverage.csv", index=False); (OUT / "x10_log.txt").write_text("\n".join(log))
