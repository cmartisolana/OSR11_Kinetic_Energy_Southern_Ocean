"""F4j: joint distribution of the lagged interannual anomalies along the ACC (28 Oct 2026), for Figure 4d.
Input: f4v_bin_series.npz (72 bins of 5 deg longitude of the common cells, monthly 2004-2025: eastward wind stress, EKE, slope).
Per bin: annual means, linear trend removed, standardised by the bin's standard deviation (z). Pairs, wind leading by one year:
(wind_{y-1}, EKE_y) and (wind_{y-1}, slope_y), y = 2005-2025, pooled over bins (72 x 21 = 1512 points).
Pearson r of the pooled z; p from 2000 phase-randomised surrogates of the wind, using the same random phases for every bin in a
surrogate (keeps the along-path correlation of the wind anomalies), each surrogate standardised and lagged identically.
Also the regression slope (z per z) and the mean EKE / slope z in wind-z terciles.
Usage: python f4j_lagjoint.py <f4v out dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

FV, OUT = (Path(a) for a in sys.argv[1:3]); OUT.mkdir(parents=True, exist_ok=True)
Z = np.load(FV / "f4v_bin_series.npz"); YY = Z["year"]; yrs = np.arange(2004, 2026)
def zan(v):
    a = np.array([[s[YY == y].mean() for y in yrs] for s in v]); t = yrs - yrs.mean()
    a = a - a.mean(1, keepdims=True); b = (a * t).sum(1) / (t ** 2).sum(); a = a - b[:, None] * t
    return a / a.std(1, ddof=1, keepdims=True)
W, E, S = zan(Z["wind"]), zan(Z["eke"]), zan(Z["slope"])
w1 = W[:, :-1].ravel(); e1 = E[:, 1:].ravel(); s1 = S[:, 1:].ravel()
rng = np.random.default_rng(20261028); n = W.shape[1]
F = np.fft.rfft(W - W.mean(1, keepdims=True), axis=1); amp = np.abs(F); ph0 = np.angle(F)
null = {"eke": [], "slope": []}
for _ in range(2000):
    ph = rng.uniform(0, 2 * np.pi, F.shape[1]); ph[0] = 0
    if n % 2 == 0: ph[-1] = 0
    Ws = np.fft.irfft(amp * np.exp(1j * (ph0 + ph)), n=n, axis=1); Ws = Ws / Ws.std(1, ddof=1, keepdims=True); ws = Ws[:, :-1].ravel()
    null["eke"].append(np.corrcoef(ws, e1)[0, 1]); null["slope"].append(np.corrcoef(ws, s1)[0, 1])
rows = []
for k, y in (("eke", e1), ("slope", s1)):
    r = np.corrcoef(w1, y)[0, 1]; p = (1 + np.sum(np.abs(null[k]) >= abs(r))) / 2001; b = np.polyfit(w1, y, 1)[0]
    q = np.quantile(w1, [1 / 3, 2 / 3]); terc = [y[w1 <= q[0]].mean(), y[(w1 > q[0]) & (w1 <= q[1])].mean(), y[w1 > q[1]].mean()]
    rows.append(dict(pair=f"wind(y-1) vs {k}(y)", r=r, p=p, slope_z=b, n=y.size, terc_low=terc[0], terc_mid=terc[1], terc_high=terc[2]))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "f4j_results.csv", index=False)
np.savez(OUT / "f4j_points.npz", wind_lag1=w1, eke=e1, slope=s1, lon=np.repeat(Z["lon"], n - 1), year=np.tile(yrs[1:], len(Z["lon"])))
