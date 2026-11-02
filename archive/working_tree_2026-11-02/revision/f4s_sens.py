"""F4s: sensitivity of KE, EKE and the sigma0 = 27.3 slope to the wind, circumpolar (28 Oct 2026), the classical eddy-saturation
diagram (response against forcing). Input: f4e_series.csv (monthly domain means on the common cells, ACC band + 100 km, ice removed).
Annual means; anomalies in % of the 2004-2025 mean; wind of year y-1 against KE, EKE, slope of year y, y = 2005-2025 (21 points;
lag 1 is the lag of the circumpolar maximum, f4e). Trends are kept (the total response is the point). Elasticity = OLS slope
(% per %), 95 % from 2000 pairs bootstraps; r; p from 2000 phase-randomised surrogates of the annual wind (2001-2025), lagged alike.
Usage: python f4s_sens.py <f4e out dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

FE, OUT = (Path(a) for a in sys.argv[1:3]); OUT.mkdir(parents=True, exist_ok=True)
S = pd.read_csv(FE / "f4e_series.csv"); yrs = np.arange(2001, 2026)
ann = lambda c: np.array([S[c][S.year == y].mean() for y in yrs])
W = ann("wind"); X = {k: ann(k) for k in ("KE", "EKE", "slope")}
c4 = yrs >= 2004; wm = W[c4].mean(); wp = 100 * (W / wm - 1)                         # wind anomaly, % of the 2004-2025 mean
wx = wp[(yrs >= 2004) & (yrs <= 2024)]                                                   # year y-1 for y = 2005-2025
rng = np.random.default_rng(20261028)
def surr(y):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size) + y.mean()
SW = [surr(wp) for _ in range(2000)]
rows = []; pts = {"year": yrs[yrs >= 2005], "wind_lag1_pct": wx}
for k, v in X.items():
    yp = 100 * (v / v[c4].mean() - 1); yy = yp[yrs >= 2005]; pts[f"{k}_pct"] = yy
    b = np.polyfit(wx, yy, 1)[0]; r = np.corrcoef(wx, yy)[0, 1]
    bb = [np.polyfit(wx[i], yy[i], 1)[0] for i in (rng.integers(0, wx.size, wx.size) for _ in range(2000))]
    null = np.array([np.corrcoef(s_[(yrs >= 2004) & (yrs <= 2024)], yy)[0, 1] for s_ in SW]); p = (1 + np.sum(np.abs(null) >= abs(r))) / 2001
    rows.append(dict(series=k, elasticity=b, el_lo=np.percentile(bb, 2.5), el_hi=np.percentile(bb, 97.5), r=r, p=p, n=yy.size))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "f4s_results.csv", index=False); pd.DataFrame(pts).to_csv(OUT / "f4s_points.csv", index=False)
