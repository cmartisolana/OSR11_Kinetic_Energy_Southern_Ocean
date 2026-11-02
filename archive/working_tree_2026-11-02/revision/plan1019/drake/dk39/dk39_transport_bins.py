"""dk39: altimetric surface transport per unit depth per 25-km bin along the mean LM Gould track, passage coordinates
(distance from the northern end, as GV23 / dk5), Product 1: sum of u_n ds over the track points in the bin (m2 s-1),
NaN in months when any point of the bin is ice-covered (ice fraction > 0.10). Trends (calibrated Theil-Sen, n_outer 150):
2001-2025 monthly where complete, else annual means (>= 6 months); Oct 2005-Apr 2019 annual means (>= 3 months), as dk5.
Also the latitude of each bin centre on the track, to plot in situ (km) bins on the latitude axis.
Usage: python dk39_transport_bins.py <revision dir> <dk24 npz> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, SEC, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
Z = np.load(SEC); lo, la, dsx = Z["lon"], Z["lat"], Z["dsx"]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
un = Z["u"] * nx + Z["v"] * ny; ice = Z["ice"] > 0.10
o = np.argsort(-la); dist = np.zeros(la.size); dist[o] = np.r_[0, np.cumsum(dsx[o][:-1] + dsx[o][1:]) / 2] / 1000   # km from the north end
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - .5) / 12
ov = ((T.x >= 2005 + 9 / 12) & (T.x <= 2019 + 4 / 12)).values
def annual(years, vals, minn):
    df = pd.DataFrame(dict(year=years, v=vals)).dropna(); g = df.groupby("year").v; a = g.mean()[g.count() >= minn]
    return a.index.values.astype(float), a.values
rows = []
for b0 in np.arange(0, 850, 25):
    c = b0 + 12.5; s = (dist >= b0) & (dist < b0 + 25)
    if not s.any(): continue
    q = (un[:, s] * dsx[s]).sum(1); q[ice[:, s].any(1)] = np.nan
    latc = float(np.interp(c, dist[o], la[o]))
    out = dict(bin=c, lat=latc, mean=np.nanmean(q), months=int(np.isfinite(q).sum()))
    r = R.trend_bcz(T.x.values, q, month=T.month.values, n_outer=150) if np.isfinite(q).all() else None
    if r is None:
        xa, ya = annual(T.year.values, q, 6); r = R.trend_bcz(xa, ya, monthly=False, n_outer=150) if len(xa) >= 8 else None
    out.update({k + "_0125": (np.nan if r is None else (10 * r[k] if k != "p" else r[k])) for k in ("slope", "ci_lo", "ci_hi", "p")})
    xa, ya = annual(T.year[ov].values, q[ov], 3); r = R.trend_bcz(xa, ya, monthly=False, n_outer=150) if len(xa) >= 8 else None
    out.update({k + "_0519": (np.nan if r is None else (10 * r[k] if k != "p" else r[k])) for k in ("slope", "ci_lo", "ci_hi", "p")})
    rows.append(out); print(c, round(latc, 2), round(out["mean"]), flush=True)
pd.DataFrame(rows).to_csv(OUT / "dk39_transport_bins.csv", index=False); print("done")
