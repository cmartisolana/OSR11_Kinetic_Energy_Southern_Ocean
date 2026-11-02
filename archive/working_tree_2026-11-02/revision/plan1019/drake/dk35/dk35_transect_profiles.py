"""dk35: EKE and surface transport profiles along the mean LM Gould track by 0.25 deg latitude (as dk25 for the velocity):
dk24 section fields (Product 1), points with ice fraction > 0.10 masked per month.
EKE (Eq. 1, cm2 s-2) bin mean; transport per unit depth accumulated from the northern end of the track to each latitude,
T(lat) = sum over points north of lat of u_n ds (m2 s-1), defined only in months when all those points are ice-free.
Trends 2001-2025: calibrated Theil-Sen (rstats.trend_bcz), monthly where the series is complete, otherwise annual means
(>= 6 months), as dk25; reported per decade and in % of the 2001-2025 mean.
Usage: python dk35_transect_profiles.py <revision dir> <dk24 npz> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, SEC, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
Z = np.load(SEC); lo, la, dsx = Z["lon"], Z["lat"], Z["dsx"]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
ice = Z["ice"] > 0.10; u, v = Z["u"], Z["v"]
un = u * nx + v * ny
EKE = (0.5 * ((u - np.nanmean(u, 0)) ** 2 + (v - np.nanmean(v, 0)) ** 2) + 0.5 * Z["vi"]) * 1e4; EKE[ice] = np.nan
q = un * dsx; q[ice] = np.nan
o = np.argsort(-la)                                    # north -> south
cum = np.cumsum(q[:, o], 1); cum[~np.isfinite(q[:, o]).cumprod(1).astype(bool)] = np.nan
TR = np.full_like(q, np.nan); TR[:, o] = cum
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - .5) / 12
def annual(years, vals, minn):
    df = pd.DataFrame(dict(year=years, v=vals)).dropna(); g = df.groupby("year").v; a = g.mean()[g.count() >= minn]
    return a.index.values.astype(float), a.values
rows = []
for c in np.arange(-63.5, -54.99, 0.25):
    s = np.abs(la - c) <= 0.125
    if not s.any(): continue
    for var, F in (("eke", EKE), ("transport", TR)):
        y = np.nanmean(F[:, s], 1) if var == "eke" else np.nanmean(F[:, s], 1)
        out = dict(latb=c, var=var, mean=np.nanmean(y), months=int(np.isfinite(y).sum()))
        if np.isfinite(y).all(): r = R.trend_bcz(T.x.values, y, month=T.month.values)
        else:
            xa, ya = annual(T.year.values, y, 6); r = R.trend_bcz(xa, ya, monthly=False) if len(xa) >= 8 else None
        for k in ("slope", "ci_lo", "ci_hi", "p"):
            out[k] = np.nan if r is None else (10 * r[k] if k != "p" else r[k])
        for k in ("slope", "ci_lo", "ci_hi"): out[k + "_pct"] = 100 * out[k] / out["mean"] if np.isfinite(out[k]) else np.nan
        rows.append(out); print(c, var, round(out["mean"], 1), round(out["slope_pct"], 1) if np.isfinite(out["slope_pct"]) else None, flush=True)
pd.DataFrame(rows).to_csv(OUT / "dk35_transect_profiles.csv", index=False); print("done")
