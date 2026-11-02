"""dk40: net eastward surface transport per unit depth across the LM Gould line (sum of u_n ds over all 121 track points,
Product 1), 2001-2025. Monthly value is NaN if any point is ice-covered (ice fraction > 0.10); annual means from >= 3
ice-free months. Trends: calibrated Theil-Sen on the annual series (2001-2025 and Oct 2005-Apr 2019), and on the complete
monthly series of the always-ice-free northern part (points ice-free in every month), for reference.
Usage: python dk40_net_transport.py <revision dir> <dk24 npz> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, SEC, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
Z = np.load(SEC); lo, la, dsx = Z["lon"], Z["lat"], Z["dsx"]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
un = Z["u"] * nx + Z["v"] * ny; ice = Z["ice"] > 0.10
q = un * dsx
full = q.sum(1); full[ice.any(1)] = np.nan
perm = ~ice.any(0); north = q[:, perm].sum(1)
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - .5) / 12
T["net_full"] = full; T["net_north"] = north; T.to_csv(OUT / "dk40_net_monthly.csv", index=False)
an = T.groupby("year").net_full.agg(["mean", "count"]); an = an[an["count"] >= 3]
an.to_csv(OUT / "dk40_net_annual.csv")
ov = ((T.x >= 2005 + 9 / 12) & (T.x <= 2019 + 4 / 12)).values
rows = []
r = R.trend_bcz(an.index.values.astype(float), an["mean"].values, monthly=False, n_outer=150)
rows.append(dict(series="net_full_annual", period="2001-2025", mean=an["mean"].mean(), n=len(an), **{k: r[k] for k in ("slope", "ci_lo", "ci_hi", "p")}))
a2 = T[ov].groupby("year").net_full.agg(["mean", "count"]); a2 = a2[a2["count"] >= 3]
r = R.trend_bcz(a2.index.values.astype(float), a2["mean"].values, monthly=False, n_outer=150)
rows.append(dict(series="net_full_annual", period="2005-2019", mean=a2["mean"].mean(), n=len(a2), **{k: r[k] for k in ("slope", "ci_lo", "ci_hi", "p")}))
r = R.trend_bcz(T.x.values, north, month=T.month.values, n_outer=150)
rows.append(dict(series="net_north_monthly", period="2001-2025", mean=north.mean(), n=300, lat_south=float(la[perm].min()), **{k: r[k] for k in ("slope", "ci_lo", "ci_hi", "p")}))
D = pd.DataFrame(rows); D["pct_dec"] = 1000 * D.slope / D["mean"]; D["pct_lo"] = 1000 * D.ci_lo / D["mean"]; D["pct_hi"] = 1000 * D.ci_hi / D["mean"]
D.to_csv(OUT / "dk40_net_trends.csv", index=False); pd.set_option("display.width", 200); print(D.round(3).to_string(index=False))
print("months per year with full line ice-free:", T.groupby("year").net_full.count().values)
