"""B3 relations (pre-specified in b3_alongstream.py): spatial correlation across 5-degree bins between the EKE trend and
the speed, sinuosity and latitude trends; significance from 2000 circular shifts of one profile. Also: number of bins
with positive trend (binomial p, descriptive since bins are correlated), and mean trends in obstacle bins
(shallowest depth < 3000 m), in the 3 bins downstream of an obstacle, and elsewhere. Excludes the Malvinas bins.
Usage: python b3_relations.py"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
D = Path(__file__).resolve().parent / "outputs" / "b3"
C = pd.read_csv(D / "b3_climatology.csv")
T = {v: pd.read_csv(D / f"b3_trends_{v}.csv").set_index("lon") for v in ("eke", "speed", "sinu", "lat")}
lon = C.lon.values; keep = ~((lon > -60) & (lon < -50))
obst = (C.depth_min.values < 3000)
down = np.zeros_like(obst)
for i in np.nonzero(obst)[0]:
    for k in (1, 2, 3):
        j = (i + k) % len(lon)
        if not obst[j]: down[j] = True
rng = np.random.default_rng(20261015); out = []
for v in ("eke", "speed", "sinu", "lat"):
    y = T[v].tr.reindex(lon).values; ok = keep & np.isfinite(y)
    out.append(dict(test=f"{v}: bins positive", value=f"{int((y[ok] > 0).sum())}/{ok.sum()}", p=stats.binomtest(int((y[ok] > 0).sum()), int(ok.sum())).pvalue))
    for name, m in (("obstacle", obst), ("downstream 1-3 bins", down), ("elsewhere", ~obst & ~down)):
        mm = ok & m; out.append(dict(test=f"{v}: mean trend {name}", value=f"{np.nanmean(y[mm]):+.2f} (n={mm.sum()})", p=np.nan))
e = T["eke"].tr.reindex(lon).values
for v in ("speed", "sinu", "lat"):
    y = T[v].tr.reindex(lon).values; ok = keep & np.isfinite(e) & np.isfinite(y)
    r0 = np.corrcoef(e[ok], y[ok])[0, 1]
    rs = []
    for _ in range(2000):
        s = rng.integers(1, len(lon)); ys = np.roll(y, s); kk = ok & np.isfinite(ys) & keep
        rs.append(np.corrcoef(e[kk], ys[kk])[0, 1])
    p = (1 + np.sum(np.abs(rs) >= abs(r0))) / 2001
    out.append(dict(test=f"spatial r(EKE trend, {v} trend)", value=f"{r0:+.2f}", p=p))
# also: EKE trend vs mean EKE (where is EKE growing, in relative terms)
y = np.log(C.eke_mean.values); ok = keep & np.isfinite(e)
r0 = np.corrcoef(e[ok], y[ok])[0, 1]; rs = [np.corrcoef(e[ok], np.roll(y, rng.integers(1, 72))[ok])[0, 1] for _ in range(2000)]
out.append(dict(test="spatial r(EKE trend %, log mean EKE)", value=f"{r0:+.2f}", p=(1 + np.sum(np.abs(rs) >= abs(r0))) / 2001))
R = pd.DataFrame(out); R.to_csv(D / "b3_relations.csv", index=False); print(R.to_string(index=False))
