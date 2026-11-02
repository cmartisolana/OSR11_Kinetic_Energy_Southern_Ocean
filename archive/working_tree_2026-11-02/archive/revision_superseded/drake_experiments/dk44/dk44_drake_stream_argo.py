"""dk44: Argo baroclinic eastward velocity across Drake Passage per pair of SSH streamlines (GV23 pairs mapped onto the mean ADT,
dk5 stream edges), 2004-2025. Cells of the RG grid at 65.5-59.5 W with complete profiles (dk42 u_surf, u_mean: relative to
1975 dbar) are assigned to pairs by the 2001-2025 mean ADT on the RG grid (dk32); pair mean per year; trend: calibrated
Theil-Sen on annual values (cm/s per decade, 95 %). Pairs with < 3 cells are skipped.
Usage: python dk44_drake_stream_argo.py <revision dir> <drake dir> <dk42 npz> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, DR, F42, OUT = (Path(a) for a in sys.argv[1:5]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
Z = np.load(F42); lon, lat, yrs, adt = Z["lon"], Z["lat"], Z["years"].astype(float), Z["adt"]
ED = pd.read_csv(DR / "out" / "dk5" / "dk5_stream_edges.csv").sort_values("psi"); A24 = pd.read_csv(DR / "out" / "dk24" / "dk24_stream_trends.csv")
ilon = (lon >= 294) & (lon <= 301); rows = []
for (p0, e0), (p1, e1) in zip(ED[["psi", "adt_edge"]].values[:-1], ED[["psi", "adt_edge"]].values[1:]):
    lo_, hi_ = sorted([e0, e1]); c = (adt >= lo_) & (adt < hi_) & ilon[None, :] & np.isfinite(Z["u_mean"]).all(0)
    if c.sum() < 3: continue
    psic = 0.5 * (p0 + p1); adtc = 0.5 * (e0 + e1)
    for k in ("u_surf", "u_mean"):
        y = Z[k][:, c].mean(1); r = R.trend_bcz(yrs, y, monthly=False, n_outer=150)
        rows.append(dict(psic=round(psic, 4), adt_c=adtc, var=k, ncells=int(c.sum()), mean=y.mean(), slope_dec=10 * r["slope"], ci_lo=10 * r["ci_lo"], ci_hi=10 * r["ci_hi"], p=r["p"]))
D = pd.DataFrame(rows); D.to_csv(OUT / "dk44_drake_stream_argo.csv", index=False); pd.set_option("display.width", 200); print(D.round(3).to_string(index=False))
