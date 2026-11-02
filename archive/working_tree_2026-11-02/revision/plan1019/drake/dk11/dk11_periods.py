"""dk11: Figure 1 quantities for two periods, 2001-2015 and 2016-2025 (Product 1).
(1) EKE trends in 25 km bins along the LM Gould transect (passage coordinates); (2) surface transport trends between pairs of
SSH streamlines (GV23 frontal zones); monthly series (deseasonalised) when complete, otherwise annual means (>= 6 months);
calibrated test (trend_bcz, 150 null series). Inputs: dk5_section_fields.npz, dk5_stream_edges.csv, dk5_geometry.py."""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[3])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, numpy as np, pandas as pd
from pathlib import Path
REV = _REV; sys.path.insert(0, REV); import rstats as R
OUT = Path(_DK_OUT); GV = Path(_GV); exec(open((_DK + "/dk5_geometry.py")).read())
Z = np.load(OUT / "dk5" / "dk5_section_fields.npz"); E = pd.read_csv(OUT / "dk5" / "dk5_stream_edges.csv"); psi2adt = dict(zip(E.psi.round(3), E.adt_edge))
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - .5) / 12
PER = {"2001-2015": (2001, 2015), "2016-2025": (2016, 2025)}
def tr(y, sel):
    x, m, yy = T.x.values[sel], T.month.values[sel], y[sel]
    if np.isfinite(yy).all(): o = R.trend_bcz(x, yy, month=m, monthly=True, n_outer=150)
    else:
        df = pd.DataFrame(dict(year=T.year.values[sel], v=yy)).dropna(); g = df.groupby("year").v; a = g.mean()[g.count() >= 6]
        if len(a) < 6: return dict(slope=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan)
        o = R.trend_bcz(a.index.values.astype(float), a.values, monthly=False, n_outer=150)
    return o
ice = Z["allsat_ice"] > 0.10; ek = np.where(ice, np.nan, Z["allsat_eke"]) * 1e4; un = Z["allsat_un"]; adt = Z["allsat_adt"]
rows = []
for b0 in np.arange(0, 850, 25):
    s = (dkm >= b0) & (dkm < b0 + 25)
    if not s.any(): continue
    y = np.nanmean(ek[:, s], 1)
    for per, (a, b) in PER.items():
        sel = ((T.year >= a) & (T.year <= b)).values; o = tr(y, sel); mu = np.nanmean(y[sel])
        rows.append(dict(kind="eke_bin", key=b0 + 12.5, period=per, mean=mu, slope_dec=10 * o["slope"], pct=1000 * o["slope"] / mu, pct_lo=1000 * o["ci_lo"] / mu, pct_hi=1000 * o["ci_hi"] / mu, p=o["p"]))
SZ = {"SAF": (-0.40, -0.60), "SAF–PF": (-0.60, -1.00), "PF": (-1.00, -1.25)}
for name, (a, b) in SZ.items():
    lo_e, hi_e = sorted((psi2adt[round(b, 3)], psi2adt[round(a, 3)])); inb = (adt >= lo_e) & (adt < hi_e); icy = (inb & ice).any(1)
    q = np.where(icy, np.nan, np.nansum(np.where(inb, un * dsx[None, :], 0), 1))
    for per, (a2, b2) in PER.items():
        sel = ((T.year >= a2) & (T.year <= b2)).values; o = tr(q, sel); mu = np.nanmean(q[sel])
        rows.append(dict(kind="transport_zone", key=name, period=per, mean=mu, slope_dec=10 * o["slope"], pct=1000 * o["slope"] / mu, pct_lo=1000 * o["ci_lo"] / mu, pct_hi=1000 * o["ci_hi"] / mu, p=o["p"]))
    print(name, flush=True)
pd.DataFrame(rows).to_csv(OUT / "dk11" / "dk11_periods.csv", index=False); print("ok")
