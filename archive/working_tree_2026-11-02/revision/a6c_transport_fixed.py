"""A6c: altimetric transport proxy with a FIXED section geometry (correction to a6_drake.py).

a6_drake.py excluded transect points month by month where ice_frac > 0.10. The resulting T_v2 is
correlated with the valid fraction of the section (r = 0.89), so its trend partly reflects the number
of ice-free points, not the flow. Here the set of points is fixed in time:
  T_perm : points with monthly ice_frac <= 0.10 in every month 2001-2025 (rule of a5_energy.py;
           v2 used an analogous permanent mask)
  T_all  : all points, no ice exclusion (under-ice values are mapped, less reliable)
for the v2 line and all variants of a6_drake.py. Calibrated trends 2001-2025 and 2005-2019 and the
date-matched agreement with the in-situ transects are recomputed.
Usage: python a6c_transport_fixed.py <revision/data dir> <repository outputs/trends dir> <out_dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R

data, trends_dir, outdir = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
R_E = 6.371e6
TR = {"v2": (-65.0, -55.0, -60.0, -63.0)}
for dl in (-1.0, -0.5, 0.5, 1.0):
    TR[f"v2_lon{dl:+.1f}"] = (-65.0 + dl, -55.0, -60.0 + dl, -63.0)
TR["v2_south+0.5"] = (-65.0, -55.0, -60.0, -62.5); TR["v2_north-0.5"] = (-65.0, -55.5, -60.0, -63.0)
TR["gould_approx"] = (-65.0, -55.0, -58.5, -62.0)


def geometry(lon0, lat0, lon1, lat1, n=121):
    tt = np.linspace(0, 1, n); lo = lon0 + tt * (lon1 - lon0); la = lat0 + tt * (lat1 - lat0)
    latm = np.radians(0.5 * (la[:-1] + la[1:]))
    seg = np.hypot(R_E * np.cos(latm) * np.radians(np.diff(lo)), R_E * np.radians(np.diff(la)))
    ds = np.r_[seg[0] / 2, 0.5 * (seg[:-1] + seg[1:]), seg[-1] / 2]
    latr = np.radians(0.5 * (lat0 + lat1))
    DX = R_E * np.cos(latr) * np.radians(lon1 - lon0); DY = R_E * np.radians(lat1 - lat0); L = np.hypot(DX, DY)
    nx, ny = DY / L, -DX / L
    if nx < 0:
        nx, ny = -nx, -ny
    return lo, la, ds, nx, ny


files = [data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
VN, ICE = {k: [] for k in TR}, {k: [] for k in TR}
G = {k: geometry(*v) for k, v in TR.items()}
for f in files:
    with xr.open_dataset(f) as ds:
        d = ds.sel(longitude=slice(-70, -55), latitude=slice(-65, -53)).load()
    for k, (lo, la, dsx, nx, ny) in G.items():
        pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
        VN[k].append(d.u_mean.interp(**pts).values * nx + d.v_mean.interp(**pts).values * ny)
        ICE[k].append(d.ice_frac.interp(**pts).values)
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
t["x"] = t.year + (t.month - 0.5) / 12
info = []
for k in TR:
    vn = np.array(VN[k]); ice = np.nan_to_num(np.array(ICE[k])); dsx = G[k][2]
    perm = (ice.max(0) <= 0.10) & np.all(np.isfinite(vn), 0)
    allp = np.all(np.isfinite(vn), 0)
    t[f"Tperm_{k}"] = np.nansum(np.where(perm[None, :], vn, 0) * dsx, 1)
    t[f"Tall_{k}"] = np.nansum(np.where(allp[None, :], vn, 0) * dsx, 1)
    info.append(dict(line=k, n_points=vn.shape[1], n_perm=int(perm.sum()), n_all=int(allp.sum()),
                     perm_southern_lat=float(G[k][1][perm].min()) if perm.any() else np.nan))
pd.DataFrame(info).to_csv(outdir / "a6c_section_points.csv", index=False)
t.to_csv(outdir / "a6c_transport_fixed_series.csv", index=False)
print(pd.DataFrame(info).to_string())
rows = []
for c in [c for c in t.columns if c.startswith("T")]:
    for per, sel in (("2001-2025", np.ones(len(t), bool)), ("2005-2019", (t.year >= 2005) & (t.year <= 2019))):
        o = R.trend_calibrated(t.x.values[sel], t[c].values[sel], month=t.month.values[sel], n_outer=100, n_boot=300)
        y = t[c].values[sel]
        o.update(series=c, period=per, mean=float(y.mean()), pct_per_decade=1000 * o["slope"] / float(y.mean()),
                 sigma_per_yr=o["slope"] / float(np.std(y)))
        rows.append(o); print(c, per, round(o["pct_per_decade"], 2), round(o["p_cal"], 3), flush=True)
pd.DataFrame(rows).to_csv(outdir / "a6c_trends.csv", index=False)
obs = pd.read_csv(f"{trends_dir}/drake_passage_total_transport_observations.csv", parse_dates=["date"]).dropna(subset=["total_transport_sv"])
obs["year"] = obs.date.dt.year; obs["month"] = obs.date.dt.month
mm = obs.merge(t, on=["year", "month"], how="left")
ag = {}
for c in ("Tperm_v2", "Tall_v2", "Tperm_gould_approx", "Tall_gould_approx"):
    g = mm.groupby("year")[[c, "total_transport_sv"]].mean()
    ag[c] = dict(n=len(mm), r_transects=float(np.corrcoef(mm[c], mm.total_transport_sv)[0, 1]),
                 r_annual=float(np.corrcoef(g[c], g.total_transport_sv)[0, 1]))
pd.DataFrame(ag).T.to_csv(outdir / "a6c_agreement.csv")
print(pd.DataFrame(ag).T.round(3))
