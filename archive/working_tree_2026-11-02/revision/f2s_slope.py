"""F2s: sigma0 = 27.3 cross-stream slope trends for Figure 2 (Drake Passage), from the monthly fields of f4m_monthly.py (28 Oct 2026).
 S1 per-cell trend map, 90-30 W, 66-46 S, within the Fig. 4 domain (ACC band + 100 km) (Theil-Sen on the deseasonalised monthly series, % of the cell mean per decade; cells
    with ice-masked altimetry pixels or with any undefined month removed). Descriptive, for contours on Fig. 2a; no test.
 S2 2 deg longitude bins of the band between the mean SAF and PF widened by 100 km (1 deg cells whose centre lies within 100 km,
    great circle, of a 1/8 deg pixel whose 2001-2025 mean ADT (clim_allsat) lies between the dk12 front levels; revised 28 Oct 2026
    at the author's request, previously the unbuffered band), centres as dk8 (85-35 W): monthly trend_bcz (n_outer 150), % per decade; cross-bin r with the dk8 EKE trend (%), circular shifts.
Usage: python f2s_slope.py <revision dir> <f4 out dir> <f4m out dir> <out dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

REV, F4, FM, OUT = (Path(a) for a in sys.argv[1:5]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402

C = np.load(FM / "f4m_p273_monthly.npz"); YY, MM, PM, lonr = C["year"], C["month"], C["p"], C["lon"]
x = YY + (MM - 0.5) / 12
M = xr.open_dataset(F4 / "f4_tilt_maps.nc"); glat, glon = M.lat.values, M.lon.values
L0 = np.where(lonr > 180, lonr - 360, lonr); o = np.argsort(L0)
KM = 111.195; adt = M["adt"].values
def grad(F):
    dy = np.gradient(F, glat, axis=-2) / KM * 100
    dx = (np.roll(F, -1, axis=-1) - np.roll(F, 1, axis=-1)) / 2.0 / (KM * np.cos(np.radians(glat))[:, None]) * 100
    return dx, dy
ax_, ay_ = grad(adt); an = np.hypot(ax_, ay_)
dx, dy = grad(PM[:, :, o].astype(float)); SL = dx * ax_ / an + dy * ay_ / an
ok = (M["ice"].values == 0) & np.isfinite(SL).all(0)
LONG, LATG = np.meshgrid(glon, glat); w = np.cos(np.radians(LATG))

# S1 map
box = ok & (M["buffer100"].values > 0) & (LONG >= -90) & (LONG <= -30) & (LATG >= -66) & (LATG <= -46)
PCT = np.full(ok.shape, np.nan)
for i, j in np.argwhere(box):
    y = R.deseasonalise(SL[:, i, j], MM); b, _ = R.theil_sen(x, y); PCT[i, j] = 1000 * b / y.mean()
xr.Dataset({"slope_pct_ts": (("lat", "lon"), PCT)}, coords={"lat": glat, "lon": glon}).to_netcdf(OUT / "f2s_slope_map.nc")

# S2 bins between the fronts
Z = np.load(REV / "plan1019" / "drake" / "out" / "dk12" / "dk12_band_series_allsat.npz"); LS, LP = float(Z["lev_saf"]), float(Z["lev_pf"])
AL = pd.read_csv(REV / "plan1019" / "drake" / "dk8" / "dk8_along_acc_fronts_allsat.csv")
from scipy.spatial import cKDTree
cl = xr.open_dataset(REV / "outputs" / "clim_allsat.nc").adt.sel(longitude=slice(-100, -20), latitude=slice(-70, -40))
CLA, CLO = np.meshgrid(cl.latitude.values, cl.longitude.values, indexing="ij"); pix = (cl.values >= LP) & (cl.values <= LS)
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
d, _ = cKDTree(xyz(CLA[pix], CLO[pix])).query(xyz(LATG.ravel(), LONG.ravel()))
near = (2 * np.arcsin(d / 2) * 6371.0 <= 100.0).reshape(LATG.shape)
frt = ok & near; rows = []
for c in AL.lon.values:
    m = frt & (LONG >= c - 1) & (LONG < c + 1)
    if m.sum() < 1: rows.append(dict(lon=c, ncell=0)); continue
    s = np.array([np.sum(SL[t][m] * w[m]) / np.sum(w[m]) for t in range(SL.shape[0])])
    r = R.trend_bcz(x, s, month=MM, n_outer=150); mu = s.mean()
    rows.append(dict(lon=c, ncell=int(m.sum()), mean=mu, pct=1000 * r["slope"] / mu, lo=1000 * r["ci_lo"] / mu, hi=1000 * r["ci_hi"] / mu, p=r["p"]))
    print(rows[-1], flush=True)
B = pd.DataFrame(rows).merge(AL[["lon", "pct"]].rename(columns={"pct": "eke_pct"}), on="lon"); B.to_csv(OUT / "f2s_bins.csv", index=False)
v = B.dropna(subset=["pct"]); a, b = v.pct.values, v.eke_pct.values; rr = np.corrcoef(a, b)[0, 1]
null = np.array([np.corrcoef(a, np.roll(b, k))[0, 1] for k in range(1, len(b))]); pp = (1 + np.sum(np.abs(null) >= abs(rr))) / (1 + null.size)
pd.DataFrame([dict(series="cross-bin r(slope trend, EKE trend), Drake 2 deg bins", r=rr, p=pp, n=len(v))]).to_csv(OUT / "f2s_results.csv", index=False)
print(f"r = {rr:+.2f}, p = {pp:.3f}, n = {len(v)}; bins p<0.05: {int((v.p < 0.05).sum())}", flush=True)
