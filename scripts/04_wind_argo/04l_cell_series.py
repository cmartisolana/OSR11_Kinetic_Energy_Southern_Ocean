"""x13 (10 Nov 2026): monthly 1-deg cell fields of wind stress, KE, EKE and the sigma0 = 27.3 cross-stream slope, as in 04g,
plus the distance of each 1-deg cell centre from the 1/8-deg ACC band (as 04c), so that domain means can be formed for
alternative buffers (sensitivity of Sect. 3.3 / Fig. 4 to the 100 km extension). Vectorised binning (bincount), same definitions.
Outputs: x13/x13_cells.npz (WC, KC, EC [300,30,360], SL [264,30,360], dmin, band, ice, lat, lon, year, month)

Origin  : x13_cell_series.py (10 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : 1-degree cell series for the domain-width test of 04m and for the Product 4 wind domain (06_product4/p4_04).
Inputs  : raw monthly files, wind file, outputs/f4/f4_tilt_maps.nc, outputs/f4m/f4m_p273_monthly.npz (04e cache)
Outputs : outputs/x13/x13_cells.npz (31 MB, distributed)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, xarray as xr, pandas as pd
from scipy.spatial import cKDTree
from osr11 import paths as PA
OUT = PA.out("x13")
M = xr.open_dataset(PA.out("f4") / "f4_tilt_maps.nc"); lat, lon = M.lat.values, M.lon.values
LON, LAT = np.meshgrid(lon, lat); KM = 111.195
la_e = np.r_[lat - 0.5, lat[-1] + 0.5]; lo_e = np.r_[lon - 0.5, lon[-1] + 0.5]; NC = lat.size * lon.size
# distance from the 1/8-deg band (04c)
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
madt = xr.open_dataset(PA.MEAN_ADT).mean_adt; A8 = madt.values; b8 = np.isfinite(A8) & (A8 >= -0.6) & (A8 <= 0.2)
L8, B8 = np.meshgrid(madt.longitude.values, madt.latitude.values)
dch, _ = cKDTree(xyz(B8[b8], L8[b8])).query(xyz(LAT.ravel(), LON.ravel()))
dmin = (2 * 6371.0 * np.arcsin(np.clip(dch / 2, 0, 1))).reshape(LAT.shape)
print("buffer100 check:", int((dmin <= 100 + 1e-6).sum()), "vs stored", int(M.buffer100.sum()), flush=True)
# slope (04g)
C = np.load(PA.out("f4m") / "f4m_p273_monthly.npz"); PM, lonr = C["p"], C["lon"]
L0 = np.where(lonr > 180, lonr - 360, lonr); o = np.argsort(L0); adt = M["adt"].values
def grad(F):
    dy = np.gradient(F, lat, axis=-2) / KM * 100
    dx = (np.roll(F, -1, axis=-1) - np.roll(F, 1, axis=-1)) / 2.0 / (KM * np.cos(np.radians(lat))[:, None]) * 100
    return dx, dy
ax_, ay_ = grad(adt); an = np.hypot(ax_, ay_); dx, dy = grad(PM[:, :, o].astype(float)); SL = (dx * ax_ / an + dy * ay_ / an).astype(np.float32)
def binner(la, lo):
    i = np.digitize(la, la_e) - 1; j = np.digitize(lo, lo_e) - 1
    I, J = np.meshgrid(i, j, indexing="ij"); ok = (I >= 0) & (I < lat.size) & (J >= 0) & (J < lon.size)
    idx = np.where(ok, I * lon.size + J, NC); return idx.ravel()
def cellmean(F, idx):
    f = F.ravel(); g = np.isfinite(f); s = np.bincount(idx[g], f[g], NC + 1)[:NC]; n = np.bincount(idx[g], None, NC + 1)[:NC]
    with np.errstate(invalid="ignore"): return (s / np.where(n > 0, n, np.nan)).reshape(lat.size, lon.size)
W = xr.open_dataset(PA.WIND).eastward_stress.sel(latitude=slice(lat[0] - 0.5, lat[-1] + 0.5))
t = pd.DatetimeIndex(W.time.values); iw = binner(W.latitude.values, W.longitude.values)
WC = np.stack([cellmean(W.values[k], iw) for k in range(t.size)]).astype(np.float32)
files = PA.monthly_files("allsat"); assert len(files) == 300
with xr.open_dataset(files[0]) as d0: la8, lo8 = d0.latitude.values, d0.longitude.values
ice8 = xr.open_dataset(PA.ICE_MASK).ice_mask.values > 0.5; i8 = binner(la8, lo8)
su = 0; sv = 0
for f in files:
    with xr.open_dataset(f) as d: su = su + d.u_mean.values.astype(float); sv = sv + d.v_mean.values.astype(float)
ub, vb = su / 300, sv / 300
EC = np.empty((300, lat.size, lon.size), np.float32); KC = np.empty_like(EC)
for k, f in enumerate(files):
    with xr.open_dataset(f) as d: u, v, vi = d.u_mean.values.astype(float), d.v_mean.values.astype(float), d.var_intra.values.astype(float)
    e8 = (0.5 * ((u - ub) ** 2 + (v - vb) ** 2) + 0.5 * vi) * 1e4; k8 = 0.5 * (u ** 2 + v ** 2) * 1e4
    e8[ice8] = np.nan; k8[ice8] = np.nan; EC[k] = cellmean(e8, i8); KC[k] = cellmean(k8, i8)
    if k % 50 == 0: print("month", k, flush=True)
np.savez_compressed(OUT / "x13_cells.npz", WC=WC, KC=KC, EC=EC, SL=SL, dmin=dmin, band=M.band.values, ice=M.ice.values,
                    slope_ok=np.isfinite(M["signed_27p3_pct"].values), lat=lat, lon=lon, year=t.year.values, month=t.month.values)
print("done", flush=True)
