"""Coverage diagnostics for the Antarctic-slope index region (same region definition as h0_antarctic_slope_index.py):
fraction of months each region point is ice-free (ice_frac <= 0.10 and finite ADT), and the index from the fixed subset of
points that are ice-free in at least F of months. Usage: python h0_coverage_check.py <monthly dir> <bathy> <out dir>"""
import sys, numpy as np, pandas as pd, xarray as xr
from pathlib import Path
from scipy import ndimage
mdir, bat_f, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
files = [mdir / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values, d0.longitude.values
dep = -xr.open_dataset(bat_f).elevation.interp(latitude=xr.DataArray(lat, dims="latitude"), longitude=xr.DataArray(lon, dims="longitude")).values
dep = np.nan_to_num(dep, nan=-1.0)
reg = (lat[:, None] <= -60) & (dep >= 1000) & (dep <= 3000)
shallow = (dep < 3000) & (lat[:, None] <= -55); lab, n = ndimage.label(shallow, structure=np.ones((3, 3), int))
parent = np.arange(n + 1)
def find(a):
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
for i in range(lab.shape[0]):
    for di in (-1, 0, 1):
        k = i + di
        if 0 <= k < lab.shape[0] and lab[i, 0] and lab[k, -1]: parent[find(lab[i, 0])] = find(lab[k, -1])
lab = np.array([find(a) for a in range(n + 1)])[lab]
keep = np.unique(lab[0][lab[0] > 0]); reg &= np.isin(lab, keep) & (lab > 0)
ii, jj = np.nonzero(reg); A = np.empty((300, ii.size)); OK = np.empty((300, ii.size), bool)
for k, f in enumerate(files):
    with xr.open_dataset(f) as d:
        a = d.adt_mean.values[ii, jj]; c = np.nan_to_num(d.ice_frac.values[ii, jj]); A[k] = a; OK[k] = np.isfinite(a) & (c <= 0.10)
fr = OK.mean(0); w = np.cos(np.radians(lat[ii]))
rows = []
for F in (1.0, 0.95, 0.9, 0.75):
    sel = fr >= F
    rows.append(dict(F=F, n_points=int(sel.sum()), share=float(sel.mean()), lon_min=float(lon[jj][sel].min()) if sel.any() else np.nan,
                     lon_max=float(lon[jj][sel].max()) if sel.any() else np.nan,
                     n_30deg_bands=int(np.unique((np.mod(lon[jj][sel], 360) // 30)).size) if sel.any() else 0))
pd.DataFrame(rows).to_csv(out / "h0_coverage_subsets.csv", index=False); print(pd.DataFrame(rows).to_string(index=False))
sel = fr >= 0.95
B = pd.DataFrame(dict(lon=lon[jj], lat=lat[ii], frac_icefree=fr)); B.to_csv(out / "h0_region_points.csv", index=False)
print("points ice-free >= 95 % by 30-deg band:", B[B.frac_icefree >= 0.95].groupby((np.mod(B.lon, 360) // 30 * 30)).size().to_dict())
X = np.where(OK[:, sel], A[:, sel], np.nan); X = X - np.nanmean(X, 0)
h = 100 * np.nansum(X * w[sel], 1) / np.sum(np.isfinite(X) * w[sel], 1)
pd.DataFrame(dict(year=np.repeat(np.arange(2001, 2026), 12), month=np.tile(np.arange(1, 13), 25), h_fixed95_cm=h,
                  n_valid=np.isfinite(X).sum(1))).to_csv(out / "h0_fixed95_allsat.csv", index=False)
