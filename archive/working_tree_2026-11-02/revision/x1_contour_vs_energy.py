"""x1 (29 Oct 2026, exploratory): contour-vs-energy test of frontal displacement on the 900 retained half-power tracks.
For each retained track (outputs/a2/tracks_baseline.npz, keep = True) and month:
  lat_hp      half-power latitude (tracked, from the npz)
  adt_hp      monthly ADT (adt_mean) interpolated to the half-power point  -> does the energy peak stay on one streamline?
  lat_c       latitude of the climatological ADT contour through the mean half-power point (ADT level = mean ADT at
              lat_hp_clim), followed each month in the same longitude within +-3 deg of the mean position
              -> what a fixed-height-contour front would do
Writes outputs/x1/x1_tracks.npz with the three (300 x N) arrays plus metadata.
Usage: python x1_contour_vs_energy.py <revision dir>"""
import sys, glob
from pathlib import Path
import numpy as np, xarray as xr
REV = Path(sys.argv[1]); OUT = REV / "outputs" / "x1"; OUT.mkdir(parents=True, exist_ok=True)
T = np.load(REV / "outputs" / "a2" / "tracks_baseline.npz")
keep = T["keep"]; j = T["j"][keep]; lat_hp = T["lat_hp"][:, keep].astype(float); lhc = T["lat_hp_clim"][keep]
print("tracks", keep.sum())
files = sorted(glob.glob(str(REV / "data" / "monthly_allsat" / "*" / "*.nc")))
files = [f for f in files if 2001 <= int(Path(f).parent.name) <= 2025]; assert len(files) == 300
d0 = xr.open_dataset(files[0]); lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float)
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc")["mean_adt"]
madt = madt.interp(latitude=lat, longitude=lon).values if madt.shape != (lat.size, lon.size) else madt.values
# level of the mean streamline through each track's climatological half-power point
lev = np.array([np.interp(lhc[k], lat, madt[:, j[k]]) for k in range(j.size)])
W = 3.0   # +- window (deg) around the mean half-power latitude for the contour search
adt_hp = np.full(lat_hp.shape, np.nan); lat_c = np.full(lat_hp.shape, np.nan)
for m, f in enumerate(files):
    with xr.open_dataset(f) as d:
        A = d.adt_mean.values.astype(float)
    for k in range(j.size):
        col = A[:, j[k]]
        if np.isfinite(lat_hp[m, k]):
            adt_hp[m, k] = np.interp(lat_hp[m, k], lat, col) if np.isfinite(np.interp(lat_hp[m, k], lat, col)) else np.nan
        # contour crossing closest to the mean position within the window
        w = (lat >= lhc[k] - W) & (lat <= lhc[k] + W); la = lat[w]; c = col[w] - lev[k]
        ok = np.isfinite(c)
        if ok.sum() < 3: continue
        la, c = la[ok], c[ok]; s = np.where(np.diff(np.sign(c)) != 0)[0]
        if s.size == 0: continue
        x = la[s] + (0 - c[s]) * (la[s + 1] - la[s]) / (c[s + 1] - c[s])
        lat_c[m, k] = x[np.argmin(np.abs(x - lhc[k]))]
    if m % 50 == 0: print(m, flush=True)
np.savez_compressed(OUT / "x1_tracks.npz", lat_hp=lat_hp, adt_hp=adt_hp, lat_c=lat_c, lon=lon[j], lat_hp_clim=lhc, lev=lev, j=j)
print("done")
