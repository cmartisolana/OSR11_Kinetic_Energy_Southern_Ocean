"""x1 (29 Oct 2026, exploratory): contour-vs-energy test of frontal displacement on the 900 retained half-power tracks.
For each retained track (outputs/a2/tracks_baseline.npz, keep = True) and month:
  lat_hp      half-power latitude (tracked, from the npz)
  adt_hp      monthly ADT (adt_mean) interpolated to the half-power point  -> does the energy peak stay on one streamline?
  lat_c       latitude of the climatological ADT contour through the mean half-power point (ADT level = mean ADT at
              lat_hp_clim), followed each month in the same longitude within +-3 deg of the mean position
              -> what a fixed-height-contour front would do
Writes outputs/x1/x1_tracks.npz with the three (300 x N) arrays plus metadata.

Purpose : contour-vs-energy test, extraction step (SI S3, Fig. S7; Sect. 3.2 and Discussion: -11.5 km per decade for the
          contour, no significant half-power displacement). Also writes the monthly ACC-band mean ADT drawn in Fig. S7b
          (added for the publication repository: cos-latitude weighted mean of adt_mean over the pixels whose 2001-2025
          mean ADT lies in -0.6..+0.2 m; the working repository held this series without a script. Identical to the
          figure agent's reconstruction with the ice mask applied, because ice-masked pixels have no mean ADT in the band).
Inputs  : outputs/a2/tracks_baseline.npz (02c), $OSR11_DATA/monthly_allsat (00a), outputs/base/mean_adt.nc
Outputs : outputs/x1/x1_tracks.npz (lat_hp, adt_hp, lat_c, ... on the 900 tracks), outputs/x1/x1_band_adt.npy (300 values, m)
Run     : python scripts/02_circumpolar/02j_contour_vs_energy.py
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, xarray as xr
from osr11 import paths as PA
OUT = PA.out("x1")
T = np.load(PA.OUT / "a2" / "tracks_baseline.npz")
keep = T["keep"]; j = T["j"][keep]; lat_hp = T["lat_hp"][:, keep].astype(float); lhc = T["lat_hp_clim"][keep]
print("tracks", keep.sum())
files = PA.monthly_files("allsat"); assert all(f.exists() for f in files) and len(files) == 300
d0 = xr.open_dataset(files[0]); lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float)
madt = xr.open_dataset(PA.MEAN_ADT)["mean_adt"]
madt = madt.interp(latitude=lat, longitude=lon).values if madt.shape != (lat.size, lon.size) else madt.values
# level of the mean streamline through each track's climatological half-power point
lev = np.array([np.interp(lhc[k], lat, madt[:, j[k]]) for k in range(j.size)])
W = 3.0   # +- window (deg) around the mean half-power latitude for the contour search
adt_hp = np.full(lat_hp.shape, np.nan); lat_c = np.full(lat_hp.shape, np.nan)
# ACC-band mean ADT (Fig. S7b): band from the 2001-2025 mean ADT, cos-lat weights, no ice mask
band = np.isfinite(madt) & (madt >= -0.6) & (madt <= 0.2); wb = np.cos(np.radians(lat))[:, None] * band; band_adt = np.full(len(files), np.nan)
for m, f in enumerate(files):
    with xr.open_dataset(f) as d:
        A = d.adt_mean.values.astype(float)
    okb = np.isfinite(A) & band; band_adt[m] = np.sum(A[okb] * wb[okb]) / np.sum(wb[okb])
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
np.save(OUT / "x1_band_adt.npy", band_adt)
np.savez_compressed(OUT / "x1_tracks.npz", lat_hp=lat_hp, adt_hp=adt_hp, lat_c=lat_c, lon=lon[j], lat_hp_clim=lhc, lev=lev, j=j)
print("done")
