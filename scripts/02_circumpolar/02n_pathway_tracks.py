"""x6 (5 Nov 2026): half-power tracking with a longitude-dependent threshold (the 'relthr' set of halfpower.SETS: 20 % of each
longitude's 2001-2025 maximum KE instead of 200 cm2 s-2; envelope to 25 % of the peak, >= 100 km), run on the monthly KE
= 0.5 (u^2 + v^2) of the monthly-mean geostrophic velocity (data/monthly_allsat). Climatological envelopes from the mean KE field
with the same rule; monthly envelopes matched as in halfpower.matched_tracks. For every climatological envelope in the ACC band the
monthly half-power latitude and the monthly latitude of the fixed ADT contour through its mean half-power point (nearest crossing
within 3 deg). Output: outputs/x6/x6_relthr_tracks.npz (lat_hp, lat_c: 300 x ntracks; meta columns).

Origin  : x6_relthr_tracks.py (5 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : tracks of the energetic pathway with the longitude-dependent threshold used in the paper (Sect. 2.2; 1826 retained tracks).
Inputs  : $OSR11_DATA/monthly_allsat (raw data), outputs/base/mean_adt.nc, clim_allsat.nc
Outputs : outputs/x6/x6_relthr_tracks.npz (distributed)
Run     : python scripts/02_circumpolar/02n_pathway_tracks.py   (raw data; about 6 GB of memory)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, xarray as xr, time
from osr11 import halfpower as HPm, paths as PA
files = PA.monthly_files("allsat"); assert all(f.exists() for f in files), "monthly files missing (scripts/00_extract/00a)"; OUT = PA.out("x6")
d0 = xr.open_dataset(files[0]); lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float); ny, nx = len(lat), len(lon)
adt = xr.open_dataset(PA.MEAN_ADT)["mean_adt"].values; acc = np.isfinite(adt) & (adt >= -0.6) & (adt <= 0.2)
ke = np.empty((300, ny, nx), np.float32); adt_m = np.empty((300, ny, nx), np.float16)
t0 = time.time()
for m, f in enumerate(files):
    ds = xr.open_dataset(f); u = ds.u_mean.values; v = ds.v_mean.values
    ke[m] = 0.5 * (u ** 2 + v ** 2); adt_m[m] = ds.adt_mean.values.astype(np.float16)
print("loaded", round(time.time() - t0), flush=True)
kemean = np.nanmean(ke, 0); clim_max = np.nanmax(np.where(acc, kemean, np.nan), 0); clim_max = np.where(np.isfinite(clim_max) & (clim_max > 0), clim_max, np.nan)
p = HPm.SETS["relthr"]
ccounts, ctab = HPm.clim_envelopes(kemean, lat, lon, p, clim_max)
counts, tab = HPm.run_set(ke, lat, lon, p, clim_max=clim_max)
print("detected", round(time.time() - t0), flush=True)
arr, meta = HPm.matched_tracks(counts, tab, ccounts, ctab, lat, lon, acc, min_cov=0.9)
meta["lon"] = lon[meta.j.values]
print("tracks", len(meta), "acc", int(meta.acc.sum()), "keep90", int(meta.keep.sum()), "cov>=0.75", int((meta.acc & (meta.coverage >= 0.75)).sum()), flush=True)
# fixed ADT contour through the mean half-power point of each ACC envelope
A = xr.open_dataset(PA.CLIM).adt.values
sel = np.nonzero(meta.acc.values)[0]; lat_c = np.full((300, len(meta)), np.nan, np.float32); lev = np.full(len(meta), np.nan)
for i in sel:
    j = int(meta.j[i]); lc = float(meta.lat_hp_clim[i])
    if not np.isfinite(A[:, j]).all(): continue
    L = np.interp(lc, lat, A[:, j]); lev[i] = L
    i0 = int(np.argmin(np.abs(lat - lc))); lo, hi = max(i0 - 24, 0), min(i0 + 25, ny)
    seg = adt_m[:, lo:hi, j].astype(np.float32); la = lat[lo:hi]; s = np.sign(seg - L); c = s[:, :-1] * s[:, 1:] < 0
    for m in range(300):
        k = np.nonzero(c[m])[0]
        if len(k) == 0: continue
        xs = la[k] + (L - seg[m, k]) * (la[k + 1] - la[k]) / (seg[m, k + 1] - seg[m, k]); lat_c[m, i] = xs[np.argmin(np.abs(xs - lc))]
np.savez_compressed(OUT / "x6_relthr_tracks.npz", lat_hp=arr["lat_hp"], peak_val=arr["peak_val"], lat_c=lat_c, level=lev,
                    j=meta.j.values, lon=meta.lon.values, lat_hp_clim=meta.lat_hp_clim.values, lat_lo=meta.lat_lo.values, lat_hi=meta.lat_hi.values,
                    acc=meta.acc.values, coverage=meta.coverage.values, keep=meta.keep.values, clim_max=clim_max)
print("saved", round(time.time() - t0))
