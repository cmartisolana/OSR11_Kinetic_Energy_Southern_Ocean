"""Step 8: half-power tracking of Product 4 with the RELATIVE threshold of the main text (20 % of each longitude's maximum 2001-2025
mean KE in the ACC band; set 'relthr' of osr11.halfpower.SETS), counterpart of the Product 1 run x6_relthr_tracks.py
(scripts/02_circumpolar/02n, 02p; 1826 retained tracks), plus the monthly latitude of the fixed ADT contour through each track's
mean half-power point. The same script with --thr baseline reruns the FIXED 200 cm2 s-2 tracking of 02e_halfpower_twosat.py, which
is the validation of the tracking stage against outputs/a2/a2h_twosat_vs_allsat.csv (400 tracks, p4_09 compares the trends).

Purpose : Product 4 relative-threshold tracks for the position comparison of Fig. S3 (Product 1 equivalent: outputs/x6, outputs/x7).
Inputs  : Product 4 monthly fields (p4_common: u_mean, v_mean, adt_mean, ice_frac; 300 x 120 x 1440, 1/4 deg);
          outputs/base/mean_adt.nc (Product 1 2001-2025 mean ADT, 1/8 deg: ACC band mask)
Outputs : outputs/p4/p4_08_tracks_<thr>_<ice>.npz  lat_hp, peak_val, lat_c (300 x ntracks), level, j, lon, lat_hp_clim, lat_lo, lat_hi, acc,
          coverage, keep, clim_max  (same keys as outputs/x6/x6_relthr_tracks.npz)
Method  : identical code path to x6 / 02e (osr11.halfpower.run_set, clim_envelopes, matched_tracks unchanged):
            KE = 0.5 (u_mean^2 + v_mean^2) of the monthly-mean geostrophic velocity, m2 s-2;
            ACC band = cells whose Product 1 mean ADT (mean_adt.nc, nearest-neighbour onto the 1/4 deg grid, exactly the 02e call)
              lies in [-0.6, +0.2] m;
            relthr: clim_max = maximum over the ACC band of the 2001-2025 mean KE at each longitude (x6 line), threshold 0.20 of it;
              baseline: 200 cm2 s-2 (02e; its clim_max = unrestricted column maximum, unused for a fixed threshold);
            envelope to 25 % of the peak, minimum width 100 km (osr11.fronts.meridional_dists: km on any grid, so the 1/8 deg rule
              carries over to 1/4 deg unchanged); monthly envelopes matched to the climatological envelopes of the mean KE field;
              retained = climatological half-power point in the ACC band and matched in >= 90 % of months;
            fixed contour: level = the product's OWN 2001-2025 mean ADT (unweighted mean of the 300 monthly adt_mean, as
              01d_climatology.py builds clim_allsat.nc for Product 1) at the climatological half-power latitude; monthly latitude =
              the crossing of that level by the monthly ADT nearest to it within +-3 deg (x6: +-24 cells of 1/8 deg; here +-12
              cells of 1/4 deg, i.e. lo = i0 - 12, hi = i0 + 13).
          Ice: x6 used the raw monthly fields of Product 1 (no ice masking; the monthly files carry no NaN under ice), so the
          primary run here does the same (--ice none). --ice perm applies the permanent ice mask of 01c (any JJA month with
          ice_frac > 0 in 2001-2025 -> NaN in every month), which is how the 02e input fields were built: used for the 02e
          validation and as a sensitivity.
Run     : python scripts/06_product4/p4_08_relthr_tracks.py --thr relthr --ice none     (primary)
          python scripts/06_product4/p4_08_relthr_tracks.py --thr baseline --ice perm   (validation against a2h)
          python scripts/06_product4/p4_08_relthr_tracks.py --thr relthr --ice perm     (sensitivity)
"""
import argparse
import time

import numpy as np
import xarray as xr

import p4_common as PC
from osr11 import halfpower as HPm

ap = argparse.ArgumentParser()
ap.add_argument("--thr", choices=("relthr", "baseline"), default="relthr")
ap.add_argument("--ice", choices=("none", "perm"), default="none")
a = ap.parse_args()
t0 = time.time()
lat, lon = PC.grid(); ny, nx = lat.size, lon.size
u = PC.stack("u_mean"); v = PC.stack("v_mean")
ke = np.empty((300, ny, nx), np.float32)
for m in range(300):
    ke[m] = 0.5 * (np.asarray(u[m], np.float32) ** 2 + np.asarray(v[m], np.float32) ** 2)
if a.ice == "perm":
    ice = PC.stack("ice_frac"); jja = np.isin(PC.MON, (6, 7, 8))
    perm = (np.nan_to_num(np.asarray(ice[jja]), nan=0.0) > 0).any(0)
    ke[:, perm] = np.nan
    print("permanent ice mask cells", int(perm.sum()), flush=True)
adt_m = PC.stack("adt_mean")
A4 = np.nanmean(np.asarray(adt_m, np.float64), 0)                     # Product 4 own 2001-2025 mean ADT
adt1 = xr.open_dataset(PC.PA.MEAN_ADT)["mean_adt"].interp(latitude=lat, longitude=lon, method="nearest").values
acc = np.isfinite(adt1) & (adt1 >= -0.6) & (adt1 <= 0.2)
print("loaded", round(time.time() - t0), "s; ACC cells", int(acc.sum()), flush=True)
kemean = np.nanmean(ke.astype(np.float64), 0)
if a.thr == "relthr":                                                  # x6
    clim_max = np.nanmax(np.where(acc, kemean, np.nan), 0)
    clim_max = np.where(np.isfinite(clim_max) & (clim_max > 0), clim_max, np.nan)
else:                                                                  # 02e
    clim_max = np.nanmax(kemean, 0)
p = HPm.SETS[a.thr]
ccounts, ctab = HPm.clim_envelopes(kemean, lat, lon, p, clim_max)
counts, tab = HPm.run_set(ke, lat, lon, p, clim_max=clim_max)
print("detected", round(time.time() - t0), "s", flush=True)
arr, meta = HPm.matched_tracks(counts, tab, ccounts, ctab, lat, lon, acc, min_cov=0.9)
meta["lon"] = lon[meta.j.values]
print("tracks", len(meta), "acc", int(meta.acc.sum()), "keep90", int(meta.keep.sum()),
      f"({meta.keep.sum() / nx:.1%} of longitudes)", flush=True)
# fixed ADT contour through the mean half-power point (x6 loop; +-3 deg = +-12 cells of 1/4 deg)
HW = int(round(3.0 / (lat[1] - lat[0])))
sel = np.nonzero(meta.acc.values)[0]; lat_c = np.full((300, len(meta)), np.nan, np.float32); lev = np.full(len(meta), np.nan)
adt_f = np.asarray(adt_m, np.float32)
for i in sel:
    j = int(meta.j[i]); lc = float(meta.lat_hp_clim[i])
    if not np.isfinite(A4[:, j]).all():
        continue
    L = np.interp(lc, lat, A4[:, j]); lev[i] = L
    i0 = int(np.argmin(np.abs(lat - lc))); lo, hi = max(i0 - HW, 0), min(i0 + HW + 1, ny)
    seg = adt_f[:, lo:hi, j]; la = lat[lo:hi]; s = np.sign(seg - L); c = s[:, :-1] * s[:, 1:] < 0
    for m in range(300):
        k = np.nonzero(c[m])[0]
        if len(k) == 0:
            continue
        xs = la[k] + (L - seg[m, k]) * (la[k + 1] - la[k]) / (seg[m, k + 1] - seg[m, k]); lat_c[m, i] = xs[np.argmin(np.abs(xs - lc))]
out = PC.HERE / f"p4_08_tracks_{a.thr}_{a.ice}.npz"
np.savez_compressed(out, lat_hp=arr["lat_hp"], peak_val=arr["peak_val"], lat_c=lat_c, level=lev,
                    j=meta.j.values, lon=meta.lon.values, lat_hp_clim=meta.lat_hp_clim.values, lat_lo=meta.lat_lo.values,
                    lat_hi=meta.lat_hi.values, acc=meta.acc.values, coverage=meta.coverage.values, keep=meta.keep.values,
                    clim_max=clim_max)
print("saved", out.name, round(time.time() - t0), "s")
