"""x7 (5 Nov 2026): products of the relative-threshold half-power tracking (x6) needed by the figures:
 - all monthly half-power points in the ACC band (month index, lon, lat) for the occurrence maps (outputs/x7/x7_points_all.npy)
 - per retained track: monthly ADT at the half-power point (adt_hp), ADT gradient at the climatological point
 - composites (monthly KE and ADT within +-28 cells of the climatological half-power latitude) for the circumpolar set and the 86-40 W subset
 - trends (calibrated test) of the mean half-power latitude, contour latitude and adt_hp for circumpolar, 86-40 W and 70-55 W; period
   shifts 2001-2017/2018-2025 and 2001-2015/2016-2025; implied contour displacement from the ADT rise over the median gradient.

Origin  : x7_relthr_products.py (5 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : everything the figures need from the relative-threshold tracks: retained tracks with the fixed ADT contour, composites across the
          pathway (Figs. 2b, 3b), position and ADT trends (Sects. 3.1-3.2; Figs. S3f, S5, S6c), implied contour displacement (SI S4.2).
Inputs  : $OSR11_DATA/monthly_allsat (raw data), outputs/x6/x6_relthr_tracks.npz (02n), outputs/base/*
Outputs : outputs/x7/x7_tracks.npz, x7_trends.csv, x7_period_shifts.json, x7_implied_shift.csv, x7_composite_annual_{circ,drake}.npz (distributed);
          x7_composite_{circ,drake}.npz (monthly composites, 255 MB) and x7_points_all.npy (28 MB) are written but not distributed.
          The annual composites were added on 17 Nov 2026: the same averages the figure scripts made from the monthly composites,
          stored so that the figures run from the distributed outputs.
Run     : python scripts/02_circumpolar/02p_pathway_products.py   (raw data; about 8 GB of memory)
Seeds   : rstats.trend_bcz defaults.
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, xarray as xr, pandas as pd, json
from osr11 import halfpower as HPm, rstats as R, paths as PA
files = PA.monthly_files("allsat"); assert all(f.exists() for f in files), "monthly files missing (scripts/00_extract/00a)"; X6, OUT = PA.out("x6"), PA.out("x7")
d0 = xr.open_dataset(files[0]); lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float); ny, nx = len(lat), len(lon)
Z = np.load(X6 / "x6_relthr_tracks.npz"); keep = Z["keep"]; j_all = Z["j"]; lhc = Z["lat_hp_clim"]; acc = Z["acc"]
adt = xr.open_dataset(PA.MEAN_ADT)["mean_adt"].values; accmask = np.isfinite(adt) & (adt >= -0.6) & (adt <= 0.2)
A = xr.open_dataset(PA.CLIM).adt.values
KM = 111.195; HW = 28; t = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25); yr = np.repeat(np.arange(2001, 2026), 12)
ke_all = np.empty((300, ny, nx), np.float32); adt_m = np.empty((300, ny, nx), np.float32)
for m, f in enumerate(files):
    ds = xr.open_dataset(f); u = ds.u_mean.values; v = ds.v_mean.values; ke_all[m] = 0.5 * (u ** 2 + v ** 2); adt_m[m] = ds.adt_mean.values
kemean = np.nanmean(ke_all, 0); clim_max = np.nanmax(np.where(accmask, kemean, np.nan), 0); clim_max = np.where(np.isfinite(clim_max) & (clim_max > 0), clim_max, np.nan)
p = HPm.SETS["relthr"]; counts, tab = HPm.run_set(ke_all, lat, lon, p, clim_max=clim_max); FI = HPm.FI
pts = []
for k in range(300):
    jj, kk = np.nonzero(np.arange(HPm.MAX_ENV)[None, :] < counts[k][:, None]); r = tab[k, jj, kk]; lh = r[:, FI["lat_hp"]]
    ok = HPm.in_acc(accmask, lat, lon, lh, jj); pts.append(np.c_[np.full(ok.sum(), k), lon[jj[ok]], lh[ok]])
pts = np.concatenate(pts); np.save(OUT / "x7_points_all.npy", pts.astype(np.float64)); print("points", len(pts), flush=True)
# per retained track: adt at the half-power point, gradient, composites
ki = np.nonzero(keep)[0]; nk = len(ki); lat_hp = Z["lat_hp"][:, ki]; lat_c = Z["lat_c"][:, ki]; jj = j_all[ki]; lc = lhc[ki]
adt_hp = np.full((300, nk), np.nan, np.float32); grad = np.full(nk, np.nan)
ke_c = np.full((300, nk, 2 * HW + 1), np.nan, np.float32); adt_c = np.full((300, nk, 2 * HW + 1), np.nan, np.float32); i0s = np.zeros(nk, int)
dlat = lat[1] - lat[0]
for n, (j, c) in enumerate(zip(jj, lc)):
    col = adt_m[:, :, j]; lh = lat_hp[:, n]; ok = np.isfinite(lh)
    fi = np.clip((lh[ok] - lat[0]) / dlat, 0, ny - 1.001); i = fi.astype(int); w = fi - i
    adt_hp[ok, n] = col[ok, i] * (1 - w) + col[ok, i + 1] * w
    i0 = int(np.argmin(np.abs(lat - c))); i0s[n] = i0
    if 1 <= i0 < ny - 1: grad[n] = 100 * (A[i0 + 1, j] - A[i0 - 1, j]) / (2 * dlat * KM) * 10   # cm per 10 km
    lo, hi = i0 - HW, i0 + HW + 1
    if lo >= 0 and hi <= ny: ke_c[:, n] = ke_all[:, lo:hi, j] * 1e4; adt_c[:, n] = adt_m[:, lo:hi, j]
lev = Z["level"][ki]; off = np.arange(-HW, HW + 1) * dlat; lonk = lon[jj]
np.savez_compressed(OUT / "x7_composite_circ.npz", ke=ke_c, adt=adt_c, offset_deg=off, lat_hp_clim=lc, lev=lev, lon=lonk, lat0=lat[i0s])
sel = (lonk >= -86) & (lonk <= -40)
np.savez_compressed(OUT / "x7_composite_drake.npz", ke=ke_c[:, sel], adt=adt_c[:, sel], offset_deg=off, lat_hp_clim=lc[sel], lev=lev[sel], lon=lonk[sel], lat0=lat[i0s[sel]])
# annual composites (the curves of Figs. 2b and 3b): mean over the months of each year, then over the tracks; ADT relative to each track's contour level (cm)
YRS_ = np.arange(2001, 2026)
for name_, s_ in (("circ", np.ones(nk, bool)), ("drake", sel)):
    ke_y = np.array([np.nanmean(np.nanmean(ke_c[yr == y][:, s_], 0), 0) for y in YRS_])
    ad_y = np.array([np.nanmean(np.nanmean((adt_c[yr == y][:, s_] - lev[s_][None, :, None]) * 100, 0), 0) for y in YRS_])
    np.savez_compressed(OUT / f"x7_composite_annual_{name_}.npz", ke_y=ke_y, ad_y=ad_y, offset_deg=off, years=YRS_, ntracks=int(s_.sum()), lon=lonk[s_], lev=lev[s_], lat_hp_clim=lc[s_])
np.savez_compressed(OUT / "x7_tracks.npz", lat_hp=lat_hp, lat_c=lat_c, adt_hp=adt_hp, lon=lonk, lat_hp_clim=lc, lev=lev, grad_cm_per_10km=grad)
print("composites saved", flush=True)
# trends and shifts
def series(a): return np.nanmean(a - np.nanmean(a, 0), 1)
rows = []; shifts = {}; implied = []
for region, m in (("circumpolar", np.ones(nk, bool)), ("Drake region 86-40W", sel), ("Drake 70-55W", (lonk >= -70) & (lonk <= -55))):
    for key, arr, fac, unit in (("lat_hp", lat_hp, KM, "km per decade"), ("lat_contour", lat_c, KM, "km per decade"), ("adt_hp", adt_hp, 100, "cm per decade")):
        y = series(arr[:, m]) * fac; r = R.trend_bcz(t, y, month=month, monthly=True)
        rows.append(dict(region=region, series=key, unit=unit, slope=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"], n_tracks=int(m.sum()))); print(rows[-1], flush=True)
    sh = {}
    for split in (2016, 2018):
        e, l = yr < split, yr >= split
        for key, arr in (("halfpower", lat_hp), ("contour", lat_c)):
            an = (arr[:, m] - np.nanmean(arr[:, m], 0)) * KM; sh[f"{key}_late_minus_early_km_split{split}"] = round(float(np.nanmean(np.nanmean(an[l], 0)) - np.nanmean(np.nanmean(an[e], 0))), 1)
    sh["n_tracks"] = int(m.sum()); shifts[region] = sh
    g = np.nanmedian(grad[m]); a_tr = [x for x in rows if x["region"] == region and x["series"] == "adt_hp"][0]["slope"]
    implied.append(dict(region=region, median_grad_cm_per_10km=g, adt_hp_trend_cm_per_decade=a_tr, implied_contour_shift_km_per_decade=-10 * a_tr / g))
pd.DataFrame(rows).to_csv(OUT / "x7_trends.csv", index=False); json.dump(shifts, open(OUT / "x7_period_shifts.json", "w"), indent=1); pd.DataFrame(implied).to_csv(OUT / "x7_implied_shift.csv", index=False)
print(json.dumps(shifts, indent=1)); print(pd.DataFrame(implied).round(2).to_string()); print("done")
