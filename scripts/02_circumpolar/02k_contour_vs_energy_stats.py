"""Contour-vs-energy test, statistics (former x1_stats.py, plus the bathymetric split that was computed interactively).

Purpose : trends (calibrated Theil-Sen, rstats.trend_bcz, 300 null series) of the circumpolar-mean anomaly of
            lat_hp       half-power latitude of the 900 retained tracks          (km per decade)
            lat_contour  latitude of the fixed ADT contour through the mean half-power point  (km per decade)
            adt_hp       ADT at the half-power point                             (cm per decade)
          and of the ACC-band mean ADT; per-track OLS trends; sector trends (printed); and the split of the half-power trend by
          bathymetric roughness and depth (Discussion: "over rough and smooth topography alike").
          Roughness = standard deviation of depth in a +-2 deg box (+-16 pixels of bathymetry_so.nc) around the nearest pixel
          to the climatological half-power point; depth = depth at that pixel. Split at the median roughness (448 m) and at
          3500 m depth. This definition reproduces the working-repository x1_track_bathy.csv exactly (checked 2 Nov 2026).
Inputs  : outputs/x1/x1_tracks.npz, outputs/x1/x1_band_adt.npy (02j), outputs/base/bathymetry_so.nc, outputs/base/mean_adt.nc
Outputs : outputs/x1/x1_circumpolar_trends.csv, x1_pertrack.npz, x1_track_bathy.csv, x1_bathy_split.csv, x1_band_adt_trend.csv,
          x1_implied_shift.csv (contour shift implied by the ADT rise / median local ADT gradient; SI Table S3)
Run     : python scripts/02_circumpolar/02k_contour_vs_energy_stats.py      (no raw data needed; ~10 min)
Seeds   : trend_bcz(seed=0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
from osr11 import rstats as R, paths as PA

O = PA.out("x1")
Z = np.load(O / "x1_tracks.npz"); KM = 111.195
lat_hp, adt_hp, lat_c, lon = Z['lat_hp'], Z['adt_hp'], Z['lat_c'], Z['lon']


def series(a):
    an = a - np.nanmean(a, 0); return np.nanmean(an, 1)


def bcz(y):
    return R.trend_bcz(t, y, month=month, monthly=True, n_outer=300, n_bias=100, seed=0)


t = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25); yr = np.repeat(np.arange(2001, 2026), 12)
rows = []
for nm, a, sc, u in (('lat_hp', lat_hp, KM, 'km'), ('lat_contour', lat_c, KM, 'km'), ('adt_hp', adt_hp, 100, 'cm')):
    y = series(a) * sc; r = bcz(y)
    rows.append(dict(series=nm, unit=u + ' per decade', slope=10 * r['slope'], lo=10 * r['ci_lo'], hi=10 * r['ci_hi'], p=r['p']))
    print(f"{nm:12s} {10*r['slope']:+.2f} {u}/decade ({10*r['ci_lo']:+.2f} to {10*r['ci_hi']:+.2f}) p={r['p']:.3f}", flush=True)
pd.DataFrame(rows).to_csv(O / 'x1_circumpolar_trends.csv', index=False)


def tr(a):
    out = []
    for k in range(a.shape[1]):
        y = a[:, k]; ok = np.isfinite(y); out.append(stats.linregress(t[ok], y[ok]).slope * 10 if ok.sum() > 100 else np.nan)
    return np.array(out)


th, tc, ta = tr(lat_hp) * KM, tr(lat_c) * KM, tr(adt_hp) * 100
print('per-track median trend: hp %.1f km/dec, contour %.1f km/dec, adt_hp %.2f cm/dec' % (np.nanmedian(th), np.nanmedian(tc), np.nanmedian(ta)))
print('fraction of tracks southward: contour %.2f, hp %.2f' % ((tc < 0).mean(), (th < 0).mean()))
ok = np.isfinite(th) & np.isfinite(tc); print('corr per-track trends hp vs contour r=%.2f' % np.corrcoef(th[ok], tc[ok])[0, 1])
np.savez(O / 'x1_pertrack.npz', th=th, tc=tc, ta=ta, lon=lon)
sec = np.where((lon >= -70) & (lon < 20), 'Atlantic', np.where((lon >= 20) & (lon < 147), 'Indian', 'Pacific'))
for s in ('Atlantic', 'Indian', 'Pacific'):
    m = sec == s; line = f'  {s:9s} n={m.sum():3d}'
    for nm, a, sc in (('hp', lat_hp, KM), ('contour', lat_c, KM), ('adt_hp', adt_hp, 100)):
        y = series(a[:, m]) * sc; ann = pd.Series(y).groupby(yr).mean(); rr = stats.linregress(ann.index, ann.values)
        line += f' | {nm} {10*rr.slope:+.2f} p={rr.pvalue:.2f}'
    print(line)

# contour displacement implied by the ADT rise at the half-power point / median local meridional ADT gradient (SI Table S3)
madt = xr.open_dataset(PA.MEAN_ADT)["mean_adt"]; Mv = madt.values; mlat = madt.latitude.values.astype(float)
g = np.array([abs(np.interp(lhc_k, mlat, np.gradient(Mv[:, k], mlat * KM))) for lhc_k, k in zip(Z['lat_hp_clim'], Z['j'])]) * 100   # cm per km, at the half-power point
imp = -(rows[2]['slope'] / (np.nanmedian(g) * 10)) * 10
pd.DataFrame([dict(region='circumpolar', median_grad_cm_per_10km=float(np.nanmedian(g) * 10), adt_hp_trend_cm_per_decade=rows[2]['slope'],
                   implied_contour_shift_km_per_decade=float(imp))]).to_csv(O / 'x1_implied_shift.csv', index=False)
print(f"median |dADT/dlat| {np.nanmedian(g)*10:.1f} cm per 10 km -> implied contour shift {imp:+.1f} km/decade")

# ACC-band mean ADT (Fig. S7b)
band = np.load(O / 'x1_band_adt.npy') * 100
r = bcz(band)
pd.DataFrame([dict(series='band_adt', unit='cm per decade', slope=10 * r['slope'], lo=10 * r['ci_lo'], hi=10 * r['ci_hi'], p=r['p'])]).to_csv(
    O / 'x1_band_adt_trend.csv', index=False)
print(f"ACC-band mean ADT {10*r['slope']:+.2f} cm/decade ({10*r['ci_lo']:+.2f} to {10*r['ci_hi']:+.2f}) p={r['p']:.3f}")

# bathymetric split of the half-power trend
b = xr.open_dataset(PA.BATHY).elevation; bla, blo, E = b.latitude.values, b.longitude.values, -b.values.astype(float)
lhc = Z['lat_hp_clim']
dep, rough = np.empty(lon.size), np.empty(lon.size)
for k in range(lon.size):
    i = int(np.argmin(np.abs(bla - lhc[k]))); j = int(np.argmin(np.abs(blo - lon[k])))
    dep[k] = E[i, j]; rough[k] = np.nanstd(E[i - 16:i + 17, j - 16:j + 17])
pd.DataFrame(dict(lon=lon, lat_hp_clim=lhc, depth=dep, rough=rough)).to_csv(O / 'x1_track_bathy.csv', index=False)
med = np.median(rough); rows = []
for nm, m in (('rough_half', rough >= med), ('smooth_half', rough < med), ('depth_ge_3500', dep >= 3500), ('depth_lt_3500', dep < 3500)):
    y = series(lat_hp[:, m]) * KM; r = bcz(y)
    rows.append(dict(subset=nm, n_tracks=int(m.sum()), split_value=med if 'half' in nm else 3500.0, slope=10 * r['slope'],
                     lo=10 * r['ci_lo'], hi=10 * r['ci_hi'], p=r['p']))
    print(f"  {nm:14s} n={m.sum():3d} {10*r['slope']:+.1f} km/decade ({10*r['ci_lo']:+.1f} to {10*r['ci_hi']:+.1f}) p={r['p']:.2f}", flush=True)
pd.DataFrame(rows).to_csv(O / 'x1_bathy_split.csv', index=False)
