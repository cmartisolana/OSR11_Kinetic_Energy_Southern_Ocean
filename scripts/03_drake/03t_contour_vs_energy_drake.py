"""Drake Passage version of the contour-vs-energy test (former x1_drake.py, repaired for the publication repository).

Purpose : on the retained circumpolar half-power tracks (02c, keep = True; 62 tracks in 70-55 W, 105 in 86-40 W; NOT the
          Drake-rule envelopes of Fig. 2b), compare the half-power latitude with the latitude of the fixed ADT contour through
          the mean half-power point: monthly trends (calibrated test) and the 2016-2025 minus 2001-2015 change of the mean
          position (Fig. 2b inset; Sect. 3.1).
          Changes from the working copy: (1) the two interleaved docstrings that made the file fail to compile are merged;
          (2) the Drake columns (90-35 W) of adt_mean are read directly from the monthly files, replacing extract_drake_cols.py,
          which deleted every monthly file it read (that script is kept, inert, only in archive/); (3) outputs go to
          outputs/x1/ and the late-minus-early shifts, previously only printed, are written to x1_drake_period_shift.json.
Inputs  : outputs/a2/tracks_baseline.npz (02c), outputs/base/mean_adt.nc, $OSR11_DATA/monthly_allsat (00a)
Outputs : outputs/x1/x1_drake_trends.csv, x1_drake_period_shift.json (Fig. 2b inset), x1_drake_tracks.npz,
          x1_drake_implied_shift.csv (contour shift implied by the ADT rise / median local ADT gradient; SI Table S3)
Run     : python scripts/03_drake/03t_contour_vs_energy_drake.py     (~1 min)
Seeds   : trend_bcz(seed=0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import json
import numpy as np, pandas as pd, xarray as xr
from osr11 import rstats as R, paths as PA

KM = 111.195
O = PA.out("x1")
T = np.load(PA.OUT / "a2" / "tracks_baseline.npz")
keep = T['keep']; j = T['j'][keep]; lat_hp = T['lat_hp'][:, keep].astype(float); lhc = T['lat_hp_clim'][keep]
madt = xr.open_dataset(PA.MEAN_ADT)['mean_adt']
lat_full = madt.latitude.values.astype(float); lon_full = madt.longitude.values.astype(float)
lon_tr = lon_full[j]
files = PA.monthly_files("allsat"); assert len(files) == 300 and all(f.exists() for f in files)
# Drake columns of the monthly ADT (90-35 W), read once
COLS = []
for f in files:
    with xr.open_dataset(f) as d:
        sub = d.sel(longitude=slice(-90, -35)); COLS.append(sub.adt_mean.values.astype('float32'))
        lat, lon_c = sub.latitude.values.astype(float), sub.longitude.values.astype(float)
assert np.allclose(lat, lat_full)


def run(lo, hi, label):
    sel = np.where((lon_tr >= lo) & (lon_tr <= hi))[0]
    jj = j[sel]; jc = np.array([np.argmin(np.abs(lon_c - lon_full[k])) for k in jj]); assert np.allclose(lon_c[jc], lon_full[jj])
    hp = lat_hp[:, sel]; lhc_s = lhc[sel]
    M = madt.values
    lev = np.array([np.interp(lhc_s[k], lat, M[:, jj[k]]) for k in range(sel.size)])
    W = 3.0; adt_hp = np.full(hp.shape, np.nan); lat_c = np.full(hp.shape, np.nan)
    for m in range(len(files)):
        A = COLS[m].astype(float)
        for k in range(sel.size):
            col = A[:, jc[k]]
            if np.isfinite(hp[m, k]):
                v = np.interp(hp[m, k], lat, col); adt_hp[m, k] = v if np.isfinite(v) else np.nan
            w = (lat >= lhc_s[k] - W) & (lat <= lhc_s[k] + W); la = lat[w]; c = col[w] - lev[k]; ok = np.isfinite(c)
            if ok.sum() < 3: continue
            la, c = la[ok], c[ok]; s = np.where(np.diff(np.sign(c)) != 0)[0]
            if s.size == 0: continue
            x = la[s] + (0 - c[s]) * (la[s + 1] - la[s]) / (c[s + 1] - c[s]); lat_c[m, k] = x[np.argmin(np.abs(x - lhc_s[k]))]
    t = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25); yr = np.repeat(np.arange(2001, 2026), 12)

    def series(a):
        an = a - np.nanmean(a, 0); return np.nanmean(an, 1)
    print(f'\n== {label}: {sel.size} tracks, lon {lon_tr[sel].min():.1f} to {lon_tr[sel].max():.1f}; contour found in '
          f'{np.isfinite(lat_c).mean()*100:.0f} % of track-months, hp in {np.isfinite(hp).mean()*100:.0f} %')
    rows = []
    for nm, a, sc, u in (('lat_hp', hp, KM, 'km'), ('lat_contour', lat_c, KM, 'km'), ('adt_hp', adt_hp, 100, 'cm')):
        y = series(a) * sc; r = R.trend_bcz(t, y, month=month, monthly=True, n_outer=300, n_bias=100, seed=0)
        rows.append(dict(region=label, series=nm, unit=u + ' per decade', slope=10 * r['slope'], lo=10 * r['ci_lo'], hi=10 * r['ci_hi'], p=r['p']))
        print(f"  {nm:12s} {10*r['slope']:+.1f} {u}/decade ({10*r['ci_lo']:+.1f} to {10*r['ci_hi']:+.1f}) p={r['p']:.3f}")
    e = (yr <= 2015); l_ = (yr >= 2016); shift = {}
    for nm, a in (('hp', hp), ('contour', lat_c)):
        an = a - np.nanmean(a, 0); d = (np.nanmean(an[l_], 0) - np.nanmean(an[e], 0)) * KM
        shift[nm] = dict(mean=float(np.nanmean(d)), median=float(np.nanmedian(d)), frac_south=float((d < 0).mean()))
        print(f"  {nm:8s} late-minus-early mean position: {np.nanmean(d):+.1f} km (median {np.nanmedian(d):+.1f}; southward at {(d<0).mean()*100:.0f} % of tracks)")
    # |dADT/dy| of the mean ADT at each track's climatological half-power latitude (cm per km); the working copy took the median
    # over the whole latitude column, which mixes in the flat flanks; SI Table S3 uses the value at the half-power point
    g = np.array([abs(np.interp(lhc_s[q], lat, np.gradient(M[:, k], lat * KM))) for q, k in enumerate(jj)]) * 100
    imp = -(rows[2]['slope'] / (np.nanmedian(g) * 10)) * 10
    print(f"  median |dADT/dlat| {np.nanmedian(g)*10:.1f} cm per 10 km -> implied contour shift from adt_hp trend: {imp:+.1f} km/decade")
    IMPL.append(dict(region=label, median_grad_cm_per_10km=float(np.nanmedian(g) * 10), adt_hp_trend_cm_per_decade=rows[2]['slope'],
                     implied_contour_shift_km_per_decade=float(imp)))
    return rows, dict(hp=hp, lat_c=lat_c, adt_hp=adt_hp, lon=lon_tr[sel], lhc=lhc_s, lev=lev), sel.size, shift


out = []; js = {}; D = {}; IMPL = []
for key, lo, hi, label in (("7055", -70, -55, 'Drake 70-55W'), ("8640", -86, -40, 'Drake region 86-40W')):
    r, D[key], n, sh = run(lo, hi, label); out += r
    js[key] = dict(n_tracks=int(n), halfpower_late_minus_early_km=round(sh['hp']['mean'], 1),
                   contour_late_minus_early_km=round(sh['contour']['mean'], 1))
js["note"] = ("2016-2025 minus 2001-2015 mean latitude anomaly (km, + north), mean over the retained circumpolar tracks "
              "(>=90 % rule) in the region; scripts/03_drake/03t_contour_vs_energy_drake.py")
pd.DataFrame(out).to_csv(O / 'x1_drake_trends.csv', index=False)
pd.DataFrame(IMPL).to_csv(O / 'x1_drake_implied_shift.csv', index=False)
np.savez_compressed(O / 'x1_drake_tracks.npz', **{k + '_7055': v for k, v in D["7055"].items()}, **{k + '_8640': v for k, v in D["8640"].items()})
json.dump(js, open(O / 'x1_drake_period_shift.json', 'w'), indent=1); print(json.dumps(js, indent=1))
