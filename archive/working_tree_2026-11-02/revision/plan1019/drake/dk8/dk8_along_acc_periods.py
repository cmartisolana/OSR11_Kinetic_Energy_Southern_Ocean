"""dk8: along-ACC EKE trend around Drake Passage. Band = between the mean SAF and PF streamlines (central streamlines of the GV23 frontal zones mapped to the 2001-2025 mean ADT, as drawn in Fig. 1a),
pixels ice-free (ice fraction <= 0.10) in every month. Monthly EKE (Eq. 1) averaged (cos-lat weights) in 2 deg longitude bins (calibrated test with 150 null series),
86-34W (Drake Passage domain); trend: Theil-Sen on the deseasonalised series with the calibrated test (trend_bcz). Depth along the band from
bathymetry_so.nc: mean and shallowest depth inside the band per bin. Product 1 (and 4 for reference)."""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[3])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, glob, numpy as np, pandas as pd, xarray as xr
REV = _REV; sys.path.insert(0, REV); import rstats
OUT = (_DK_OUT + "/dk8"); BOX = dict(longitude=slice(-112, -18), latitude=slice(-66, -44))
adt = xr.open_dataset(f"{REV}/outputs/clim_allsat.nc").adt.sel(**BOX)
_e = pd.read_csv((_DK_OUT + "/dk5/dk5_stream_edges.csv")).set_index("psi").adt_edge
LEV_SAF, LEV_PF = (float(np.interp(p, _e.index[::-1], _e.values[::-1])) for p in (-0.5, -1.125))
bat = xr.open_dataset(f"{REV}/data/bathymetry_so.nc").elevation.sel(**BOX)
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
for prod in ("allsat",):
    fs = sorted(glob.glob(f"{REV}/data/monthly_{prod}/*/*.nc")); fs = [f for f in fs if 2001 <= int(f.split('/')[-2]) <= 2025]
    U, V, VI, IC = [], [], [], []
    for f in fs:
        with xr.open_dataset(f) as d:
            s = d.sel(**BOX); U.append(s.u_mean.values); V.append(s.v_mean.values); VI.append(s.var_intra.values); IC.append(s.ice_frac.values)
            lat, lon = s.latitude.values, s.longitude.values
    U, V, VI, IC = map(lambda a: np.array(a, np.float32), (U, V, VI, IC))
    ub, vb = np.nanmean(U, 0), np.nanmean(V, 0)
    E = (0.5 * ((U - ub) ** 2 + (V - vb) ** 2) + 0.5 * VI) * 1e4
    a = adt.interp(latitude=lat, longitude=lon).values
    band = np.isfinite(a) & (a >= LEV_PF) & (a <= LEV_SAF) & (np.nanmax(IC, 0) <= 0.10)
    w = np.cos(np.radians(lat))[:, None] * band
    b = bat.interp(latitude=lat, longitude=lon).values
    rows = []
    for l0 in np.arange(-86, -34, 2.0):
        m = (lon >= l0) & (lon < l0 + 2)
        ww = w[:, m]
        if ww.sum() == 0: continue
        y = np.nansum(E[:, :, m] * ww, axis=(1, 2)) / np.sum(ww * np.isfinite(E[0][:, m]))
        dep = -b[:, m][band[:, m]]
        for per, (a0, a1) in (("2001-2015", (2001, 2015)), ("2016-2025", (2016, 2025))):
            sel = (T >= a0) & (T < a1 + 1); r = rstats.trend_bcz(T[sel], y[sel], month=MON[sel], n_outer=150); mu = np.nanmean(y[sel])
            rows.append(dict(lon=l0 + 1.0, period=per, mean=mu, slope_dec=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"],
                             depth_mean=np.nanmean(dep), depth_min=np.nanmin(dep)))
        print(prod, l0, flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/dk8_along_acc_fronts_periods_{prod}.csv", index=False)
