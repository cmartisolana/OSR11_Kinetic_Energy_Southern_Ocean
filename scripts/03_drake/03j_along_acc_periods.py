"""dk8: along-ACC EKE trend around Drake Passage. Band = between the mean SAF and PF streamlines (central streamlines of the GV23 frontal zones mapped to the 2001-2025 mean ADT, as drawn in Fig. 1a),
pixels ice-free (ice fraction <= 0.10) in every month. Monthly EKE (Eq. 1) averaged (cos-lat weights) in 2 deg longitude bins (calibrated test with 150 null series),
86-34W (Drake Passage domain); trend: Theil-Sen on the deseasonalised series with the calibrated test (trend_bcz). Depth along the band from
bathymetry_so.nc: mean and shallowest depth inside the band per bin. Product 1 (and 4 for reference).

Purpose : as 03i, separately for 2001-2015 and 2016-2025 (input of 05d, Sect. 3.1 percentages).
Inputs  : $OSR11_DATA/monthly_allsat (00a), outputs/base/clim_allsat.nc, outputs/base/bathymetry_so.nc, outputs/drake/dk5/dk5_stream_edges.csv (03b)
Outputs : outputs/drake/dk8/dk8_along_acc_fronts_periods_allsat.csv (distributed). The working repository wrote here (drake/out/dk8) after the 26 Oct path fix but
          its outputs and the figure input sat in drake/dk8/; both now live in outputs/drake/dk8/.
Run     : python scripts/03_drake/03j_along_acc_periods.py
Seeds   : rstats.trend_bcz default seed (0); 150 null series.
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, glob, numpy as np, pandas as pd, xarray as xr
from osr11 import rstats, paths as PA
OUT = PA.dk(8); BOX = dict(longitude=slice(-112, -18), latitude=slice(-66, -44))
adt = xr.open_dataset(PA.CLIM).adt.sel(**BOX)
_e = pd.read_csv(PA.dk(5) / "dk5_stream_edges.csv").set_index("psi").adt_edge
LEV_SAF, LEV_PF = (float(np.interp(p, _e.index[::-1], _e.values[::-1])) for p in (-0.5, -1.125))
bat = xr.open_dataset(PA.BATHY).elevation.sel(**BOX)
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
for prod in ("allsat",):
    fs = [str(f) for f in PA.monthly_files(prod)]
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
