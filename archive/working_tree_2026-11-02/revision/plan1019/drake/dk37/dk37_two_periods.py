"""dk37: trends in two periods, 2001-2017 and 2018-2025 (2001-2025 is in dk8 / dk25), Product 1, calibrated Theil-Sen on the
deseasonalised monthly series (rstats.trend_bcz, n_outer = 150), as dk8 and dk25.
(eke) EKE between the SAF and PF per 2 deg longitude, 86-34 W (band, ice and weighting exactly as dk8).
(velocity) Cross-track surface geostrophic velocity along the mean LM Gould track per 0.25 deg latitude (as dk25): monthly where
    the bin series is complete in the period, otherwise annual means (>= 6 months; skipped if < 6 years).
Usage: python dk37_two_periods.py <revision dir> <dk5 stream edges csv> <dk24 npz> <out dir> <eke|velocity>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, EDG, SEC, OUT = (Path(a) for a in sys.argv[1:5]); PART = sys.argv[5]; OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25); YR = np.repeat(np.arange(2001, 2026), 12)
PER = {"2001-2017": (2001, 2017), "2018-2025": (2018, 2025)}
rows = []
def eke_part():
    BOX = dict(longitude=slice(-87, -33), latitude=slice(-66, -44))
    _e = pd.read_csv(EDG).set_index("psi").adt_edge
    LEV_SAF, LEV_PF = (float(np.interp(p, _e.index[::-1], _e.values[::-1])) for p in (-0.5, -1.125))
    fs = sorted((REV / "data" / "monthly_allsat").glob("*/*.nc")); fs = [f for f in fs if 2001 <= int(f.parent.name) <= 2025]
    U, V, VI, IC = [], [], [], []
    for f in fs:
        with xr.open_dataset(f) as d:
            s = d.sel(**BOX); U.append(s.u_mean.values); V.append(s.v_mean.values); VI.append(s.var_intra.values); IC.append(s.ice_frac.values)
            lat, lon = s.latitude.values, s.longitude.values
    U, V, VI, IC = map(lambda a: np.array(a, np.float32), (U, V, VI, IC))
    E = (0.5 * ((U - np.nanmean(U, 0)) ** 2 + (V - np.nanmean(V, 0)) ** 2) + 0.5 * VI) * 1e4
    a = xr.open_dataset(REV / "outputs" / "clim_allsat.nc").adt.interp(latitude=lat, longitude=lon).values
    band = np.isfinite(a) & (a >= LEV_PF) & (a <= LEV_SAF) & (np.nanmax(IC, 0) <= 0.10); w = np.cos(np.radians(lat))[:, None] * band
    for l0 in np.arange(-86, -34, 2.0):
        m = (lon >= l0) & (lon < l0 + 2); ww = w[:, m]
        if ww.sum() == 0: continue
        y = np.nansum(E[:, :, m] * ww, axis=(1, 2)) / np.sum(ww * np.isfinite(E[0][:, m]))
        for pn, (y0, y1) in PER.items():
            k = (YR >= y0) & (YR <= y1); r = rstats.trend_bcz(T[k], y[k], month=MON[k], n_outer=150)
            rows.append(dict(kind="eke_band", x=l0 + 1.0, period=pn, mean=np.nanmean(y), slope=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"]))
        print("eke", l0, flush=True)
def velocity_part():
    Z = np.load(SEC); lo, la = Z["lon"], Z["lat"]
    tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
    nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
    un = (Z["u"] * nx + Z["v"] * ny) * 100; un[Z["ice"] > 0.10] = np.nan
    for c in np.arange(-63.5, -54.99, 0.25):
        s = np.abs(la - c) <= 0.125
        if not s.any(): continue
        y = np.nanmean(un[:, s], 1)
        for pn, (y0, y1) in PER.items():
            k = (YR >= y0) & (YR <= y1); yy = y[k]; r = None
            if np.isfinite(yy).all(): r = rstats.trend_bcz(T[k], yy, month=MON[k], n_outer=150)
            else:
                df = pd.DataFrame(dict(year=YR[k], v=yy)).dropna(); g = df.groupby("year").v; an = g.mean()[g.count() >= 6]
                if len(an) >= 6: r = rstats.trend_bcz(an.index.values.astype(float), an.values, monthly=False, n_outer=150)
            rows.append(dict(kind="velocity", x=c, period=pn, mean=np.nanmean(y), slope=np.nan if r is None else 10 * r["slope"],
                             lo=np.nan if r is None else 10 * r["ci_lo"], hi=np.nan if r is None else 10 * r["ci_hi"], p=np.nan if r is None else r["p"]))
        print("vel", c, flush=True)
eke_part() if PART == "eke" else velocity_part()
pd.DataFrame(rows).to_csv(OUT / f"dk37_two_periods_{PART}.csv", index=False); print("done")
