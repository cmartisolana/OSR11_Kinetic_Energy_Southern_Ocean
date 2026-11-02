"""M2: our data processed with the procedure of Martinez-Moreno et al. (2021) (github.com/josuemtzmo/EKE_SST_trends):
EKE area-integrated with rho0 = 1025 kg m-3 (PJ m-1) over 35-65S, NO ice mask (all pixels finite in every month),
12-month running mean (their 365-day rolling mean), Theil-Sen slope of the smoothed series, % of its mean.
Differences that remain: monthly instead of daily fields; our EKE = 0.5|U_m - ubar_2001-2025|^2 + 0.5 var_intra
(theirs: 0.5(ugosa^2+vgosa^2), anomalies w.r.t. the 1993-2012 DUACS reference); DT2024 1/8 (all-sat) and 1/4 (two-sat)
instead of vDT2018 (v2019) 1/4 deg coarsened to 1 deg; their ACC region = all ocean south of 35S in a 65S-65N domain.
Their values (repository notebooks, 1993-2020): ACC dynamical region +0.0368 PJ m-1 per decade (+5.1 %), Southern basin
+0.05 PJ m-1 per decade (+4.0 %), global +0.09 (+1.2 %). Their significance: modified Mann-Kendall (effective n).
The calibrated test (trend_bcz) is also given on the unsmoothed monthly series.
Usage: python m2_mm21_procedure.py <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
rows = []
for p in ("allsat", "twosat"):
    fs = [REV/"data"/f"monthly_{p}"/f"{y}"/f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
    with xr.open_dataset(fs[0]) as d0:
        lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float)
    sh = (lat.size, lon.size); su = np.zeros(sh); sv = np.zeros(sh); nv = np.zeros(sh)
    for f in fs:
        with xr.open_dataset(f) as d:
            u = d.u_mean.values.astype(float); v = d.v_mean.values.astype(float); ok = np.isfinite(u) & np.isfinite(v)
            su += np.where(ok, u, 0); sv += np.where(ok, v, 0); nv += ok
    valid = nv == len(fs); ub = su/np.maximum(nv, 1); vb = sv/np.maximum(nv, 1)
    area = (6.371e6**2)*np.radians(abs(lat[1]-lat[0]))*np.radians(abs(lon[1]-lon[0]))*np.cos(np.radians(lat))[:, None]*np.ones(sh)
    ser = []
    for f in fs:
        with xr.open_dataset(f) as d:
            u = d.u_mean.values.astype(float); v = d.v_mean.values.astype(float); s2 = d.var_intra.values.astype(float)
        eke = 0.5*((u-ub)**2+(v-vb)**2) + 0.5*s2
        ser.append(1025*np.nansum(np.where(valid, eke*area, 0))/1e15)
    t = pd.DataFrame(dict(year=np.repeat(range(2001, 2026), 12), month=np.tile(range(1, 13), 25), EKE_PJm=ser))
    t["x"] = t.year + (t.month-0.5)/12
    t["EKE_12mo"] = t.EKE_PJm.rolling(12, center=True).mean()
    t.to_csv(OUT/f"m2_series_{p}.csv", index=False)
    for per, (a, b) in {"2001-2020": (2001, 2020), "2001-2025": (2001, 2025)}.items():
        s = (t.year >= a) & (t.year <= b)
        sm = t[s].dropna(subset=["EKE_12mo"])
        ts = stats.theilslopes(sm.EKE_12mo.values, sm.x.values)
        o = R.trend_bcz(t.x.values[s], t.EKE_PJm.values[s], month=t.month.values[s])
        rows.append(dict(product=p, period=per, mean_PJm=float(t.EKE_PJm[s].mean()), n_valid_pix=int(valid.sum()),
                         mm21proc_slope_PJm_dec=10*ts[0], mm21proc_pct_dec=1000*ts[0]/sm.EKE_12mo.mean(),
                         bcz_slope_PJm_dec=10*o["slope"], bcz_ci=(10*o["ci_lo"], 10*o["ci_hi"]), bcz_p=o["p"],
                         bcz_pct_dec=1000*o["slope"]/t.EKE_PJm[s].mean()))
df = pd.DataFrame(rows); df.to_csv(OUT/"m2_trends.csv", index=False); print(df.round(4).to_string())
