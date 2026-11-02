"""DK1: Drake Passage KE, EKE, along-stream speed and wind-stress trends (Drake-first structure).
Regions (fixed earlier, a6/a5d):  box 80-52W, 63-54S (Fig. 3 box; pixels with ice fraction > 0.10 in any month excluded;
'nomask' variant keeps them); upstream 100-80W and downstream (Scotia) 55-35W inside the ACC band (2001-2025 mean ADT
-0.6..+0.2 m), ice-masked as in a5_energy.py. Wind: eastward stress (WIND_GLO_PHY_L4_MY_012_006) averaged over the same
longitudes and 63-54S. Definitions as a5: KE = 0.5|U_m|^2, EKE = 0.5|U_m - ubar|^2 + 0.5 var_intra, SPD = U_m.ubar/|ubar|.
Contrast: trend of (EKE/mean - KE/mean) in % points per decade. Test: rstats.trend_bcz. Periods 2001-2025, 2005-2019.

Purpose : KE, EKE, speed and wind-stress trends in the Drake box, upstream and downstream regions, both products (Fig. S2).
Inputs  : $OSR11_DATA/monthly_allsat, monthly_twosat (00a), $OSR11_DATA/wind_taux.nc (00b)
Outputs : outputs/drake/dk1/dk1_series.csv, dk1_trends.csv (distributed)
Run     : python scripts/03_drake/03r_drake_regions.py
Seeds   : rstats.trend_bcz default seed (0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from multiprocessing import Pool
from osr11 import rstats as R, paths as PA
OUT = PA.dk(1)
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month-0.5)/12
REG = {"box": (-80, -52, -63, -54, False), "box_nomask": (-80, -52, -63, -54, None), "up": (-100, -80, None, None, True), "down": (-55, -35, None, None, True)}
def series(p):
    fs = PA.monthly_files(p)
    U, V, S2, A, I = [], [], [], [], []
    for f in fs:
        with xr.open_dataset(f) as d:
            s = d.sel(longitude=slice(-102, -33), latitude=slice(-66, -45)).load()
        U.append(s.u_mean.values.astype(float)); V.append(s.v_mean.values.astype(float)); S2.append(s.var_intra.values.astype(float))
        A.append(s.adt_mean.values.astype(float)); I.append(np.nan_to_num(s.ice_frac.values))
    lat, lon = s.latitude.values, s.longitude.values
    U, V, S2, A, I = map(np.array, (U, V, S2, A, I))
    ub, vb = np.nanmean(U, 0), np.nanmean(V, 0); am = np.nanmean(A, 0); sp = np.hypot(ub, vb)
    allfin = np.all(np.isfinite(U) & np.isfinite(V), 0); ice_any = I.max(0) > 0.10
    acc = (am >= -0.6) & (am <= 0.2)
    L, Y = np.meshgrid(lon, lat); w = np.cos(np.radians(Y))
    KE = 0.5*(U**2+V**2)*1e4; EKE = (0.5*((U-ub)**2+(V-vb)**2) + 0.5*S2)*1e4; SPD = (U*ub+V*vb)/np.where(sp > 0, sp, np.nan)*100
    out = {}
    for r, (l0, l1, a0, a1, band) in REG.items():
        m = (L >= l0) & (L < l1) & allfin
        m &= ((Y >= a0) & (Y <= a1)) if band in (False, None) else acc
        if band is not None: m &= ~ice_any
        for vn, F in (("KE", KE), ("EKE", EKE), ("SPD", SPD)):
            out[f"{p}_{r}_{vn}"] = np.array([np.nansum(F[k][m]*w[m])/np.sum(w[m]) for k in range(F.shape[0])])
        out[f"{p}_{r}_npix"] = np.full(F.shape[0], m.sum())
    return out
def one(a):
    lab, x, y, mo = a
    o = R.trend_bcz(x, y, month=mo); mu = float(np.nanmean(y))
    return dict(series=lab, slope_per_decade=10*o["slope"], ci_lo=10*o["ci_lo"], ci_hi=10*o["ci_hi"], p=o["p"], mean=mu,
                pct=1000*o["slope"]/mu, pct_lo=1000*o["ci_lo"]/mu, pct_hi=1000*o["ci_hi"]/mu)
if __name__ == "__main__":
    S = {}
    for p in ("allsat", "twosat"): S.update(series(p))
    W = xr.open_dataset(PA.WIND).eastward_stress
    for r, (l0, l1, *_ ) in REG.items():
        if r == "box_nomask": continue
        ww = W.sel(longitude=slice(l0, l1), latitude=slice(-63, -54))
        S[f"wind_{r}_TAUX"] = ww.weighted(np.cos(np.radians(ww.latitude))).mean(("latitude", "longitude")).values[:300]
    D = pd.concat([T, pd.DataFrame(S)], axis=1); D.to_csv(OUT/"dk1_series.csv", index=False)
    jobs = []
    for per, (a, b) in {"2001-2025": (2001, 2025), "2005-2019": (2005, 2019)}.items():
        s = (D.year >= a) & (D.year <= b)
        for c in [c for c in D.columns if c.endswith(("_KE", "_EKE", "_SPD", "_TAUX"))]:
            jobs.append((f"{c}|{per}", D.x[s].values, D[c][s].values, D.month[s].values))
        for p in ("allsat", "twosat"):
            for r in ("box", "box_nomask", "up", "down"):
                e, k = D[f"{p}_{r}_EKE"][s], D[f"{p}_{r}_KE"][s]
                jobs.append((f"{p}_{r}_EKEminusKE|{per}", D.x[s].values, (100*(e/e.mean()-k/k.mean())+100).values, D.month[s].values))
            e1, e2 = D[f"allsat_box_EKE"][s], D[f"twosat_box_EKE"][s]
            jobs.append((f"box_EKE_all-minus-two|{per}", D.x[s].values, (100*(e1/e1.mean()-e2/e2.mean())+100).values, D.month[s].values))
    with Pool(2) as pool: res = pool.map(one, jobs)
    df = pd.DataFrame(res); df[["name", "period"]] = df.series.str.split("|", expand=True)
    dm = df.name.str.contains("minus"); df.loc[dm, ["pct", "pct_lo", "pct_hi"]] = df.loc[dm, ["slope_per_decade", "ci_lo", "ci_hi"]].values
    df.to_csv(OUT/"dk1_trends.csv", index=False)
    print(df[["name", "period", "pct", "pct_lo", "pct_hi", "p"]].round(3).to_string())
