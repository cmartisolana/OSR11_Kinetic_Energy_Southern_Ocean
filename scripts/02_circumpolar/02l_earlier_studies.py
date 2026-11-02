"""M1: matched comparison with earlier altimetric EKE estimates (plan section 2).
Regions (cos-lat area means; area integral also stored for the 35-65S box, rho0 = 1025 kg m-3, PJ m-1):
  Hogg et al. (2015) sector boxes: Indian 40-150E 57-44S; Pacific 150E-72W (288E) 62-48S; Atlantic 35W-10E 56-46S
  box35_65 (domain of the monthly fields; Martinez-Moreno et al. 2021-type Southern Ocean average, NOT their definition)
  ACC band (2001-2025 mean ADT -0.6..+0.2 m of each product)
Mask: pixels finite in all months and ice fraction <= 0.10 in every month (as a5_energy.py).
EKE_m = 0.5|U_m - ubar|^2 + 0.5 var_intra (ubar = 2001-2025 mean of that product), as a5/v2.
LF4 = 0.5|U_m - ubar|^2 computed on a common 1/4 deg grid (all-sat 2x2 block means of U_m) -> scale-matched
      comparison of the low-frequency part between products (var_intra cannot be block-averaged).
Periods 2001-2012 (overlap with Hogg's record), 2001-2020 (MM21 end year), 2001-2025.
Trend test: rstats.trend_bcz (deseasonalised Theil-Sen; calibrated AR noise; max p of AIC/AR(12)).
Note: with M1_Y0=1993 the mean flow ubar is the 1993-2025 mean and the ice mask covers 1993-2025, so 2001-2025
values will differ slightly from the default run.

Purpose : EKE trends over the regions of earlier altimetric studies, both products (Fig. S1d).
Inputs  : $OSR11_DATA/monthly_allsat and monthly_twosat (00a)
Outputs : outputs/m1/m1_series_<product>.csv, m1_masks_<product>.csv, m1_trends.csv (m1_trends.csv distributed)
Run     : python scripts/02_circumpolar/02l_earlier_studies.py      (env M1_Y0, default 2001)
Seeds   : rstats.trend_bcz default seed (0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, os
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from multiprocessing import Pool
from osr11 import rstats as R, paths as PA
OUT = PA.out("m1")
Y0 = int(os.environ.get("M1_Y0", "2001"))   # set M1_Y0=1993 once 1993-2000 is extracted
YRS = range(Y0, 2026)
def files(p): return PA.monthly_files(p, Y0, 2025)
def blk(a):  # 2x2 block mean (all-sat 1/8 -> 1/4)
    n0, n1 = a.shape[0]//2*2, a.shape[1]//2*2
    return a[:n0, :n1].reshape(n0//2, 2, n1//2, 2).mean((1, 3))
def regions(lat, lon, adt_mean, valid):
    L = np.mod(lon, 360)[None, :]; A = lat[:, None]
    m = {"Hogg_Indian": (L >= 40) & (L < 150) & (A >= -57) & (A <= -44),
         "Hogg_Pacific": (L >= 150) & (L < 288) & (A >= -62) & (A <= -48),
         "Hogg_Atlantic": ((L >= 325) | (L < 10)) & (A >= -56) & (A <= -46),
         "box35_65": np.ones_like(valid, bool),
         "ACC": (adt_mean >= -0.6) & (adt_mean <= 0.2)}
    return {k: v & valid for k, v in m.items()}
def series(p):
    fs = files(p)
    with xr.open_dataset(fs[0]) as d0:
        lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float)
    sh = (lat.size, lon.size); su = np.zeros(sh); sv = np.zeros(sh); sa = np.zeros(sh); nv = np.zeros(sh); imax = np.zeros(sh)
    for f in fs:
        with xr.open_dataset(f) as d:
            u = d.u_mean.values.astype(float); v = d.v_mean.values.astype(float)
            ok = np.isfinite(u) & np.isfinite(v); su += np.where(ok, u, 0); sv += np.where(ok, v, 0); nv += ok
            sa += np.nan_to_num(d.adt_mean.values); imax = np.fmax(imax, np.nan_to_num(d.ice_frac.values))
    nt = len(fs); valid = (nv == nt) & (imax <= 0.10)
    ub, vb, am = su/np.maximum(nv, 1), sv/np.maximum(nv, 1), sa/nt
    M = regions(lat, lon, am, valid)
    w = np.broadcast_to(np.cos(np.radians(lat))[:, None], sh)
    area = (6.371e6**2)*np.radians(np.abs(np.diff(lat)).mean())*np.radians(np.abs(np.diff(lon)).mean())*w
    # 1/4-deg grid for LF4
    if p == "allsat":
        lat4, lon4 = blk(lat[:, None]*np.ones(sh))[:, 0], blk(np.ones(sh)*lon[None, :])[0]
        ub4, vb4 = blk(ub), blk(vb); valid4 = blk(valid.astype(float)) == 1; am4 = blk(am)
    else:
        lat4, lon4, ub4, vb4, valid4, am4 = lat, lon, ub, vb, valid, am
    M4 = regions(lat4, lon4, am4, valid4); w4 = np.broadcast_to(np.cos(np.radians(lat4))[:, None], ub4.shape)
    rows = []
    for k, f in enumerate(fs):
        with xr.open_dataset(f) as d:
            u = d.u_mean.values.astype(float); v = d.v_mean.values.astype(float); s2 = d.var_intra.values.astype(float)
        eke = 0.5*((u-ub)**2+(v-vb)**2) + 0.5*s2
        if p == "allsat": u4, v4 = blk(u), blk(v)
        else: u4, v4 = u, v
        lf4 = 0.5*((u4-ub4)**2+(v4-vb4)**2)
        r = {"i": k}
        for n, msk in M.items():
            r[f"{n}_EKE"] = 1e4*np.sum(eke[msk]*w[msk])/np.sum(w[msk])
            r[f"{n}_LF4"] = 1e4*np.sum(lf4[M4[n]]*w4[M4[n]])/np.sum(w4[M4[n]])
        r["box35_65_EKEint_PJm"] = 1025*np.sum(eke[M["box35_65"]]*area[M["box35_65"]])/1e15
        rows.append(r)
    t = pd.DataFrame(rows); t["year"] = np.repeat(list(YRS), 12); t["month"] = np.tile(np.arange(1, 13), len(YRS))
    t["x"] = t.year + (t.month-0.5)/12
    info = {n: dict(n_pix=int(msk.sum()), frac_valid=float(msk.sum()/max(1, regions(lat, lon, am, np.ones(sh, bool))[n].sum()))) for n, msk in M.items()}
    return t, info
def one(args):
    lab, x, y, mo, ref = args
    o = R.trend_bcz(x, y, month=mo)
    return dict(series=lab, slope_per_decade=10*o["slope"], ci_lo=10*o["ci_lo"], ci_hi=10*o["ci_hi"], p=o["p"],
                mean=float(np.mean(y)), pct_per_decade=1000*o["slope"]/float(np.mean(y)),
                pct_ci_lo=1000*o["ci_lo"]/float(np.mean(y)), pct_ci_hi=1000*o["ci_hi"]/float(np.mean(y)),
                pct_of_ref=(1000*o["slope"]/ref if ref else np.nan), n=o["n"])
if __name__ == "__main__":
    S = {}
    for p in ("allsat", "twosat"):
        S[p], info = series(p); S[p].to_csv(OUT/f"m1_series_{p}.csv", index=False)
        pd.DataFrame(info).T.to_csv(OUT/f"m1_masks_{p}.csv"); print(p, info, flush=True)
    jobs = []
    PER = {"2001-2012": (2001, 2012), "2001-2020": (2001, 2020), "2001-2025": (2001, 2025)}
    if Y0 <= 1993:
        PER.update({"1993-2012": (1993, 2012), "1993-2020": (1993, 2020), "1993-2025": (1993, 2025)})
    for p in S:
        t = S[p]
        for c in [c for c in t.columns if c.endswith(("_EKE", "_LF4", "_PJm"))]:
            ref = float(t.loc[t.year <= Y0 + 8, c].mean())   # analogue of Hogg's first-9-year reference
            for per, (a, b) in PER.items():
                s = (t.year >= a) & (t.year <= b)
                jobs.append((f"{p}|{c}|{per}", t.x.values[s], t[c].values[s], t.month.values[s], ref))
    # product difference of mean-normalised series (percentage points per decade)
    A, B = S["allsat"], S["twosat"]
    for c in [c for c in A.columns if c.endswith(("_EKE", "_LF4"))]:
        for per, (a, b) in PER.items():
            s = (A.year >= a) & (A.year <= b)
            dnorm = 100*(A[c][s]/A[c][s].mean() - B[c][s]/B[c][s].mean())
            jobs.append((f"diff_all-two|{c}|{per}", A.x.values[s], dnorm.values + 100.0, A.month.values[s], None))
    with Pool(2) as pool:
        res = pool.map(one, jobs)
    df = pd.DataFrame(res)
    df[["product", "variable", "period"]] = df.series.str.split("|", expand=True)
    dmask = df["product"].str.startswith("diff")
    df.loc[dmask, ["pct_per_decade", "pct_ci_lo", "pct_ci_hi", "mean", "pct_of_ref"]] = np.nan  # slope is already in pct points
    df.to_csv(OUT/"m1_trends.csv", index=False); print(df.round(3).to_string())
