"""E7: trends of tracked-eddy properties on 2 deg cells south of 35 S (candidate Figure 4, map version; 26 Oct 2026).

Atlas META4.0 DT all-sat subset, tracked observations (track > 0), 2001-2022 (full years; 2023 is partial).
Per 2 deg cell and year: eddies present per day, mean amplitude, mean rotational speed, sum(amplitude^2) per day
(cells with >= 15 valid years, a year valid if >= 10 observations). Trend = Theil-Sen on the annual series, % of the
cell mean per decade; p from a t test on the OLS slope with the effective sample size n (1 - r1) / (1 + r1) (AR(1)
of the residuals; draft statistic, NOT the manuscript's calibrated test, see note); BH-FDR 10 % over cells per variable.
Usage: python e7_maps.py <subset_dir> <out_dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, netCDF4 as nc
from scipy import stats

sub, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
Y0, Y1 = 2001, 2022; RES = 3.0


def load(pol):
    ds = nc.Dataset(sub / f"eddies_so_{pol}_all.nc")
    cols = ("track", "time", "latitude", "longitude", "amplitude", "speed_average")
    df = pd.DataFrame({v: np.asarray(ds.variables[v][:]) for v in cols})
    df["time"] = pd.Timestamp("1950-01-01") + pd.to_timedelta(df.time.astype(int), unit="D")
    df["lon"] = ((df.longitude + 180) % 360) - 180
    return df[(df.track > 0) & (df.time.dt.year >= Y0) & (df.time.dt.year <= Y1)]


E = pd.concat([load("anticyclonic"), load("cyclonic")], ignore_index=True)
E["year"] = E.time.dt.year
E["i"] = np.floor((E.latitude + 80) / RES).astype(int); E["j"] = np.floor((E.lon + 180) / RES).astype(int)
E = E[(E.i >= 0) & (E.i < int(45 / RES))]
ny = Y1 - Y0 + 1; ndays = np.array([pd.Timestamp(f"{y}-12-31").dayofyear for y in range(Y0, Y1 + 1)])
g = E.groupby(["i", "j", "year"])
agg = pd.DataFrame({"n": g.size(), "amp": g.amplitude.mean(), "spd": g.speed_average.mean(), "amp2": g.amplitude.apply(lambda a: np.sum(a ** 2))}).reset_index()
agg["present_per_day"] = agg.n / ndays[agg.year - Y0]; agg["amp2_per_day"] = agg.amp2 / ndays[agg.year - Y0]
agg.loc[agg.n < 30, ["amp", "spd"]] = np.nan
ni, nj = int(45 / RES), int(360 / RES)
lat_c = -80 + RES * (np.arange(ni) + 0.5); lon_c = -180 + RES * (np.arange(nj) + 0.5)
VARS = ["present_per_day", "amp", "spd", "amp2_per_day"]
x = np.arange(Y0, Y1 + 1).astype(float)
res = {v: {k: np.full((ni, nj), np.nan) for k in ("pct", "p", "mean")} for v in VARS}
for (i, j), d in agg.groupby(["i", "j"]):
    for v in VARS:
        s = pd.Series(d[v].values, index=d.year).reindex(range(Y0, Y1 + 1)).values.astype(float)
        if v in ("present_per_day", "amp2_per_day"):
            s = np.nan_to_num(s)       # no eddies = 0
        ok = np.isfinite(s)
        if ok.sum() < 15 or np.nanmean(s) <= 0: continue
        if v in ("present_per_day", "amp2_per_day") and np.mean(s) < (1.0 if v == "present_per_day" else 0): continue
        xs, ys = x[ok], s[ok]
        slope = stats.theilslopes(ys, xs)[0]
        b, a = np.polyfit(xs, ys, 1); r = ys - (a + b * xs)
        r1 = np.corrcoef(r[:-1], r[1:])[0, 1] if r.std() > 0 else 0.0
        r1 = max(min(r1, 0.95), -0.95)
        neff = max(len(ys) * (1 - r1) / (1 + r1), 3.0)
        se = np.sqrt(np.sum(r ** 2) / (len(ys) - 2) / np.sum((xs - xs.mean()) ** 2)) * np.sqrt((len(ys) - 2) / (neff - 2))
        t = b / se if se > 0 else 0.0
        p = 2 * stats.t.sf(abs(t), neff - 2)
        m = ys.mean()
        res[v]["pct"][i, j] = 1000 * slope / m; res[v]["p"][i, j] = p; res[v]["mean"][i, j] = m
import xarray as xr
ds = xr.Dataset({f"{v}_{k}": (("lat", "lon"), res[v][k]) for v in VARS for k in ("pct", "p", "mean")}, coords={"lat": lat_c, "lon": lon_c})
for v in VARS:
    p = res[v]["p"].ravel(); ok = np.isfinite(p); f = np.zeros_like(p, bool)
    if ok.sum():
        q = p[ok]; o = np.argsort(q); m = ok.sum(); thr = 0.10 * (np.arange(1, m + 1) / m)
        below = q[o] <= thr; k = np.max(np.where(below)[0]) + 1 if below.any() else 0
        sel = np.zeros(m, bool); sel[o[:k]] = True; f[np.where(ok)[0]] = sel
    ds[f"{v}_fdr"] = (("lat", "lon"), f.reshape(ni, nj))
    print(v, "cells", int(ok.sum()), "p<0.05", int((p[ok] < 0.05).sum()), "FDR", int(f.sum()), "median %/dec", round(float(np.nanmedian(res[v]["pct"])), 1), flush=True)
ds.to_netcdf(out / "e7_cell_trends.nc")
print("done")
