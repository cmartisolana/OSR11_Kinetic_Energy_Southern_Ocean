"""B1 pilot: Pacific-Antarctic Ridge (PAR) standing meander, 2001-2025.

Regions are fixed in advance from the published, bathymetry-based segmentation
of Liu et al. (2026, Table 1), latitude 60-48S:
   upstream   = Upper-Ridge section   150.00-135.13 W
   meander    = Lower-Ridge section   135.13-109.38 W
   downstream = Downstream section    109.38- 80.00 W
Pathway: ADT contour (relative to the box-mean ADT) at the climatological jet-
core level of the meander box, and +/-0.10 m (b1_meander.py).
Upstream jet speed: per longitude in the upstream box, the maximum monthly-mean
surface geostrophic speed between 60S and 48S, averaged over longitudes.
Downstream EKE: box mean, fixed 2001-2025 reference (as v2) and 36-month
running reference.
Gate (decided before trends are looked at): pathway found in >= 90 % of months
at the core level, and the minimum detectable 25-yr change in the curvature
metric (2 x calibrated SE x 25 yr) smaller than its interannual SD.
Usage: python b1_pilot.py <revision/data dir> <out_dir>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b1_meander as B  # noqa: E402
import rstats as R      # noqa: E402

data, outdir = Path(sys.argv[1]), Path(sys.argv[2])
outdir.mkdir(parents=True, exist_ok=True)
REG = dict(lon0=-152.0, lon1=-78.0, lat0=-62.0, lat1=-45.0)
UP = (-150.0, -135.13); MEA = (-135.13, -109.38); DOWN = (-109.38, -80.0); LATB = (-60.0, -48.0)
files = [data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
t["x"] = t.year + (t.month - 0.5) / 12


def sub(ds):
    return ds.sel(longitude=slice(REG["lon0"], REG["lon1"]), latitude=slice(REG["lat0"], REG["lat1"]))


stack = {v: [] for v in ("u_mean", "v_mean", "adt_mean", "var_intra")}
for f in files:
    with xr.open_dataset(f) as ds:
        s = sub(ds)
        for v in stack:
            stack[v].append(s[v].values.astype(np.float64))
lat = s.latitude.values.astype(float); lon = s.longitude.values.astype(float)
U, V, A, S2 = (np.array(stack[v]) for v in ("u_mean", "v_mean", "adt_mean", "var_intra"))
ub, vb, Ab = U.mean(0), V.mean(0), A.mean(0)

mbox = dict(lon0=MEA[0], lon1=MEA[1], lat0=-62.0, lat1=-45.0)
core, core_spd = B.core_level(Ab, ub, vb, lat, lon, mbox)
levels = {"core": core, "core-0.10": core - 0.10, "core+0.10": core + 0.10}
print("climatological core level (relative ADT, m):", round(core, 3), "speed", round(core_spd, 3))

la_m = (lat >= LATB[0]) & (lat <= LATB[1])
def lonm(a, b):
    return (lon >= a) & (lon <= b)

# running 36-month reference
def runmean(X, half=18):
    out = np.empty_like(X)
    for k in range(X.shape[0]):
        out[k] = X[max(0, k - half):min(X.shape[0], k + half)].mean(0)
    return out
U3, V3 = runmean(U), runmean(V)
w = np.cos(np.radians(lat))[:, None]

rows = []
for k in range(len(files)):
    rec = {"i": k}
    rel, la, lo = B.box_relative(A[k], lat, lon, mbox)
    for name, lev in levels.items():
        line = B.spanning_contour(rel, la, lo, lev)
        for sm in (50.0, 100.0, 150.0):
            m = None if line is None else B.path_metrics(line, smooth=sm)
            for key in ("length_km", "sinuosity", "amp_km", "rms_kappa", "int_abs_kappa"):
                rec[f"{name}_s{int(sm)}_{key}"] = np.nan if m is None else m[key]
        rec[f"{name}_meanlat"] = np.nan if line is None else float(np.mean(line[:, 1]))
    spd = np.hypot(U[k], V[k])
    cols = lonm(*UP)
    rec["up_core_speed"] = float(np.nanmean(np.nanmax(spd[np.ix_(la_m, cols)], axis=0)))
    cols_m = lonm(*MEA)
    rec["mea_core_speed"] = float(np.nanmean(np.nanmax(spd[np.ix_(la_m, cols_m)], axis=0)))
    dm = la_m[:, None] & lonm(*DOWN)[None, :]
    eke = 0.5 * ((U[k] - ub) ** 2 + (V[k] - vb) ** 2) + 0.5 * S2[k]
    eke3 = 0.5 * ((U[k] - U3[k]) ** 2 + (V[k] - V3[k]) ** 2) + 0.5 * S2[k]
    for nm, E in (("down_eke", eke), ("down_eke3", eke3)):
        ok = dm & np.isfinite(E)
        rec[nm] = float(1e4 * np.sum(E[ok] * np.broadcast_to(w, E.shape)[ok]) / np.sum(np.broadcast_to(w, E.shape)[ok]))
    mm = la_m[:, None] & cols_m[None, :]
    eke_m = eke[mm]; rec["mea_eke"] = float(1e4 * np.nanmean(eke_m))
    rows.append(rec)
P = pd.concat([t, pd.DataFrame(rows)], axis=1)
P.to_csv(outdir / "b1_pilot_series.csv", index=False)

cov = {c: float(np.isfinite(P[c]).mean()) for c in P.columns if c.endswith("_s100_int_abs_kappa")}
print("pathway coverage (s100):", cov)

out = []
for c in [c for c in P.columns if c not in ("year", "month", "x", "i")]:
    y = P[c].values
    if np.isfinite(y).mean() < 0.5:
        continue
    o = R.trend_calibrated(P.x.values, y, month=P.month.values, n_outer=100, n_boot=300)
    yy = R.deseasonalise(y, P.month.values)
    ann_sd = float(pd.Series(yy).groupby(P.year.values).mean().std())
    se = (o["ci_cal_hi"] - o["ci_cal_lo"]) / (2 * o["z_crit_cal"]) if np.isfinite(o["ci_cal_hi"]) else np.nan
    o.update(series=c, coverage=float(np.isfinite(y).mean()), mean=float(np.nanmean(y)),
             interannual_sd=ann_sd, mdc_25yr=float(2 * se * 25) if np.isfinite(se) else np.nan)
    out.append(o)
pd.DataFrame(out).to_csv(outdir / "b1_pilot_trends.csv", index=False)

# pre-specified lagged associations (0-24 months), deseasonalised, detrended, 12-month low-passed
def prep(y):
    y = R.deseasonalise(y, P.month.values)
    ok = np.isfinite(y)
    b, a = R.theil_sen(P.x.values[ok], y[ok])
    r = y - (a + b * P.x.values)
    return pd.Series(r).interpolate(limit=2).rolling(12, center=True, min_periods=9).mean().values

chain = [("up_core_speed", "core_s100_int_abs_kappa"), ("core_s100_int_abs_kappa", "down_eke"),
         ("up_core_speed", "down_eke"), ("core_s100_length_km", "down_eke")]
lr = []
for a_, b_ in chain:
    xa, yb = prep(P[a_].values), prep(P[b_].values)
    for lag in range(0, 25, 3):
        xx = xa[:len(xa) - lag] if lag else xa
        yy = yb[lag:]
        ok = np.isfinite(xx) & np.isfinite(yy)
        r = float(np.corrcoef(xx[ok], yy[ok])[0, 1])
        r1 = lambda z: np.corrcoef(z[:-1], z[1:])[0, 1]
        neff = ok.sum() * (1 - r1(xx[ok]) * r1(yy[ok])) / (1 + r1(xx[ok]) * r1(yy[ok]))
        from scipy import stats
        tt = r * np.sqrt(max(neff - 2, 1) / max(1 - r * r, 1e-9))
        lr.append(dict(x=a_, y=b_, lag_months=lag, r=r, n_eff=float(neff), p=float(2 * stats.t.sf(abs(tt), max(neff - 2, 1)))))
pd.DataFrame(lr).to_csv(outdir / "b1_pilot_lags.csv", index=False)
print("done")
