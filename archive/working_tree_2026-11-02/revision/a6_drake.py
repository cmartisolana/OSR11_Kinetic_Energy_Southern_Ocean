"""A6: Drake Passage transport proxies, energy and wind stress.

Usage: python a6_drake.py <revision/data dir> <repository outputs/trends dir> <out_dir>

Altimetric surface-transport proxy T(t) = integral of the transect-normal
monthly-mean geostrophic velocity along the transect [m2 s-1]. Ice handling:
a transect point is excluded in a given month only if ice_frac > 0.10 that
month (v2 excluded every point ever ice-covered: 69 of 121 points).
Transect variants (sensitivity): v2 line 65W55S-60W63S; same line with
endpoints moved by +/-0.5 and +/-1 deg; an approximate LM Gould line
(65W55S-58.5W62S; approximate endpoints); and the ADT difference between the
v2 endpoints (north minus south).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

data, trends_dir, outdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
outdir.mkdir(parents=True, exist_ok=True)
files = [data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
t["x"] = t.year + (t.month - 0.5) / 12
R_E = 6.371e6

TRANSECTS = {"v2": (-65.0, -55.0, -60.0, -63.0)}
for dl in (-1.0, -0.5, 0.5, 1.0):
    TRANSECTS[f"v2_lon{dl:+.1f}"] = (-65.0 + dl, -55.0, -60.0 + dl, -63.0)
TRANSECTS["v2_south+0.5"] = (-65.0, -55.0, -60.0, -62.5)
TRANSECTS["v2_north-0.5"] = (-65.0, -55.5, -60.0, -63.0)
TRANSECTS["gould_approx"] = (-65.0, -55.0, -58.5, -62.0)


def geometry(lon0, lat0, lon1, lat1, n=121):
    tt = np.linspace(0, 1, n)
    lo = lon0 + tt * (lon1 - lon0); la = lat0 + tt * (lat1 - lat0)
    latm = np.radians(0.5 * (la[:-1] + la[1:]))
    dx = R_E * np.cos(latm) * np.radians(np.diff(lo)); dy = R_E * np.radians(np.diff(la))
    seg = np.hypot(dx, dy)
    ds = np.r_[seg[0] / 2, 0.5 * (seg[:-1] + seg[1:]), seg[-1] / 2]
    latr = np.radians(0.5 * (lat0 + lat1))
    DX = R_E * np.cos(latr) * np.radians(lon1 - lon0); DY = R_E * np.radians(lat1 - lat0)
    L = np.hypot(DX, DY)
    nx, ny = DY / L, -DX / L            # unit normal (eastward-positive for a N-S line)
    if nx < 0:
        nx, ny = -nx, -ny
    return lo, la, ds, nx, ny


rows = []
box = dict(lon0=-80, lon1=-52, lat0=-63, lat1=-54)
ub_sum = None
for k, f in enumerate(files):
    ds = xr.open_dataset(f).sel(longitude=slice(-82, -50), latitude=slice(-65, -52))
    if ub_sum is None:
        ub_sum = np.zeros(ds.u_mean.shape); vb_sum = np.zeros(ds.u_mean.shape)
    ub_sum += np.nan_to_num(ds.u_mean.values); vb_sum += np.nan_to_num(ds.v_mean.values)
    rec = {"i": k}
    for name, (a, b, c, d) in TRANSECTS.items():
        lo, la, dsx, nx, ny = geometry(a, b, c, d)
        pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
        u = ds.u_mean.interp(**pts).values; v = ds.v_mean.interp(**pts).values
        ice = ds.ice_frac.interp(**pts).values
        vn = u * nx + v * ny
        ok = np.isfinite(vn) & ~(np.nan_to_num(ice) > 0.10)
        rec[f"T_{name}"] = float(np.sum(vn[ok] * dsx[ok]))
        rec[f"fvalid_{name}"] = float(ok.mean())
        if name in ("v2", "gould_approx"):               # sensitivity: no ice masking at all
            ok2 = np.isfinite(vn)
            rec[f"Tnoice_{name}"] = float(np.sum(vn[ok2] * dsx[ok2]))
    e = ds.adt_mean.interp(longitude=xr.DataArray([-65.0, -60.0]), latitude=xr.DataArray([-55.0, -63.0])).values
    rec["dADT"] = float(e[0] - e[1])
    rows.append(rec)
    ds.close()
T = pd.concat([t, pd.DataFrame(rows)], axis=1)
T.to_csv(outdir / "a6_transport_series.csv", index=False)

out = []
def add(name, x, y, monthly=True, month=None, unit=""):
    o = R.trend_calibrated(x, y, month=month, monthly=monthly, n_outer=100, n_boot=300)
    m = float(np.nanmean(y)); s = float(np.nanstd(y))
    o.update(series=name, unit=unit, mean=m, pct_per_decade=1000 * o["slope"] / m if m else np.nan,
             sigma_per_yr=o["slope"] / s, n=int(np.isfinite(y).sum()), x0=float(np.nanmin(x)), x1=float(np.nanmax(x)))
    out.append(o); print(name, round(o["slope"], 4), round(o["p_cal"], 3), flush=True)

for c in [c for c in T.columns if c.startswith("T_") or c.startswith("Tnoice_")] + ["dADT"]:
    add(f"{c}_2001-2025", T.x.values, T[c].values, month=T.month.values, unit="m2 s-1" if c != "dADT" else "m")
    sel = (T.year >= 2005) & (T.year <= 2019)
    add(f"{c}_2005-2019", T.x.values[sel], T[c].values[sel], month=T.month.values[sel])

# in-situ (LM Gould ADCP, 0-760 m total transport): annual means, 2005-2019
obs = pd.read_csv(trends_dir / "drake_passage_total_transport_observations.csv", parse_dates=["date"])
obs = obs.dropna(subset=["total_transport_sv"])
obs["year"] = obs.date.dt.year; obs["month"] = obs.date.dt.month
ann = obs.groupby("year").total_transport_sv.mean()
add("insitu_annual_2005-2019", ann.index.values.astype(float) + 0.5, ann.values, monthly=False, unit="Sv")
# date-matched comparison (altimetry monthly value of the transect month)
mm = obs.merge(T[["year", "month", "T_v2", "T_gould_approx"]], on=["year", "month"], how="left")
agree = {}
for c in ["T_v2", "T_gould_approx"]:
    a = mm[c].values; b = mm.total_transport_sv.values; ok = np.isfinite(a) & np.isfinite(b)
    za = (a[ok] - a[ok].mean()) / a[ok].std(); zb = (b[ok] - b[ok].mean()) / b[ok].std()
    r = float(np.corrcoef(za, zb)[0, 1])
    # annual-mean correlation (fewer, more independent samples)
    g = mm[ok].groupby("year")[[c, "total_transport_sv"]].mean()
    ra = float(np.corrcoef(g[c], g.total_transport_sv)[0, 1])
    agree[c] = dict(n_transects=int(ok.sum()), r_transects=r, rmsd_z=float(np.sqrt(np.mean((za - zb) ** 2))),
                    n_years=len(g), r_annual=ra)
pd.DataFrame(agree).T.to_csv(outdir / "a6_agreement.csv")

# Drake-box energy (fixed 2001-2025 reference), periods
ubar = ub_sum / len(files); vbar = vb_sum / len(files)
E = []
for k, f in enumerate(files):
    ds = xr.open_dataset(f).sel(longitude=slice(-82, -50), latitude=slice(-65, -52))
    la = ds.latitude.values; lo = ds.longitude.values
    m = (lo[None, :] >= box["lon0"]) & (lo[None, :] <= box["lon1"]) & (la[:, None] >= box["lat0"]) & (la[:, None] <= box["lat1"])
    w = np.cos(np.radians(la))[:, None] * np.ones_like(m, float)
    U, V, S2 = ds.u_mean.values, ds.v_mean.values, ds.var_intra.values
    ke = 0.5 * (U ** 2 + V ** 2); eke = 0.5 * ((U - ubar) ** 2 + (V - vbar) ** 2) + 0.5 * S2
    ok = m & np.isfinite(ke)
    E.append(dict(KE=1e4 * np.sum(ke[ok] * w[ok]) / np.sum(w[ok]), EKE=1e4 * np.sum(eke[ok] * w[ok]) / np.sum(w[ok])))
    ds.close()
E = pd.concat([t, pd.DataFrame(E)], axis=1)
E.to_csv(outdir / "a6_drake_energy.csv", index=False)
for v in ["KE", "EKE"]:
    for per, sel in [("2001-2025", E.year <= 2025), ("2001-2024", E.year <= 2024), ("2005-2019", (E.year >= 2005) & (E.year <= 2019))]:
        add(f"Drake_{v}_{per}", E.x.values[sel], E[v].values[sel], month=E.month.values[sel], unit="cm2 s-2")

# wind stress (optional until extract_wind.py has been run)
wf = data / "wind_monthly.nc"
if not wf.exists():
    wf = data / "wind_taux.nc"          # eastward_stress-only float32 subset of wind_monthly.nc
if wf.exists():
    W = xr.open_dataset(wf)
    tau = W["eastward_stress"] if "eastward_stress" in W else None
    if tau is not None:
        wt = pd.to_datetime(W.time.values)
        latw = W.latitude.values; lonw = W.longitude.values
        adt_mean = None
        # ACC band from the all-sat mean ADT (regridded to the wind grid)
        acc_src = xr.open_dataset(files[0]).adt_mean * 0
        S = []
        for name, (lo0, lo1, la0, la1) in {"Drake_box": (-80, -52, -63, -54), "Pacific_45_60S": (147 - 360, -70, -60, -45),
                                           "Circumpolar_45_60S": (-180, 180, -60, -45)}.items():
            sub = tau.sel(latitude=slice(la0, la1))
            if name == "Pacific_45_60S":
                sub = sub.where((sub.longitude >= 147) | (sub.longitude <= -70), drop=True)
            else:
                sub = sub.sel(longitude=slice(lo0, lo1))
            wts = np.cos(np.radians(sub.latitude))
            s = sub.weighted(wts).mean(("latitude", "longitude")).values
            S.append(pd.Series(s, index=wt, name=name))
        S = pd.concat(S, axis=1)
        S["year"] = S.index.year; S["month"] = S.index.month; S["x"] = S.year + (S.month - 0.5) / 12
        S.to_csv(outdir / "a6_wind_stress_series.csv")
        for c in ["Drake_box", "Pacific_45_60S", "Circumpolar_45_60S"]:
            for per, sel in [("2001-2025", (S.year >= 2001) & (S.year <= 2025)), ("2008-2025", (S.year >= 2008) & (S.year <= 2025))]:
                add(f"tau_{c}_{per}", S.x.values[sel], S[c].values[sel], month=S.month.values[sel], unit="N m-2")
        # pre-specified lags 0-3 yr, annual means, detrended
        A = S[(S.year >= 2001) & (S.year <= 2025)].groupby("year").mean()
        EA = E.groupby("year")[["EKE"]].mean(); TA = T.groupby("year")[["T_v2"]].mean()
        J = A.join(EA).join(TA)
        def detr(y):
            x = np.arange(y.size); b, a = R.theil_sen(x, y); return y - (a + b * x)
        lagrows = []
        for c in ["Circumpolar_45_60S", "Pacific_45_60S", "Drake_box"]:
            for tgt in ["EKE", "T_v2"]:
                for lag in range(4):
                    xa = detr(J[c].values); ya = detr(J[tgt].values)
                    if lag:
                        xa, ya = xa[:-lag], ya[lag:]
                    r = float(np.corrcoef(xa, ya)[0, 1])
                    r1x = np.corrcoef(xa[:-1], xa[1:])[0, 1]; r1y = np.corrcoef(ya[:-1], ya[1:])[0, 1]
                    neff = xa.size * (1 - r1x * r1y) / (1 + r1x * r1y)
                    from scipy import stats
                    tstat = r * np.sqrt(max(neff - 2, 1) / max(1 - r * r, 1e-9))
                    p = 2 * stats.t.sf(abs(tstat), max(neff - 2, 1))
                    lagrows.append(dict(wind=c, target=tgt, lag_yr=lag, r=r, n=xa.size, n_eff=neff, p=p))
        pd.DataFrame(lagrows).to_csv(outdir / "a6_lags.csv", index=False)
pd.DataFrame(out).to_csv(outdir / "a6_trends.csv", index=False)
print("done")
