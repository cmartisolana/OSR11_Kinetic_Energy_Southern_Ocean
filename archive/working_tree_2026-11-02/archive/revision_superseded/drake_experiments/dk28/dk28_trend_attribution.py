"""dk28: how much of the 2001-2025 EKE increase can the wind, the SAM and ENSO account for? (Product 1 EKE, Product 2 wind
stress; 3 deg cells from dk22; ACC band; annual means 2001-2025.)
Congruent-trend method (e.g. Thompson et al. 2000 for SAM-congruent trends): for each driver X, regress detrended annual
EKE (% of mean) on the detrended driver at lag 1 yr (EKE lagging; lag 0 also computed) -> sensitivity b; the X-congruent
EKE trend = b * trend(X). Done (i) for the ACC band and sectors, with 95 % intervals from the regression standard error
(effective sample size) combined with the driver-trend interval, and (ii) per 3 deg cell for maps. Drivers: local eastward
wind stress (% of mean), SAM, ONI; also multiple regression on all three. An equilibrium upper bound for the wind uses
d ln EKE / d ln tau = 1 (EKE proportional to tau in eddy-saturated models; Meredith et al. 2012, Marshall et al. 2017).
Observed trends: least squares on annual means (as the maps), % of the mean per decade.
Usage: python dk28_trend_attribution.py <revision dir> <drake out dir> <dk22 out dir> <oni csv> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, DO, D22, ONI, OUT = (Path(a) for a in sys.argv[1:6]); OUT.mkdir(parents=True, exist_ok=True)
Z = np.load(D22 / "dk22_fields.npz"); EK, TAU, lat, lon = Z["EK"], Z["TAU"], Z["lat"], Z["lon"]
m1 = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.coarsen(latitude=8, longitude=8, boundary="trim").mean().values
acc1 = (m1 >= -0.6) & (m1 <= 0.2) & np.isfinite(EK).all(0) & np.isfinite(TAU).all(0)
w1 = np.cos(np.radians(lat))[:, None] * np.ones(lon.size)[None, :]
yrs = np.arange(2001, 2026); t = yrs - yrs.mean()
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv").merge(pd.read_csv(ONI, comment="#"), on=["year", "month"], how="left")
an = A.groupby("year").mean(numeric_only=True)
def ann(a): return a.reshape(25, 12, *a.shape[1:]).mean(1)
def lsq(y):
    y = np.asarray(y, float); tt = t.reshape((-1,) + (1,) * (y.ndim - 1)); m = y.mean(0)
    b = (tt * (y - m)).sum(0) / (t ** 2).sum(); res = y - m - tt * b
    return b, res
def trend_ci(y):
    b, res = lsq(y); n = len(y); r1 = np.corrcoef(res[1:], res[:-1])[0, 1]; ne = max(3, n * (1 - r1) / (1 + r1))
    se = np.sqrt((res ** 2).sum() / (ne - 2) / (t ** 2).sum()); return b, se, ne
def reg(x, y, L=1):
    xa, ya = (x[:25 - L], y[L:]) if L else (x, y); xa = xa - xa.mean(); ya = ya - ya.mean(); n = len(xa)
    b = (xa * ya).sum() / (xa ** 2).sum(); res = ya - b * xa
    a1 = np.corrcoef(xa[1:], xa[:-1])[0, 1]; b1 = np.corrcoef(ya[1:], ya[:-1])[0, 1]; ne = max(3, n * (1 - a1 * b1) / (1 + a1 * b1))
    se = np.sqrt((res ** 2).sum() / (ne - 2) / (xa ** 2).sum()); r = np.corrcoef(xa, ya)[0, 1]
    return b, se, r, 2 * stats.t.sf(abs(r) * np.sqrt((ne - 2) / max(1 - r ** 2, 1e-9)), ne - 2)
rows = []
SECT = {"ACC band": (-180, 180), "Pacific": (147, 290), "Indian": (20, 147), "Atlantic": (290, 380)}
EA, TA = ann(EK), ann(TAU)
for sec, (a_, b_) in SECT.items():
    l = np.mod(lon, 360); lm = ((l >= a_) & (l < b_)) | ((l + 360 >= a_) & (l + 360 < b_)) if sec != "ACC band" else np.ones(lon.size, bool)
    m = acc1 & lm[None, :]
    e = np.nansum(EA[:, m] * w1[m], 1) / w1[m].sum(); tau = np.nansum(TA[:, m] * w1[m], 1) / w1[m].sum()
    ep, tp = 100 * (e / e.mean() - 1), 100 * (tau / tau.mean() - 1)
    be, see, _ = trend_ci(ep); obs = 10 * be
    drivers = {"wind stress": tp, "SAM": an.SAM.values, "Nino 3.4": an.oni.values}
    for dn, x in drivers.items():
        bx, sex, _ = trend_ci(x); _, rx = lsq(x); _, ry = lsq(ep)
        for L in (0, 1):
            b, se, r, p = reg(rx, ry, L); cong = 10 * b * bx
            se_c = 10 * np.sqrt((b * sex) ** 2 + (se * bx) ** 2)
            rows.append(dict(region=sec, driver=dn, lag=L, obs_trend=obs, obs_lo=obs - 1.96 * 10 * see, obs_hi=obs + 1.96 * 10 * see,
                             driver_trend_dec=10 * bx, sens=b, sens_se=se, r=r, p=p, congruent=cong, cong_lo=cong - 1.96 * se_c, cong_hi=cong + 1.96 * se_c,
                             frac=cong / obs))
    rows.append(dict(region=sec, driver="wind stress, equilibrium EKE ~ tau", lag=np.nan, obs_trend=obs, driver_trend_dec=10 * trend_ci(tp)[0],
                     sens=1.0, congruent=10 * trend_ci(tp)[0], frac=10 * trend_ci(tp)[0] / obs))
    # multiple regression (lag 1) on the three drivers
    X = np.column_stack([lsq(v)[1][:24] for v in drivers.values()]); y = lsq(ep)[1][1:]
    bb, *_ = np.linalg.lstsq(X - X.mean(0), y - y.mean(), rcond=None); tr = np.array([10 * trend_ci(v)[0] for v in drivers.values()])
    R2 = 1 - ((y - y.mean() - (X - X.mean(0)) @ bb) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    rows.append(dict(region=sec, driver="all three (lag 1, multiple regression)", lag=1, obs_trend=obs, congruent=float(bb @ tr), frac=float(bb @ tr) / obs, r=np.sqrt(max(R2, 0))))
D = pd.DataFrame(rows); D.to_csv(OUT / "dk28_attribution.csv", index=False)
pd.set_option("display.width", 220); print(D[["region", "driver", "lag", "obs_trend", "driver_trend_dec", "sens", "r", "p", "congruent", "cong_lo", "cong_hi", "frac"]].round(2).to_string(index=False))
# maps (3 deg): observed EKE trend (%/dec), wind-congruent (local wind, lag 1), SAM- and ENSO-congruent, residual
def cg(a): return a.reshape(25, 10, 3, 120, 3).mean(axis=(2, 4))
E3, T3 = cg(EA), cg(TA); Ep = 100 * (E3 / E3.mean(0) - 1); Tp = 100 * (T3 / T3.mean(0) - 1)
bE, rE = lsq(Ep); bT, rT = lsq(Tp)
def regmap(xres, yres, L=1):
    xa, ya = xres[:25 - L] - xres[:25 - L].mean(0), yres[L:] - yres[L:].mean(0)
    return (xa * ya).sum(0) / (xa ** 2).sum(0)
out = dict(lat3=lat.reshape(-1, 3).mean(1), lon3=lon.reshape(-1, 3).mean(1), obs=10 * bE, wind_trend=10 * bT)
out["cong_wind"] = 10 * regmap(rT, rE) * bT
for dn, x in (("sam", an.SAM.values), ("oni", an.oni.values)):
    bx, rx = lsq(x); out[f"cong_{dn}"] = 10 * regmap(rx[:, None, None] * np.ones(Ep.shape[1:]), rE) * bx; out[f"{dn}_trend"] = 10 * bx
out["resid"] = out["obs"] - out["cong_wind"] - out["cong_sam"] - out["cong_oni"]
m3 = cg(np.broadcast_to(acc1, EA.shape).astype(float))[0] > 0.5
for k in ("obs", "cong_wind", "cong_sam", "cong_oni", "resid"):
    v = out[k][m3 & np.isfinite(out[k])]; print(f"map {k}: ACC-band median {np.median(v):+.1f} %/dec")
ok = m3 & np.isfinite(out["obs"]) & np.isfinite(out["wind_trend"])
print("pattern r(obs EKE trend, wind trend) =", np.corrcoef(out["obs"][ok], out["wind_trend"][ok])[0, 1].round(2))
out["accmask3"] = m3
np.savez(OUT / "dk28_maps.npz", **out)
