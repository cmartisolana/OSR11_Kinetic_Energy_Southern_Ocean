"""dk38: differences between 2018-2025 and 2001-2017 (Argo: 2004-2017), Product 1, Drake Passage domain.
(band) per 2 deg longitude between the SAF and PF (band, ice and weights as dk8): EKE and mean-flow KE (0.5|<u>|^2 of the
    annual-mean velocity), annual means; difference of the period means in % of the 2001-2025 mean; 95 % from Welch's t on the
    annual values.
(line) per 0.25 deg latitude along the mean LM Gould track (as dk25/dk35): cross-track velocity (cm/s), EKE (% of the mean),
    depth of sigma0 = 27.5 (dbar, Argo annual, dk33); annual means (velocity/EKE: years with >= 6 months); Welch 95 %.
Usage: python dk38_period_differences.py <revision dir> <dk5 stream edges csv> <dk24 npz> <dk33 npz> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, EDG, SEC, D33, OUT = (Path(a) for a in sys.argv[1:6]); OUT.mkdir(parents=True, exist_ok=True)
YR = np.repeat(np.arange(2001, 2026), 12); YY = np.arange(2001, 2026)
def welch(a, b):
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3: return np.nan, np.nan, np.nan
    d = b.mean() - a.mean(); se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    df = se ** 4 / ((a.var(ddof=1) / len(a)) ** 2 / (len(a) - 1) + (b.var(ddof=1) / len(b)) ** 2 / (len(b) - 1))
    return d, stats.t.ppf(0.975, df) * se, 2 * stats.t.sf(abs(d / se), df)
rows = []
# band
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
EA = np.array([np.nanmean(E[YR == y], 0) for y in YY])
UA, VA = (np.array([np.nanmean(X[YR == y], 0) for y in YY]) for X in (U, V)); KA = 0.5 * (UA ** 2 + VA ** 2) * 1e4
a = xr.open_dataset(REV / "outputs" / "clim_allsat.nc").adt.interp(latitude=lat, longitude=lon).values
band = np.isfinite(a) & (a >= LEV_PF) & (a <= LEV_SAF) & (np.nanmax(IC, 0) <= 0.10); w = np.cos(np.radians(lat))[:, None] * band
e_, l_ = YY <= 2017, YY >= 2018
for l0 in np.arange(-86, -34, 2.0):
    m = (lon >= l0) & (lon < l0 + 2); ww = w[:, m]
    if ww.sum() == 0: continue
    for nm, F in (("eke", EA), ("mke", KA)):
        y = np.nansum(F[:, :, m] * ww, axis=(1, 2)) / np.sum(ww * np.isfinite(F[0][:, m])); mu = y.mean()
        d, ci, p = welch(y[e_], y[l_]); rows.append(dict(kind=f"band_{nm}", x=l0 + 1, mean=mu, early=y[e_].mean(), late=y[l_].mean(), diff=d, diff_pct=100 * d / mu, ci_pct=100 * ci / mu, p=p))
# line
Z = np.load(SEC); lo, la = Z["lon"], Z["lat"]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
ice = Z["ice"] > 0.10; u, v = Z["u"], Z["v"]
un = (u * nx + v * ny) * 100; un[ice] = np.nan
ek = (0.5 * ((u - np.nanmean(u, 0)) ** 2 + (v - np.nanmean(v, 0)) ** 2) + 0.5 * Z["vi"]) * 1e4; ek[ice] = np.nan
def annual(y):
    df = pd.DataFrame(dict(year=YR, v=y)).dropna(); g = df.groupby("year").v; an = g.mean()[g.count() >= 6]
    out = np.full(YY.size, np.nan); out[np.searchsorted(YY, an.index.values)] = an.values; return out
S33 = np.load(D33); yA = S33["years"]; ZA = S33["p27.5"]; laA = S33["tlat"]
for c in np.arange(-63.5, -54.99, 0.25):
    s = np.abs(la - c) <= 0.125
    if not s.any(): continue
    for nm, F in (("vel", un), ("eke", ek)):
        y = annual(np.nanmean(F[:, s], 1)); mu = np.nanmean(y); d, ci, p = welch(y[e_], y[l_])
        rows.append(dict(kind=f"line_{nm}", x=c, mean=mu, early=np.nanmean(y[e_]), late=np.nanmean(y[l_]), diff=d, diff_pct=100 * d / mu, ci_pct=100 * ci / mu, ci=ci, p=p))
    sA = np.abs(laA - c) <= 0.125
    if sA.any():
        y = np.nanmean(ZA[:, sA], 1); ea, lt = (yA >= 2004) & (yA <= 2017), yA >= 2018
        if np.isfinite(y).all():
            d, ci, p = welch(y[ea], y[lt]); rows.append(dict(kind="line_iso", x=c, mean=y.mean(), early=y[ea].mean(), late=y[lt].mean(), diff=d, ci=ci, p=p))
D = pd.DataFrame(rows); D.to_csv(OUT / "dk38_period_differences.csv", index=False)
pd.set_option("display.width", 200); print(D.groupby("kind")[["diff", "diff_pct", "p"]].describe().round(2).to_string())
