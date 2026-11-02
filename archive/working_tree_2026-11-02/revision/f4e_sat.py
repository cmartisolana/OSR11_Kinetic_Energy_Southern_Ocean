"""F4e: eddy-saturation diagnostics for Figure 4d (28 Oct 2026). Specification fixed before results.
Common cells: f4m domain (ACC band + 100 km, ice removed, sigma0 = 27.3 slope defined in all 264 months) with EKE in every month.
Domain means (cos-lat) of monthly series:
  wind  Product 2 eastward stress (2001-2025)
  KE    Product 1, 0.5 |U_m|^2 of the monthly-mean surface geostrophic velocity (2001-2025)
  EKE   Product 1, 0.5 |U_m - U_bar|^2 + 0.5 var_intra, U_bar = 2001-2025 mean (2001-2025), ice pixels dropped
  slope sigma0 = 27.3 signed cross-stream slope (f4m_p273_monthly.npz, 2004-2025)
D1 response per unit wind, 2004-2025 (common period): relative trends (rstats.trend_bcz monthly, % of the mean per decade) and the
   ratio trend_X / trend_wind; 95 % interval of the ratio from a joint moving-block residual bootstrap (blocks of 24 months, the same
   block draws for X and the wind, 2000 replicates; Theil-Sen slopes of the deseasonalised series).
D2 lagged response: detrended annual means, r(X_y, wind_{y-L}) for L = 0..3 yr (wind leading); EKE and KE use 2004-2025 for X so
   that all three share the same years and the wind lags stay within 2001-2025. Significance from 2000 phase-randomised surrogates
   of the detrended annual wind (2001-2025), then the same lag alignment: pointwise 95 % envelope and family-wise p over the lags.
Usage: python f4e_sat.py <revision dir> <f4 out dir> <f4m out dir> <out dir>
"""
import sys, glob
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

REV, F4, FM, OUT = (Path(a) for a in sys.argv[1:5]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402

M = xr.open_dataset(F4 / "f4_tilt_maps.nc"); lat, lon = M.lat.values, M.lon.values
LON, LAT = np.meshgrid(lon, lat); w = np.cos(np.radians(LAT)); KM = 111.195
la_e = np.r_[lat - 0.5, lat[-1] + 0.5]; lo_e = np.r_[lon - 0.5, lon[-1] + 0.5]
C = np.load(FM / "f4m_p273_monthly.npz"); YYs, MMs, PM, lonr = C["year"], C["month"], C["p"], C["lon"]
L0 = np.where(lonr > 180, lonr - 360, lonr); o = np.argsort(L0); adt = M["adt"].values
def grad(F):
    dy = np.gradient(F, lat, axis=-2) / KM * 100
    dx = (np.roll(F, -1, axis=-1) - np.roll(F, 1, axis=-1)) / 2.0 / (KM * np.cos(np.radians(lat))[:, None]) * 100
    return dx, dy
ax_, ay_ = grad(adt); an = np.hypot(ax_, ay_)
dx, dy = grad(PM[:, :, o].astype(float)); SL = dx * ax_ / an + dy * ay_ / an
dom = (M["buffer100"].values > 0) & (M["ice"].values == 0) & np.isfinite(M["signed_27p3_pct"].values) & np.isfinite(SL).all(0)

W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress.sel(latitude=slice(lat[0] - 0.5, lat[-1] + 0.5))
t = pd.DatetimeIndex(W.time.values); YY, MMo = t.year.values, t.month.values
il = np.digitize(W.latitude.values, la_e) - 1; jl = np.digitize(W.longitude.values, lo_e) - 1; A = W.values
WC = np.full((A.shape[0], lat.size, lon.size), np.nan)
for i in range(lat.size):
    sub = A[:, il == i, :]
    if sub.shape[1] == 0: continue
    for j in range(lon.size): WC[:, i, j] = np.nanmean(sub[:, :, jl == j], axis=(1, 2))
files = sorted(glob.glob(str(REV / "data" / "monthly_allsat" / "*" / "*.nc"))); assert len(files) == 300 and WC.shape[0] == 300
with xr.open_dataset(files[0]) as d0: la8, lo8 = d0.latitude.values, d0.longitude.values
ice8 = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
su = np.zeros((la8.size, lo8.size)); sv = np.zeros_like(su)
for f in files:
    with xr.open_dataset(f) as d: su += d.u_mean.values; sv += d.v_mean.values
ub, vb = su / 300, sv / 300; i8 = np.digitize(la8, la_e) - 1; j8 = np.digitize(lo8, lo_e) - 1
EC = np.full((300, lat.size, lon.size), np.nan); KC = np.full_like(EC, np.nan)
def tocell(e8, out, k):
    for i in range(lat.size):
        ri = i8 == i
        if not ri.any(): continue
        sub = e8[ri, :]
        out[k, i, :] = [np.nanmean(sub[:, j8 == j]) if np.isfinite(sub[:, j8 == j]).any() else np.nan for j in range(lon.size)]
for k, f in enumerate(files):
    with xr.open_dataset(f) as d:
        u, v, vi = d.u_mean.values, d.v_mean.values, d.var_intra.values
    e8 = (0.5 * ((u - ub) ** 2 + (v - vb) ** 2) + 0.5 * vi) * 1e4; k8 = 0.5 * (u ** 2 + v ** 2) * 1e4
    e8[ice8] = np.nan; k8[ice8] = np.nan; tocell(e8, EC, k); tocell(k8, KC, k)
dom &= np.isfinite(EC).all(0) & np.isfinite(KC).all(0) & np.isfinite(WC).all(0)
print("common cells", int(dom.sum()), flush=True)
def dm(F): return np.array([np.sum(F[q][dom] * w[dom]) / np.sum(w[dom]) for q in range(F.shape[0])])
S = {"wind": dm(WC), "KE": dm(KC), "EKE": dm(EC)}; s_sl = dm(SL)
x = YY + (MMo - 0.5) / 12
pd.DataFrame({"year": YY, "month": MMo, "x": x, **S, "slope": np.r_[np.full(36, np.nan), s_sl]}).to_csv(OUT / "f4e_series.csv", index=False)

# D1
c4 = YY >= 2004; X = {k: v[c4] for k, v in S.items()}; X["slope"] = s_sl; x4, m4 = x[c4], MMo[c4]
rows = []; TR = {}
for k, y in X.items():
    r = R.trend_bcz(x4, y, month=m4, n_outer=300); mu = y.mean(); TR[k] = 1000 * r["slope"] / mu
    rows.append(dict(test="D1", series=k, pct=TR[k], pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"]))
    print(rows[-1], flush=True)
rng = np.random.default_rng(20261028); n = x4.size; Lb = 24; nb = int(np.ceil(n / Lb))
D = {k: R.deseasonalise(y, m4) for k, y in X.items()}
fit = {}; res = {}
for k, y in D.items():
    b, a = R.theil_sen(x4, y); fit[k] = a + b * x4; res[k] = y - fit[k]
BR = {k: [] for k in X if k != "wind"}; BW = []
for _ in range(2000):
    st = rng.integers(0, n - Lb + 1, nb); idx = np.concatenate([np.arange(s_, s_ + Lb) for s_ in st])[:n]
    bw = R.theil_sen(x4, fit["wind"] + res["wind"][idx])[0] / D["wind"].mean(); BW.append(bw)
    for k in BR:
        bk = R.theil_sen(x4, fit[k] + res[k][idx])[0] / D[k].mean(); BR[k].append(bk / bw)
for k in BR:
    v = np.array(BR[k]); lo, med, hi = np.percentile(v, [2.5, 50, 97.5])
    rows.append(dict(test="D1 ratio", series=f"{k}/wind", ratio=TR[k] / TR["wind"], ratio_med=med, ratio_lo=lo, ratio_hi=hi,
                     frac_wind_nonpos=float(np.mean(np.array(BW) <= 0))))
    print(rows[-1], flush=True)

# D2
yrs = np.arange(2001, 2026); yr4 = np.arange(2004, 2026)
def ann(v, yy): return np.array([v[yy == q].mean() for q in np.unique(yy)])
def detr(y):
    tt = np.arange(y.size); b = np.polyfit(tt, y, 1); return y - np.polyval(b, tt)
wa = detr(ann(S["wind"], YY)); XA = {"KE": detr(ann(S["KE"][c4], YY[c4])), "EKE": detr(ann(S["EKE"][c4], YY[c4])), "slope": detr(ann(s_sl, YYs))}
LAGS = (0, 1, 2, 3)
def lagged(wfull, L): return wfull[(yrs >= 2004 - L) & (yrs <= 2025 - L)]
def surr(y, rng):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size)
NS = 2000; SW = np.array([surr(wa, rng) for _ in range(NS)])
for k, xa in XA.items():
    robs = np.array([np.corrcoef(xa, lagged(wa, L))[0, 1] for L in LAGS])
    null = np.array([[np.corrcoef(xa, lagged(s_, L))[0, 1] for L in LAGS] for s_ in SW])
    env = np.quantile(np.abs(null), 0.95, axis=0); mx = np.abs(null).max(1)
    pfw = (1 + np.sum(mx >= np.abs(robs).max())) / (1 + NS)
    for q, L in enumerate(LAGS):
        p = (1 + np.sum(np.abs(null[:, q]) >= abs(robs[q]))) / (1 + NS)
        rows.append(dict(test="D2", series=k, lag=L, r=robs[q], p=p, env95=env[q], p_fw=pfw))
    print(k, np.round(robs, 2), "env", np.round(env, 2), "family-wise p", round(pfw, 3), flush=True)
pd.DataFrame(rows).to_csv(OUT / "f4e_results.csv", index=False)
print("done", flush=True)
