"""F4v: along-path (5 deg longitude) relations between eastward wind stress, EKE and the sigma0 = 27.3 cross-stream slope, using
values as well as trends (28 Oct 2026, for Figure 4d). Domain: the f4m cells (ACC band + 100 km, ice removed, slope defined in all
264 months) that also have EKE in every month; 72 bins of 5 deg; cos-lat weighted bin means. Common period 2004-2025 (monthly).
Products: wind Product 2 eastward stress (wind_taux.nc), EKE Product 1 (0.5|U_m - U_bar|^2 + 0.5 var_intra, U_bar 2001-2025, ice
pixels dropped; as f4d), slope f4m_p273_monthly.npz.
For each pair (wind-EKE, wind-slope, EKE-slope), Pearson r across bins and p from all circular shifts of the along-path profile:
 V1 time-mean values (2004-2025)
 V2 trends (Theil-Sen of the deseasonalised monthly bin series, % of the bin mean per decade; descriptive slope, same period)
 V3 interannual values: annual bin anomalies (bin mean removed, linearly detrended per bin), pooled over bins and years;
    p from 2000 random circular shifts in time of the second variable (same shift for all bins), which keeps the spatial and
    temporal autocorrelation.
Usage: python f4v_values.py <revision dir> <f4 out dir> <f4m out dir> <out dir>
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
# slope, monthly 2004-2025
C = np.load(FM / "f4m_p273_monthly.npz"); YY, MM, PM, lonr = C["year"], C["month"], C["p"], C["lon"]
L0 = np.where(lonr > 180, lonr - 360, lonr); o = np.argsort(L0); adt = M["adt"].values
def grad(F):
    dy = np.gradient(F, lat, axis=-2) / KM * 100
    dx = (np.roll(F, -1, axis=-1) - np.roll(F, 1, axis=-1)) / 2.0 / (KM * np.cos(np.radians(lat))[:, None]) * 100
    return dx, dy
ax_, ay_ = grad(adt); an = np.hypot(ax_, ay_)
dx, dy = grad(PM[:, :, o].astype(float)); SL = dx * ax_ / an + dy * ay_ / an
dom = (M["buffer100"].values > 0) & (M["ice"].values == 0) & np.isfinite(M["signed_27p3_pct"].values) & np.isfinite(SL).all(0)
# wind, monthly 2004-2025 on cells
W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress.sel(latitude=slice(lat[0] - 0.5, lat[-1] + 0.5), time=slice("2004-01-01", "2025-12-31"))
il = np.digitize(W.latitude.values, la_e) - 1; jl = np.digitize(W.longitude.values, lo_e) - 1; A = W.values
WC = np.full((A.shape[0], lat.size, lon.size), np.nan)
for i in range(lat.size):
    sub = A[:, il == i, :]
    if sub.shape[1] == 0: continue
    for j in range(lon.size): WC[:, i, j] = np.nanmean(sub[:, :, jl == j], axis=(1, 2))
assert WC.shape[0] == 264
# EKE, monthly on cells (U_bar over 2001-2025 as in the main text), keep 2004-2025
files = sorted(glob.glob(str(REV / "data" / "monthly_allsat" / "*" / "*.nc"))); assert len(files) == 300
with xr.open_dataset(files[0]) as d0: la8, lo8 = d0.latitude.values, d0.longitude.values
ice8 = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
su = np.zeros((la8.size, lo8.size)); sv = np.zeros_like(su)
for f in files:
    with xr.open_dataset(f) as d: su += d.u_mean.values; sv += d.v_mean.values
ub, vb = su / 300, sv / 300; i8 = np.digitize(la8, la_e) - 1; j8 = np.digitize(lo8, lo_e) - 1
EC = np.full((264, lat.size, lon.size), np.nan); k = 0
for f in files[36:]:
    with xr.open_dataset(f) as d:
        e8 = (0.5 * ((d.u_mean.values - ub) ** 2 + (d.v_mean.values - vb) ** 2) + 0.5 * d.var_intra.values) * 1e4
    e8[ice8] = np.nan
    for i in range(lat.size):
        ri = i8 == i
        if not ri.any(): continue
        sub = e8[ri, :]
        EC[k, i, :] = [np.nanmean(sub[:, j8 == j]) if np.isfinite(sub[:, j8 == j]).any() else np.nan for j in range(lon.size)]
    k += 1
print("EKE done", flush=True)
dom &= np.isfinite(EC).all(0) & np.isfinite(WC).all(0)
print("common cells", int(dom.sum()), flush=True)
x = YY + (MM - 0.5) / 12; edges = np.arange(-180, 181, 5); cen = edges[:-1] + 2.5
V = {"wind": WC, "eke": EC, "slope": SL}; BS = {k_: [] for k_ in V}; keep = []
for c in cen:
    m = dom & (LON >= c - 2.5) & (LON < c + 2.5)
    if m.sum() < 2: continue
    keep.append(c)
    for k_, F in V.items(): BS[k_].append(np.array([np.sum(F[t][m] * w[m]) / np.sum(w[m]) for t in range(264)]))
cen = np.array(keep); BS = {k_: np.array(v) for k_, v in BS.items()}            # (nbin, 264)
np.savez(OUT / "f4v_bin_series.npz", lon=cen, x=x, year=YY, month=MM, **BS)
mean = {k_: v.mean(1) for k_, v in BS.items()}
trend = {k_: np.array([1000 * R.theil_sen(x, R.deseasonalise(s, MM))[0] / s.mean() for s in v]) for k_, v in BS.items()}
yrs = np.arange(2004, 2026)
def ann_anom(v):
    a = np.array([[s[YY == y].mean() for y in yrs] for s in v]); a = a - a.mean(1, keepdims=True)
    t = yrs - yrs.mean(); b = (a * t).sum(1) / (t ** 2).sum(); return a - b[:, None] * t
anom = {k_: ann_anom(v) for k_, v in BS.items()}
def circ(a, b):
    r = np.corrcoef(a, b)[0, 1]; null = np.array([np.corrcoef(a, np.roll(b, s))[0, 1] for s in range(1, len(b))])
    return r, (1 + np.sum(np.abs(null) >= abs(r))) / (1 + null.size)
rng = np.random.default_rng(20261028); rows = []
for a, b in (("wind", "eke"), ("wind", "slope"), ("eke", "slope")):
    r1, p1 = circ(mean[a], mean[b]); r2, p2 = circ(trend[a], trend[b])
    A_, B_ = anom[a], anom[b]; r3 = np.corrcoef(A_.ravel(), B_.ravel())[0, 1]
    null = np.array([np.corrcoef(A_.ravel(), np.roll(B_, rng.integers(1, len(yrs)), axis=1).ravel())[0, 1] for _ in range(2000)])
    p3 = (1 + np.sum(np.abs(null) >= abs(r3))) / 2001
    rows += [dict(x=a, y=b, kind="mean values", r=r1, p=p1, n=len(cen)), dict(x=a, y=b, kind="trends", r=r2, p=p2, n=len(cen)),
             dict(x=a, y=b, kind="interannual anomalies (pooled)", r=r3, p=p3, n=A_.size)]
    print(f"{a}-{b}: mean r={r1:+.2f} p={p1:.3f} | trend r={r2:+.2f} p={p2:.3f} | anomalies r={r3:+.2f} p={p3:.3f}", flush=True)
pd.DataFrame(rows).to_csv(OUT / "f4v_results.csv", index=False)
pd.DataFrame({"lon": cen, **{f"{k_}_mean": mean[k_] for k_ in V}, **{f"{k_}_trend": trend[k_] for k_ in V}}).to_csv(OUT / "f4v_bins.csv", index=False)
print("done", flush=True)
