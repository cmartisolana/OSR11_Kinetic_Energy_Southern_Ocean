"""F3w: eastward wind stress on the Figure 3 domains (28 Oct 2026). Specification fixed before any result.
Product 2 monthly eastward stress, 0.25 deg, 2001-2025. Masks from the altimetric 2001-2025 mean (1/8 deg, block-averaged 2x2 to
0.25 deg): ACC band = mean ADT -0.6 to +0.2 m (outputs/mean_adt.nc, as Fig. 3c); band between the mean SAF and PF = clim_allsat
ADT between the dk12 levels (as Fig. 3d); 0.25 deg cells containing any ice-masked 1/8 deg pixel are dropped.
 S1 ACC-band monthly series (cos-lat weighted); trend with rstats.trend_bcz (monthly test, both AR orders) and annual test.
 S2 2 deg longitude bins between the SAF and PF (centres as dk12): annual trend 2001-2025 (trend_bcz annual, 150 nulls), % per decade.
 S3 cross-bin correlation of the EKE trend (dk12, % per decade, 2001-2025) with the S2 wind trend; p from all circular shifts.
 S4 per-cell wind trend map on the 1 deg grid south of 35 S (cells with no ice-masked pixel), calibrated as a1_maps_bc / f4_tilt,
    BH-FDR 10 % (run with argument "map"; slow).
Usage: python f3w_wind.py <revision dir> <out dir> [map]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy.stats import norm

REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True); DO_MAP = len(sys.argv) > 3 and sys.argv[3] == "map"
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402

W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress.sel(latitude=slice(-65, -35))
la, lo = W.latitude.values, W.longitude.values
t = pd.DatetimeIndex(W.time.values); yy, mm = t.year.values, t.month.values; x = yy + (mm - 0.5) / 12
A = W.values                                                     # (300, 120, 1440)
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc").adt
ice8 = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask
def to025(da, how="mean"):
    v = da.values; ny, nx = v.shape; v = v[: ny // 2 * 2, : nx // 2 * 2].reshape(ny // 2, 2, nx // 2, 2)
    return np.nanmean(v, (1, 3)) if how == "mean" else np.nanmax(v, (1, 3))
lat25 = madt.latitude.values.reshape(-1, 2).mean(1); lon25 = madt.longitude.values.reshape(-1, 2).mean(1)
assert np.allclose(lat25, la, atol=0.01) and np.allclose(lon25, lo, atol=0.01), (lat25[:3], la[:3], lon25[:3], lo[:3])
adt25, cl25 = to025(madt), to025(clim); ice25 = to025(ice8.astype(float), "max") > 0.5
Z = np.load(REV / "plan1019" / "drake" / "out" / "dk12" / "dk12_band_series_allsat.npz"); LS, LP = float(Z["lev_saf"]), float(Z["lev_pf"])
ok = np.isfinite(A).all(0) & ~ice25
acc = ok & (adt25 >= -0.6) & (adt25 <= 0.2); frt = ok & (cl25 >= LP) & (cl25 <= LS)
LAT = np.broadcast_to(la[:, None], acc.shape); LON = np.broadcast_to(lo[None, :], acc.shape); w = np.cos(np.radians(LAT))
def rm(m): return np.array([np.sum(A[k][m] * w[m]) / np.sum(w[m]) for k in range(A.shape[0])])

rows = []
s = rm(acc); pd.DataFrame({"year": yy, "month": mm, "x": x, "tau_acc": s}).to_csv(OUT / "f3w_acc_monthly.csv", index=False)
r = R.trend_bcz(x, s, month=mm, n_outer=300); mu = s.mean()
rows.append(dict(series="tau ACC band, monthly test", slope_dec=10 * r["slope"], lo_dec=10 * r["ci_lo"], hi_dec=10 * r["ci_hi"],
                 pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"], mean=mu))
ya = np.array([s[yy == k].mean() for k in range(2001, 2026)]); r = R.trend_bcz(np.arange(2001, 2026) + 0.5, ya, monthly=False, n_outer=300)
rows.append(dict(series="tau ACC band, annual test", pct=1000 * r["slope"] / ya.mean(), pct_lo=1000 * r["ci_lo"] / ya.mean(),
                 pct_hi=1000 * r["ci_hi"] / ya.mean(), p=r["p"], mean=ya.mean()))
for q in rows: print(q, flush=True)

T2 = pd.read_csv(REV / "plan1019" / "drake" / "out" / "dk12" / "dk12_band_trends_allsat.csv").sort_values("lon")
brow = []
for c in T2.lon.values:
    m = frt & (LON >= c - 1) & (LON < c + 1)
    if m.sum() < 2: brow.append(dict(lon=c, n=int(m.sum()))); continue
    sa = rm(m); ya = np.array([sa[yy == k].mean() for k in range(2001, 2026)])
    r = R.trend_bcz(np.arange(2001, 2026) + 0.5, ya, monthly=False, n_outer=150)
    brow.append(dict(lon=c, n=int(m.sum()), tau_mean=ya.mean(), tau_pct=1000 * r["slope"] / ya.mean(), tau_lo=1000 * r["ci_lo"] / ya.mean(),
                     tau_hi=1000 * r["ci_hi"] / ya.mean(), tau_p=r["p"]))
B = pd.DataFrame(brow).merge(T2[["lon", "pct", "pct_lo", "pct_hi", "p"]].rename(columns={"pct": "eke_pct", "pct_lo": "eke_lo", "pct_hi": "eke_hi", "p": "eke_p"}), on="lon")
B.to_csv(OUT / "f3w_bins.csv", index=False)
v = B.dropna(subset=["tau_pct", "eke_pct"]); a, b = v.eke_pct.values, v.tau_pct.values
rr = np.corrcoef(a, b)[0, 1]; null = np.array([np.corrcoef(a, np.roll(b, k))[0, 1] for k in range(1, len(b))])
pp = (1 + np.sum(np.abs(null) >= abs(rr))) / (1 + null.size)
rows.append(dict(series="cross-bin r(EKE trend, wind trend), 2 deg bins between SAF and PF", r=rr, p=pp, n=len(v)))
print(f"S3 r = {rr:+.2f}, p = {pp:.3f}, n = {len(v)}; wind bins positive {int((v.tau_pct > 0).sum())}, p<0.05 {int((v.tau_p < 0.05).sum())}", flush=True)
pd.DataFrame(rows).to_csv(OUT / "f3w_results.csv", index=False)

if DO_MAP:
    lat1 = np.arange(-64.5, -35, 1.0); lon1 = np.arange(-179.5, 180, 1.0)
    il = np.digitize(la, np.r_[lat1 - 0.5, lat1[-1] + 0.5]) - 1; jl = np.digitize(lo, np.r_[lon1 - 0.5, lon1[-1] + 0.5]) - 1
    Ay = np.array([A[yy == k].mean(0) for k in range(2001, 2026)])
    C1 = np.full((25, lat1.size, lon1.size), np.nan); bad = np.zeros((lat1.size, lon1.size), bool)
    for i in range(lat1.size):
        ri = il == i
        for j in range(lon1.size):
            cj = jl == j; blk = Ay[:, ri][:, :, cj]
            bad[i, j] = ice25[np.ix_(ri, cj)].any() or not np.isfinite(blk).all()
            C1[:, i, j] = np.nanmean(blk, axis=(1, 2))
    m = ~bad & np.isfinite(C1).all(0); idx = np.argwhere(m); xs = np.arange(2001, 2026) + 0.5
    rng = np.random.default_rng(0); bb = []; zz_ = []; ph = []; sg = []; mu_ = []; mods = []
    for (i, j) in idx:
        b_, a_, sd, phis, phi_bc, sig = R._z_stat(xs, C1[:, i, j], False, 200, rng)
        bb.append(b_); zz_.append(b_ / sd); ph.append(float(phi_bc.sum())); sg.append(sig); mu_.append(C1[:, i, j].mean()); mods.append(phi_bc)
    bb, zz_, ph, sg, mu_ = map(np.array, (bb, zz_, ph, sg, mu_))
    e = np.quantile(ph, np.linspace(0, 1, 11)); e[-1] += 1e-9; cls = np.clip(np.digitize(ph, e) - 1, 0, 9); kappa = np.zeros(10)
    rng = np.random.default_rng(1)
    for k in range(10):
        pool = np.nonzero(cls == k)[0]; pick = rng.choice(pool, min(25, pool.size), replace=False); nz = []
        for q in pick:
            for s_ in R._simulate_ar(mods[q], sg[q], xs.size, 40, rng, burn=500):
                bi, _, sdi, *_ = R._z_stat(xs, s_, False, 200, rng); nz.append(abs(bi / sdi))
        nz = np.array(nz); kappa[k] = max(np.quantile(nz, 0.95) / 1.96, np.quantile(nz, 0.99) / 2.576)
    p = 2 * norm.sf(np.abs(zz_) / kappa[cls]); fdr = R.bh_fdr(p, 0.10).astype(bool)
    PCT = np.full(m.shape, np.nan); P = np.full(m.shape, np.nan); F = np.zeros(m.shape, np.int8)
    MU = np.full(m.shape, np.nan); TD = np.full(m.shape, np.nan)
    for q, (i, j) in enumerate(idx):
        PCT[i, j] = 1000 * bb[q] / mu_[q]; P[i, j] = p[q]; F[i, j] = fdr[q]; MU[i, j] = mu_[q]; TD[i, j] = 10 * bb[q]
    # wind_trend_dec (N m-2 per decade) is drawn in Fig. 3b: the percentage is unstable where the mean stress is near zero
    xr.Dataset({"wind_pct": (("lat", "lon"), PCT), "wind_p": (("lat", "lon"), P), "wind_fdr": (("lat", "lon"), F),
                "wind_mean": (("lat", "lon"), MU), "wind_trend_dec": (("lat", "lon"), TD)},
               coords={"lat": lat1, "lon": lon1}).to_netcdf(OUT / "f3w_wind_map.nc")
    print(f"S4 map cells {len(idx)}, positive {int((bb > 0).sum())}, p<0.05 {int((p < 0.05).sum())}, FDR {int(fdr.sum())}", flush=True)
print("done", flush=True)
