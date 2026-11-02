"""F4d: inputs for the revised Figure 4 (28 Oct 2026): wind-stress trend map, monthly wind series, EKE on the slope cells and
the slope-EKE comparison. Specification fixed before any result.

Domain: the sigma0 = 27.3 slope cells of f4_tilt.py (ACC band + 100 km, ice removed, slope defined in all years).
EKE (Product 1, as in the main text): monthly EKE = 0.5 |U_m - U_bar|^2 + 0.5 var_intra (U_bar = 2001-2025 mean), cm2 s-2, on
the 1/8 deg pixels; pixels with ice in any month dropped; averaged onto the 1 deg RG cells (8 x 8 pixels); annual means.
Wind (Product 2 eastward stress): as f4b/f4c (1 deg cells, annual means); monthly domain-mean series for panel c.
 W-map: per-cell trend of annual stress 2001-2025 (% per decade), z statistic of rstats.trend_bcz calibrated per persistence
        class as in a1_maps_bc / f4_tilt (10 classes, 25 cells x 40 nulls), BH-FDR 10 %.
 E1 trends 2004-2025 (Argo period) of domain-mean EKE (rstats.trend_bcz, annual) and 2001-2025 for reference.
 E2 interannual: detrended annual domain means, r(slope_y, EKE_{y+L}) for L = -2..+2 (positive L: EKE lagging the slope);
    phase-randomised surrogates of the EKE series (2000), per lag and family-wise over the five lags.
 E3 spatial: 5 deg bins, EKE trend 2004-2025 vs slope trend; cross-bin r with circular-shift p (as f4b T3).
Usage: python f4d_eke_wind.py <revision dir> <f4 out dir> <out dir>
"""
import sys, glob
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy.stats import norm

REV, F4, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402

M = xr.open_dataset(F4 / "f4_tilt_maps.nc"); lat, lon = M.lat.values, M.lon.values
dom = (M["buffer100"].values > 0) & (M["ice"].values == 0) & np.isfinite(M["signed_27p3_pct"].values)
LON, LAT = np.meshgrid(lon, lat); w = np.cos(np.radians(LAT))
SER = np.load(F4 / "f4_regional_series.npz"); yS = SER["years"].astype(int); s_dom = SER["27.3__signed__buffer100__ACC"]

def regmean(F, m):
    return np.array([np.nansum(F[t][m] * w[m]) / np.sum(w[m] * np.isfinite(F[t][m])) for t in range(F.shape[0])])

# ---------------- wind (annual cells + monthly domain mean)
W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress.sel(latitude=slice(lat[0] - 0.5, lat[-1] + 0.5))
la_e = np.r_[lat - 0.5, lat[-1] + 0.5]; lo_e = np.r_[lon - 0.5, lon[-1] + 0.5]
il = np.digitize(W.latitude.values, la_e) - 1; jl = np.digitize(W.longitude.values, lo_e) - 1
def to_cells(A):
    T = np.full((A.shape[0], lat.size, lon.size), np.nan)
    for i in range(lat.size):
        sub = A[:, il == i, :]
        if sub.shape[1] == 0: continue
        for j in range(lon.size):
            T[:, i, j] = np.nanmean(sub[:, :, jl == j], axis=(1, 2))
    return T
Wm = to_cells(W.values)                                            # monthly, 300
tW = W.time.values; yW = pd.DatetimeIndex(tW).year.values; mW = pd.DatetimeIndex(tW).month.values
Wy_years = np.unique(yW); Wy = np.array([np.nanmean(Wm[yW == y], 0) for y in Wy_years])
tau_month = regmean(Wm, dom)
pd.DataFrame({"year": yW, "month": mW, "x": yW + (mW - 0.5) / 12, "tau": tau_month}).to_csv(OUT / "f4d_wind_monthly.csv", index=False)
r = R.trend_bcz(yW + (mW - 0.5) / 12, tau_month, month=mW, n_outer=300); mu = tau_month.mean()
rows = [dict(test="W", series="tau domain, monthly test, 2001-2025", slope_dec=10 * r["slope"], lo_dec=10 * r["ci_lo"], hi_dec=10 * r["ci_hi"],
             pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"])]
print("tau domain monthly test:", rows[-1], flush=True)

# ---------------- per-cell calibrated trend map (annual)
def cellmap(F, m, x):
    idx = np.argwhere(m); n = len(idx)
    b = np.full(n, np.nan); z = np.full(n, np.nan); ph = np.full(n, np.nan); sg = np.full(n, np.nan); mu = np.full(n, np.nan); models = []
    rng = np.random.default_rng(0)
    for q, (i, j) in enumerate(idx):
        y = F[:, i, j]
        bb, a, sd, phis, phi_bc, sig = R._z_stat(x, y, False, 200, rng)
        b[q], z[q], ph[q], sg[q], mu[q] = bb, bb / sd, float(phi_bc.sum()), sig, y.mean(); models.append(phi_bc)
    e = np.quantile(ph, np.linspace(0, 1, 11)); e[-1] += 1e-9
    cls = np.clip(np.digitize(ph, e) - 1, 0, 9); kappa = np.zeros(10); rng = np.random.default_rng(1)
    for k in range(10):
        pool = np.nonzero(cls == k)[0]; pick = rng.choice(pool, min(25, pool.size), replace=False); zz = []
        for q in pick:
            for s_ in R._simulate_ar(models[q], sg[q], x.size, 40, rng, burn=500):
                bi, _, sdi, *_ = R._z_stat(x, s_, False, 200, rng); zz.append(abs(bi / sdi))
        zz = np.array(zz); kappa[k] = max(np.quantile(zz, 0.95) / 1.96, np.quantile(zz, 0.99) / 2.576)
    p = 2 * norm.sf(np.abs(z) / kappa[cls]); fdr = R.bh_fdr(p, 0.10).astype(bool)
    out = {k_: np.full(m.shape, np.nan) for k_ in ("pct", "p")}; out["fdr"] = np.zeros(m.shape, bool)
    for q, (i, j) in enumerate(idx):
        out["pct"][i, j] = 1000 * b[q] / mu[q]; out["p"][i, j] = p[q]; out["fdr"][i, j] = fdr[q]
    print(f"  map cells {n}, positive {int((b > 0).sum())}, p<0.05 {int((p < 0.05).sum())}, FDR {int(fdr.sum())}", flush=True)
    return out
print("wind map", flush=True)
WM = cellmap(Wy, dom & np.isfinite(Wy).all(0), Wy_years + 0.5)

# ---------------- EKE on the cells
files = sorted(glob.glob(str(REV / "data" / "monthly_allsat" / "*" / "*.nc")))
assert len(files) == 300, len(files)
with xr.open_dataset(files[0]) as d0:
    la8, lo8 = d0.latitude.values, d0.longitude.values
ice8 = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
su = np.zeros((la8.size, lo8.size)); sv = np.zeros_like(su)
for f in files:                                   # pass 1: 2001-2025 mean velocity
    with xr.open_dataset(f) as d:
        su += d.u_mean.values; sv += d.v_mean.values
ub, vb = su / 300, sv / 300
i8 = np.digitize(la8, la_e) - 1; j8 = np.digitize(lo8, lo_e) - 1
E1m = np.full((300, lat.size, lon.size), np.nan)
for k, f in enumerate(files):                     # pass 2: monthly EKE, ice pixels dropped, onto 1 deg cells
    with xr.open_dataset(f) as d:
        e8 = (0.5 * ((d.u_mean.values - ub) ** 2 + (d.v_mean.values - vb) ** 2) + 0.5 * d.var_intra.values) * 1e4
    e8[ice8] = np.nan
    for i in range(lat.size):
        ri = i8 == i
        if not ri.any(): continue
        sub = e8[ri, :]
        E1m[k, i, :] = [np.nanmean(sub[:, j8 == j]) if np.isfinite(sub[:, j8 == j]).any() else np.nan for j in range(lon.size)]
yE = np.repeat(np.arange(2001, 2026), 12)
Ey = np.array([np.nanmean(E1m[yE == y], 0) for y in range(2001, 2026)]); Ey_years = np.arange(2001, 2026)
domE = dom & np.isfinite(Ey).all(0)
print("EKE cells", int(domE.sum()), "of", int(dom.sum()), flush=True)
e_dom = regmean(Ey, domE); e_month = regmean(E1m, domE)
pd.DataFrame({"x": yE + (np.tile(np.arange(12), 25) + 0.5) / 12, "eke": e_month}).to_csv(OUT / "f4d_eke_monthly.csv", index=False)

def trA(name, x, y):
    r = R.trend_bcz(x + 0.5, y, monthly=False, n_outer=300); mu = y.mean()
    rows.append(dict(test="E1", series=name, slope_dec=10 * r["slope"], lo_dec=10 * r["ci_lo"], hi_dec=10 * r["ci_hi"],
                     pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"]))
    print(name, f"{rows[-1]['pct']:+.2f} [{rows[-1]['pct_lo']:+.2f}, {rows[-1]['pct_hi']:+.2f}] p={r['p']:.3f}", flush=True)
s4 = Ey_years >= 2004
trA("EKE domain 2004-2025", Ey_years[s4], e_dom[s4]); trA("EKE domain 2001-2025", Ey_years, e_dom)

# E2 interannual
def detr(y):
    x = np.arange(y.size); return y - np.polyval(np.polyfit(x, y, 1), x)
def phase_surr(y, rng):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size)
rng = np.random.default_rng(20261028); LAGS = (-2, -1, 0, 1, 2)
def pairs(sv, ev, L):
    # slope years 2004-2025 (index 0..21); EKE years 2001-2025 (index 0..24); EKE_{y+L}
    ys = np.arange(2004, 2026); ye = ys + L; ok = (ye >= 2001) & (ye <= 2025)
    return sv[ok], ev[ye[ok] - 2001]
r_obs = {}; null = np.zeros((2000, len(LAGS)))
for q, L in enumerate(LAGS):
    a, b = pairs(s_dom, e_dom, L); r_obs[L] = np.corrcoef(detr(a), detr(b))[0, 1]
ed = detr(e_dom)
for k in range(2000):
    es = phase_surr(ed, rng)
    for q, L in enumerate(LAGS):
        a, b = pairs(s_dom, es, L); null[k, q] = np.corrcoef(detr(a), detr(b))[0, 1]
for q, L in enumerate(LAGS):
    p = (1 + np.sum(np.abs(null[:, q]) >= abs(r_obs[L]))) / 2001
    rows.append(dict(test="E2", series=f"slope vs EKE, lag {L:+d} yr (EKE lagging if +)", r=r_obs[L], p=p))
    print(f"E2 lag {L:+d}: r = {r_obs[L]:+.2f}, p = {p:.3f}", flush=True)
robs = max(abs(v) for v in r_obs.values()); pmax = (1 + np.sum(np.abs(null).max(1) >= robs)) / 2001
rows.append(dict(test="E2", series="slope vs EKE, max over lags -2..+2", r=robs, p=pmax))
print(f"E2 max |r| = {robs:.2f}, family-wise p = {pmax:.3f}", flush=True)

# E3 spatial
L5 = pd.read_csv(F4 / "f4_lon5_trends.csv"); L5 = L5[L5.surface.astype(str) == "27.3"].set_index("lon")
edges = np.arange(-180, 181, 5); brow = []
for k in range(edges.size - 1):
    c = edges[k] + 2.5; m = domE & (LON >= edges[k]) & (LON < edges[k + 1])
    if m.sum() < 2 or c not in L5.index: continue
    e = regmean(Ey, m)[s4]; r = R.trend_bcz(Ey_years[s4] + 0.5, e, monthly=False, n_outer=150)
    brow.append(dict(lon=c, eke_pct=1000 * r["slope"] / e.mean(), eke_lo=1000 * r["ci_lo"] / e.mean(), eke_hi=1000 * r["ci_hi"] / e.mean(),
                     eke_p=r["p"], eke_mean=e.mean(), slope_pct=L5.loc[c, "pct"]))
B = pd.DataFrame(brow); B.to_csv(OUT / "f4d_eke_bins.csv", index=False)
x, y = B.slope_pct.values, B.eke_pct.values; r_sp = np.corrcoef(x, y)[0, 1]
null = np.array([np.corrcoef(x, np.roll(y, k))[0, 1] for k in range(1, len(y))])
p_sp = (1 + np.sum(np.abs(null) >= abs(r_sp))) / (1 + null.size)
rows.append(dict(test="E3", series="cross-bin r(slope trend, EKE trend), 2004-2025", r=r_sp, p=p_sp, n=len(B)))
print(f"E3 cross-bin r = {r_sp:+.2f}, circular-shift p = {p_sp:.3f} (all {len(y) - 1} shifts)", flush=True)
pd.DataFrame(rows).to_csv(OUT / "f4d_results.csv", index=False)
xr.Dataset({"wind_pct": (("lat", "lon"), WM["pct"]), "wind_p": (("lat", "lon"), WM["p"]), "wind_fdr": (("lat", "lon"), WM["fdr"].astype(np.int8)),
            "eke_mean": (("lat", "lon"), np.where(domE, np.nanmean(Ey[s4], 0), np.nan))},
           coords={"lat": lat, "lon": lon}).to_netcdf(OUT / "f4d_maps.nc")
np.savez(OUT / "f4d_series.npz", years=Ey_years, eke_dom=e_dom, tau_dom=regmean(Wy, dom), slope_years=yS, slope=s_dom)
print("done", flush=True)
