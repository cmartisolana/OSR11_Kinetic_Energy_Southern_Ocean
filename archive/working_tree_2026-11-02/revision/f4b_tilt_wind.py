"""F4b: isopycnal slope (sigma0 = 27.3) against eastward wind stress (28 Oct 2026). Specification fixed before any result.

Inputs: f4_tilt.py outputs (f4_tilt_maps.nc: domain = ACC band + 100 km, ice removed; f4_regional_series.npz) and the
annual sigma0 = 27.3 cross-stream slope recomputed from dk32_iso.npz exactly as in f4_tilt.py; Product 2 monthly eastward
stress (data/wind_taux.nc, 0.25 deg, 2001-2025), block-averaged onto the 1 deg RG cells, annual means.
Wind regions: (W1) the same cells as the slope (domain, cos-lat weights); (W2) 45-60 S circumpolar zonal band (as in SI).
Tests:
 T1 trends 2004-2025 of the domain-mean slope and of W1, W2 (rstats.trend_bcz, annual); descriptive ratio of relative trends.
 T2 interannual: linearly detrended annual series, r(slope_y, tau_{y-L}) for L = 0, 1, 2, 3 yr (wind leading; wind from
    2001 so the slope keeps 2004-2025); significance from 2000 phase-randomised surrogates of the wind series, both for each
    lag and for the maximum |r| over the four lags (family-wise p). W1 and W2.
 T3 spatial (5 deg longitude bins of the domain): cross-bin correlation of the slope trend with the W1 trend, significance from
    2000 circular shifts of the wind-trend profile; and per-bin detrended r at L = 0 and 1 (counts with p < 0.05, BH-FDR 10 %).
Revision (28 Oct 2026, author request): wind trends use the full wind record (2001-2025 by default, WIND_Y0 = 4th argument);
the slope trends stay 2004-2025 (Argo). T2 is unchanged (it already uses wind from 2001 for the lags).
Usage: python f4b_tilt_wind.py <revision dir> <f4 out dir> <out dir> [WIND_Y0, default 2001]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

REV, F4, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
WIND_Y0 = int(sys.argv[4]) if len(sys.argv) > 4 else 2001
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402

M = xr.open_dataset(F4 / "f4_tilt_maps.nc")
lat, lon = M.lat.values, M.lon.values
I = np.load(REV / "plan1019" / "drake" / "out" / "dk32" / "dk32_iso.npz")
yrs = I["years"].astype(int); lon0 = I["lon"]; L0 = np.where(lon0 > 180, lon0 - 360, lon0); o = np.argsort(L0)
assert np.allclose(L0[o], lon)
adt = M["adt"].values; KM = 111.195
def grad(F):
    dy = np.gradient(F, lat, axis=-2) / KM * 100
    dx = (np.roll(F, -1, axis=-1) - np.roll(F, 1, axis=-1)) / 2.0 / (KM * np.cos(np.radians(lat))[:, None]) * 100
    return dx, dy
ax_, ay_ = grad(adt); an = np.hypot(ax_, ay_)
P = I["p27.3"][:, :, o]; dx, dy = grad(P); SL = dx * ax_ / an + dy * ay_ / an          # (22, lat, lon)
dom = (M["buffer100"].values > 0) & (M["ice"].values == 0) & np.isfinite(SL).all(0)
LON, LAT = np.meshgrid(lon, lat); w = np.cos(np.radians(LAT))

# wind on the 1 deg cells, annual 2001-2025
W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress
W = W.sel(latitude=slice(lat[0] - 0.5, lat[-1] + 0.5))
Wy = W.groupby("time.year").mean("time")                                                  # (25, 0.25 deg)
ny_w = Wy.year.values
la_e = np.r_[lat - 0.5, lat[-1] + 0.5]; lo_e = np.r_[lon - 0.5, lon[-1] + 0.5]
il = np.digitize(Wy.latitude.values, la_e) - 1; jl = np.digitize(Wy.longitude.values, lo_e) - 1
T1 = np.full((ny_w.size, lat.size, lon.size), np.nan)
A = Wy.values
for i in range(lat.size):
    ri = il == i
    if not ri.any(): continue
    sub = A[:, ri, :]
    for j in range(lon.size):
        cj = jl == j
        T1[:, i, j] = np.nanmean(sub[:, :, cj], axis=(1, 2))
print("wind years", ny_w[0], ny_w[-1], "NaN fraction in domain", float(np.isnan(T1[:, dom]).mean()), flush=True)

def regmean(F, m):
    return np.array([np.nansum(F[t][m] * w[m]) / np.sum(w[m] * np.isfinite(F[t][m])) for t in range(F.shape[0])])

s_dom = regmean(SL, dom)                                    # 2004-2025
tau_dom = regmean(T1, dom)                                  # 2001-2025
z45 = (LAT >= -60) & (LAT <= -45) & np.isfinite(T1).all(0)
tau_45 = regmean(T1, z45)

rows = []
def tr(name, x, y):
    r = R.trend_bcz(x + 0.5, y, monthly=False, n_outer=300); mu = y.mean()
    rows.append(dict(test="T1", series=name, pct=1000 * r["slope"] / mu, lo=1000 * r["ci_lo"] / mu, hi=1000 * r["ci_hi"] / mu, p=r["p"]))
    print(f"T1 {name}: {rows[-1]['pct']:+.2f} %/dec [{rows[-1]['lo']:+.2f}, {rows[-1]['hi']:+.2f}] p={r['p']:.3f}", flush=True)
sel = (ny_w >= WIND_Y0); s4 = ny_w >= 2004
tr("slope 27.3, domain, 2004-2025", yrs, s_dom)
tr(f"tau, domain (W1) {WIND_Y0}-2025", ny_w[sel], tau_dom[sel])
tr(f"tau, 45-60S (W2) {WIND_Y0}-2025", ny_w[sel], tau_45[sel])
tr("tau, domain (W1) 2004-2025", ny_w[s4], tau_dom[s4])
tr("tau, 45-60S (W2) 2004-2025", ny_w[s4], tau_45[s4])

def detr(y):
    x = np.arange(y.size); b = np.polyfit(x, y, 1); return y - np.polyval(b, x)
def phase_surr(y, rng):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size)
rng = np.random.default_rng(20261028); NS = 2000; LAGS = (0, 1, 2, 3)
def lagset(tau_full):
    # tau_{y-L} aligned with slope years 2004-2025
    return {L: tau_full[(ny_w >= 2004 - L) & (ny_w <= 2025 - L)] for L in LAGS}
for wname, tau_full in (("W1 domain", tau_dom), ("W2 45-60S", tau_45)):
    sd = detr(s_dom); TL = {L: detr(v) for L, v in lagset(tau_full).items()}
    r_obs = {L: np.corrcoef(sd, TL[L])[0, 1] for L in LAGS}
    # surrogates of the full detrended wind series, then the same lag alignment
    td = detr(tau_full); null = np.zeros((NS, len(LAGS)))
    for k in range(NS):
        ts = phase_surr(td, rng)
        for q, L in enumerate(LAGS):
            v = ts[(ny_w >= 2004 - L) & (ny_w <= 2025 - L)]; null[k, q] = np.corrcoef(sd, detr(v))[0, 1]
    mx = np.abs(null).max(1); robs = max(abs(v) for v in r_obs.values())
    for q, L in enumerate(LAGS):
        p = (1 + np.sum(np.abs(null[:, q]) >= abs(r_obs[L]))) / (1 + NS)
        rows.append(dict(test="T2", series=f"{wname}, lag {L} yr", r=r_obs[L], p=p))
        print(f"T2 {wname} lag {L}: r = {r_obs[L]:+.2f}, p = {p:.3f}", flush=True)
    pmax = (1 + np.sum(mx >= robs)) / (1 + NS)
    rows.append(dict(test="T2", series=f"{wname}, max over lags 0-3", r=robs, p=pmax))
    print(f"T2 {wname} max|r| over lags = {robs:.2f}, family-wise p = {pmax:.3f}", flush=True)

# T3 spatial, 5 deg bins
edges = np.arange(-180, 181, 5); brow = []
for k in range(edges.size - 1):
    m = dom & (LON >= edges[k]) & (LON < edges[k + 1])
    if m.sum() < 2: continue
    s = regmean(SL, m); t = regmean(T1, m)
    rs = R.trend_bcz(yrs + 0.5, s, monthly=False, n_outer=150); rt = R.trend_bcz(ny_w[sel] + 0.5, t[sel], monthly=False, n_outer=150)
    tl = lagset(t); sd = detr(s); rr = {}
    for L in (0, 1):
        v = detr(tl[L]); r_ = np.corrcoef(sd, v)[0, 1]
        nl = np.array([np.corrcoef(sd, detr(phase_surr(detr(t), rng)[(ny_w >= 2004 - L) & (ny_w <= 2025 - L)]))[0, 1] for _ in range(1000)])
        rr[L] = (r_, (1 + np.sum(np.abs(nl) >= abs(r_))) / 1001)
    brow.append(dict(lon=edges[k] + 2.5, ncell=int(m.sum()), slope_pct=1000 * rs["slope"] / s.mean(), slope_p=rs["p"],
                     tau_pct=1000 * rt["slope"] / t[sel].mean(), tau_p=rt["p"], tau_mean=t[sel].mean(),
                     r_lag0=rr[0][0], p_lag0=rr[0][1], r_lag1=rr[1][0], p_lag1=rr[1][1]))
B = pd.DataFrame(brow)
for L in (0, 1):
    B[f"fdr_lag{L}"] = R.bh_fdr(B[f"p_lag{L}"].values, 0.10)
B.to_csv(OUT / "f4b_bins.csv", index=False)
x, y = B.slope_pct.values, B.tau_pct.values
r_sp = np.corrcoef(x, y)[0, 1]
null = np.array([np.corrcoef(x, np.roll(y, k))[0, 1] for k in range(1, len(y))])
null = np.r_[null, [np.corrcoef(x, np.roll(y, rng.integers(1, len(y))))[0, 1] for _ in range(2000 - null.size)]]
p_sp = (1 + np.sum(np.abs(null) >= abs(r_sp))) / (1 + null.size)
rows.append(dict(test="T3", series="cross-bin r(slope trend, tau trend)", r=r_sp, p=p_sp, n=len(B)))
for L in (0, 1):
    rows.append(dict(test="T3", series=f"bins with p<0.05, lag {L}", r=float(np.median(B[f'r_lag{L}'])),
                     n=int((B[f'p_lag{L}'] < 0.05).sum()), p=np.nan, fdr=int(B[f'fdr_lag{L}'].sum()), npos=int((B[f'r_lag{L}'] > 0).sum())))
print(f"T3 cross-bin r(trends) = {r_sp:+.2f} (n = {len(B)} bins), circular-shift p = {p_sp:.3f}", flush=True)
for L in (0, 1):
    print(f"T3 per-bin lag {L}: median r {np.median(B[f'r_lag{L}']):+.2f}, positive {int((B[f'r_lag{L}'] > 0).sum())}/{len(B)}, p<0.05 {int((B[f'p_lag{L}'] < 0.05).sum())}, FDR {int(B[f'fdr_lag{L}'].sum())}", flush=True)
pd.DataFrame(rows).to_csv(OUT / "f4b_results.csv", index=False)
np.savez(OUT / "f4b_series.npz", years_slope=yrs, slope=s_dom, years_tau=ny_w, tau_dom=tau_dom, tau_45=tau_45)
print("done", flush=True)
