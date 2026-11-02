"""B5: two further measures of meander shape at the five meanders of dk13 (Referee Gutierrez-Villanueva, comment 3).
(1) Meridional-flow fraction of the annual-mean flow: f = sum(w v^2) / sum(w (u^2 + v^2)) over the SAF-PF band of dk12
    (mean SAF and PF streamlines, ice-free in every month) in the ridge + downstream longitudes of each meander; w = cos(lat).
    Annual means suppress most transient eddies, so f measures how far the mean jet turns meridionally (flexing raises it).
    Trend: trend_bcz on the 25 annual values (AR(1), calibrated), % of the mean per decade.
(2) Meander amplitude of the mean path: +0.18 m contour of ADT relative to its ACC-band mean (as B3/B4, 150 km smoothing)
    from 5-year mean fields (2001-05 ... 2021-25) and from the two halves (2001-2012, 2013-2025); amplitude = maximum minus
    minimum path latitude between the start of the upstream and the end of the downstream region (km); crest and trough
    longitudes. Difference late - early with a bootstrap over years (resampling years within each half, 500 draws).
Usage: python b5_meander_shape.py <revision dir> <dk12 npz> <dk5 out dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, contourpy
from scipy.ndimage import gaussian_filter1d
REV, Z12, D5, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
MEANDERS = {"Kerguelen Plateau": ((55, 65), (65, 80), (80, 95)), "Southeast Indian Ridge": ((120, 130), (130, 140), (140, 150)),
            "Macquarie Ridge": ((150, 155), (155, 160), (160, 175)), "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
            "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35))}
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
af = OUT / "b5_annual_fields.npz"
if not af.exists():
    U = np.zeros((25, lat.size, lon.size), np.float32); V = np.zeros_like(U); A = np.zeros_like(U); N = np.zeros_like(U)
    for k, f in enumerate(files):
        with xr.open_dataset(f) as d:
            u, v, a = d.u_mean.values, d.v_mean.values, d.adt_mean.values; ok = np.isfinite(u) & np.isfinite(v) & np.isfinite(a)
        y = k // 12; U[y][ok] += u[ok]; V[y][ok] += v[ok]; A[y][ok] += a[ok]; N[y][ok] += 1
    N[N == 0] = np.nan; U /= N; V /= N; A /= N
    np.savez(af, U=U, V=V, A=A, N=N); print("annual fields saved", flush=True)
Zf = np.load(af); U, V, A = Zf["U"], Zf["V"], Zf["A"]
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
ACCM = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice
EXCL = (mean_adt > 0.45) | (lat[:, None] > -40) | ~np.isfinite(mean_adt)
# (1) band as dk12
z = np.load(Z12); LEV_SAF, LEV_PF = float(z["lev_saf"]), float(z["lev_pf"])
adtc = xr.open_dataset(REV / "outputs" / "clim_allsat.nc").adt.interp(latitude=lat, longitude=lon).values
band = np.isfinite(adtc) & (adtc >= LEV_PF) & (adtc <= LEV_SAF) & ~ice & np.all(np.isfinite(U), 0)
W = np.cos(np.radians(lat))[:, None] * band
X = np.arange(2001, 2026) + 0.5
rows = []
for nm, (u_, r_, d_) in MEANDERS.items():
    m = ((lon > r_[0]) & (lon < r_[1])) | ((lon > d_[0]) & (lon < d_[1]))
    if nm.startswith("Drake"): m &= ~((lon > -60) & (lon < -50))
    w = W[:, m]
    f = np.array([np.nansum(w * V[y][:, m] ** 2) / np.nansum(w * (U[y][:, m] ** 2 + V[y][:, m] ** 2)) for y in range(25)])
    r = R.trend_bcz(X, f, monthly=False, n_outer=300); mu = f.mean()
    rows.append(dict(meander=nm, metric="meridional-flow fraction, annual-mean flow", mean=mu, pct=1000 * r["slope"] / mu,
                     lo=1000 * r["ci_lo"] / mu, hi=1000 * r["ci_hi"] / mu, p=r["p"]))
    print(nm, "fraction", round(mu, 3), round(rows[-1]["pct"], 2), round(r["p"], 3), flush=True)
pd.DataFrame(rows).to_csv(OUT / "b5_meridional_fraction.csv", index=False)
# (2) mean-path amplitude
PAD = int(30 / np.diff(lon).mean()); lonw = np.r_[lon[-PAD:] - 360, lon, lon[:PAD] + 360]
wrap = lambda a: np.c_[a[:, -PAD:], a, a[:, :PAD]]
RE, DS, SM = 6371.0, 10.0, 150.0


def path(adt_mean):
    rel = adt_mean - np.nanmean(adt_mean[ACCM]); rel[EXCL] = np.nan
    gen = contourpy.contour_generator(lonw, lat, np.ma.masked_invalid(wrap(rel)), name="serial"); P = []
    for L in gen.lines(0.18):
        if np.ptp(L[:, 0]) < 5: continue
        if L[0, 0] > L[-1, 0]: L = L[::-1]
        x = RE * np.cos(np.radians(L[:, 1])) * np.radians(L[:, 0]); y = RE * np.radians(L[:, 1])
        s = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
        if s[-1] < 5 * DS: continue
        si = np.arange(0, s[-1], DS)
        P.append((gaussian_filter1d(np.interp(si, s, L[:, 0]), SM / DS, mode="nearest"), gaussian_filter1d(np.interp(si, s, L[:, 1]), SM / DS, mode="nearest")))
    return P


def amp(P, a, b):
    best = None
    for lo_s, la_s in P:
        sel = (lo_s >= a) & (lo_s <= b)
        if sel.sum() < 3: continue
        cov = np.ptp(lo_s[sel])
        if best is None or cov > best[0]: best = (cov, lo_s[sel], la_s[sel])
    if best is None or best[0] < 0.8 * (b - a): return np.nan, np.nan, np.nan
    _, lo_b, la_b = best
    return 111.195 * (la_b.max() - la_b.min()), lo_b[np.argmax(la_b)], lo_b[np.argmin(la_b)]


rows = []
for nm, (u_, r_, d_) in MEANDERS.items():
    a, b = u_[0], d_[1]
    for lab, ys in [(f"{y0}-{y0 + 4}", list(range(y0 - 2001, y0 - 2001 + 5))) for y0 in range(2001, 2026, 5)] + [("2001-2012", list(range(12))), ("2013-2025", list(range(12, 25)))]:
        am, lc, lt = amp(path(np.nanmean(A[ys], 0)), a, b)
        rows.append(dict(meander=nm, period=lab, amplitude_km=am, crest_lon=lc, trough_lon=lt))
rng = np.random.default_rng(20261020); NBOOT = 500
DD = {nm: [] for nm in MEANDERS}
for it in range(NBOOT):                                   # one draw of years per half, shared by all meanders
    e = rng.integers(0, 12, 12); l = 12 + rng.integers(0, 13, 13)
    Pe, Pl = path(np.nanmean(A[e], 0)), path(np.nanmean(A[l], 0))
    for nm, (u_, r_, d_) in MEANDERS.items():
        DD[nm].append(amp(Pl, u_[0], d_[1])[0] - amp(Pe, u_[0], d_[1])[0])
    if it % 100 == 0: print("boot", it, flush=True)
for nm in MEANDERS:
    dd = np.array(DD[nm]); dd = dd[np.isfinite(dd)]
    if dd.size < 50:                                      # path does not span the region in most draws (e.g. Malvinas loop)
        rows.append(dict(meander=nm, period="late-early", amplitude_km=np.nan, n_boot=dd.size)); print(nm, "not measurable", dd.size); continue
    e0 = [r for r in rows if r["meander"] == nm and r["period"] == "2001-2012"][0]["amplitude_km"]
    l0 = [r for r in rows if r["meander"] == nm and r["period"] == "2013-2025"][0]["amplitude_km"]
    rows.append(dict(meander=nm, period="late-early", amplitude_km=l0 - e0, boot_lo=np.percentile(dd, 2.5), boot_hi=np.percentile(dd, 97.5),
                     boot_p=2 * min((dd <= 0).mean(), (dd >= 0).mean()), n_boot=dd.size))
    print(nm, "amplitude early/late", round(e0, 1), round(l0, 1), "diff CI", np.round(np.percentile(dd, [2.5, 97.5]), 1), flush=True)
D = pd.DataFrame(rows); D.to_csv(OUT / "b5_meander_amplitude.csv", index=False); print(D.to_string(index=False)); print("done")
