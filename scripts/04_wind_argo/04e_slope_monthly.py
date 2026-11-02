"""F4m: monthly cross-stream slope of sigma0 = 27.3 (28 Oct 2026; author request: monthly points for the isopycnal trends).
Same definitions as f4_tilt.py, but from the monthly Roemmich-Gilson fields instead of annual means:
 - RG 2019 release, 2004-2018 mean + monthly anomalies (dk32a extracts so_T, so_S), extended 2019-2025 with the monthly RG files;
   T, S -> TEOS-10 (SA, CT, sigma0) each month; pressure of sigma0 = 27.3 by linear interpolation at the first crossing
   (NaN if it outcrops or lies below 1975 dbar), exactly as dk32b but per month.
 - Signed slope s = grad(p) . n, n along grad(mean ADT) (f4_tilt_maps.nc), dbar per 100 km, centred differences, periodic.
 - Domain: f4 domain (ACC band + 100 km, ice removed); a cell is kept only if s is defined in all 264 months.
Tests (specification fixed before results):
 M1 domain mean (cos-lat), monthly series 2004-2025: rstats.trend_bcz monthly (deseasonalised; max p of AIC and AR(12)).
    For comparison, the annual means of the monthly series with the annual test.
 M2 sectors (as f4_tilt.py), monthly test.
 M3 5 deg longitude bins: monthly test (n_outer 150); cross-bin r with the EKE trend (f4d_eke_bins.csv, 2004-2025), p from all
    circular shifts.
 M4 consistency: r between the annual means of the monthly domain series and the annual f4 series.

Purpose : monthly sigma0 = 27.3 cross-stream slope on the fixed domain; domain, sector and 5-deg trends (Sect. 3.3, SI S5).
Inputs  : $ARGO_RG_DIR (00d), outputs/f4/f4_tilt_maps.nc, f4_regional_series.npz (04c), outputs/f4/f4d_eke_bins.csv (04d)
Outputs : outputs/f4m/f4m_p273_monthly.npz (8.8 MB cache; input of 04f, 04g), f4m_results.csv (sector values, SI S5), f4m_series.npz, f4m_lon5_trends.csv
Run     : python scripts/04_wind_argo/04e_slope_monthly.py
Seeds   : trend_bcz seed 0.
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, glob, os, struct
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, gsw

from osr11 import rstats as R, paths as PA  # noqa: E402
ARG, F4, FD, OUT = PA.ARGO_RG, PA.out("f4"), PA.out("f4"), PA.out("f4m")

L = lambda d, k: np.load(ARG / d / f"{k}.npy")
lat, lonr, P = L("so_T", "LATITUDE").astype(float), L("so_T", "LONGITUDE").astype(float), L("so_T", "PRESSURE").astype(float)
Tm, Sm = L("so_T", "ARGO_TEMPERATURE_MEAN").astype(float), L("so_S", "ARGO_SALINITY_MEAN").astype(float)
bad = (Tm < -90) | (Sm < -90); Tm[bad] = np.nan; Sm[bad] = np.nan
Pz = P[:, None, None] * np.ones_like(Tm); LO = np.broadcast_to(lonr[None, None, :], Tm.shape); LA = np.broadcast_to(lat[None, :, None], Tm.shape)
SIG = 27.3
CACHE = OUT / "f4m_p273_monthly.npz"

TY = {1: ("b", 1), 2: ("c", 1), 3: (">i2", 2), 4: (">i4", 4), 5: (">f4", 4), 6: (">f8", 8)}
def nc3(F, keys):
    raw = open(F, "rb").read(); pos = [0]
    def rd(n): b = raw[pos[0]:pos[0] + n]; pos[0] += n; return b
    def i4(): return struct.unpack(">i", rd(4))[0]
    def name(): n = i4(); s = rd(n).decode(); rd((4 - n % 4) % 4); return s
    def atts():
        i4(); n = i4()
        for _ in range(n): name(); t, m = i4(), i4(); sz = TY[t][1] * m; rd(sz); rd((4 - sz % 4) % 4)
    mg = rd(4); off8 = mg[3] == 2; i4(); i4(); nd = i4(); dims = [(name(), i4()) for _ in range(nd)]; atts(); i4(); nv = i4(); V = {}
    for _ in range(nv):
        nm = name(); nid = i4(); did = [i4() for _ in range(nid)]; atts(); t = i4(); vs = i4(); beg = struct.unpack(">q", rd(8))[0] if off8 else i4()
        V[nm] = (did, t, vs, beg)
    out = {}
    for k in keys:
        k2 = [v for v in V if v.upper() == k.upper()][0]; did, t, vs, beg = V[k2]; shp = [max(dims[d][1], 1) for d in did]
        out[k] = np.frombuffer(raw[beg:beg + int(np.prod(shp)) * TY[t][1]], TY[t][0]).reshape(shp).astype(float)
    return out


def iso_press(s0, s):
    below = s0 >= s; first = np.argmax(below, 0); has = below.any(0) & (first > 0) & ~below[0]
    k1 = np.clip(first, 1, len(P) - 1); k0 = k1 - 1
    g0 = np.take_along_axis(s0, k0[None], 0)[0]; g1 = np.take_along_axis(s0, k1[None], 0)[0]
    p = P[k0] + (s - g0) / (g1 - g0) * (P[k1] - P[k0]); return np.where(has & np.isfinite(p), p, np.nan)


def months():
    TA = sorted(glob.glob(str(ARG / "so_T" / "ARGO_TEMPERATURE_ANOMALY_m*.npy"))); SA_ = sorted(glob.glob(str(ARG / "so_S" / "ARGO_SALINITY_ANOMALY_m*.npy")))
    for i, (ft, fs) in enumerate(zip(TA, SA_)):
        a, b = np.load(ft).astype(float), np.load(fs).astype(float); a[np.abs(a) > 900] = np.nan; b[np.abs(b) > 900] = np.nan
        for k in range(a.shape[0]):
            yield 2004 + i, k + 1, Tm + a[k], Sm + b[k]
    for f in sorted(glob.glob(str(ARG / "RG_ArgoClim_20*.nc"))):
        y, m = int(os.path.basename(f)[12:16]), int(os.path.basename(f)[16:18])
        d = nc3(f, ["LATITUDE", "ARGO_TEMPERATURE_ANOMALY", "ARGO_SALINITY_ANOMALY"]); j = np.isin(d["LATITUDE"], lat); nl = d["LATITUDE"].size
        a = d["ARGO_TEMPERATURE_ANOMALY"].reshape(-1, len(P), nl, 360)[0][:, j]; b = d["ARGO_SALINITY_ANOMALY"].reshape(-1, len(P), nl, 360)[0][:, j]
        a[np.abs(a) > 900] = np.nan; b[np.abs(b) > 900] = np.nan
        yield y, m, Tm + a, Sm + b


if CACHE.exists():
    C = np.load(CACHE); YY, MM, PM = C["year"], C["month"], C["p"]
else:
    YY, MM, PM = [], [], []
    for y, m, T_, S_ in months():
        SA = gsw.SA_from_SP(S_, Pz, LO, LA); CT = gsw.CT_from_t(SA, T_, Pz); s0 = gsw.sigma0(SA, CT)
        PM.append(iso_press(s0, SIG).astype(np.float32)); YY.append(y); MM.append(m)
        if m == 1: print(y, flush=True)
    YY, MM, PM = np.array(YY), np.array(MM), np.array(PM)
    np.savez_compressed(CACHE, year=YY, month=MM, p=PM, lat=lat, lon=lonr)
assert len(YY) == 264 and YY[0] == 2004 and YY[-1] == 2025, (len(YY), YY[[0, -1]])
x = YY + (MM - 0.5) / 12

# slope on the f4 grid (lon -180..180)
M = xr.open_dataset(F4 / "f4_tilt_maps.nc"); glat, glon = M.lat.values, M.lon.values
L0 = np.where(lonr > 180, lonr - 360, lonr); o = np.argsort(L0); assert np.allclose(L0[o], glon) and np.allclose(lat, glat)
KM = 111.195; adt = M["adt"].values
def grad(F):
    dy = np.gradient(F, glat, axis=-2) / KM * 100
    dx = (np.roll(F, -1, axis=-1) - np.roll(F, 1, axis=-1)) / 2.0 / (KM * np.cos(np.radians(glat))[:, None]) * 100
    return dx, dy
ax_, ay_ = grad(adt); an = np.hypot(ax_, ay_)
Pm = PM[:, :, o].astype(float); dx, dy = grad(Pm); SL = dx * ax_ / an + dy * ay_ / an
dom0 = (M["buffer100"].values > 0) & (M["ice"].values == 0) & np.isfinite(M["signed_27p3_pct"].values)
dom = dom0 & np.isfinite(SL).all(0)
print(f"cells: annual domain {int(dom0.sum())}, defined in all 264 months {int(dom.sum())}", flush=True)
LONG, LATG = np.meshgrid(glon, glat); w = np.cos(np.radians(LATG))
def regmean(m): return np.array([np.sum(SL[t][m] * w[m]) / np.sum(w[m]) for t in range(SL.shape[0])])

rows = []
def mtest(name, y, n_outer=300):
    r = R.trend_bcz(x, y, month=MM, n_outer=n_outer); mu = y.mean()
    rows.append(dict(series=name, slope_dec=10 * r["slope"], lo_dec=10 * r["ci_lo"], hi_dec=10 * r["ci_hi"], pct=1000 * r["slope"] / mu,
                     pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"], p_aic=r["p_aic"], p_ar12=r["p_ar12"], mean=mu, n=len(y)))
    print(f"{name}: {rows[-1]['pct']:+.2f} % per decade [{rows[-1]['pct_lo']:+.2f}, {rows[-1]['pct_hi']:+.2f}] p={r['p']:.3f}", flush=True)
    return r
s_dom = regmean(dom); mtest("27.3 signed, domain, monthly test", s_dom)
ya = np.array([s_dom[YY == k].mean() for k in range(2004, 2026)])
r = R.trend_bcz(np.arange(2004, 2026) + 0.5, ya, monthly=False, n_outer=300)
rows.append(dict(series="27.3 signed, domain, annual means of monthly, annual test", slope_dec=10 * r["slope"], pct=1000 * r["slope"] / ya.mean(),
                 pct_lo=1000 * r["ci_lo"] / ya.mean(), pct_hi=1000 * r["ci_hi"] / ya.mean(), p=r["p"], mean=ya.mean(), n=22))
SER = np.load(F4 / "f4_regional_series.npz"); yA = SER["27.3__signed__buffer100__ACC"]
rows.append(dict(series="r(annual means of monthly, f4 annual series)", r=float(np.corrcoef(ya, yA)[0, 1]), mean=float(yA.mean())))
print(rows[-2], rows[-1], flush=True)
SECT = {"Atlantic": (LONG >= -70) & (LONG < 20), "Indian": (LONG >= 20) & (LONG < 147), "Pacific": (LONG >= 147) | (LONG < -70)}
sec = {}
for k, ms in SECT.items():
    sec[k] = regmean(dom & ms); mtest(f"27.3 signed, {k}, monthly test", sec[k])

EB = pd.read_csv(FD / "f4d_eke_bins.csv").sort_values("lon"); brow = []
for c in EB.lon.values:
    m = dom & (LONG >= c - 2.5) & (LONG < c + 2.5)
    if m.sum() < 2: continue
    s = regmean(m); r = R.trend_bcz(x, s, month=MM, n_outer=150)
    brow.append(dict(lon=c, ncell=int(m.sum()), pct=1000 * r["slope"] / s.mean(), lo=1000 * r["ci_lo"] / s.mean(), hi=1000 * r["ci_hi"] / s.mean(), p=r["p"]))
B = pd.DataFrame(brow).merge(EB[["lon", "eke_pct"]], on="lon"); B["fdr10"] = R.bh_fdr(B.p.values, 0.10)
B.to_csv(OUT / "f4m_lon5_trends.csv", index=False)
a, b = B.pct.values, B.eke_pct.values; rr = np.corrcoef(a, b)[0, 1]
null = np.array([np.corrcoef(a, np.roll(b, k))[0, 1] for k in range(1, len(b))]); pp = (1 + np.sum(np.abs(null) >= abs(rr))) / (1 + null.size)
rows.append(dict(series="cross-bin r(monthly-test slope trend, EKE trend 2004-2025)", r=rr, p=pp, n=len(B)))
print(f"bins {len(B)}: p<0.05 {int((B.p < 0.05).sum())} ({int(((B.p < 0.05) & (B.pct > 0)).sum())} positive), FDR {int(B.fdr10.sum())}; r with EKE {rr:+.2f} p={pp:.3f}", flush=True)
pd.DataFrame(rows).to_csv(OUT / "f4m_results.csv", index=False)
np.savez(OUT / "f4m_series.npz", x=x, year=YY, month=MM, domain=s_dom, **sec)
print("done", flush=True)
