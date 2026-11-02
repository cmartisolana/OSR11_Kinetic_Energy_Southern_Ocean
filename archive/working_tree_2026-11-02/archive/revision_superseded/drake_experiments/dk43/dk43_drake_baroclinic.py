"""dk43 (derived from dk42): baroclinic transport across Drake Passage from Argo (RG), 2004-2025, by thermal wind between a fixed
northern and southern profile per longitude (65.5-59.5 W): T = (1/|f|) int (Phi_N - Phi_S) dz, Phi = dynamic height anomaly
relative to 1975 dbar; the pair is the northernmost and southernmost cells of the 55-64 S range with complete profiles in all
years. Passage mean over the longitudes; trends 2004-2025 and 2005-2019 (calibrated Theil-Sen, annual).
Original dk42 docstring: baroclinic geostrophic flow of the ACC from the Roemmich-Gilson Argo climatology (1 deg, 2004-2025, annual means
as dk32). Dynamic height anomaly relative to 1975 dbar (TEOS-10 specific volume anomaly integrated in pressure),
geostrophic velocity relative to 1975 dbar (centred differences on the 1 deg grid), eastward component:
  u_surf  surface baroclinic velocity (cm/s)
  u_mean  0-1975 dbar depth mean (cm/s)
  T_bc    baroclinic transport across the ACC band (mean ADT -0.6..+0.2 m on the RG grid, as dk32) per longitude, Sv,
          summed over band cells with complete profiles to 1975 dbar (RG BATHYMETRY_MASK = 1 at every level: values that RG
          maps below the seafloor are discarded); longitudes kept only if every band cell is complete in every year (fixed set).
Series: ACC band and sectors (longitude means of T_bc; band-cell means of u_surf and u_mean), annual 2004-2025.
Streamline bins: 0.1 m bins of mean ADT (-0.9..+0.5), band-cell means of u_mean and u_surf.
Trends: calibrated Theil-Sen on annual values (rstats.trend_bcz, monthly=False, n_outer 150), % of the mean per decade.
Maps: least-squares trend of u_mean and u_surf per cell with n_eff from lag-1 autocorrelation (p-values).
Usage: python dk43_drake_baroclinic.py <revision dir> <argo_rg dir> <dk32 out dir> <out dir>"""
import sys, glob, os, struct
from pathlib import Path
import numpy as np, pandas as pd, gsw
from scipy import stats
REV, ARG, D32, OUT = (Path(a) for a in sys.argv[1:5]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
L = lambda d, k: np.load(ARG / d / f"{k}.npy")
lat, lon, P = L("so_T", "LATITUDE").astype(float), L("so_T", "LONGITUDE").astype(float), L("so_T", "PRESSURE").astype(float)
Tm, Sm = L("so_T", "ARGO_TEMPERATURE_MEAN"), L("so_S", "ARGO_SALINITY_MEAN")
BM = L("so_T", "BATHYMETRY_MASK")                   # 1 = ocean; 0 = below the seafloor (RG maps values there); -9 = none
bad = (Tm < -90) | (Sm < -90) | (BM != 1); Tm = np.where(bad, np.nan, Tm); Sm = np.where(bad, np.nan, Sm)
def years_from_npy():
    TA = sorted(glob.glob(str(ARG / "so_T" / "ARGO_TEMPERATURE_ANOMALY_m*.npy"))); SA_ = sorted(glob.glob(str(ARG / "so_S" / "ARGO_SALINITY_ANOMALY_m*.npy")))
    for i, (ft, fs) in enumerate(zip(TA, SA_)):
        a, b = np.load(ft), np.load(fs); a[np.abs(a) > 900] = np.nan; b[np.abs(b) > 900] = np.nan
        yield 2004 + i, Tm + np.nanmean(a, 0), Sm + np.nanmean(b, 0)
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
def years_from_ext():
    fl = sorted(glob.glob(str(ARG / "RG_ArgoClim_20*.nc"))); byy = {}
    for f in fl: byy.setdefault(int(os.path.basename(f)[12:16]), []).append(f)
    for y, ff in sorted(byy.items()):
        TT, SS = [], []
        for f in ff:
            d = nc3(f, ["LATITUDE", "ARGO_TEMPERATURE_ANOMALY", "ARGO_SALINITY_ANOMALY"]); j = np.isin(d["LATITUDE"], lat)
            a = d["ARGO_TEMPERATURE_ANOMALY"].reshape(-1, Tm.shape[0], d["LATITUDE"].size, 360)[0][:, j]; b = d["ARGO_SALINITY_ANOMALY"].reshape(-1, Tm.shape[0], d["LATITUDE"].size, 360)[0][:, j]
            a[np.abs(a) > 900] = np.nan; b[np.abs(b) > 900] = np.nan; TT.append(a); SS.append(b)
        yield y, Tm + np.nanmean(TT, 0), Sm + np.nanmean(SS, 0)
Pz = P[:, None, None] * np.ones_like(Tm); LO = np.broadcast_to(lon[None, None, :], Tm.shape); LA = np.broadcast_to(lat[None, :, None], Tm.shape)
Z = -gsw.z_from_p(P, -55.0)                        # depth (m) of the pressure levels, positive down
dz = np.diff(np.r_[0, Z]); dz[0] = Z[0]            # layer thickness for depth integrals (surface layer to first level)
f = 2 * 7.2921e-5 * np.sin(np.radians(lat)); DY = 111195.0
def dyn_height(SA, CT):
    """Dynamic height anomaly (m2 s-2) relative to P[-1] at each level: integral of specific volume anomaly in pressure (Pa)."""
    sva = gsw.specvol_anom_standard(SA, CT, Pz)
    dp = np.diff(P) * 1e4; mid = 0.5 * (sva[1:] + sva[:-1]) * dp[:, None, None]
    phi = np.zeros_like(sva); phi[:-1] = np.cumsum(mid[::-1], 0)[::-1]
    phi[~np.isfinite(sva).all(0)[None].repeat(len(P), 0)] = np.nan   # incomplete profiles -> NaN
    return phi
YR, PHI = [], []
ILON = (lon >= 294) & (lon <= 301); JLAT = (lat <= -54.5) & (lat >= -64)
for gen in (years_from_npy(), years_from_ext()):
    for y, T_, S_ in gen:
        SA = gsw.SA_from_SP(S_, Pz, LO, LA); CT = gsw.CT_from_t(SA, T_, Pz); phi = dyn_height(SA, CT)
        YR.append(y); PHI.append(phi[:, JLAT][:, :, ILON]); print(y)
YR = np.array(YR); PHI = np.array(PHI); la, lo = lat[JLAT], lon[ILON]; t = YR.astype(float)
ok = np.isfinite(PHI).all((0, 1))                      # (lat, lon) complete in every year
rows, TL = [], np.full((len(YR), lo.size), np.nan)
for i in range(lo.size):
    jj = np.nonzero(ok[:, i])[0]
    if jj.size < 2: continue
    jS, jN = jj.min(), jj.max(); fm = 2 * 7.2921e-5 * np.sin(np.radians(0.5 * (la[jS] + la[jN])))
    TL[:, i] = 1e-6 * np.nansum((PHI[:, :, jN, i] - PHI[:, :, jS, i]) * dz[None, :], 1) / abs(fm)
    rows.append(dict(lon=lo[i] - 360, lat_S=la[jS], lat_N=la[jN], mean_Sv=TL[:, i].mean()))
print(pd.DataFrame(rows).round(2).to_string(index=False))
y = np.nanmean(TL, 1); out = []
for per, m in (("2004-2025", np.ones(len(YR), bool)), ("2005-2019", (YR >= 2005) & (YR <= 2019))):
    r = R.trend_bcz(t[m], y[m], monthly=False, n_outer=150); mu = y[m].mean()
    out.append(dict(series="drake_Tbc_annual", period=per, mean=mu, n=int(m.sum()), slope=r["slope"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"], p=r["p"], pct_dec=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu))
D = pd.DataFrame(out); D.to_csv(OUT / "dk43_drake_trends.csv", index=False); pd.set_option("display.width", 220); print(D.round(3).to_string(index=False))
pd.DataFrame({"year": YR, "Tbc_Sv": y}).to_csv(OUT / "dk43_drake_annual.csv", index=False); np.savez(OUT / "dk43_drake_lon.npz", years=YR, lon=lo, Tlon=TL)
