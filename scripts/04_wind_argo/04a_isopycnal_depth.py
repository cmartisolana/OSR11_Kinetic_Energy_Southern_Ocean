"""dk32b: isopycnal depths and cross-ACC isopycnal slopes from the Roemmich-Gilson Argo climatology (RG 2019 release,
2004-2018 mean + monthly anomalies, extended with the RG monthly files when present), south of 35 S, 1 deg, 0-2000 dbar.
T (in situ, ITS-90) and practical salinity -> annual means -> TEOS-10 (gsw): SA, CT, sigma0. Pressure of sigma0 surfaces
by linear interpolation (first crossing from the surface; NaN if the surface outcrops or lies below 1975 dbar).
ACC band: mean ADT (altimetry, 2001-2025 mean, -0.6 to +0.2 m) interpolated to the RG grid (nearest 1/8 deg pixel average).
Cross-ACC slope per 1 deg longitude: least-squares slope of isopycnal pressure against latitude across the band cells
(dbar per degree latitude; positive = deepening northward), also the depth difference between the northern and southern
thirds of the band. Sector and ACC-band averages; trends (least squares on annual values; 95 % with n_eff from lag-1
autocorrelation of residuals).

Purpose : annual-mean pressure of sigma0 = 27.3, 27.5, 27.7 on the RG 1-deg grid, 2004-2025 (input of 04c).
Inputs  : $ARGO_RG_DIR (00d; plus RG_ArgoClim_20*.nc monthly extensions), outputs/base/mean_adt.nc
Outputs : outputs/drake/dk32/dk32_iso.npz (5.8 MB; input of 04c), dk32_slopes.csv, dk32_slope_lon.npz, dk32_drake_sigma0.npz
Run     : python scripts/04_wind_argo/04a_isopycnal_depth.py
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, glob, os, struct
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, gsw
from scipy import stats
from osr11 import paths as PA
ARG, OUT = PA.ARGO_RG, PA.dk(32)
L = lambda d, k: np.load(ARG / d / f"{k}.npy")
lat, lon, P = L("so_T", "LATITUDE").astype(float), L("so_T", "LONGITUDE").astype(float), L("so_T", "PRESSURE").astype(float)
Tm, Sm = L("so_T", "ARGO_TEMPERATURE_MEAN"), L("so_S", "ARGO_SALINITY_MEAN")
bad = (Tm < -90) | (Sm < -90); Tm = np.where(bad, np.nan, Tm); Sm = np.where(bad, np.nan, Sm)
def years_from_npy():
    TA = sorted(glob.glob(str(ARG / "so_T" / "ARGO_TEMPERATURE_ANOMALY_m*.npy"))); SA_ = sorted(glob.glob(str(ARG / "so_S" / "ARGO_SALINITY_ANOMALY_m*.npy")))
    for i, (ft, fs) in enumerate(zip(TA, SA_)):
        a, b = np.load(ft), np.load(fs); a[np.abs(a) > 900] = np.nan; b[np.abs(b) > 900] = np.nan
        yield 2004 + i, Tm + a, Sm + b, a.shape[0]
# --- RG monthly extension files (netCDF-3): minimal reader
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
    fl = sorted(glob.glob(str(ARG / "RG_ArgoClim_20*.nc")))
    if not fl: return
    js = None; byy = {}
    for f in fl:
        y = int(os.path.basename(f)[12:16]); byy.setdefault(y, []).append(f)
    for y, ff in sorted(byy.items()):
        TT, SS = [], []
        for f in ff:
            d = nc3(f, ["LATITUDE", "ARGO_TEMPERATURE_ANOMALY", "ARGO_SALINITY_ANOMALY"])
            j = np.isin(d["LATITUDE"], lat); a, b = d["ARGO_TEMPERATURE_ANOMALY"].reshape(-1, *Tm.shape[0:1], d["LATITUDE"].size, 360)[0][:, j], d["ARGO_SALINITY_ANOMALY"].reshape(-1, Tm.shape[0], d["LATITUDE"].size, 360)[0][:, j]
            a[np.abs(a) > 900] = np.nan; b[np.abs(b) > 900] = np.nan; TT.append(a); SS.append(b)
        yield y, Tm + np.nanmean(TT, 0)[None], Sm + np.nanmean(SS, 0)[None], len(ff)
SIG = [27.3, 27.5, 27.7]
Pz = P[:, None, None] * np.ones_like(Tm); LO = np.broadcast_to(lon[None, None, :], Tm.shape); LA = np.broadcast_to(lat[None, :, None], Tm.shape)
def iso_press(s0):
    out = {}
    for s in SIG:
        below = s0 >= s; first = np.argmax(below, 0); has = below.any(0) & (first > 0) & ~below[0]
        k1 = np.clip(first, 1, len(P) - 1); k0 = k1 - 1
        g0 = np.take_along_axis(s0, k0[None], 0)[0]; g1 = np.take_along_axis(s0, k1[None], 0)[0]
        p = P[k0] + (s - g0) / (g1 - g0) * (P[k1] - P[k0]); out[s] = np.where(has & np.isfinite(p), p, np.nan)
    return out
YR, ISO, NM = [], {s: [] for s in SIG}, []
DBOX = (lon >= 288) & (lon <= 302); DLAT = lat <= -53; S0D = []
for gen in (years_from_npy(), years_from_ext() or []):
    for y, T_, S_, nm in gen:
        Ta, Sa = np.nanmean(T_, 0) if T_.ndim == 4 else T_, np.nanmean(S_, 0) if S_.ndim == 4 else S_
        SA = gsw.SA_from_SP(Sa, Pz, LO, LA); CT = gsw.CT_from_t(SA, Ta, Pz); s0 = gsw.sigma0(SA, CT)
        ip = iso_press(s0); S0D.append(s0[:, DLAT][:, :, DBOX]); YR.append(y); NM.append(nm); [ISO[s].append(ip[s]) for s in SIG]
        print(y, nm, "months", {s: round(float(np.nanmean(ip[s])), 0) for s in SIG})
YR = np.array(YR); ISO = {s: np.array(v) for s, v in ISO.items()}
# ACC band on the RG grid from the altimetric mean ADT
adt = xr.open_dataset(PA.MEAN_ADT).mean_adt
lon180 = np.where(lon > 180, lon - 360, lon)
A1 = np.array([[float(adt.sel(latitude=slice(la - 0.5, la + 0.5), longitude=slice(lo - 0.5, lo + 0.5)).mean()) for lo in lon180] for la in lat])
band = (A1 >= -0.6) & (A1 <= 0.2)
np.savez(OUT / "dk32_drake_sigma0.npz", years=np.array(YR), P=P, lat=lat[DLAT], lon=lon[DBOX], s0=np.array(S0D))
np.savez(OUT / "dk32_iso.npz", years=YR, nmonths=np.array(NM), lat=lat, lon=lon, adt=A1, band=band, **{f"p{s}": v for s, v in ISO.items()})
# cross-band slope per longitude
def slope_lon(Z):
    out = np.full((Z.shape[0], lon.size), np.nan)
    for i in range(lon.size):
        j = band[:, i]
        if j.sum() < 3: continue
        z = Z[:, j, i]; ok = np.isfinite(z).all(0)
        if ok.sum() < 3: continue
        x = lat[j][ok]; zz = z[:, ok]; xc = x - x.mean()
        out[:, i] = (zz - zz.mean(1, keepdims=True)) @ xc / (xc @ xc)
    return out
SECT = {"ACC band": None, "Pacific": (147, 290), "Indian": (20, 147), "Atlantic": (290, 380)}
t = YR - YR.mean()
def trend(y):
    b, a = np.polyfit(t, y, 1); r = y - a - b * t; r1 = np.corrcoef(r[1:], r[:-1])[0, 1]; ne = max(3, len(y) * (1 - r1) / (1 + r1))
    se = np.sqrt((r ** 2).sum() / (ne - 2) / (t ** 2).sum()); return 10 * b, 10 * stats.t.ppf(0.975, ne - 2) * se, 2 * stats.t.sf(abs(b / se), ne - 2)
rows = []; SL = {}
for s in SIG:
    sl = slope_lon(ISO[s]); SL[s] = sl
    for sec, rg in SECT.items():
        l = np.mod(lon, 360); lm = np.ones(lon.size, bool) if rg is None else (((l >= rg[0]) & (l < rg[1])) | ((l + 360 >= rg[0]) & (l + 360 < rg[1])))
        cols = lm & np.isfinite(sl).all(0)
        y = sl[:, cols].mean(1); b, ci, p = trend(y)
        zb = np.nanmean(np.where(band[None] & lm[None, None, :], ISO[s], np.nan), axis=(1, 2)); bz, ciz, pz = trend(zb)
        rows.append(dict(sigma0=s, region=sec, ncols=int(cols.sum()), slope_mean=y.mean(), slope_trend_pct=100 * b / y.mean(), slope_ci_pct=100 * ci / abs(y.mean()), slope_p=p,
                         depth_mean=zb.mean(), depth_trend_m=bz, depth_ci_m=ciz, depth_p=pz))
D = pd.DataFrame(rows); D.to_csv(OUT / "dk32_slopes.csv", index=False); pd.set_option("display.width", 200); print(D.round(3).to_string(index=False))
np.savez(OUT / "dk32_slope_lon.npz", years=YR, lon=lon, **{f"sl{s}": v for s, v in SL.items()})
