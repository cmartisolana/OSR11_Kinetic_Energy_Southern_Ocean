"""dk42: baroclinic geostrophic flow of the ACC from the Roemmich-Gilson Argo climatology (1 deg, 2004-2025, annual means
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
Usage: python dk42_argo_baroclinic.py <revision dir> <argo_rg dir> <dk32 out dir> <out dir>"""
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
YR, US, UM, TB = [], [], [], []
A32 = np.load(D32 / "dk32_iso.npz"); band = A32["band"]; ADT1 = A32["adt"]
for gen in (years_from_npy(), years_from_ext()):
    for y, T_, S_ in gen:
        SA = gsw.SA_from_SP(S_, Pz, LO, LA); CT = gsw.CT_from_t(SA, T_, Pz); phi = dyn_height(SA, CT)
        dphidy = np.full_like(phi, np.nan); dphidy[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) / (2 * DY)
        u = -dphidy / f[None, :, None]                               # m/s, eastward, relative to 1975 dbar
        us = u[0]; um = np.nansum(u * dz[:, None, None], 0) / Z[-1]; um[~np.isfinite(u).all(0)] = np.nan
        tb = np.nansum(u * dz[:, None, None], 0) * DY * 1e-6          # Sv per cell (full-depth 0-1975, 1 deg of latitude)
        tb[~np.isfinite(u).all(0)] = np.nan
        YR.append(y); US.append(100 * us); UM.append(100 * um); TB.append(tb); print(y, "band mean u_surf %.1f cm/s" % np.nanmean(us[band] * 100))
YR = np.array(YR); US, UM, TB = (np.array(v) for v in (US, UM, TB)); t = YR.astype(float)
# longitudes with every band cell complete in every year
okcol = np.array([band[:, i].sum() >= 2 and np.isfinite(TB[:, band[:, i], i]).all() for i in range(lon.size)])
Tlon = np.where(okcol[None], np.where(band[None], TB, 0).sum(1), np.nan)      # Sv per longitude
print("longitudes kept", okcol.sum(), "of", (band.any(0)).sum())
SECT = {"ACC band": None, "Pacific": (147, 290), "Indian": (20, 147), "Atlantic": (290, 380)}
def secmask(rg):
    l = np.mod(lon, 360); return np.ones(lon.size, bool) if rg is None else (((l >= rg[0]) & (l < rg[1])) | ((l + 360 >= rg[0]) & (l + 360 < rg[1])))
def tr(y):
    r = R.trend_bcz(t, y, monthly=False, n_outer=150); mu = y.mean()
    return dict(mean=mu, slope=r["slope"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"], p=r["p"], pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu)
rows, SER = [], {"year": YR}
for sec, rg in SECT.items():
    m = secmask(rg); cols = m & okcol; cells = band & m[None, :] & np.isfinite(UM).all(0)
    series = {"T_bc": Tlon[:, cols].mean(1), "u_mean": UM[:, cells].mean(1), "u_surf": US[:, cells].mean(1)}
    for k, y in series.items():
        rows.append(dict(region=sec, var=k, n_cols=int(cols.sum()) if k == "T_bc" else int(cells.sum()), **tr(y))); SER[f"{sec}|{k}"] = y
D = pd.DataFrame(rows); D.to_csv(OUT / "dk42_trends.csv", index=False); pd.DataFrame(SER).to_csv(OUT / "dk42_series.csv", index=False)
pd.set_option("display.width", 220); print(D.round(3).to_string(index=False))
# streamline bins (0.1 m of mean ADT), band of valid cells
edges = np.arange(-0.9, 0.51, 0.1); rows = []
for lo_, hi_ in zip(edges[:-1], edges[1:]):
    c = (ADT1 >= lo_) & (ADT1 < hi_) & np.isfinite(UM).all(0) & (lat[:, None] > -64)
    if c.sum() < 20: continue
    for k, F in (("u_mean", UM), ("u_surf", US)):
        rows.append(dict(lo=lo_, hi=hi_, mid=0.5 * (lo_ + hi_), var=k, ncells=int(c.sum()), **tr(F[:, c].mean(1))))
B = pd.DataFrame(rows); B.to_csv(OUT / "dk42_stream_bins.csv", index=False); print(B.round(3).to_string(index=False))
# maps: LS trend with n_eff
def lstrend(F):
    tc = t - t.mean(); Fm = F - np.nanmean(F, 0); b = np.einsum("t,tij->ij", tc, np.nan_to_num(Fm)) / (tc @ tc)
    res = Fm - b[None] * tc[:, None, None]; r1 = np.nansum(res[1:] * res[:-1], 0) / np.nansum(res ** 2, 0)
    ne = np.clip(len(t) * (1 - r1) / (1 + r1), 3, len(t)); se = np.sqrt(np.nansum(res ** 2, 0) / (ne - 2) / (tc @ tc))
    p = 2 * stats.t.sf(np.abs(b / se), ne - 2); b[~np.isfinite(F).all(0)] = np.nan; p[~np.isfinite(F).all(0)] = np.nan
    return 10 * b, p
bm, pm = lstrend(UM); bs, ps = lstrend(US)
np.savez(OUT / "dk42_maps.npz", lat=lat, lon=lon, years=YR, u_mean_trend=bm, u_mean_p=pm, u_surf_trend=bs, u_surf_p=ps,
         u_mean_mean=np.nanmean(UM, 0), u_surf_mean=np.nanmean(US, 0), Tlon=Tlon, band=band, adt=ADT1, okcol=okcol, u_mean=UM, u_surf=US)
print("saved", OUT)
