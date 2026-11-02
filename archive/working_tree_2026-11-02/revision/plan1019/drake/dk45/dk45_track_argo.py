"""dk45 (derived from dk42): cross-track baroclinic velocity along the LM Gould line from Argo (RG), 2004-2025. Dynamic
height anomaly relative to 1975 dbar (surface value and 0-1975 dbar depth mean) on the 1 deg grid, bilinearly interpolated onto
the 121 track points (NaN where any of the four cells is incomplete); cross-track velocity u_n = (1/f) dPhi/ds along the track
(eastward-positive normal, as dk40), cm/s. Trends (calibrated Theil-Sen, annual, cm/s per decade): per 25 km passage bin and
per pair of SSH streamlines (dk5 stream edges, points assigned by the 2001-2025 mean ADT at the track points).
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
Usage: python dk45_track_argo.py <revision dir> <drake dir> <argo_rg dir> <dk32 out dir> <dk24 npz> <out dir>"""
import sys, glob, os, struct
from pathlib import Path
import numpy as np, pandas as pd, gsw
from scipy import stats
REV, DR, ARG, D32, SECF, OUT = (Path(a) for a in sys.argv[1:7]); OUT.mkdir(parents=True, exist_ok=True)
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
YR, PS, PM = [], [], []
for gen in (years_from_npy(), years_from_ext()):
    for y, T_, S_ in gen:
        SA = gsw.SA_from_SP(S_, Pz, LO, LA); CT = gsw.CT_from_t(SA, T_, Pz); phi = dyn_height(SA, CT)
        YR.append(y); PS.append(phi[0]); PM.append(np.nansum(phi * dz[:, None, None], 0) / Z[-1] * np.where(np.isfinite(phi).all(0), 1, np.nan)); print(y)
YR = np.array(YR); PS, PM = np.array(PS), np.array(PM); t = YR.astype(float)
# track
SEC = np.load(SECF); tlo, tla, dsx = SEC["lon"], SEC["lat"], SEC["dsx"]; tlo360 = np.mod(tlo, 360)
import xarray as xr
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt; tadt = np.array([float(madt.interp(longitude=a, latitude=b)) for a, b in zip(tlo, tla)])
def interp_track(F):
    """bilinear from the 1 deg grid (cell centres lon, lat) to the track; NaN if any corner is NaN"""
    out = np.full((F.shape[0], tlo.size), np.nan)
    for k, (a, b) in enumerate(zip(tlo360, tla)):
        i0 = int(np.floor(a - lon[0])); j0 = int(np.floor(b - lat[0]))
        if j0 < 0 or j0 + 1 >= lat.size: continue
        i1 = (i0 + 1) % lon.size; wa = a - lon[i0]; wb = b - lat[j0]
        out[:, k] = (1 - wa) * (1 - wb) * F[:, j0, i0] + wa * (1 - wb) * F[:, j0, i1] + (1 - wa) * wb * F[:, j0 + 1, i0] + wa * wb * F[:, j0 + 1, i1]
    return out
tx = np.gradient(tlo) * np.cos(np.radians(tla)); ty = np.gradient(tla); nn = np.hypot(tx, ty); nx, ny = -ty / nn, tx / nn
flip = nx < 0; sgn = np.where(flip, -1.0, 1.0)                         # eastward-positive normal (dk40)
s_m = np.r_[0, np.cumsum(np.hypot(np.diff(tlo) * np.cos(np.radians(0.5 * (tla[1:] + tla[:-1]))), np.diff(tla)))] * DY   # along-track distance (m)
ft = 2 * 7.2921e-5 * np.sin(np.radians(tla))
U = {}
for nm, F in (("u_surf", PS), ("u_mean", PM)):
    ph = interp_track(F); U[nm] = 100 * sgn[None, :] * np.gradient(ph, s_m, axis=1) / ft[None, :]
    print(nm, "mean along track (cm/s): %.1f" % np.nanmean(U[nm]), "valid points", int(np.isfinite(ph).all(0).sum()))
np.savez(OUT / "dk45_track_fields.npz", years=YR, lon=tlo, lat=tla, adt=tadt, dist_km=s_m / 1000, **U)
def tr(y):
    r = R.trend_bcz(t, y, monthly=False, n_outer=150); return dict(mean=y.mean(), slope_dec=10 * r["slope"], ci_lo=10 * r["ci_lo"], ci_hi=10 * r["ci_hi"], p=r["p"])
rows = []
ED = pd.read_csv(DR / "out" / "dk5" / "dk5_stream_edges.csv").sort_values("psi")
for (p0, e0), (p1, e1) in zip(ED[["psi", "adt_edge"]].values[:-1], ED[["psi", "adt_edge"]].values[1:]):
    lo_, hi_ = sorted([e0, e1]); c = (tadt >= lo_) & (tadt < hi_)
    for nm in U:
        ok = c & np.isfinite(U[nm]).all(0)
        if ok.sum() < 2: continue
        rows.append(dict(coord="stream", bin=round(0.5 * (p0 + p1), 4), adt_c=0.5 * (e0 + e1), var=nm, npts=int(ok.sum()), **tr(U[nm][:, ok].mean(1))))
kmb = np.floor(s_m / 25000) * 25 + 12.5
for b in np.unique(kmb):
    for nm in U:
        ok = (kmb == b) & np.isfinite(U[nm]).all(0)
        if ok.sum() < 2: continue
        rows.append(dict(coord="passage", bin=b, adt_c=np.nan, var=nm, npts=int(ok.sum()), **tr(U[nm][:, ok].mean(1))))
D = pd.DataFrame(rows); D.to_csv(OUT / "dk45_track_trends.csv", index=False); pd.set_option("display.width", 200); print(D.round(3).to_string(index=False))
