"""DK5: changes along the LM Gould transect from altimetry, compared with Gutierrez-Villanueva et al. (2023, GV23).
Two coordinate systems, as in GV23:
 (1) passage coordinates: 25-km distance bins from the north (GV23 Figs. 3-4);
 (2) streamwise coordinates: pairs of sea-surface-height streamlines (GV23 Figs. 5, 7), with GV23's own bin edges
     (psi = -0.2 ... -1.7 m; SAF -0.40..-0.60, PF -1.00..-1.30, SACCF -1.55..-1.65 m). GV23's SSH = Niiler-Maximenko MDT
     - 0.60 m + SLA, whereas our ADT uses another MDT: the offset is estimated empirically by sampling our monthly ADT at
     the positions where GV23 found each streamline bin on each transect (stream.lon/lat), and a linear fit of our ADT on
     psic gives the mapping of their bin edges onto our ADT.
Altimetry (each product): cross-track surface geostrophic velocity u_n (eastward-positive normal) and EKE
(0.5|U_m - Ubar|^2 + 0.5 var_intra) at 121 points along the mean track; points with ice fraction > 0.10 in a month are
set to NaN for that month. Per bin: monthly mean of u_n, of EKE, and (streamwise) of the surface transport per unit
depth (sum of u_n * ds over the points falling in the bin, m2 s-1).
In situ (GV23 files): (1) velocity at 70-118 m per 25-km bin (lmgvelos38nb_50_300.mat, u) and 0-760 m transport per bin
(TotalTransport); (2) streamwise 0-760 m transport and depth-integrated EKE per streamline pair (Streamwise_binning_*.mat).
Trends: Theil-Sen + calibrated test (rstats.trend_bcz). 2005-2019: annual means (Oct 2005-Apr 2019 data, years with
>= 3 transects / >= 6 months), AR(1) test, for altimetry and in situ alike; 2001-2025: monthly altimetric series
(monthly test) where complete, else annual means.
Usage: python dk5_transect_vs_gv23.py <revision dir> <GV repo dir> <out dir>"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, scipy.io as sio
warnings.filterwarnings("ignore")
REV, GV, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
R_E = 6.371e6
L = lambda f, k: sio.loadmat(GV / "Datasets" / f, squeeze_me=True, struct_as_record=False)[k]
TT = L("TotalTransport_os38nbmaxz_760_2.mat", "totaltransp")
blat = np.array(TT.lat, float); blon = np.array(TT.lon, float) - 360; dist = np.array(TT.dist, float)
keep = dist <= 840
mlat = np.nanmean(blat[keep], 1); mlon = np.nanmean(blon[keep], 1)
s = np.r_[0, np.cumsum(np.hypot(np.diff(mlon) * np.cos(np.radians(mlat[:-1])), np.diff(mlat)))]
ss = np.linspace(0, s[-1], 121); lo = np.interp(ss, s, mlon); la = np.interp(ss, s, mlat)
dkm = np.interp(ss, s, dist[keep])                              # distance from north (km), GV23 convention
latm = np.radians(0.5 * (la[:-1] + la[1:]))
seg = np.hypot(R_E * np.cos(latm) * np.radians(np.diff(lo)), R_E * np.radians(np.diff(la)))
dsx = np.r_[seg[0] / 2, 0.5 * (seg[:-1] + seg[1:]), seg[-1] / 2]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
YRS = range(2001, 2026)
T = pd.DataFrame({"year": np.repeat(list(YRS), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - 0.5) / 12

ST = L("Streamwise_binning_os38nb.mat", "stream"); SE = L("Streamwise_binning_eke.mat", "stream")
psi = np.array(ST.psi, float); psic = np.array(ST.psic, float)
st_time = pd.to_datetime(np.array(ST.time, float) - 719529, unit="D")

# ------------------------------------------------------------------ altimetry along the section
cache = OUT / "dk5_section_fields.npz"
if cache.exists():
    Z = np.load(cache); A = {k: Z[k] for k in Z.files}
else:
    A = {}
    for p in ("allsat", "twosat"):
        UN, EK, AD, IC = [], [], [], []
        ub = vb = None
        rows = []
        for y in YRS:
            for mo in range(1, 13):
                with xr.open_dataset(REV / "data" / f"monthly_{p}" / f"{y}" / f"{mo:02d}.nc") as d:
                    dd = d.sel(longitude=slice(-70, -55), latitude=slice(-65, -53)).load()
                rows.append((dd.u_mean.interp(**pts).values, dd.v_mean.interp(**pts).values, dd.var_intra.interp(**pts).values,
                             dd.adt_mean.interp(**pts).values, np.nan_to_num(dd.ice_frac.interp(**pts).values)))
        u = np.array([r[0] for r in rows]); v = np.array([r[1] for r in rows]); vi = np.array([r[2] for r in rows])
        A[f"{p}_un"] = u * nx + v * ny
        A[f"{p}_eke"] = 0.5 * ((u - np.nanmean(u, 0)) ** 2 + (v - np.nanmean(v, 0)) ** 2) + 0.5 * vi
        A[f"{p}_adt"] = np.array([r[3] for r in rows]); A[f"{p}_ice"] = np.array([r[4] for r in rows])
        print(p, "read", flush=True)
    np.savez_compressed(cache, **A)
if "adt_at_gvbins" not in A:
    # ADT at the GV23 streamline-bin positions (for the psi -> ADT mapping), all-satellite monthly fields;
    # each cell holds the positions of the profiles in that streamline pair: use their mean position
    def cellmean(c):
        a = np.atleast_1d(np.asarray(c, float)); return np.nanmean(a) if a.size else np.nan
    slat = np.vectorize(cellmean)(ST.lat); slon = np.vectorize(cellmean)(ST.lon); slon = np.where(slon > 180, slon - 360, slon)
    adt_at = np.full(slat.shape, np.nan)
    for k, tt in enumerate(st_time):
        if not (2001 <= tt.year <= 2025): continue
        with xr.open_dataset(REV / "data" / "monthly_allsat" / f"{tt.year}" / f"{tt.month:02d}.nc") as d:
            dd = d.sel(longitude=slice(-70, -55), latitude=slice(-65, -53)).adt_mean.load()
        ok = np.isfinite(slat[k]) & np.isfinite(slon[k])
        if ok.any():
            adt_at[k, ok] = dd.interp(longitude=xr.DataArray(slon[k, ok], dims="q"), latitude=xr.DataArray(slat[k, ok], dims="q")).values
    A["adt_at_gvbins"] = adt_at
    np.savez_compressed(cache, **A)

# psi -> ADT mapping (robust linear fit over bins with data)
y_ = A["adt_at_gvbins"]; X = np.broadcast_to(psic[None, :], y_.shape)
ok = np.isfinite(y_)
b1, b0 = np.polyfit(X[ok], y_[ok], 1)
off_med = np.nanmedian(y_ - X)
edges_adt = b0 + b1 * psi
pd.Series(dict(slope=b1, intercept=b0, median_offset=off_med, r=np.corrcoef(X[ok], y_[ok])[0, 1], n=int(ok.sum()))).to_csv(OUT / "dk5_psi_to_adt.csv")
print("psi->ADT: slope %.3f intercept %.3f median offset %.3f" % (b1, b0, off_med), flush=True)


def trend_row(x, yv, monthly, month=None):
    ok = np.isfinite(yv)
    if ok.sum() < (60 if monthly else 8): return dict(slope=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan, n=int(ok.sum()))
    o = R.trend_bcz(x[ok], yv[ok], month=None if month is None else month[ok], monthly=monthly)
    return dict(slope=o["slope"], ci_lo=o["ci_lo"], ci_hi=o["ci_hi"], p=o["p"], n=int(ok.sum()))


def annual(years, vals, minn):
    df = pd.DataFrame(dict(year=years, v=vals)).dropna()
    g = df.groupby("year").v; a = g.mean()[g.count() >= minn]
    return a.index.values.astype(float), a.values


ov = (T.x >= 2005 + 9 / 12) & (T.x <= 2019 + 4 / 12)            # Oct 2005 - Apr 2019
rows = []
# ------------------------------------------------------------------ (1) passage coordinates, 25-km bins
DB = np.arange(0, 850, 25)
for p in ("allsat", "twosat"):
    un = np.where(A[f"{p}_ice"] > 0.10, np.nan, A[f"{p}_un"]); ek = np.where(A[f"{p}_ice"] > 0.10, np.nan, A[f"{p}_eke"])
    for b0_ in DB:
        sel = (dkm >= b0_) & (dkm < b0_ + 25)
        if not sel.any(): continue
        for var, F, scale in (("u", un, 100.0), ("eke", ek, 1e4)):
            y = np.nanmean(F[:, sel], 1) * scale
            xa, ya = annual(T.year[ov].values, y[ov.values], 3)
            r1 = trend_row(xa, ya, False)
            if np.isfinite(y).all(): r2 = trend_row(T.x.values, y, True, T.month.values)
            else:
                xa2, ya2 = annual(T.year.values, y, 6); r2 = trend_row(xa2, ya2, False)
            for per, r in (("2005-2019", r1), ("2001-2025", r2)):
                rows.append(dict(coord="passage", src=p, var=var, bin=b0_ + 12.5, period=per, mean=np.nanmean(y), ice_months=int(np.isnan(y).sum()), **r))
# in situ, passage coordinates: velocity 70-118 m and transport per bin
LV = L("lmgvelos38nb_50_300.mat", "lmgdata")
z = np.array(LV.z, float); zs = (z >= 70) & (z <= 118)
U = np.array(LV.u, float)                                    # (z, bin, transect), m/s; gaps filled with OM
vt = pd.to_datetime(np.array(LV.time, float) - 719529, unit="D")
vd = np.array(LV.dist, float)                                 # (bin, transect)
tr = np.array(TT.transp, float)
tt_ = pd.to_datetime(np.array(TT.time, float) - 719529, unit="D")
for b0_ in DB:
    c = b0_ + 12.5
    uu = np.full(U.shape[2], np.nan)
    for k in range(U.shape[2]):
        j = np.where(np.abs(vd[:, k] - c) < 12.6)[0]
        if j.size: uu[k] = np.nanmean(U[zs][:, j[0], k]) * 100
    j2 = np.where(np.abs(dist - c) < 12.6)[0]
    q = tr[j2[0], :] / 1e6 if j2.size else np.full(tr.shape[1], np.nan)
    for var, yv, yrs in (("u", uu, vt.year), ("transport", q, tt_.year)):
        inov = (pd.Series(yrs) >= 2005).values & (pd.Series(yrs) <= 2019).values
        xa, ya = annual(np.array(yrs)[inov], yv[inov], 3)
        rows.append(dict(coord="passage", src="insitu", var=var, bin=c, period="2005-2019", mean=np.nanmean(yv), **trend_row(xa, ya, False)))

# ------------------------------------------------------------------ (2) streamwise coordinates
for p in ("allsat", "twosat"):
    adt = A[f"{p}_adt"]; ice = A[f"{p}_ice"] > 0.10
    un = np.where(ice, np.nan, A[f"{p}_un"]); ek = np.where(ice, np.nan, A[f"{p}_eke"])
    for j in range(len(psic)):
        lo_e, hi_e = sorted((edges_adt[j], edges_adt[j + 1]))
        inb = (adt >= lo_e) & (adt < hi_e)
        icy = (inb & ice).any(1)                                 # month in which part of this streamline band is ice-covered
        q = np.where(icy, np.nan, np.nansum(np.where(inb, A[f"{p}_un"] * dsx[None, :], 0), 1))   # m2 s-1
        uu = np.array([np.nanmean(un[k, inb[k]]) * 100 if inb[k].any() else np.nan for k in range(len(adt))])
        ee = np.array([np.nanmean(ek[k, inb[k]]) * 1e4 if inb[k].any() else np.nan for k in range(len(adt))])
        for var, y in (("transport", q), ("u", uu), ("eke", ee)):
            xa, ya = annual(T.year[ov].values, y[ov.values], 3); r1 = trend_row(xa, ya, False)
            xa2, ya2 = annual(T.year.values, y, 6); r2 = trend_row(xa2, ya2, False)
            for per, r in (("2005-2019", r1), ("2001-2025", r2)):
                rows.append(dict(coord="stream", src=p, var=var, bin=psic[j], period=per, mean=np.nanmean(y), ice_months=int(np.isnan(y).sum()), **r))
# in situ streamwise: transport (0-760 m) and depth-integrated EKE per streamline pair
for var, arr, scale, tm in (("transport", np.array(ST.transp, float), 1e-6, st_time),
                            ("eke", np.array(SE.eke, float), 1.0, pd.to_datetime(np.array(SE.time, float) - 719529, unit="D"))):
    for j in range(len(psic)):
        y = arr[:, j] * scale
        xa, ya = annual(tm.year.values, y, 3)
        sel = (xa >= 2005) & (xa <= 2019)
        rows.append(dict(coord="stream", src="insitu", var=var, bin=psic[j], period="2005-2019", mean=np.nanmean(y), **trend_row(xa[sel], ya[sel], False)))

D = pd.DataFrame(rows)
D.to_csv(OUT / "dk5_trends.csv", index=False)
pd.DataFrame(dict(psi=psi, adt_edge=edges_adt)).to_csv(OUT / "dk5_stream_edges.csv", index=False)
print(D.groupby(["coord", "src", "var", "period"]).apply(lambda g: f"{(g.p < 0.05).sum()} sig / {g.p.notna().sum()}").to_string())
