"""dk26: is the southeast Pacific wind pattern of dk25 (Fig. 4d) the SAM and/or the ENSO (Pacific-South American) pattern?
Indices: SAM (Marshall station index, dk17/s3 'SAM'); ENSO = NOAA CPC Oceanic Nino Index (Nino 3.4, 3-month running mean,
ERSSTv5; file oni_cpc_2001_2025.csv); SE Pacific wind = eastward wind stress over the ACC band at 100-70W (dk25).
All annual means 2001-2025 (also austral winter JJA and summer DJF for the indices), linearly detrended.
(1) Correlations: SE Pacific wind vs SAM and vs ONI; multiple regression of the SE Pacific wind on SAM and ONI (standardised
    coefficients, R^2). p from t tests with effective sample size (lag-1 autocorrelations).
(2) Maps (3 deg): correlation of the wind stress with SAM and with ONI (same year); spatial (pattern) correlation of each with
    the SE Pacific one-point map of dk25 over ocean cells south of 35 S.
(3) Regional EKE (dk25 regions) vs SAM and ONI at lags 0-2 yr.
Usage: python dk26_wind_modes.py <revision dir> <drake out dir> <dk22 out dir> <dk25 out dir> <oni csv> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, DO, D22, D25, ONI, OUT = (Path(a) for a in sys.argv[1:7]); OUT.mkdir(parents=True, exist_ok=True)
Z = np.load(D22 / "dk22_fields.npz"); EK, TAU, lat, lon = Z["EK"], Z["TAU"], Z["lat"], Z["lon"]
m1 = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.coarsen(latitude=8, longitude=8, boundary="trim").mean().values
acc = (m1 >= -0.6) & (m1 <= 0.2)
w = np.cos(np.radians(lat))[:, None] * np.ones(lon.size)[None, :]
t = np.arange(25) - 12.0
def detr(a):
    a = np.asarray(a, float); tt = t.reshape((-1,) + (1,) * (a.ndim - 1)); m = a.mean(0)
    return a - m - tt * ((tt * (a - m)).sum(0) / (t ** 2).sum())
def corr(x, y, L=0):
    xa, ya = (x[:25 - L], y[L:]) if L >= 0 else (x[-L:], y[:25 + L])
    xa = xa - xa.mean(0); ya = ya - ya.mean(0); n = len(xa)
    r = (xa * ya).sum(0) / np.sqrt((xa ** 2).sum(0) * (ya ** 2).sum(0))
    a1 = (xa[1:] * xa[:-1]).sum(0) / (xa ** 2).sum(0); b1 = (ya[1:] * ya[:-1]).sum(0) / (ya ** 2).sum(0)
    ne = np.clip(n * (1 - a1 * b1) / (1 + a1 * b1), 3, n)
    return r, 2 * stats.t.sf(np.abs(r) * np.sqrt((ne - 2) / np.maximum(1 - r ** 2, 1e-9)), ne - 2)
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); O = pd.read_csv(ONI, comment="#")
A = A.merge(O, on=["year", "month"], how="left")
TAUm = TAU  # monthly 1 deg
mT = acc & np.isfinite(TAUm).all(0) & ((lon >= -100) & (lon < -70))[None, :]
A["tau_sep"] = np.nansum(TAUm[:, mT] * w[mT], 1) / w[mT].sum()
res = []
for season, months in (("annual", range(1, 13)), ("JJA", (6, 7, 8)), ("DJF", (12, 1, 2))):
    B = A[A.month.isin(months)].copy()
    if season == "DJF": B["year"] = np.where(B.month == 12, B.year + 1, B.year); B = B[(B.year >= 2002) & (B.year <= 2025)]
    g = B.groupby("year")[["tau_sep", "SAM", "oni"]].mean(); n = len(g)
    tt = np.arange(n) - (n - 1) / 2
    d = {c: g[c].values - np.polyval(np.polyfit(tt, g[c].values, 1), tt) for c in g}
    for ix in ("SAM", "oni"):
        r = np.corrcoef(d[ix], d["tau_sep"])[0, 1]; p = stats.pearsonr(d[ix], d["tau_sep"]).pvalue
        res.append(dict(season=season, pair=f"SE Pacific wind vs {ix}", r=r, p_naive=p, n=n))
    Xs = np.column_stack([(d["SAM"] - d["SAM"].mean()) / d["SAM"].std(), (d["oni"] - d["oni"].mean()) / d["oni"].std()])
    y = (d["tau_sep"] - d["tau_sep"].mean()) / d["tau_sep"].std(); b, *_ = np.linalg.lstsq(Xs, y, rcond=None)
    R2 = 1 - ((y - Xs @ b) ** 2).sum() / (y ** 2).sum()
    res.append(dict(season=season, pair="regression on SAM, ONI (std coeffs)", r=np.nan, p_naive=np.nan, n=n, b_SAM=b[0], b_ONI=b[1], R2=R2))
    res.append(dict(season=season, pair="SAM vs ONI", r=np.corrcoef(d["SAM"], d["oni"])[0, 1], p_naive=stats.pearsonr(d["SAM"], d["oni"]).pvalue, n=n))
R = pd.DataFrame(res); R.to_csv(OUT / "dk26_index_corr.csv", index=False); print(R.round(3).to_string(index=False))
an = A.groupby("year").mean(numeric_only=True); sam, oni = detr(an.SAM.values), detr(an.oni.values)
TA = TAU.reshape(25, 12, *TAU.shape[1:]).mean(1); T3 = detr(TA.reshape(25, lat.size // 3, 3, lon.size // 3, 3).mean(axis=(2, 4)))
WC = np.load(D25 / "dk25_windcoherence.npz"); out = dict(lat3=WC["lat3"], lon3=WC["lon3"])
for nm, x in (("sam", sam), ("oni", oni)):
    r, p = corr(x[:, None, None] * np.ones(T3.shape), T3); out[f"r_{nm}"] = r; out[f"p_{nm}"] = p
    ok = np.isfinite(r) & np.isfinite(WC["r"]); out[f"pattern_r_{nm}"] = np.corrcoef(r[ok], WC["r"][ok])[0, 1]
    print(f"pattern correlation of wind-{nm} map with SE Pacific one-point map: {out[f'pattern_r_{nm}']:.2f}")
np.savez(OUT / "dk26_maps.npz", **out)
# (3) regional EKE vs SAM, ONI
REG = {"Kerguelen Plateau": (65, 95), "Southeast Indian Ridge": (130, 150), "Macquarie Ridge": (155, 175),
       "Pacific–Antarctic Ridge": (-145, -100), "Southeast Pacific": (-100, -70), "Drake Passage–Scotia Sea": (-70, -35)}
accE = acc & np.isfinite(EK).all(0); EA = EK.reshape(25, 12, *EK.shape[1:]).mean(1); rows = []
for nm, (a_, b_) in REG.items():
    mE = accE & ((lon >= a_) & (lon < b_))[None, :]; e = detr(np.nansum(EA[:, mE] * w[mE], 1) / w[mE].sum())
    for ix, x in (("SAM", sam), ("ONI", oni), ("SE Pacific wind", detr(an.tau_sep.values))):
        for L in range(-1, 4):
            r, p = corr(x, e, L); rows.append(dict(region=nm, index=ix, lag=L, r=float(r), p=float(p)))
E = pd.DataFrame(rows); E.to_csv(OUT / "dk26_regional_eke.csv", index=False)
pd.set_option("display.width", 200); print(E.pivot_table(index=["index", "region"], columns="lag", values="r").round(2).to_string())
