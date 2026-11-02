"""DK24: EKE, KE and surface transport per streamline pair along the LM Gould transect (Sect. 2.5 method, as dk5),
adding KE. Altimetry Product 1 at the 121 points of the mean track (d1b_track.csv, same construction as dk5), monthly
fields; points with ice fraction > 0.10 masked per month; streamline pairs of GV23 mapped onto the monthly ADT with the
dk5 edges (dk5_stream_edges.csv). Per pair and month: mean EKE and KE (cm2 s-2) over the points in the pair; surface
transport per unit depth = sum u_n ds (m2 s-1), NaN if any point of the pair is ice-covered.
Trends (rstats.trend_bcz): 2001-2025 annual means (>= 6 months), 2005-2019 annual means of Oct 2005-Apr 2019 (>= 3 months).
Usage: python dk24_section_streamwise.py <revision dir> <monthly dir> <track csv> <dk5 dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, MON, TRK, DK5, OUT = (Path(a) for a in sys.argv[1:6]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
R_E = 6.371e6
tk = pd.read_csv(TRK); lo, la = tk.lon.values, tk.lat.values
latm = np.radians(0.5 * (la[:-1] + la[1:]))
seg = np.hypot(R_E * np.cos(latm) * np.radians(np.diff(lo)), R_E * np.radians(np.diff(la)))
dsx = np.r_[seg[0] / 2, 0.5 * (seg[:-1] + seg[1:]), seg[-1] / 2]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - .5) / 12
cache = OUT / "dk24_section_fields_allsat.npz"
if cache.exists():
    Z = dict(np.load(cache))
else:
    rows = []
    for y, mo in zip(T.year, T.month):
        with xr.open_dataset(MON / f"{y}" / f"{mo:02d}.nc") as d:
            dd = d.sel(longitude=slice(-70, -55), latitude=slice(-65, -53)).load()
        rows.append([dd[v].interp(**pts).values for v in ("u_mean", "v_mean", "var_intra", "adt_mean", "ice_frac")])
    A = np.array(rows)                                    # (300, 5, 121)
    Z = dict(u=A[:, 0], v=A[:, 1], vi=A[:, 2], adt=A[:, 3], ice=np.nan_to_num(A[:, 4]), lon=lo, lat=la, dsx=dsx)
    np.savez_compressed(cache, **Z)
u, v, vi, adt, ice = Z["u"], Z["v"], Z["vi"], Z["adt"], Z["ice"] > 0.10
un = u * nx + v * ny
EKE = (0.5 * ((u - np.nanmean(u, 0)) ** 2 + (v - np.nanmean(v, 0)) ** 2) + 0.5 * vi) * 1e4
KE = 0.5 * (u ** 2 + v ** 2) * 1e4
EKE[ice] = np.nan; KE[ice] = np.nan
E = pd.read_csv(DK5 / "dk5_stream_edges.csv"); psi = E.psi.values; edg = E.adt_edge.values
ov = ((T.x >= 2005 + 9 / 12) & (T.x <= 2019 + 4 / 12)).values
def annual(years, vals, minn):
    df = pd.DataFrame(dict(year=years, v=vals)).dropna(); g = df.groupby("year").v; a = g.mean()[g.count() >= minn]
    return a.index.values.astype(float), a.values
def trow(x, y):
    ok = np.isfinite(y)
    if ok.sum() < 8: return dict(slope=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan, n=int(ok.sum()))
    o = R.trend_bcz(x[ok], y[ok], monthly=False); return dict(slope=o["slope"], ci_lo=o["ci_lo"], ci_hi=o["ci_hi"], p=o["p"], n=int(ok.sum()))
rows, ser = [], {}
for j in range(len(psi) - 1):
    a, b = sorted((edg[j], edg[j + 1])); inb = (adt >= a) & (adt < b); icy = (inb & ice).any(1)
    q = np.where(icy, np.nan, np.nansum(np.where(inb, un * dsx[None, :], 0), 1)); q[~inb.any(1)] = np.nan
    e = np.array([np.nanmean(EKE[k, inb[k]]) if inb[k].any() else np.nan for k in range(300)])
    kk = np.array([np.nanmean(KE[k, inb[k]]) if inb[k].any() else np.nan for k in range(300)])
    psic = 0.5 * (psi[j] + psi[j + 1]); adtc = 0.5 * (edg[j] + edg[j + 1])
    for var, y in (("transport", q), ("eke", e), ("ke", kk)):
        ser[f"{var}|{psic:.3f}"] = y
        xa, ya = annual(T.year[ov].values, y[ov], 3); r1 = trow(xa, ya)
        xa2, ya2 = annual(T.year.values, y, 6); r2 = trow(xa2, ya2)
        for per, r in (("2005-2019", r1), ("2001-2025", r2)):
            rows.append(dict(src="allsat", var=var, psic=psic, adt_c=adtc, adt_lo=a, adt_hi=b, period=per, mean=np.nanmean(y),
                             months=int(np.isfinite(y).sum()), **r))
D = pd.DataFrame(rows); D["pct"] = 1000 * D.slope / D["mean"]; D["pct_lo"] = 1000 * D.ci_lo / D["mean"]; D["pct_hi"] = 1000 * D.ci_hi / D["mean"]
D.to_csv(OUT / "dk24_stream_trends.csv", index=False)
pd.DataFrame(ser).assign(year=T.year, month=T.month).to_csv(OUT / "dk24_stream_series.csv", index=False)
print(D[D.period == "2001-2025"].round(3).to_string())
