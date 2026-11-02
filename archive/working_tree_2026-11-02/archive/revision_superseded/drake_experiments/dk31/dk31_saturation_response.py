"""dk31b: response of the eddies and of the mean flow to the past wind stress (eddy-saturation test), 2001-2025.
Annual EKE (dk22), MKE (dk31a) and eastward wind stress (dk22), ACC band (mean ADT -0.6..0.2 m), each as % of its
2001-2025 mean, linearly detrended (interannual response only). Past-wind predictors:
  lag L   : wind of year y-L (L = 0..4)
  window N: wind averaged over years y-N+1..y (N = 1..6; N = 1 is the same year)
Regression slope = sensitivity (% EKE or MKE per % wind stress); r and p with effective sample size (Bretherton 1999).
Done for the ACC band and sectors (area-weighted regional means of each field) and per 3 deg cell (local wind).
Usage: python dk31_saturation_response.py <revision dir> <dk22 out> <dk31 out>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, D22, OUT = (Path(a) for a in sys.argv[1:4])
Z = np.load(D22 / "dk22_fields.npz"); EK, TAU, lat, lon = Z["EK"], Z["TAU"], Z["lat"], Z["lon"]
MK = np.load(OUT / "dk31_mke_fields.npz")["MKE"]
ann = lambda a: a.reshape(25, 12, *a.shape[1:]).mean(1)
EA, TA = ann(EK), ann(TAU)
m1 = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.coarsen(latitude=8, longitude=8, boundary="trim").mean().values
acc = (m1 >= -0.6) & (m1 <= 0.2) & np.isfinite(EA).all(0) & np.isfinite(MK).all(0) & np.isfinite(TA).all(0)
w = np.cos(np.radians(lat))[:, None] * np.ones(lon.size)[None, :]
t = np.arange(25) - 12.0
def detr(y):
    tt = t.reshape((-1,) + (1,) * (y.ndim - 1)); m = y.mean(0); b = (tt * (y - m)).sum(0) / (t ** 2).sum(); return y - m - tt * b
def pct(a): return 100 * (a / a.mean(0) - 1)
def pred(x, kind, k):          # returns (x_aligned, slice of years for y)
    if kind == "lag": return x[:25 - k], slice(k, 25)
    xs = np.array([x[i - k + 1:i + 1].mean(0) for i in range(k - 1, 25)]); return xs, slice(k - 1, 25)
def reg(x, y):
    x = x - x.mean(0); y = y - y.mean(0); n = x.shape[0]
    b = (x * y).sum(0) / (x * x).sum(0); r = (x * y).sum(0) / np.sqrt((x * x).sum(0) * (y * y).sum(0))
    a1 = (x[1:] * x[:-1]).sum(0) / (x * x).sum(0); a2 = (y[1:] * y[:-1]).sum(0) / (y * y).sum(0)
    ne = np.clip(n * (1 - a1 * a2) / (1 + a1 * a2), 3, n); res = y - b * x
    se = np.sqrt((res ** 2).sum(0) / (ne - 2) / (x * x).sum(0))
    p = 2 * stats.t.sf(np.abs(r) * np.sqrt((ne - 2) / np.maximum(1 - r ** 2, 1e-9)), ne - 2)
    return b, se, r, p, ne
SECT = {"ACC band": None, "Pacific": (147, 290), "Indian": (20, 147), "Atlantic": (290, 380)}
rows, series = [], {}
for sec, rg in SECT.items():
    l = np.mod(lon, 360); lm = np.ones(lon.size, bool) if rg is None else (((l >= rg[0]) & (l < rg[1])) | ((l + 360 >= rg[0]) & (l + 360 < rg[1])))
    m = acc & lm[None, :]
    reg_mean = lambda F: (F[:, m] * w[m]).sum(1) / w[m].sum()
    e, k_, tau = (detr(pct(reg_mean(F))) for F in (EA, MK, TA))
    series[sec] = dict(eke=e, mke=k_, tau=tau, eke_raw=pct(reg_mean(EA)), mke_raw=pct(reg_mean(MK)), tau_raw=pct(reg_mean(TA)))
    for kind, ks in (("lag", range(0, 5)), ("window", range(1, 7))):
        for k in ks:
            x, sl = pred(tau, kind, k)
            for nm, y in (("EKE", e), ("MKE", k_)):
                b, se, r, p, ne = reg(x, y[sl])
                rows.append(dict(region=sec, kind=kind, k=k, var=nm, sens=b, lo=b - 1.96 * se, hi=b + 1.96 * se, r=r, p=p, neff=ne))
D = pd.DataFrame(rows); D.to_csv(OUT / "dk31_response.csv", index=False)
pd.set_option("display.width", 200)
print(D.pivot_table(index=["region", "kind", "k"], columns="var", values=["sens", "r", "p"]).round(2).to_string())
# maps, 3 deg cells, local wind
def cg(a): return a.reshape(a.shape[0], 10, 3, 120, 3).mean(axis=(2, 4))
m3 = acc.reshape(10, 3, 120, 3).mean(axis=(1, 3)) > 0.5
E3, K3, T3 = (detr(pct(cg(np.where(acc[None], F, np.nan)))) for F in (EA, MK, TA))
mp = dict(lat3=lat.reshape(-1, 3).mean(1), lon3=lon.reshape(-1, 3).mean(1), mask=m3)
for kind, k in (("lag", 0), ("lag", 1), ("window", 2), ("window", 3)):
    x, sl = pred(T3, kind, k)
    for nm, Y in (("EKE", E3), ("MKE", K3)):
        b, se, r, p, ne = reg(x, Y[sl])
        for q, v in (("sens", b), ("r", r), ("p", p)): mp[f"{q}_{nm}_{kind}{k}"] = np.where(m3, v, np.nan)
        ok = m3 & np.isfinite(p); print(f"map {nm} {kind}{k}: median sens {np.nanmedian(b[ok]):+.2f}, median r {np.nanmedian(r[ok]):+.2f}, cells p<0.05 {100*np.mean(p[ok]<0.05):.0f} % (pos {100*np.mean((p[ok]<0.05)&(r[ok]>0)):.0f} %)")
np.savez(OUT / "dk31_maps.npz", **mp); np.savez(OUT / "dk31_series.npz", **{f"{s}|{k}": v for s, d in series.items() for k, v in d.items()})
