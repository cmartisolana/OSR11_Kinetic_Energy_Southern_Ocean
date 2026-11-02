"""dk25: regional EKE and the wind stress that precedes it (after Hogg et al. 2022, Sect. 4: regional EKE correlated with
coarse-grained wind stress at a range of lags, to separate local from upstream/remote forcing). Product 1 EKE, Product 2
eastward wind stress, 1 deg fields of dk22, annual means 2001-2025, both linearly detrended.
Regions (ACC band, mean ADT -0.6..+0.2 m): Kerguelen Plateau 65-95E, Southeast Indian Ridge 130-150E, Macquarie Ridge
155-175E, Pacific-Antarctic Ridge 145-100W, southeast Pacific 100-70W (upstream of Drake Passage), Drake Passage-Scotia
Sea 70-35W (the ridge/downstream ranges of dk13 for the first four; the southeast Pacific and Drake ranges as in dk1).
(1) Maps: correlation of regional EKE with the wind stress on 3 deg cells at lags 0, 1, 2 yr; p from a t test with effective
    sample size (lag-1 autocorrelations).
(2) Lag curves, lags -2..+4 yr: regional EKE against the wind stress averaged over the ACC band of the same longitudes
    ('local'), and of the 30 deg immediately upstream (west; 'upstream').
Usage: python dk25_regional_wind.py <revision dir> <dk22 out dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, D22, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
Z = np.load(D22 / "dk22_fields.npz"); EK, TAU, lat, lon = Z["EK"], Z["TAU"], Z["lat"], Z["lon"]
m1 = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.coarsen(latitude=8, longitude=8, boundary="trim").mean().values
accE = (m1 >= -0.6) & (m1 <= 0.2) & np.isfinite(EK).all(0); accT = (m1 >= -0.6) & (m1 <= 0.2) & np.isfinite(TAU).all(0)
w = np.cos(np.radians(lat))[:, None] * np.ones(lon.size)[None, :]
ann = lambda a: a.reshape(25, 12, *a.shape[1:]).mean(1)
t = np.arange(25) - 12.0


def detr(a):
    a = np.asarray(a, float); tt = t.reshape((-1,) + (1,) * (a.ndim - 1)); m = a.mean(0)
    return a - m - tt * ((tt * (a - m)).sum(0) / (t ** 2).sum())


def lagc(x, y, L):
    xa, ya = (x[:25 - L], y[L:]) if L >= 0 else (x[-L:], y[:25 + L])
    xa = xa - xa.mean(0); ya = ya - ya.mean(0); n = len(ya)
    r = (xa * ya).sum(0) / np.sqrt((xa ** 2).sum(0) * (ya ** 2).sum(0))
    a1 = (xa[1:] * xa[:-1]).sum(0) / (xa ** 2).sum(0); b1 = (ya[1:] * ya[:-1]).sum(0) / (ya ** 2).sum(0)
    ne = np.clip(n * (1 - a1 * b1) / (1 + a1 * b1), 3, n)
    p = 2 * stats.t.sf(np.abs(r) * np.sqrt((ne - 2) / np.maximum(1 - r ** 2, 1e-9)), ne - 2)
    return r, p


def lonmask(a, b):
    return ((lon >= a) & (lon < b)) if a < b else ((lon >= a) | (lon < b))


REG = {"Kerguelen Plateau": (65, 95), "Southeast Indian Ridge": (130, 150), "Macquarie Ridge": (155, 175),
       "Pacific–Antarctic Ridge": (-145, -100), "Southeast Pacific": (-100, -70), "Drake Passage–Scotia Sea": (-70, -35)}
EA, TA = ann(EK), ann(TAU)
T3 = detr(TA.reshape(25, lat.size // 3, 3, lon.size // 3, 3).mean(axis=(2, 4)))
out = dict(lat3=lat.reshape(-1, 3).mean(1), lon3=lon.reshape(-1, 3).mean(1), regnames=np.array(list(REG)))
rows = []
for nm, (a_, b_) in REG.items():
    mE = accE & lonmask(a_, b_)[None, :]; e = detr(np.nansum(EA[:, mE] * w[mE], 1) / w[mE].sum())
    for L in (0, 1, 2):
        r, p = lagc(T3, e[:, None, None] * np.ones(T3.shape), L); out[f"r|{nm}|{L}"] = r; out[f"p|{nm}|{L}"] = p
    ua, ub = (a_ - 30 + 180) % 360 - 180, a_
    for kind, (c0, c1) in (("local", (a_, b_)), ("upstream", (ua, ub))):
        mT = accT & lonmask(c0, c1)[None, :]; x = detr(np.nansum(TA[:, mT] * w[mT], 1) / w[mT].sum())
        for L in range(-2, 5):
            r, p = lagc(x, e, L); rows.append(dict(region=nm, wind=kind, lag=L, r=float(r), p=float(p)))
np.savez(OUT / "dk25_maps.npz", **out)
D = pd.DataFrame(rows); D.to_csv(OUT / "dk25_lagcurves.csv", index=False)
print(D.pivot_table(index=["region", "wind"], columns="lag", values="r").round(2).to_string())
for nm in REG:
    for L in (0, 1, 2):
        r, p = out[f"r|{nm}|{L}"], out[f"p|{nm}|{L}"]; ok = np.isfinite(r)
        print(nm, L, "frac p<0.05", np.round(np.mean(p[ok] < 0.05), 3), "max r", np.round(np.nanmax(r), 2))
# (3) control: one-point correlation of the wind stress (3 deg, detrended annual) with the southeast Pacific band wind,
#     to show how far wind anomalies are coherent (remote EKE-wind correlations can reflect this coherence)
mT = accT & lonmask(-100, -70)[None, :]; xsp = detr(np.nansum(TA[:, mT] * w[mT], 1) / w[mT].sum())
rw, pw = lagc(xsp[:, None, None] * np.ones(T3.shape), T3, 0)
np.savez(OUT / "dk25_windcoherence.npz", r=rw, p=pw, lat3=out["lat3"], lon3=out["lon3"])
print("wind coherence with SE Pacific wind: frac |r|>0.5", np.round(np.nanmean(np.abs(rw) > 0.5), 3))
