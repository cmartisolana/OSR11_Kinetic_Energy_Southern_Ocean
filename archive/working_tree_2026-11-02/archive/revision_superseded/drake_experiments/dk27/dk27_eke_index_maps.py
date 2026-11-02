"""dk27: correlation maps of EKE with the local eastward wind stress, the SAM and Nino 3.4 (ONI), 3 deg cells, annual means
2001-2025, both detrended, at lags 0 and 1 yr (EKE lagging). p: t test with effective sample size; FDR (BH, 10 %) count.
Inputs: dk22 fields (EKE, wind stress at 1 deg), dk17/s3 (SAM), ONI csv (NOAA CPC). Same conventions as dk22/dk23.
Usage: python dk27_eke_index_maps.py <drake out dir> <dk22 out dir> <oni csv> <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
DO, D22, ONI, REV, OUT = (Path(a) for a in sys.argv[1:6]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
Z = np.load(D22 / "dk22_fields.npz"); EK, TAU, lat, lon = Z["EK"], Z["TAU"], Z["lat"], Z["lon"]
t = np.arange(25) - 12.0
def detr(a):
    a = np.asarray(a, float); tt = t.reshape((-1,) + (1,) * (a.ndim - 1)); m = a.mean(0)
    return a - m - tt * ((tt * (a - m)).sum(0) / (t ** 2).sum())
def cg(a): return a.reshape(25, 12, 10, 3, 120, 3).mean(axis=(1, 3, 5))
E3, T3 = detr(cg(EK)), detr(cg(TAU))
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); O = pd.read_csv(ONI, comment="#"); A = A.merge(O, on=["year", "month"], how="left")
an = A.groupby("year").mean(numeric_only=True); sam, oni = detr(an.SAM.values), detr(an.oni.values)
def corr(x, y, L):
    xa, ya = (x[:25 - L], y[L:]) if L else (x, y)
    xa = xa - xa.mean(0); ya = ya - ya.mean(0); n = len(ya)
    r = (xa * ya).sum(0) / np.sqrt((xa ** 2).sum(0) * (ya ** 2).sum(0))
    a1 = (xa[1:] * xa[:-1]).sum(0) / (xa ** 2).sum(0); b1 = (ya[1:] * ya[:-1]).sum(0) / (ya ** 2).sum(0)
    ne = np.clip(n * (1 - a1 * b1) / (1 + a1 * b1), 3, n)
    return r, 2 * stats.t.sf(np.abs(r) * np.sqrt((ne - 2) / np.maximum(1 - r ** 2, 1e-9)), ne - 2)
out = dict(lat3=lat.reshape(-1, 3).mean(1), lon3=lon.reshape(-1, 3).mean(1))
one = np.ones(E3.shape[1:])
for nm, X in (("wind", T3), ("sam", sam[:, None, None] * one), ("oni", oni[:, None, None] * one)):
    for L in (0, 1):
        r, p = corr(X, E3, L); ok = np.isfinite(p); f = np.zeros_like(p, bool); f[ok] = R.bh_fdr(p[ok], 0.1)
        out[f"r_{nm}_{L}"], out[f"p_{nm}_{L}"], out[f"fdr_{nm}_{L}"] = r, p, f
        print(f"{nm} lag {L}: median r {np.nanmedian(r):+.2f}, p<0.05 in {np.mean(p[ok] < 0.05):.3f} of {ok.sum()} cells, FDR {f.sum()}")
out["sam_annual"], out["oni_annual"], out["years"] = an.SAM.values, an.oni.values, an.index.values
np.savez(OUT / "dk27_maps.npz", **out)
