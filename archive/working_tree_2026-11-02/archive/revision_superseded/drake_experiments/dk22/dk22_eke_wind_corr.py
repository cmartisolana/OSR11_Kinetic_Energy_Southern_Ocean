"""dk22: correlation maps between EKE and the eastward wind stress, annual means 2001-2025, both linearly detrended.
EKE (Product 1) = 0.5|U - Ubar|^2 + 0.5 var_intra (as dk17/s2), block-averaged from 1/8 deg to 1 deg (cells with all
64 pixels valid in every month and outside the sea-ice mask). Wind stress (Product 2, eastward) block-averaged from
1/4 deg to 1 deg. Correlations at lag 0 and lag 1 yr (EKE lagging) with (i) the local wind stress in the same 1 deg cell
and (ii) the ACC-band wind-stress index of Fig. 4 (dk17/s1 'tau_accmask_ACC').
Significance: two-sided t test with an effective sample size n_eff = n (1 - r1 r2)/(1 + r1 r2) (lag-1 autocorrelations of
the two detrended series; Bretherton et al. 1999), and Benjamini-Hochberg FDR at 10 % over ocean cells.
Usage: python dk22_eke_wind_corr.py <revision dir> <drake out dir> <out dir> [CG]\n(CG: average the 1 deg fields over CG x CG deg cells before correlating; a cell is kept only if all sub-cells are valid)"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, DO, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
ny, nx = lat.size // 8, lon.size // 8
blk = lambda a: a.reshape(ny, 8, nx, 8)
sf = OUT / "dk22_fields.npz"
if not sf.exists():
    nv = np.zeros(lat.size * lon.size).reshape(lat.size, lon.size); su = np.zeros_like(nv); sv = np.zeros_like(nv)
    for f in files:
        with xr.open_dataset(f) as d: u, v = d.u_mean.values.astype(float), d.v_mean.values.astype(float)
        ok = np.isfinite(u) & np.isfinite(v); nv += ok; su += np.where(ok, u, 0); sv += np.where(ok, v, 0)
    valid = (nv == len(files)) & ~ice; ub, vb = su / np.maximum(nv, 1), sv / np.maximum(nv, 1)
    cellok = blk(valid).all(axis=(1, 3))
    EK = np.zeros((len(files), ny, nx))
    for t, f in enumerate(files):
        with xr.open_dataset(f) as d: u, v, s2 = (d[k].values.astype(float) for k in ("u_mean", "v_mean", "var_intra"))
        e = (0.5 * ((u - ub) ** 2 + (v - vb) ** 2) + 0.5 * s2) * 1e4
        EK[t] = np.where(cellok, np.nanmean(blk(np.where(valid, e, np.nan)), axis=(1, 3)), np.nan)
    W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress
    W = W.sel(latitude=slice(lat.min() - 0.2, lat.max() + 0.2))
    Wc = W.coarsen(latitude=4, longitude=4, boundary="trim").mean()
    Wc = Wc.interp(latitude=lat.reshape(ny, 8).mean(1), longitude=lon.reshape(nx, 8).mean(1), method="nearest").values
    np.savez(sf, EK=EK, TAU=Wc, lat=lat.reshape(ny, 8).mean(1), lon=lon.reshape(nx, 8).mean(1))
Z = np.load(sf); EK, TAU = Z["EK"], Z["TAU"]; LAT, LON = Z["lat"], Z["lon"]
CG = int(sys.argv[4]) if len(sys.argv) > 4 else 1          # optional coarse-graining factor (cells of CG x CG degrees)
if CG > 1:
    cg = lambda a: a.reshape(a.shape[0], a.shape[1] // CG, CG, a.shape[2] // CG, CG).mean(axis=(2, 4))
    EK, TAU = cg(EK), cg(TAU); LAT, LON = LAT.reshape(-1, CG).mean(1), LON.reshape(-1, CG).mean(1)
ann = lambda a: a.reshape(25, 12, *a.shape[1:]).mean(1)
def detr(a):
    t = np.arange(a.shape[0]) - (a.shape[0] - 1) / 2; m = a.mean(0)
    b = (t[:, None, None] * (a - m)).sum(0) / (t ** 2).sum() if a.ndim == 3 else (t * (a - m)).sum() / (t ** 2).sum()
    return a - m - (t[:, None, None] * b if a.ndim == 3 else t * b)
E = detr(ann(EK)); T = detr(ann(TAU))
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); ti = detr(A.groupby("year").tau_accmask_ACC.mean().values)
def r1(x): return np.sum(x[1:] * x[:-1], 0) / np.sum(x * x, 0)
def corr(x, y):
    x = x - x.mean(0); y = y - y.mean(0)
    return np.sum(x * y, 0) / np.sqrt(np.sum(x * x, 0) * np.sum(y * y, 0))
out = {}
for lag in (0, 1):
    e = E[lag:]; n = e.shape[0]
    for nm, x in (("local", T[:25 - lag] if lag else T), ("index", (ti[:25 - lag] if lag else ti)[:, None, None] * np.ones(E.shape[1:]))):
        r = corr(x, e); a1, a2 = r1(x - x.mean(0)), r1(e - e.mean(0)); ne = np.clip(n * (1 - a1 * a2) / (1 + a1 * a2), 3, n)
        tt = r * np.sqrt((ne - 2) / np.maximum(1 - r ** 2, 1e-9)); p = 2 * stats.t.sf(np.abs(tt), ne - 2)
        ok = np.isfinite(p); fdr = np.zeros_like(p, bool)
        fdr[ok] = R.bh_fdr(p[ok], 0.1) if hasattr(R, "bh_fdr") else False
        out[f"r_{nm}_lag{lag}"] = r; out[f"p_{nm}_lag{lag}"] = p; out[f"fdr_{nm}_lag{lag}"] = fdr
        print(nm, lag, "median r", np.nanmedian(r).round(2), "frac p<0.05", np.round(np.nanmean(p < 0.05), 3), "FDR cells", int(fdr.sum()), "of", int(ok.sum()))
np.savez(OUT / (f"dk22_corr_{CG}deg.npz" if CG > 1 else "dk22_corr.npz"), lat=LAT, lon=LON, **out)
