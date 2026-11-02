"""dk30: memory test of dk29 per 3 deg cell (maps). Annual EKE per 3 deg cell (dk22 fields, Product 1), % of the cell mean,
NOT detrended; regressed on the SAM history X_T (same kernels as dk29: exponential, T = 0..30 yr, and Hogg et al. 2015).
Per cell: r, the EKE trend reproduced by the fit (% per decade) and the fraction of the observed cell trend, and the
surrogate p (2000 phase-randomised SAM series, same as dk29: one-sided, r_sur >= r). Also the best T per cell.
Usage: python dk30_memory_maps.py <revision dir> <drake out dir> <dk22 out dir> <sam 1957-2000 txt> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, DO, D22, SAMOLD, OUT = (Path(a) for a in sys.argv[1:6]); OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(29)
old = pd.read_csv(SAMOLD, sep=r"\s+", comment="#", header=None); old.columns = ["year"] + list(range(1, 13))
old = old.melt(id_vars="year", var_name="month", value_name="SAM")
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv")
S = pd.concat([old, A[["year", "month", "SAM"]]]).sort_values(["year", "month"]).reset_index(drop=True); assert len(S) == 69 * 12
def annual(m): return m.reshape(-1, 12).mean(1)
sam_a = annual(S.SAM.values); i0 = 44
def kern(T, n=69):
    k = np.arange(n)
    if T == "Hogg": h = (2 / 3) * np.exp(-k / 2.0) / 2.0 + (1 / 3) * np.exp(-k / 15.0) / 15.0
    elif T == 0: h = (k == 0).astype(float)
    else: h = np.exp(-k / T)
    return h / h.sum()
def filt(sa, h): return np.array([np.dot(h[:i + 1], sa[i::-1]) for i in range(i0, len(sa))])
def phase_rand(x):
    X = np.fft.rfft(x - x.mean()); ph = np.exp(1j * rng.uniform(0, 2 * np.pi, X.size)); ph[0] = 1
    if x.size % 2 == 0: ph[-1] = np.sign(np.cos(np.angle(ph[-1])))
    return np.fft.irfft(X * ph, n=x.size) + x.mean()
NS = 2000; SUR = np.array([annual(phase_rand(S.SAM.values)) for _ in range(NS)])
# EKE per 3 deg cell, annual, % of mean
Z = np.load(D22 / "dk22_fields.npz"); EK, lat, lon = Z["EK"], Z["lat"], Z["lon"]
EA = EK.reshape(25, 12, *EK.shape[1:]).mean(1)
E3 = EA.reshape(25, 10, 3, 120, 3).mean(axis=(2, 4)); Ep = 100 * (E3 / E3.mean(0) - 1)
m1 = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.coarsen(latitude=8, longitude=8, boundary="trim").mean().values
acc1 = (m1 >= -0.6) & (m1 <= 0.2) & np.isfinite(EK).all(0)
m3 = acc1.reshape(10, 3, 120, 3).mean(axis=(1, 3)) > 0.5
ok = m3 & np.isfinite(Ep).all(0)
yrs = np.arange(2001, 2026); t = yrs - yrs.mean()
Y = Ep[:, ok]; Yc = Y - Y.mean(0); obs_tr = 10 * (t @ Yc) / (t @ t)
TS = [0, 1, 2, 5, 10, 15, 20, 30]
out = dict(lat3=lat.reshape(-1, 3).mean(1), lon3=lon.reshape(-1, 3).mean(1), mask=ok, obs_trend=np.full(ok.shape, np.nan))
out["obs_trend"][ok] = obs_tr
RR = []
for T in TS + ["Hogg"]:
    h = kern(T); x = filt(sam_a, h); xc = x - x.mean()
    r = (xc @ Yc) / np.sqrt((xc @ xc) * (Yc ** 2).sum(0)); b = (xc @ Yc) / (xc @ xc)
    rep = b * 10 * (t @ xc) / (t @ t)
    xs = np.array([filt(s, h) for s in SUR]); xs -= xs.mean(1, keepdims=True)
    rs = (xs @ Yc) / np.sqrt((xs ** 2).sum(1)[:, None] * (Yc ** 2).sum(0)[None, :])
    p = (rs >= r[None, :]).mean(0)
    for nm, v in (("r", r), ("rep", rep), ("p", p), ("r95", np.quantile(rs, 0.95, axis=0))):
        f = np.full(ok.shape, np.nan); f[ok] = v; out[f"{nm}_{T}"] = f
    w = np.cos(np.radians(out["lat3"]))[:, None] * np.ones(120)[None, :]; ww = w[ok]
    RR.append(dict(T=str(T), median_r=np.median(r), frac_cells_p05=np.mean(p < 0.05), frac_cells_r_gt_0p5=np.mean(r > 0.5),
                   weighted_rep=np.sum(rep * ww) / ww.sum(), weighted_obs=np.sum(obs_tr * ww) / ww.sum()))
RT = pd.DataFrame(RR); RT["rep_over_obs"] = RT.weighted_rep / RT.weighted_obs
print(RT.round(3).to_string(index=False))
# best T per cell (exponential kernels only)
rstack = np.array([out[f"r_{T}"] for T in TS]); bi = np.nanargmax(np.where(np.isfinite(rstack), rstack, -9), 0)
out["bestT"] = np.where(ok, np.array(TS, float)[bi], np.nan)
print("best T distribution (cells):", pd.Series(out["bestT"][ok]).value_counts().sort_index().to_dict())
RT.to_csv(OUT / "dk30_summary.csv", index=False); np.savez(OUT / "dk30_maps.npz", **out)
