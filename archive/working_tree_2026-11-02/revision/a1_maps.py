"""A1 maps: KE and EKE trend maps with field significance (Fig. 1c,d).

Monthly fields are block-averaged to 1 deg x 1 deg (8 x 8 cells of 0.125 deg;
cells with < 50 % valid pixels dropped), deseasonalised, Theil-Sen slope per
cell, AR-sieve null (n_boot = 200), then Benjamini-Hochberg FDR at
alpha_FDR = 0.10 (Wilks 2016). Period 2001-2025.
The FDR uses p_z = 2*(1 - Phi(|b| / sd_null)) (Gaussian tail of the surrogate null):
the empirical p (kept as *_p_emp) has a floor of 1/201, which makes BH over ~8000
cells unreachable (needs >400 cells at the floor). Same z-approach as rstats.trend_calibrated.
Output: revision/outputs/a1_maps.nc
"""
import os, sys
from pathlib import Path
import numpy as np, xarray as xr
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import a2_run as A, rstats as R

def coarsen(cube):
    nt, ny, nx = cube.shape
    c = cube.reshape(nt, ny // 8, 8, nx // 8, 8)
    valid = np.isfinite(c).sum(axis=(2, 4))
    m = np.nanmean(c, axis=(2, 4))
    m[valid < 32] = np.nan
    return m

def _levinson_aic(r, maxp=12):
    """Yule-Walker AR fits of order 1..maxp (Levinson-Durbin); AIC-selected order."""
    n = r.size
    r = r - r.mean()
    acov = np.array([np.dot(r[: n - k], r[k:]) / n for k in range(maxp + 1)])
    if acov[0] <= 0:
        return np.zeros(0), 0.0
    phi = np.zeros(0); err = acov[0]; best = (n * np.log(err), np.zeros(0), err)
    for k in range(1, maxp + 1):
        kk = (acov[k] - np.dot(phi, acov[1:k][::-1])) / err
        phi = np.r_[phi - kk * phi[::-1], kk]
        err = err * (1 - kk * kk)
        aic = n * np.log(err) + 2 * k
        if aic < best[0]:
            best = (aic, phi.copy(), err)
    return best[1], float(np.sqrt(best[2]))


def cell(args):
    """Fast per-cell version of rstats.trend (Yule-Walker sieve, no Hamed-Rao)."""
    x, month, y = args
    if np.isfinite(y).sum() < 240:
        return np.nan, np.nan, np.nan
    y = R.deseasonalise(y, month)
    m = np.isfinite(y); x, y = x[m], y[m]
    b, a = R.theil_sen(x, y)
    phis, sig = _levinson_aic(y - (a + b * x))
    rng = np.random.default_rng(0)
    sims = R._simulate_ar(phis, sig, x.size, 200, rng)
    null = R._ts_many(x, sims)
    p_emp = (1 + np.sum(np.abs(null) >= abs(b))) / 201
    sd = np.std(null, ddof=1)
    p_z = float(2 * norm.sf(abs(b) / sd)) if sd > 0 else np.nan
    return b, p_emp, p_z

if __name__ == "__main__":
    ke, lat, lon, acc, t = A.load()
    eke = np.load(A.CACHE / "eke.npy", mmap_mode="r")
    out = {}
    latc = lat.reshape(-1, 8).mean(1); lonc = lon.reshape(-1, 8).mean(1)
    for name, cube in [("KE", ke), ("EKE", eke)]:
        C = coarsen(np.asarray(cube, np.float32)) * 1e4
        ny, nx = C.shape[1:]
        jobs = [(t.x.values, t.month.values, C[:, i, j]) for i in range(ny) for j in range(nx)]
        with ProcessPoolExecutor(max_workers=int(os.environ.get("NPROC", "2"))) as ex:
            res = list(ex.map(cell, jobs, chunksize=50))
        sl = np.array([r[0] for r in res]).reshape(ny, nx)
        pe = np.array([r[1] for r in res]).reshape(ny, nx)
        p = np.array([r[2] for r in res]).reshape(ny, nx)
        out[f"{name}_slope"] = (("lat", "lon"), sl)
        out[f"{name}_p_emp"] = (("lat", "lon"), pe)
        out[f"{name}_p"] = (("lat", "lon"), p)
        fdr = R.bh_fdr(p, 0.10).astype(np.int8)
        out[f"{name}_fdr"] = (("lat", "lon"), fdr)
        out[f"{name}_mean"] = (("lat", "lon"), np.nanmean(C, 0))
        print(name, "cells", np.isfinite(p).sum(), "p_emp<0.05", int(np.nansum(pe < 0.05)),
              "p_z<0.05", int(np.nansum(p < 0.05)), "FDR", int(fdr.sum()), flush=True)
    xr.Dataset(out, coords={"lat": latc, "lon": lonc}).to_netcdf(HERE / "outputs" / "a1_maps.nc")
