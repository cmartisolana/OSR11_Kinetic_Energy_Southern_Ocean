"""A1 maps with the bias-corrected, calibrated trend statistic (replaces the p of a1_maps.py).

Per 1° cell (same coarsening, deseasonalisation and Theil-Sen slope as a1_maps.py):
  z = slope / sd, sd = exact OLS-slope sd under the bias-corrected AR model (rstats._z_stat),
  i.e. the statistic of rstats.trend_bcz, computed with both noise models (AIC order and AR(12));
  as in trend_bcz, the reported p is the larger of the two.
Calibration: a per-cell double bootstrap is too costly for 2 x 8,184 cells, so the null distribution
of |z| is built per persistence class: cells are grouped into 10 quantile classes of the
bias-corrected AR coefficient sum; from 25 random cells per class, 40 null series each are simulated
from that cell's bias-corrected model and processed identically (1,000 null z per class). The p-value
of a cell is the Gaussian tail of |z| / kappa, kappa = max(q95 / 1.96, q99 / 2.576) of the class's null |z|
(the q99 term guards against the heavier-than-Gaussian tails found in the most persistent classes)
(continuous p, needed for BH over ~8,000 cells); the empirical class p is kept as *_p_emp.
BH-FDR alpha 0.10 as before.
Output: outputs/a1_maps_bc.nc

Purpose : KE/EKE trend maps with the bias-corrected, calibrated statistic and BH-FDR (Fig. 3a).
Inputs  : stored monthly KE/EKE fields via $OSR11_CACHE (01b); outputs/a1/a1_maps.nc (02f)
Outputs : outputs/a1/a1_maps_bc.nc (distributed; read by Fig. 3)
Run     : NPROC=8 python scripts/02_circumpolar/02g_eke_maps_bc.py
Seeds   : default_rng(0) / (1) for the class samples and null series.
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import os, sys
from pathlib import Path
import numpy as np, xarray as xr
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import norm

from osr11 import halfpower as A, rstats as R, paths as PA
from osr11.halfpower import coarsen


def cell(args):
    x, month, y, order = args
    R.MONTHLY_ORDER = order
    if np.isfinite(y).sum() < 240:
        return (np.nan,) * 5 + (None,)
    y = R.deseasonalise(y, month)
    m = np.isfinite(y); x, y = x[m], y[m]
    rng = np.random.default_rng(0)
    b, a, sd, phis, phi_bc, sig = R._z_stat(x, y, True, 100, rng)
    return b, b / sd, float(phi_bc.sum()), sig, sd, (phi_bc, m)


def null_z(args):
    x, month, phi_bc, sig, m, seed, order = args
    R.MONTHLY_ORDER = order
    rng = np.random.default_rng(seed)
    xs, ms = x[m], month[m]
    out = []
    for s_ in R._simulate_ar(phi_bc, sig, xs.size, 40, rng, burn=500):
        s_ = R.deseasonalise(s_, ms)
        bi, _, sdi, *_ = R._z_stat(xs, s_, True, 100, rng)
        out.append(bi / sdi)
    return out


if __name__ == "__main__":
    ke, lat, lon, acc, t = A.load()
    eke = np.load(A.CACHE / "eke.npy", mmap_mode="r")
    latc = lat.reshape(-1, 8).mean(1); lonc = lon.reshape(-1, 8).mean(1)
    old = xr.open_dataset(PA.OUT / "a1" / "a1_maps.nc")
    out = {}
    npr = int(os.environ.get("NPROC", "2"))
    for name, cube in [("KE", ke), ("EKE", eke)]:
        C = coarsen(np.asarray(cube, np.float32)) * 1e4
        ny, nx = C.shape[1:]
        sh = (ny, nx)
        P = {}
        for order in ("aic", "12"):
            jobs = [(t.x.values, t.month.values, C[:, i, j], order) for i in range(ny) for j in range(nx)]
            with ProcessPoolExecutor(max_workers=npr) as ex:
                res = list(ex.map(cell, jobs, chunksize=50))
            b = np.array([r[0] for r in res]); z = np.array([r[1] for r in res]); ps = np.array([r[2] for r in res])
            ok = np.isfinite(z)
            edges = np.quantile(ps[ok], np.linspace(0, 1, 11)); edges[-1] += 1e-9
            cls = np.clip(np.digitize(ps, edges) - 1, 0, 9)
            rng = np.random.default_rng(1)
            kappa = np.full(10, np.nan); nullz = {}
            for k in range(10):
                idx = rng.choice(np.nonzero(ok & (cls == k))[0], 25, replace=False)
                njobs = [(t.x.values, t.month.values.copy(), res[i][5][0], res[i][3], res[i][5][1], 100 * k + q, order)
                         for q, i in enumerate(idx)]
                with ProcessPoolExecutor(max_workers=npr) as ex:
                    zz = np.abs(np.concatenate([np.array(v) for v in ex.map(null_z, njobs)]))
                zz = zz[np.isfinite(zz)]; nullz[k] = zz
                kappa[k] = max(np.quantile(zz, 0.95) / 1.96, np.quantile(zz, 0.99) / 2.576)   # heavier-than-Gaussian tails
                print(name, order, "class", k, "phi_sum", edges[k].round(3), "-", edges[k + 1].round(3), "kappa", kappa[k].round(3),
                      "q99 ratio", (np.quantile(zz, 0.99) / 2.576).round(3), flush=True)
            p = np.full(b.size, np.nan); pe = np.full(b.size, np.nan)
            p[ok] = 2 * norm.sf(np.abs(z[ok]) / kappa[cls[ok]])
            for i in np.nonzero(ok)[0]:
                zz = nullz[cls[i]]; pe[i] = (1 + np.sum(zz >= abs(z[i]))) / (1 + zz.size)
            P[order] = p
            out[f"{name}_z_{order}"] = (("lat", "lon"), z.reshape(sh))
            out[f"{name}_phi_sum_bc_{order}"] = (("lat", "lon"), ps.reshape(sh))
            out[f"{name}_p_{order}"] = (("lat", "lon"), p.reshape(sh))
            out[f"{name}_p_emp_{order}"] = (("lat", "lon"), pe.reshape(sh))
            out[f"{name}_kappa_{order}"] = (("class",), kappa)
        p = np.fmax(P["aic"], P["12"])
        fdr = R.bh_fdr(p.reshape(sh), 0.10).astype(np.int8)
        out[f"{name}_slope"] = (("lat", "lon"), b.reshape(sh))
        out[f"{name}_p"] = (("lat", "lon"), p.reshape(sh))
        out[f"{name}_fdr"] = (("lat", "lon"), fdr)
        out[f"{name}_mean"] = (("lat", "lon"), np.nanmean(C, 0))
        print(name, "cells", int(np.isfinite(p).sum()), "p<0.05", int(np.nansum(p < 0.05)), "FDR", int(fdr.sum()),
              "| old: p<0.05", int(np.nansum(old[f"{name}_p"].values < 0.05)), "FDR", int(old[f"{name}_fdr"].values.sum()), flush=True)
    xr.Dataset(out, coords={"lat": latc, "lon": lonc, "class": np.arange(10)}).to_netcdf(PA.out("a1", "a1_maps_bc.nc"))
