"""F4m map: per-cell trend of the monthly sigma0 = 27.3 cross-stream slope, 2004-2025 (28 Oct 2026), calibrated exactly as
a1_maps_bc.py: deseasonalised monthly series, z = Theil-Sen slope / sd under the bias-corrected AR model (rstats._z_stat), for
both noise models (AIC order and AR(12)); per noise model, kappa from 10 persistence classes x 25 cells x 40 null series;
p = max of the two; BH-FDR 10 % over the domain (f4 domain cells with the slope defined in all 264 months; f4m_monthly.py).
Usage: python f4m_map.py <revision dir> <f4 out dir> <f4m out dir>     (NPROC workers, default 2)
"""
import os, sys
from pathlib import Path
import numpy as np, xarray as xr
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import norm

REV, F4, FM = (Path(a) for a in sys.argv[1:4])
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402


def cell(args):
    x, month, y, order = args
    R.MONTHLY_ORDER = order
    y = R.deseasonalise(y, month); rng = np.random.default_rng(0)
    b, a, sd, phis, phi_bc, sig = R._z_stat(x, y, True, 100, rng)
    return b, b / sd, float(phi_bc.sum()), phi_bc, sig, y.mean()


def null_z(args):
    x, month, phi_bc, sig, seed, order = args
    R.MONTHLY_ORDER = order; rng = np.random.default_rng(seed); out = []
    for s_ in R._simulate_ar(phi_bc, sig, x.size, 40, rng, burn=500):
        s_ = R.deseasonalise(s_, month); bi, _, sdi, *_ = R._z_stat(x, s_, True, 100, rng); out.append(bi / sdi)
    return out


if __name__ == "__main__":
    C = np.load(FM / "f4m_p273_monthly.npz"); YY, MM, PM, lonr = C["year"], C["month"], C["p"], C["lon"]
    x = YY + (MM - 0.5) / 12
    M = xr.open_dataset(F4 / "f4_tilt_maps.nc"); glat, glon = M.lat.values, M.lon.values
    L0 = np.where(lonr > 180, lonr - 360, lonr); o = np.argsort(L0); KM = 111.195; adt = M["adt"].values
    def grad(F):
        dy = np.gradient(F, glat, axis=-2) / KM * 100
        dx = (np.roll(F, -1, axis=-1) - np.roll(F, 1, axis=-1)) / 2.0 / (KM * np.cos(np.radians(glat))[:, None]) * 100
        return dx, dy
    ax_, ay_ = grad(adt); an = np.hypot(ax_, ay_)
    dx, dy = grad(PM[:, :, o].astype(float)); SL = dx * ax_ / an + dy * ay_ / an
    dom = (M["buffer100"].values > 0) & (M["ice"].values == 0) & np.isfinite(M["signed_27p3_pct"].values) & np.isfinite(SL).all(0)
    idx = np.argwhere(dom); npr = int(os.environ.get("NPROC", "2")); P = {}; B = None
    print("cells", len(idx), flush=True)
    for order in ("aic", "12"):
        jobs = [(x, MM, SL[:, i, j], order) for i, j in idx]
        with ProcessPoolExecutor(max_workers=npr) as ex:
            res = list(ex.map(cell, jobs, chunksize=40))
        b = np.array([r[0] for r in res]); z = np.array([r[1] for r in res]); ps = np.array([r[2] for r in res]); mu = np.array([r[5] for r in res])
        ok = np.isfinite(z); edges = np.quantile(ps[ok], np.linspace(0, 1, 11)); edges[-1] += 1e-9
        cls = np.clip(np.digitize(ps, edges) - 1, 0, 9); rng = np.random.default_rng(1); kappa = np.full(10, np.nan)
        for k in range(10):
            pick = rng.choice(np.nonzero(ok & (cls == k))[0], 25, replace=False)
            njobs = [(x, MM, res[i][3], res[i][4], 100 * k + q, order) for q, i in enumerate(pick)]
            with ProcessPoolExecutor(max_workers=npr) as ex:
                zz = np.abs(np.concatenate([np.array(v) for v in ex.map(null_z, njobs)]))
            zz = zz[np.isfinite(zz)]; kappa[k] = max(np.quantile(zz, 0.95) / 1.96, np.quantile(zz, 0.99) / 2.576)
            print(order, "class", k, "kappa", round(float(kappa[k]), 3), flush=True)
        p = np.full(b.size, np.nan); p[ok] = 2 * norm.sf(np.abs(z[ok]) / kappa[cls[ok]]); P[order] = p; B = b
    p = np.fmax(P["aic"], P["12"]); fdr = R.bh_fdr(p, 0.10).astype(np.int8)
    PCT = np.full(dom.shape, np.nan); PP = np.full(dom.shape, np.nan); F = np.zeros(dom.shape, np.int8)
    for q, (i, j) in enumerate(idx):
        PCT[i, j] = 1000 * B[q] / mu[q]; PP[i, j] = p[q]; F[i, j] = fdr[q]
    xr.Dataset({"slope_pct": (("lat", "lon"), PCT), "slope_p": (("lat", "lon"), PP), "slope_fdr": (("lat", "lon"), F),
                "domain": (("lat", "lon"), dom.astype(np.int8))}, coords={"lat": glat, "lon": glon}).to_netcdf(FM / "f4m_map.nc")
    print(f"cells {len(idx)}, positive {int((B > 0).sum())}, p<0.05 {int((p < 0.05).sum())}, FDR {int(fdr.sum())} ({int(((B > 0) & (fdr > 0)).sum())} positive)", flush=True)
