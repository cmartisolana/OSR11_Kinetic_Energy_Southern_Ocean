"""Step 3: KE and EKE trend maps on 1-degree cells, Product 4, with the statistic of 02g_eke_maps_bc.py (P1 source of main Fig. 3a:
outputs/a1/a1_maps_bc.nc), written with the same variable names and grid so that both products plot with the same code.

Purpose : Product 4 counterpart of outputs/a1/a1_maps_bc.nc.
Inputs  : Product 4 monthly fields (p4_common); outputs/a1/a1_maps_bc.nc (grid and variable names of Product 1, for the check)
Outputs : outputs/p4/p4_eke_trend_map.nc  variables as a1_maps_bc.nc: {KE,EKE}_slope (cm2 s-2 per YEAR, as in a1_maps_bc.nc; Fig. 3a plots
          10 x slope = per decade), _p (max of the two noise models), _fdr (BH 10 %), _mean, _z_{aic,12}, _phi_sum_bc_{aic,12},
          _p_{aic,12}, _p_emp_{aic,12}, _kappa_{aic,12}; coords lat (-64.5..-35.5), lon (-179.5..179.5), class (0..9)
Method  : monthly fields as 01b_halfpower_fields.py, for Product 4:
              ke  = 0.5 |U_m|^2,  eke = 0.5 |U_m - U_lt|^2 + 0.5 var_intra   (x 1e4 -> cm2 s-2)
          U_lt = 2001-2025 mean of the DAILY velocities = n_days-weighted mean of the monthly means (01a computes it from the daily
          fields; the monthly files reproduce it exactly up to float32 rounding); NaN on the permanent ice mask, built as in
          01a / 01c (any JJA day with flag_ice > 0 in 2001-2025) = any JJA month with ice_frac > 0 (ice_frac is the fraction of days
          with flag_ice == 1). Coarsening as osr11.halfpower.coarsen: block mean of the native pixels in each 1-degree cell (4 x 4 at
          1/4 deg instead of 8 x 8 at 1/8 deg), NaN where fewer than 50 % of the pixels are valid. Statistic, calibration
          (10 persistence classes, 25 cells x 40 null series per class, seeds 0 and 1, kappa = max(q95/1.96, q99/2.576)), p = max of
          AIC-order and AR(12) p, BH-FDR 10 %: code of 02g_eke_maps_bc.py unchanged (functions copied verbatim).
Run     : NPROC=2 python scripts/06_product4/p4_03_trend_map.py
Seeds   : default_rng(0) per cell, default_rng(1) for the class samples, 100*k+q for the null series (as 02g).
"""
import os
from concurrent.futures import ProcessPoolExecutor
import p4_common as C
import numpy as np, xarray as xr
from scipy.stats import norm
R = C.R


def cell(args):                                   # verbatim from 02g
    x, month, y, order = args
    R.MONTHLY_ORDER = order
    if np.isfinite(y).sum() < 240:
        return (np.nan,) * 5 + (None,)
    y = R.deseasonalise(y, month)
    m = np.isfinite(y); x, y = x[m], y[m]
    rng = np.random.default_rng(0)
    b, a, sd, phis, phi_bc, sig = R._z_stat(x, y, True, 100, rng)
    return b, b / sd, float(phi_bc.sum()), sig, sd, (phi_bc, m)


def null_z(args):                                 # verbatim from 02g
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


def coarsen4(cube):
    """osr11.halfpower.coarsen for the 1/4-deg grid: 4 x 4 pixels per 1-deg cell, NaN where < 50 % (8 of 16) are valid."""
    nt, ny, nx = cube.shape
    c = cube.reshape(nt, ny // 4, 4, nx // 4, 4)
    valid = np.isfinite(c).sum(axis=(2, 4))
    with np.errstate(invalid="ignore"):
        m = np.nanmean(c, axis=(2, 4))
    m[valid < 8] = np.nan
    return m


def fields():
    lat, lon = C.grid()
    U, V, S2, I, N = (C.stack(v) for v in ("u_mean", "v_mean", "var_intra", "ice_frac", "n_days"))
    ice = np.zeros(U.shape[1:], bool)
    for k in np.nonzero(np.isin(C.MON, (6, 7, 8)))[0]:
        ice |= np.nan_to_num(np.asarray(I[k], float)) > 0
    su = np.zeros(U.shape[1:]); sv = np.zeros_like(su); sn = np.zeros_like(su)
    for k in range(U.shape[0]):
        u, v, n = np.asarray(U[k], float), np.asarray(V[k], float), np.asarray(N[k], float)
        ok = np.isfinite(u) & np.isfinite(v) & (n > 0)
        su[ok] += (n * u)[ok]; sv[ok] += (n * v)[ok]; sn[ok] += n[ok]
    ult, vlt = su / np.where(sn > 0, sn, np.nan), sv / np.where(sn > 0, sn, np.nan)
    ke = np.empty(U.shape, np.float32); eke = np.empty(U.shape, np.float32)
    for k in range(U.shape[0]):
        u, v, s2 = np.asarray(U[k], float), np.asarray(V[k], float), np.asarray(S2[k], float)
        a = 0.5 * (u ** 2 + v ** 2); b = 0.5 * ((u - ult) ** 2 + (v - vlt) ** 2) + 0.5 * s2
        a[ice] = np.nan; b[ice] = np.nan; ke[k] = a; eke[k] = b
    return ke, eke, lat, lon, ice


if __name__ == "__main__":
    ke, eke, lat, lon, ice = fields()
    latc = lat.reshape(-1, 4).mean(1); lonc = lon.reshape(-1, 4).mean(1)
    ref = xr.open_dataset(C.OUT / "a1" / "a1_maps_bc.nc")
    assert np.allclose(latc, ref.lat.values) and np.allclose(lonc, ref.lon.values), "1-deg grid differs from a1_maps_bc.nc"
    print("permanent ice pixels", int(ice.sum()), flush=True)
    t = C.T
    out = {}
    npr = int(os.environ.get("NPROC", "2"))
    for name, cube in [("KE", ke), ("EKE", eke)]:
        Cc = coarsen4(np.asarray(cube, np.float32)) * 1e4
        ny, nx = Cc.shape[1:]
        sh = (ny, nx)
        P = {}
        for order in ("aic", "12"):
            jobs = [(t.x.values, t.month.values, Cc[:, i, j], order) for i in range(ny) for j in range(nx)]
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
                kappa[k] = max(np.quantile(zz, 0.95) / 1.96, np.quantile(zz, 0.99) / 2.576)
                print(name, order, "class", k, "kappa", kappa[k].round(3), flush=True)
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
        with np.errstate(invalid="ignore"):
            out[f"{name}_mean"] = (("lat", "lon"), np.nanmean(Cc, 0).astype(np.float32))
        print(name, "cells", int(np.isfinite(p).sum()), "p<0.05", int(np.nansum(p < 0.05)), "FDR", int(fdr.sum()),
              "| P1: cells", int(np.isfinite(ref[f"{name}_p"].values).sum()), "p<0.05", int(np.nansum(ref[f"{name}_p"].values < 0.05)),
              "FDR", int(ref[f"{name}_fdr"].values.sum()), flush=True)
    ds = xr.Dataset(out, coords={"lat": latc, "lon": lonc, "class": np.arange(10)})
    ds.attrs.update(product="Product 4 (two-satellite climate altimetry, 1/4 deg), 2001-2025",
                    slope_units="cm2 s-2 per year (x10 for per decade, as Fig. 3a)", script="scripts/06_product4/p4_03_trend_map.py",
                    method="02g_eke_maps_bc.py applied to Product 4; 4x4 coarsening; permanent JJA ice mask from ice_frac")
    ds.to_netcdf(C.HERE / "p4_eke_trend_map.nc"); print("done", flush=True)
