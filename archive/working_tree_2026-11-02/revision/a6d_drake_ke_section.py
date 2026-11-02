"""A6d: distribution of KE and EKE across Drake Passage along the altimetric transect (v2 line of a6c, 121 points,
(65 W, 55 S) -> (60 W, 63 S)), and how it changes, 2001-2025, Product 1.
Per month and transect point: cross-track surface geostrophic velocity vn (positive eastward across the section),
KE = 0.5 |U_m|^2 of the monthly-mean flow, EKE = 0.5 |U_m - Ubar|^2 + 0.5 var_intra (fixed 2001-2025 reference, as a5).
Outputs (outputs/a6/):
  a6d_section_mean.csv  per point: lon, lat, distance, mean vn, mean KE, mean EKE (all months and 2001-2012 / 2014-2025),
                         fraction of months with ice_frac > 0.10
  a6d_section_series.csv monthly KE-weighted centroid latitude and KE-weighted rms width (ice-free points of every month)
  a6d_section_trends.csv trends per 50 km distance bin (from the northern end, as Gutierrez-Villanueva et al. 2023) of
                         vn (cm s-1 per decade), KE and EKE (% per decade), for 2001-2025 and for October 2005-April 2019
                         (the in situ period), and 2001-2025 trends of the KE centroid / width; statistic of trend_bcz
                         with 300 null series per series (single bootstrap CI, as b3_trends_fast); bins need 80 % of months
Usage: python a6d_drake_ke_section.py <revision/data dir> <out dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from concurrent.futures import ProcessPoolExecutor

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import a1_maps_bc as M

R_E = 6.371e6; KM = 111.195
LON0, LAT0, LON1, LAT1, NP = -65.0, -55.0, -60.0, -63.0, 121
x = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25)


def geometry():
    tt = np.linspace(0, 1, NP); lo = LON0 + tt * (LON1 - LON0); la = LAT0 + tt * (LAT1 - LAT0)
    latr = np.radians(0.5 * (LAT0 + LAT1))
    DX = R_E * np.cos(latr) * np.radians(LON1 - LON0); DY = R_E * np.radians(LAT1 - LAT0); L = np.hypot(DX, DY)
    nx, ny = DY / L, -DX / L
    if nx < 0: nx, ny = -nx, -ny
    dist = tt * L / 1e3
    return lo, la, dist, nx, ny


def cell(xx, mm, y, order):
    M.R.MONTHLY_ORDER = order
    y = M.R.deseasonalise(y, mm); m = np.isfinite(y)
    rng = np.random.default_rng(0)
    b, a, sd, phis, phi_bc, sig = M.R._z_stat(xx[m], y[m], True, 100, rng)
    return b, b / sd, sd, phi_bc, sig, m


def trend(args):
    key, y, scale, sl = args
    xx, mm = x[sl], month[sl]; y = y[sl]
    if np.isfinite(y).mean() < 0.8:
        return dict(key=key, cov=float(np.isfinite(y).mean()))
    res = {}
    for order in ("aic", "12"):
        b, z, sd, phi, sig, m = cell(xx, mm, y.copy(), order)
        zz = []
        for s in range(300 // 40 + 1):
            zz += M.null_z((xx, mm, phi, sig, m, 7000 + s, order))
        zz = np.abs(np.array(zz[:300])); zz = zz[np.isfinite(zz)]
        res[order] = dict(b=b, z=z, sd=sd, p=(1 + np.sum(zz >= abs(z))) / (1 + zz.size), q95=np.quantile(zz, 0.95))
    b = res["aic"]["b"]; p = max(res["aic"]["p"], res["12"]["p"])
    lo = min(b - res[o]["q95"] * res[o]["sd"] for o in res); hi = max(b + res[o]["q95"] * res[o]["sd"] for o in res)
    f = scale if scale else 1000 / np.nanmean(y)
    return dict(key=key, cov=float(np.isfinite(y).mean()), mean=float(np.nanmean(y)), slope=b, tr=b * f, lo=lo * f, hi=hi * f, p=p)


if __name__ == "__main__":
    data, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
    lo, la, dist, nx, ny = geometry()
    pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
    U, V, VI, ICE = [], [], [], []
    for yr in range(2001, 2026):
        for mo in range(1, 13):
            with xr.open_dataset(data / "monthly_allsat" / f"{yr}" / f"{mo:02d}.nc") as ds:
                d = ds.sel(longitude=slice(-68, -57), latitude=slice(-65, -53)).load()
            U.append(d.u_mean.interp(**pts).values); V.append(d.v_mean.interp(**pts).values)
            VI.append(d.var_intra.interp(**pts).values); ICE.append(d.ice_frac.interp(**pts).values)
    U, V, VI, ICE = map(np.array, (U, V, VI, ICE)); ICE = np.nan_to_num(ICE)
    ok = (ICE <= 0.10) & np.isfinite(U) & np.isfinite(V)
    U = np.where(ok, U, np.nan); V = np.where(ok, V, np.nan); VI = np.where(ok, VI, np.nan)
    vn = U * nx + V * ny
    KE = 0.5 * (U ** 2 + V ** 2) * 1e4                      # cm2 s-2
    ub, vb = np.nanmean(U, 0), np.nanmean(V, 0)
    EKE = (0.5 * ((U - ub) ** 2 + (V - vb) ** 2) + 0.5 * VI) * 1e4
    early, late = slice(0, 144), slice(156, 300)            # 2001-2012, 2014-2025
    S = pd.DataFrame(dict(lon=lo, lat=la, dist_km=dist, vn=np.nanmean(vn, 0), KE=np.nanmean(KE, 0), EKE=np.nanmean(EKE, 0),
                          KE_early=np.nanmean(KE[early], 0), KE_late=np.nanmean(KE[late], 0),
                          EKE_early=np.nanmean(EKE[early], 0), EKE_late=np.nanmean(EKE[late], 0),
                          ice_months=(ICE > 0.10).mean(0), valid=ok.mean(0)))
    S.to_csv(out / "a6d_section_mean.csv", index=False)
    # KE-weighted centroid and rms width over the points free of ice in every month
    perm = ok.all(0); w = np.where(perm[None, :], KE, 0.0)
    cen = (w * la).sum(1) / w.sum(1); wid = np.sqrt((w * (la - cen[:, None]) ** 2).sum(1) / w.sum(1)) * KM
    T = pd.DataFrame(dict(year=np.repeat(np.arange(2001, 2026), 12), month=month, x=x, centroid_lat=cen, rms_width_km=wid,
                          KE_section=np.nanmean(np.where(perm[None, :], KE, np.nan), 1)))
    T.to_csv(out / "a6d_section_series.csv", index=False)
    ALL, GV = slice(0, 300), slice(57, 220)                  # 2001-01..2025-12; 2005-10..2019-04
    jobs = [("centroid_lat|2001-2025", cen, KM * 10, ALL), ("rms_width_km|2001-2025", wid, 10, ALL),
            ("KE_section|2001-2025", T.KE_section.values, None, ALL)]
    edges = np.arange(0, 951, 50.0)
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (dist >= a) & (dist < b)
        if not sel.any(): continue
        with np.errstate(all="ignore"):
            ser = {"vn": np.nanmean(vn[:, sel], 1) * 100, "KE": np.nanmean(KE[:, sel], 1), "EKE": np.nanmean(EKE[:, sel], 1)}
        for v, yv in ser.items():
            for pn, sl in (("2001-2025", ALL), ("2005-2019", GV)):
                jobs.append((f"{v}|{a + 25:.0f}|{pn}", yv, 10 if v == "vn" else None, sl))
    with ProcessPoolExecutor(max_workers=2) as ex:
        rows = list(ex.map(trend, jobs))
    D = pd.DataFrame(rows); D.to_csv(out / "a6d_section_trends.csv", index=False)
    print(D.to_string()); print("perm points", int(perm.sum()), "southernmost perm lat", la[perm].min())
