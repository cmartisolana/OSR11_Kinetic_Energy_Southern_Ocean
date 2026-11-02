#!/usr/bin/env python3
"""Daily-field extraction for reproducing Liu et al. (2026, JGR Oceans, doi:10.1029/2025JC023527).

Liu et al. method (their Sect. 2.2.1-2.2.2), applied here to SEALEVEL_GLO_PHY_L4_MY_008_047:
  1. |grad ADT| (m per 100 km) on daily fields over the gradient domain 150W-60W, 70S-30S
  2. a cell is "meander" on a day if |grad ADT| > 20 % of that day's maximum over the domain
  3. occurrence = number of detection days, summed over 4-month windows (done later from
     the monthly counts written here)
  4. dominant axis, width, speed, peaks/troughs: computed later (liu_reproduce.py)
This script streams daily fields month by month and writes, for the study region
(150W-80W, 60S-48S, with a 2-degree margin), per month:
  count      number of days the cell is detected (int16)
  gradsum    sum over days of |grad ADT| (m/100 km), for Liu's tie-break rule
  speedsum   sum over days of geostrophic speed sqrt(ugos^2+vgos^2) (m/s)
  ndays      number of valid days per cell
  gmax       daily domain-maximum gradient (time series), and the alternative threshold
             based on the study-region maximum (count_reg, sensitivity)
Two grids are processed: the native 0.125 deg grid, and a 0.25 deg grid obtained by 2x2
block-averaging ADT/ugos/vgos (Liu et al. used 0.25 deg fields; DT2024 008_047 is 0.125 deg).

Usage (same environment and login as extract_monthly.py):
    python revision/extract_liu.py                      # 1993-2025, PAR
    python revision/extract_liu.py --start-year 1993 --end-year 2025
Output: revision/data/liu_par/{g0125,g025}/YYYY_MM.nc ; re-running skips existing months.
"""
from __future__ import annotations

import argparse
import calendar
import datetime as dt
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

try:
    import copernicusmarine
except ImportError:  # pragma: no cover
    sys.exit("copernicusmarine is not installed in this environment.")

DATASET = "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D"
GRAD_DOMAIN = dict(lon=(-150.0, -60.0), lat=(-70.0, -30.0))     # Liu et al. step 1
STUDY = dict(lon=(-152.0, -78.0), lat=(-62.0, -46.0))           # 150W-80W, 60S-48S plus 2 deg margin
STUDY_CORE = dict(lon=(-150.0, -80.0), lat=(-60.0, -48.0))
FRAC = 0.20
R_E = 6.371e6
REPO_ROOT = Path(__file__).resolve().parents[1]


def grad_mag(adt, lat, lon):
    """|grad ADT| in m per 100 km (centred differences on the sphere)."""
    phi = np.radians(lat)[:, None]
    dy = R_E * np.radians(np.gradient(lat))[:, None]
    dx = R_E * np.cos(phi) * np.radians(np.gradient(lon))[None, :]
    gy = np.gradient(adt, axis=-2) / dy
    gx = np.gradient(adt, axis=-1) / dx
    return np.hypot(gx, gy) * 1e5


def coarsen2(a):
    t, ny, nx = a.shape
    a = a[:, : ny // 2 * 2, : nx // 2 * 2]
    with np.errstate(invalid="ignore"):
        return np.nanmean(a.reshape(t, ny // 2, 2, nx // 2, 2), axis=(2, 4))


def process(adt, u, v, lat, lon):
    g = grad_mag(adt, lat, lon)                                  # (t, y, x)
    gmax = np.nanmax(g.reshape(g.shape[0], -1), axis=1)          # per day, gradient domain
    sl = (lat >= STUDY["lat"][0]) & (lat <= STUDY["lat"][1])
    so = (lon >= STUDY["lon"][0]) & (lon <= STUDY["lon"][1])
    cl = (lat >= STUDY_CORE["lat"][0]) & (lat <= STUDY_CORE["lat"][1])
    co = (lon >= STUDY_CORE["lon"][0]) & (lon <= STUDY_CORE["lon"][1])
    gs = g[:, sl][:, :, so]
    gmax_reg = np.nanmax(g[:, cl][:, :, co].reshape(g.shape[0], -1), axis=1)
    det = gs > FRAC * gmax[:, None, None]
    det_reg = gs > FRAC * gmax_reg[:, None, None]
    spd = np.hypot(u, v)[:, sl][:, :, so]
    ok = np.isfinite(gs)
    out = {
        "count": det.sum(0).astype(np.int16),
        "count_reg": det_reg.sum(0).astype(np.int16),
        "gradsum": np.nansum(gs, 0).astype(np.float32),
        "speedsum": np.nansum(spd, 0).astype(np.float32),
        "ndays": ok.sum(0).astype(np.int16),
    }
    return out, gmax.astype(np.float32), gmax_reg.astype(np.float32), lat[sl], lon[so]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--start-year", type=int, default=1993)
    p.add_argument("--end-year", type=int, default=2025)
    p.add_argument("--retries", type=int, default=3)
    a = p.parse_args()
    root = REPO_ROOT / "revision" / "data" / "liu_par"
    for gname in ("g0125", "g025"):
        (root / gname).mkdir(parents=True, exist_ok=True)
    months = [(y, m) for y in range(a.start_year, a.end_year + 1) for m in range(1, 13)]
    t0 = time.time(); failed = []
    for i, (y, m) in enumerate(months, 1):
        if all((root / g / f"{y}_{m:02d}.nc").exists() for g in ("g0125", "g025")):
            continue
        last = calendar.monthrange(y, m)[1]
        for attempt in range(1, a.retries + 1):
            try:
                ds = copernicusmarine.open_dataset(
                    dataset_id=DATASET, variables=["adt", "ugos", "vgos"],
                    minimum_longitude=GRAD_DOMAIN["lon"][0], maximum_longitude=GRAD_DOMAIN["lon"][1],
                    minimum_latitude=GRAD_DOMAIN["lat"][0], maximum_latitude=GRAD_DOMAIN["lat"][1],
                    start_datetime=f"{y}-{m:02d}-01T00:00:00", end_datetime=f"{y}-{m:02d}-{last:02d}T23:59:59",
                ).load()
                lat = ds.latitude.values.astype(float); lon = ds.longitude.values.astype(float)
                A = ds.adt.values.astype(np.float64); U = ds.ugos.values.astype(np.float64); V = ds.vgos.values.astype(np.float64)
                grids = {"g0125": (A, U, V, lat, lon),
                         "g025": (coarsen2(A), coarsen2(U), coarsen2(V),
                                  lat[: lat.size // 2 * 2].reshape(-1, 2).mean(1), lon[: lon.size // 2 * 2].reshape(-1, 2).mean(1))}
                for gname, (aa, uu, vv, la, lo) in grids.items():
                    out, gmax, gmax_reg, las, los = process(aa, uu, vv, la, lo)
                    dims = ("latitude", "longitude")
                    D = xr.Dataset({k: (dims, v_) for k, v_ in out.items()},
                                   coords=dict(latitude=las, longitude=los, day=np.arange(1, aa.shape[0] + 1)))
                    D["gmax"] = (("day",), gmax); D["gmax_reg"] = (("day",), gmax_reg)
                    D.attrs.update(dataset_id=DATASET, year=y, month=m, frac=FRAC, grid=gname,
                                   grad_domain=str(GRAD_DOMAIN), study=str(STUDY), script="revision/extract_liu.py",
                                   created=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
                    f = root / gname / f"{y}_{m:02d}.nc"; tmp = f.with_suffix(".tmp.nc")
                    D.to_netcdf(tmp, encoding={k: {"zlib": True, "complevel": 4} for k in out}); tmp.rename(f)
                break
            except Exception as exc:
                print(f"  {y}-{m:02d} attempt {attempt} failed: {exc}", flush=True)
                time.sleep(10 * attempt)
        else:
            failed.append(f"{y}-{m:02d}")
        print(f"[{i:3d}/{len(months)}] {y}-{m:02d} done ({(time.time() - t0) / 60:.1f} min)", flush=True)
    print("FAILED: " + ", ".join(failed) if failed else f"All months present in {root}")


if __name__ == "__main__":
    main()
