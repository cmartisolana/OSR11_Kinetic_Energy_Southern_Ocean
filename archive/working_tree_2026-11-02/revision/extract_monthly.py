#!/usr/bin/env python3
"""One-off monthly extraction pass for the reviewer-2 revision (plan item A0).

Streams daily DUACS L4 fields from Copernicus Marine (same access route as
notebook 01) and writes one small NetCDF per month containing everything the
later analyses need, so no further downloads are required:

    u_mean, v_mean   monthly-mean surface geostrophic velocity          [m s-1]
    adt_mean         monthly-mean absolute dynamic topography           [m]
    var_intra        within-month velocity variance
                     mean_t(|u - U_m|^2) = mean_t((u-U)^2 + (v-V)^2)      [m2 s-2]
    ice_frac         fraction of days in the month with flag_ice == 1   [0-1]
    n_days           number of valid daily velocity values              [-]

From these, for any reference mean flow ubar (computed later):
    KE_m  = 0.5 |U_m|^2
    EKE_m = 0.5 |U_m - ubar|^2 + 0.5 var_intra          (exact)

Usage (run in the environment you used for notebook 01, logged in with
`copernicusmarine login`):

    python revision/extract_monthly.py --product allsat
    python revision/extract_monthly.py --product twosat

Re-running skips months that already exist, so an interrupted run can simply
be restarted.
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

PRODUCTS = {
    # SEALEVEL_GLO_PHY_L4_MY_008_047 (all-satellite, 0.125 deg)
    "allsat": "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D",
    # SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057 (two-satellite, 0.25 deg)
    "twosat": "c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg_P1D",
}

REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--product", choices=PRODUCTS, required=True)
    p.add_argument("--start-year", type=int, default=2001)
    p.add_argument("--end-year", type=int, default=2025)
    p.add_argument("--lat-min", type=float, default=-65.0)
    p.add_argument("--lat-max", type=float, default=-35.0)
    p.add_argument("--lon-min", type=float, default=-180.0)
    p.add_argument("--lon-max", type=float, default=180.0)
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Default: revision/data/monthly_<product>/")
    p.add_argument("--retries", type=int, default=3)
    return p.parse_args()


def open_month(dataset_id: str, year: int, month: int, a: argparse.Namespace) -> xr.Dataset:
    last = calendar.monthrange(year, month)[1]
    ds = copernicusmarine.open_dataset(
        dataset_id=dataset_id,
        minimum_longitude=a.lon_min,
        maximum_longitude=a.lon_max,
        minimum_latitude=a.lat_min,
        maximum_latitude=a.lat_max,
        start_datetime=f"{year}-{month:02d}-01T00:00:00",
        end_datetime=f"{year}-{month:02d}-{last:02d}T23:59:59",
    )
    keep = [v for v in ("ugos", "vgos", "adt", "flag_ice") if v in ds]
    missing = {"ugos", "vgos", "adt"} - set(keep)
    if missing:
        raise RuntimeError(f"{dataset_id} lacks variables {sorted(missing)}")
    return ds[keep]


def monthly_stats(ds: xr.Dataset) -> xr.Dataset:
    ds = ds.load()
    u = ds["ugos"].values.astype(np.float64)   # (time, lat, lon)
    v = ds["vgos"].values.astype(np.float64)
    adt = ds["adt"].values.astype(np.float64)

    valid = np.isfinite(u) & np.isfinite(v)
    n_days = valid.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        U = np.nanmean(np.where(valid, u, np.nan), axis=0)
        V = np.nanmean(np.where(valid, v, np.nan), axis=0)
        var_intra = np.nanmean(np.where(valid, (u - U) ** 2 + (v - V) ** 2, np.nan), axis=0)
        A = np.nanmean(adt, axis=0)
        if "flag_ice" in ds:
            ice = ds["flag_ice"].values.astype(np.float64)
            ice_frac = np.nanmean(np.where(np.isfinite(ice), ice > 0.5, np.nan), axis=0)
        else:
            ice_frac = np.full(U.shape, np.nan)

    coords = {"latitude": ds["latitude"].values, "longitude": ds["longitude"].values}
    dims = ("latitude", "longitude")
    out = xr.Dataset(
        {
            "u_mean": (dims, U.astype(np.float32), {"units": "m s-1", "long_name": "monthly-mean ugos"}),
            "v_mean": (dims, V.astype(np.float32), {"units": "m s-1", "long_name": "monthly-mean vgos"}),
            "adt_mean": (dims, A.astype(np.float32), {"units": "m", "long_name": "monthly-mean adt"}),
            "var_intra": (dims, var_intra.astype(np.float32),
                          {"units": "m2 s-2", "long_name": "within-month velocity variance mean_t(|u-U_m|^2)"}),
            "ice_frac": (dims, ice_frac.astype(np.float32),
                         {"units": "1", "long_name": "fraction of days with flag_ice == 1"}),
            "n_days": (dims, n_days.astype(np.int16), {"long_name": "valid daily velocity values"}),
        },
        coords=coords,
    )
    out.attrs["n_time_steps"] = int(ds.sizes["time"])
    return out


def main() -> None:
    a = parse_args()
    dataset_id = PRODUCTS[a.product]
    out_dir = a.out_dir or (REPO_ROOT / "revision" / "data" / f"monthly_{a.product}")
    out_dir.mkdir(parents=True, exist_ok=True)
    enc = {v: {"zlib": True, "complevel": 4} for v in
           ("u_mean", "v_mean", "adt_mean", "var_intra", "ice_frac", "n_days")}

    months = [(y, m) for y in range(a.start_year, a.end_year + 1) for m in range(1, 13)]
    t0 = time.time()
    failed = []
    for i, (y, m) in enumerate(months, 1):
        f = out_dir / f"{y}" / f"{m:02d}.nc"
        if f.exists():
            continue
        f.parent.mkdir(parents=True, exist_ok=True)
        for attempt in range(1, a.retries + 1):
            try:
                ds = open_month(dataset_id, y, m, a)
                out = monthly_stats(ds)
                out.attrs.update({
                    "dataset_id": dataset_id,
                    "year": y, "month": m,
                    "created": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                    "script": "revision/extract_monthly.py",
                    "source_attrs": str({k: ds.attrs.get(k) for k in ("product_version", "dataset_version", "history") if k in ds.attrs}),
                })
                tmp = f.with_suffix(".tmp.nc")
                out.to_netcdf(tmp, encoding=enc)
                tmp.rename(f)
                break
            except Exception as exc:  # network hiccups are common on long runs
                print(f"  {y}-{m:02d} attempt {attempt} failed: {exc}", flush=True)
                time.sleep(10 * attempt)
        else:
            failed.append(f"{y}-{m:02d}")
        el = (time.time() - t0) / 60
        print(f"[{i:3d}/{len(months)}] {y}-{m:02d} done  ({el:.1f} min elapsed)", flush=True)

    if failed:
        print("FAILED months (re-run the script to retry):", ", ".join(failed))
    else:
        print(f"All {len(months)} months present in {out_dir}")


if __name__ == "__main__":
    main()
