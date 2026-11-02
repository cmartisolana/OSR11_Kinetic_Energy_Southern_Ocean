#!/usr/bin/env python3
"""Monthly wind-stress extraction for the reviewer-2 revision (plan item A6e).

Uses the same product as notebook 05 (WIND_GLO_PHY_L4_MY_012_006, monthly
dataset) but keeps the wind-STRESS variables, which the product provides, as
well as the stress-equivalent wind used so far. Output: one NetCDF,
revision/data/wind_monthly.nc, for 2001-2025 over 70S-30S.

Usage (same environment as notebook 05):
    python revision/extract_wind.py
    python revision/extract_wind.py --dataset-id <other P1M id>   # if the default is not found
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import xarray as xr

try:
    import copernicusmarine
except ImportError:  # pragma: no cover
    sys.exit("copernicusmarine is not installed in this environment.")

REPO_ROOT = Path(__file__).resolve().parents[1]
WANTED = ["eastward_stress", "northward_stress", "eastward_wind", "northward_wind"]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset-id", default="cmems_obs-wind_glo_phy_my_l4_P1M")
    p.add_argument("--start", default="2001-01-01")
    p.add_argument("--end", default="2025-12-31")
    p.add_argument("--lat-min", type=float, default=-70.0)
    p.add_argument("--lat-max", type=float, default=-30.0)
    p.add_argument("--out", type=Path, default=REPO_ROOT / "revision" / "data" / "wind_monthly.nc")
    a = p.parse_args()

    ds = copernicusmarine.open_dataset(
        dataset_id=a.dataset_id,
        minimum_longitude=-180, maximum_longitude=180,
        minimum_latitude=a.lat_min, maximum_latitude=a.lat_max,
        start_datetime=a.start, end_datetime=a.end,
    )
    have = [v for v in WANTED if v in ds]
    print("Variables in dataset:", list(ds.data_vars))
    if "eastward_stress" not in have:
        print("WARNING: no eastward_stress in this dataset; try another --dataset-id.")
    ds = ds[have].load()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    ds.attrs["dataset_id"] = a.dataset_id
    ds.to_netcdf(a.out, encoding={v: {"zlib": True, "complevel": 4} for v in have})
    print(f"Saved {a.out}  time steps: {ds.sizes.get('time')}")


if __name__ == "__main__":
    main()
