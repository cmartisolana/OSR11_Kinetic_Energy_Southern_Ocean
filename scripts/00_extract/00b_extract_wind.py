"""Monthly wind-stress extraction for the reviewer-2 revision (plan item A6e).

Uses the same product as notebook 05 (WIND_GLO_PHY_L4_MY_012_006, monthly
dataset) but keeps the wind-STRESS variables, which the product provides, as
well as the stress-equivalent wind used so far. Output: one NetCDF,
revision/data/wind_monthly.nc, for 2001-2025 over 70S-30S.

    python revision/extract_wind.py
    python revision/extract_wind.py --dataset-id <other P1M id>   # if the default is not found

Purpose : download monthly wind stress (Product 2) for 2001-2025, 70-30 S.
Inputs  : Copernicus Marine WIND_GLO_PHY_L4_MY_012_006 (doi:10.48670/moi-00185); needs `copernicusmarine login`.
Outputs : $OSR11_DATA/wind_taux.nc (eastward_stress, northward_stress, eastward_wind, northward_wind).
          (The working repository wrote data/wind_monthly.nc and the analyses read the same file renamed wind_taux.nc;
          the default output name is now the name the analyses read.)
Run     : python scripts/00_extract/00b_extract_wind.py
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import argparse
import sys
from pathlib import Path

import xarray as xr

try:
    import copernicusmarine
except ImportError:  # pragma: no cover
    sys.exit("copernicusmarine is not installed in this environment.")

from osr11 import paths as PA  # noqa: E402
WANTED = ["eastward_stress", "northward_stress", "eastward_wind", "northward_wind"]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset-id", default="cmems_obs-wind_glo_phy_my_l4_P1M")
    p.add_argument("--start", default="2001-01-01")
    p.add_argument("--end", default="2025-12-31")
    p.add_argument("--lat-min", type=float, default=-70.0)
    p.add_argument("--lat-max", type=float, default=-30.0)
    p.add_argument("--out", type=Path, default=PA.WIND)
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
