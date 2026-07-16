"""Run the two-sat pipeline directly from Copernicus Marine remote access.

This wrapper streams monthly two-sat altimetry via the Copernicus Marine
Toolbox, computes a study-period permanent ice mask from winter months,
computes long-term mean velocities from all months, produces monthly
half-power point outputs under ``outputs/monthly_half-power_points_two_sat``,
then runs the trend scripts included in this repo.

Usage example
-------------
python scripts/run_all_two_sat.py --years 2001 2025

"""

from __future__ import annotations

import argparse
import calendar as _cal
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import xarray as xr

from utils.chambers_fronts import detect_fronts_from_dataset
from utils.regions_and_masks import preprocessing_region

try:
    import copernicusmarine

    _CMEMS_AVAILABLE = True
except ImportError:
    copernicusmarine = None
    _CMEMS_AVAILABLE = False


DEFAULT_TWO_SAT_DATASET_ID = "c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg_P1D"
WINTER_MONTHS = {6, 7, 8}


def _require_cmems() -> None:
    if not _CMEMS_AVAILABLE:
        raise ImportError(
            "copernicusmarine is not installed. Install it with: pip install copernicusmarine"
        )


def iter_year_months(years: Iterable[int]) -> list[tuple[int, int]]:
    return [(year, month) for year in years for month in range(1, 13)]


def open_remote_month(
    dataset_id: str,
    year: int,
    month: int,
    lonmin: float,
    lonmax: float,
    latmin: float,
    latmax: float,
    variables: list[str],
) -> xr.Dataset:
    last_day = _cal.monthrange(year, month)[1]
    return copernicusmarine.open_dataset(
        dataset_id=dataset_id,
        minimum_longitude=lonmin,
        maximum_longitude=lonmax,
        minimum_latitude=latmin,
        maximum_latitude=latmax,
        start_datetime=f"{year}-{month:02d}-01",
        end_datetime=f"{year}-{month:02d}-{last_day:02d}",
        variables=variables,
    )


def _slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


def _run_cache_dir(
    out_root: Path,
    dataset_id: str,
    years: list[int],
    lonmin: float,
    lonmax: float,
    latmin: float,
    latmax: float,
) -> Path:
    run_key = (
        f"{_slugify(dataset_id)}__"
        f"yrs{years[0]}_{years[-1]}__"
        f"lon{lonmin:g}_{lonmax:g}__"
        f"lat{latmin:g}_{latmax:g}"
    )
    return out_root / "_cache" / run_key


def _load_cached_dataarray(path: Path) -> Optional[xr.DataArray]:
    if not path.exists():
        return None
    return xr.load_dataarray(path)


def _save_dataarray(dataarray: xr.DataArray, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dataarray.to_netcdf(path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run full pipeline using remote two_sat data")
    parser.add_argument(
        "--dataset-id",
        default=DEFAULT_TWO_SAT_DATASET_ID,
        help="Copernicus Marine dataset ID for the two-sat product.",
    )
    parser.add_argument("--years", type=int, nargs=2, required=True, metavar=("START", "END"), help="Start and end year (inclusive)")
    parser.add_argument("--lonmin", type=float, default=-180.0)
    parser.add_argument("--lonmax", type=float, default=180.0)
    parser.add_argument("--latmin", type=float, default=-65.0)
    parser.add_argument("--latmax", type=float, default=-35.0)
    parser.add_argument("--out-dir", default="outputs/monthly_half-power_points_two_sat")
    parser.add_argument("--test-years", type=int, default=0, help="If >0, process only the last N years for testing")
    parser.add_argument("--overwrite", action="store_true", help="Recompute monthly outputs and trends even if files already exist.")
    args = parser.parse_args()

    start_year, end_year = args.years
    years = list(range(start_year, end_year + 1))
    if args.test_years > 0:
        years = years[-args.test_years :]

    if not years:
        raise ValueError("No years selected for processing.")

    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    repo_root = Path(__file__).resolve().parent.parent
    trend_script = repo_root / "scripts" / "ke_eke_2d_trends_two_sat.py"
    trend_output = repo_root / "outputs" / "trends" / "ke_eke_2d_trends_two_sat.nc"

    expected_month_files = [out_root / str(y) / f"{m:02d}.nc" for y, m in iter_year_months(years)]
    all_month_outputs_exist = bool(expected_month_files) and all(path.exists() for path in expected_month_files)

    if not args.overwrite and all_month_outputs_exist:
        if trend_output.exists():
            print("All requested monthly outputs and trend output already exist; nothing to do.")
            return
        print("All requested monthly outputs already exist; rebuilding trends only.")
        trend_only = True
    else:
        trend_only = False

    all_year_months = iter_year_months(years)
    winter_year_months = [(year, month) for year, month in all_year_months if month in WINTER_MONTHS]

    ice_mask: Optional[xr.DataArray] = None
    u_ltm: Optional[xr.DataArray] = None
    v_ltm: Optional[xr.DataArray] = None
    made_any_month = False

    if not trend_only or args.overwrite:
        _require_cmems()
        print(f"Fetching remote files for years {years[0]}-{years[-1]} from dataset: {args.dataset_id}")
        print("Authenticating with Copernicus Marine Toolbox...")
        copernicusmarine.login()

        print(f"Total months (all)  : {len(all_year_months)}")
        print(f"Winter months (JJA) : {len(winter_year_months)}")
        print("Opening remote datasets to compute ice mask and long-term means (this may take a while)...")

        cache_dir = _run_cache_dir(out_root, args.dataset_id, years, args.lonmin, args.lonmax, args.latmin, args.latmax)
        ice_mask_file = cache_dir / "ice_mask.nc"
        u_ltm_file = cache_dir / "u_ltm.nc"
        v_ltm_file = cache_dir / "v_ltm.nc"

        ice_mask = _load_cached_dataarray(ice_mask_file)
        u_ltm = _load_cached_dataarray(u_ltm_file)
        v_ltm = _load_cached_dataarray(v_ltm_file)

        if ice_mask is None or u_ltm is None or v_ltm is None:
            ice_ever = None
            lats_ref = None
            lons_ref = None

            for idx, (year, month) in enumerate(winter_year_months, start=1):
                ds_mo = open_remote_month(
                    args.dataset_id,
                    year,
                    month,
                    args.lonmin,
                    args.lonmax,
                    args.latmin,
                    args.latmax,
                    ["flag_ice"],
                )
                ice_vals = ds_mo["flag_ice"].values
                if ice_ever is None:
                    lats_ref = ds_mo["latitude"].values
                    lons_ref = ds_mo["longitude"].values
                    ice_ever = np.zeros(ice_vals.shape[1:], dtype=bool)
                ice_ever |= (ice_vals > 0).any(axis=0)
                ds_mo.close()
                print(f"  ice mask {idx}/{len(winter_year_months)} - {year}/{month:02d}")

            if ice_ever is None or lats_ref is None or lons_ref is None:
                raise RuntimeError("Unable to derive the study-period ice mask from remote data.")

            shape = ice_ever.shape
            u_sum = np.zeros(shape, dtype=np.float64)
            v_sum = np.zeros(shape, dtype=np.float64)
            n_valid_u = np.zeros(shape, dtype=np.int32)
            n_valid_v = np.zeros(shape, dtype=np.int32)

            for idx, (year, month) in enumerate(all_year_months, start=1):
                ds_mo = open_remote_month(
                    args.dataset_id,
                    year,
                    month,
                    args.lonmin,
                    args.lonmax,
                    args.latmin,
                    args.latmax,
                    ["ugos", "vgos", "flag_ice"],
                )
                ds_mo.load()

                for t_idx in range(ds_mo.sizes["time"]):
                    ds_t = ds_mo.isel(time=t_idx)
                    u = ds_t["ugos"].values.astype(np.float64)
                    v = ds_t["vgos"].values.astype(np.float64)

                    mask_out = ice_ever | ~np.isfinite(u)
                    valid_u = ~mask_out
                    valid_v = ~(ice_ever | ~np.isfinite(v))

                    u[~valid_u] = 0.0
                    v[~valid_v] = 0.0
                    u_sum += u
                    v_sum += v
                    n_valid_u += valid_u.astype(np.int32)
                    n_valid_v += valid_v.astype(np.int32)

                ds_mo.close()
                print(f"  means {idx}/{len(all_year_months)} - {year}/{month:02d}")

            print("Passes A & B complete.")

            coords = {"latitude": lats_ref, "longitude": lons_ref}
            ice_mask = xr.DataArray(ice_ever, dims=["latitude", "longitude"], coords=coords, name="ice_mask")
            with np.errstate(invalid="ignore"):
                u_mean_arr = np.where(n_valid_u > 0, u_sum / n_valid_u, np.nan)
                v_mean_arr = np.where(n_valid_v > 0, v_sum / n_valid_v, np.nan)

            u_ltm = xr.DataArray(u_mean_arr, dims=["latitude", "longitude"], coords=coords, name="u_ltm")
            v_ltm = xr.DataArray(v_mean_arr, dims=["latitude", "longitude"], coords=coords, name="v_ltm")

            _save_dataarray(ice_mask, ice_mask_file)
            _save_dataarray(u_ltm, u_ltm_file)
            _save_dataarray(v_ltm, v_ltm_file)
        else:
            print(f"Reusing cached ice mask and long-term means from {cache_dir}")

        for y, m in all_year_months:
            out_dir = out_root / str(y)
            out_file = out_dir / f"{m:02d}.nc"

            if out_file.exists() and not args.overwrite:
                print(f"Skipping existing monthly output {out_file}")
                continue

            try:
                ds = open_remote_month(
                    args.dataset_id,
                    y,
                    m,
                    args.lonmin,
                    args.lonmax,
                    args.latmin,
                    args.latmax,
                    ["ugos", "vgos", "flag_ice"],
                )
            except Exception:
                continue

            ds = preprocessing_region(ds, args.lonmin, args.lonmax, args.latmin, args.latmax, ice_mask=ice_mask)

            if "time" in ds.dims:
                ds_mean = ds.mean(dim="time")
            else:
                ds_mean = ds.copy()

            ds_mean = ds_mean.load()

            u_anom = ds_mean["ugos"] - u_ltm
            v_anom = ds_mean["vgos"] - v_ltm
            eke = 0.5 * (u_anom ** 2 + v_anom ** 2)
            ds_mean["eke"] = eke
            ds_mean["ke"] = 0.5 * (ds_mean["ugos"] ** 2 + ds_mean["vgos"] ** 2)

            print(f"Detecting fronts for {y}/{m:02d}...")
            results = detect_fronts_from_dataset(ds_mean, use_numba=True, sample_every=1)

            lat_min, lat_max = args.latmin, args.latmax
            filtered_results = {}
            for lon, envs in results.items():
                valid = [e for e in envs if lat_min <= e["lat_mid"] <= lat_max]
                if valid:
                    filtered_results[lon] = valid

            n_fronts = sum(len(envs) for envs in filtered_results.values())
            if n_fronts == 0:
                print(f"No fronts detected for {y}/{m:02d}")
                ds.close()
                continue

            front_lons = np.zeros(n_fronts)
            front_lats = np.zeros(n_fronts)
            envelope_lat_south = np.zeros(n_fronts)
            envelope_lat_north = np.zeros(n_fronts)
            envelope_lon = np.zeros(n_fronts)
            peak_lat = np.zeros(n_fronts)
            peak_val = np.zeros(n_fronts)
            width_km = np.zeros(n_fronts)

            lat_vals = ds_mean.latitude.values
            idx = 0
            for lon, envs in filtered_results.items():
                for e in envs:
                    front_lons[idx] = e["lon_mid"]
                    front_lats[idx] = e["lat_mid"]
                    envelope_lat_south[idx] = lat_vals[e["isouth"]]
                    envelope_lat_north[idx] = lat_vals[e["inorth"]]
                    envelope_lon[idx] = lon
                    peak_lat[idx] = e["peak_lat"]
                    peak_val[idx] = e["peak_val"]
                    width_km[idx] = e["width_km"]
                    idx += 1

            ds_output = xr.Dataset(
                {
                    "ke": (("latitude", "longitude"), ds_mean["ke"].values),
                    "eke": (("latitude", "longitude"), ds_mean["eke"].values),
                    "front_lon": (("front",), front_lons),
                    "front_lat": (("front",), front_lats),
                    "envelope_lon": (("front",), envelope_lon),
                    "envelope_lat_south": (("front",), envelope_lat_south),
                    "envelope_lat_north": (("front",), envelope_lat_north),
                    "peak_lat": (("front",), peak_lat),
                    "peak_val": (("front",), peak_val),
                    "width_km": (("front",), width_km),
                },
                coords={
                    "latitude": ds_mean.latitude.values,
                    "longitude": ds_mean.longitude.values,
                    "front": np.arange(n_fronts),
                },
                attrs={
                    "year": int(y),
                    "month": f"{m:02d}",
                    "n_fronts_detected": int(n_fronts),
                    "source": "two-sat Copernicus Marine remote access",
                },
            )

            out_dir.mkdir(parents=True, exist_ok=True)
            ds_output.to_netcdf(out_file)
            print(f"Saved {n_fronts} fronts -> {out_file}")
            made_any_month = True
            ds.close()

    if trend_only and not args.overwrite:
        candidate_scripts = [trend_script]
    elif not made_any_month and not args.overwrite:
        candidate_scripts = [trend_script] if not trend_output.exists() else []
    else:
        candidate_scripts = [
            trend_script,
            repo_root / "scripts" / "ke_eke_2d_trends_since2010.py",
            repo_root / "scripts" / "ke_trends_analysis.py",
            repo_root / "scripts" / "plot_ke_trends.py",
        ]

    cmds = []
    trend_overwrite = args.overwrite or made_any_month or not trend_output.exists()
    for script_path in candidate_scripts:
        if not script_path.exists():
            continue
        if script_path == trend_script and not trend_overwrite and trend_output.exists():
            print(f"Skipping trend recomputation because {trend_output} already exists")
            continue
        if script_path == trend_script and trend_overwrite:
            cmds.append([sys.executable, str(script_path), "--overwrite"])
        else:
            cmds.append([sys.executable, str(script_path)])

    if not cmds:
        print("No downstream scripts need to run.")
        return

    print("Running trend and plotting scripts...")
    for cmd in cmds:
        print("->", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command failed: {e}")

    missing = [script_path.name for script_path in candidate_scripts if not script_path.exists()]
    if missing:
        print(f"Skipped missing scripts: {', '.join(missing)}")


if __name__ == "__main__":
    main()
