#!/usr/bin/env python3
"""Daily DUACS all-satellite fields (ugos, vgos, adt) for one longitude-latitude window, one NetCDF per year.

Test of the path-distortion -> EKE relation at daily resolution (b2 is monthly). Default window: the
Pacific-Antarctic Ridge window centred at 150W (135-165W) plus its downstream band to 115W, 62-46S, 2015-2019.
About 240 MB per year at 1/8 deg.

Run in the environment used for extract_monthly.py (copernicusmarine login done):
    python revision/extract_daily_window.py                      # defaults below
    python revision/extract_daily_window.py --lon-min 65 --lon-max 115 --lat-min -58 --lat-max -40 --tag kerg
Output: revision/data/daily_<tag>/<year>.nc with variables ugos, vgos, adt (time, latitude, longitude), float32.
Re-running skips years that already exist.
"""
from __future__ import annotations
import argparse, datetime as dt, sys, time
from pathlib import Path
import numpy as np, xarray as xr
try:
    import copernicusmarine
except ImportError:
    sys.exit("copernicusmarine is not installed in this environment.")

DATASET = "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D"
REPO_ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--start-year", type=int, default=2015); p.add_argument("--end-year", type=int, default=2019)
p.add_argument("--lon-min", type=float, default=-165.0); p.add_argument("--lon-max", type=float, default=-115.0)
p.add_argument("--lat-min", type=float, default=-62.0); p.add_argument("--lat-max", type=float, default=-46.0)
p.add_argument("--tag", default="par150w"); p.add_argument("--retries", type=int, default=3)
a = p.parse_args()
out_dir = REPO_ROOT / "revision" / "data" / f"daily_{a.tag}"; out_dir.mkdir(parents=True, exist_ok=True)
enc = {v: {"zlib": True, "complevel": 4, "dtype": "float32"} for v in ("ugos", "vgos", "adt")}
t0 = time.time()
for y in range(a.start_year, a.end_year + 1):
    f = out_dir / f"{y}.nc"
    if f.exists():
        print(y, "exists"); continue
    for attempt in range(1, a.retries + 1):
        try:
            ds = copernicusmarine.open_dataset(dataset_id=DATASET, minimum_longitude=a.lon_min, maximum_longitude=a.lon_max,
                                               minimum_latitude=a.lat_min, maximum_latitude=a.lat_max,
                                               start_datetime=f"{y}-01-01T00:00:00", end_datetime=f"{y}-12-31T23:59:59")
            ds = ds[["ugos", "vgos", "adt"]].load()
            ds.attrs.update(dataset_id=DATASET, created=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                            script="revision/extract_daily_window.py", window=f"{a.lon_min},{a.lon_max},{a.lat_min},{a.lat_max}")
            tmp = f.with_suffix(".tmp.nc"); ds.to_netcdf(tmp, encoding=enc); tmp.rename(f)
            print(f"{y} done, {ds.sizes['time']} days, {(time.time() - t0) / 60:.1f} min elapsed", flush=True)
            break
        except Exception as exc:
            print(f"  {y} attempt {attempt} failed: {exc}", flush=True); time.sleep(10 * attempt)
    else:
        print("FAILED", y)
print("finished")
