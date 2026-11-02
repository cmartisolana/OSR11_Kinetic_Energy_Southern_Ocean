#!/usr/bin/env python3
"""Monthly fields from the C3S two-satellite sea-level product in any version (vdt2018 / vdt2021 / vdt2024), in the
same format as revision/extract_monthly.py (u_mean, v_mean, adt_mean, var_intra, ice_frac, n_days), for the
reprocessing comparison of plan section 2 (row 4: earlier vs current product version, same period, domain, method).
Martinez-Moreno et al. (2021) used the all-sat vDT2018, which is no longer distributed; the two-sat vDT2018 is
(CDS dataset 'satellite-sea-level-global', version 'vdt2018'; checked 19 Oct 2026 via the CDS API).

Needs a CDS account and ~/.cdsapirc (https://cds.climate.copernicus.eu/how-to-api), and `pip install cdsapi`.
Daily files are downloaded one month at a time, reduced, and deleted. Restartable (finished months skipped).
Usage:
  python extract_c3s_twosat.py --version vdt2018 --start-year 2001 --end-year 2020
  -> revision/data/monthly_twosat_vdt2018/YYYY/MM.nc
Then run m1_comparison.py on a copy whose `files("twosat")` points to that folder (or symlink it), same period.
"""
from __future__ import annotations
import argparse, calendar, tempfile, zipfile, shutil
from pathlib import Path
import numpy as np, xarray as xr

REPO = Path(__file__).resolve().parents[2]


def monthly_stats(ds: xr.Dataset, a) -> xr.Dataset:
    ds = ds.sel(latitude=slice(a.lat_min, a.lat_max))
    if float(ds.longitude.max()) > 180:     # C3S grids are 0..360
        ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180).sortby("longitude")
    ds = ds.load()
    u = ds["ugos"].values.astype(np.float64); v = ds["vgos"].values.astype(np.float64); adt = ds["adt"].values.astype(np.float64)
    valid = np.isfinite(u) & np.isfinite(v); n_days = valid.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        U = np.nanmean(np.where(valid, u, np.nan), 0); V = np.nanmean(np.where(valid, v, np.nan), 0)
        var_intra = np.nanmean(np.where(valid, (u-U)**2 + (v-V)**2, np.nan), 0); A = np.nanmean(adt, 0)
        if "flag_ice" in ds:
            ice = ds["flag_ice"].values.astype(np.float64)
            ice_frac = np.nanmean(np.where(np.isfinite(ice), ice > 0.5, np.nan), 0)
        else:
            ice_frac = np.full(U.shape, np.nan)
    dims = ("latitude", "longitude")
    return xr.Dataset({"u_mean": (dims, U.astype(np.float32)), "v_mean": (dims, V.astype(np.float32)),
                       "adt_mean": (dims, A.astype(np.float32)), "var_intra": (dims, var_intra.astype(np.float32)),
                       "ice_frac": (dims, ice_frac.astype(np.float32)), "n_days": (dims, n_days.astype(np.int16))},
                      coords={"latitude": ds.latitude.values, "longitude": ds.longitude.values})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--version", default="vdt2018", choices=["vdt2018", "vdt2021", "vdt2024"])
    p.add_argument("--start-year", type=int, default=2001); p.add_argument("--end-year", type=int, default=2020)
    p.add_argument("--lat-min", type=float, default=-65.0); p.add_argument("--lat-max", type=float, default=-35.0)
    p.add_argument("--out-dir", type=Path, default=None)
    a = p.parse_args()
    import cdsapi
    c = cdsapi.Client()
    out = a.out_dir or REPO / "revision" / "data" / f"monthly_twosat_{a.version}"
    for y in range(a.start_year, a.end_year + 1):
        for m in range(1, 13):
            f = out / f"{y}" / f"{m:02d}.nc"
            if f.exists():
                continue
            f.parent.mkdir(parents=True, exist_ok=True)
            tmp = Path(tempfile.mkdtemp())
            try:
                z = tmp / "d.zip"
                c.retrieve("satellite-sea-level-global", {"variable": "daily", "year": [f"{y}"], "month": [f"{m:02d}"],
                           "day": [f"{d:02d}" for d in range(1, calendar.monthrange(y, m)[1] + 1)], "version": a.version}, str(z))
                if zipfile.is_zipfile(z):
                    zipfile.ZipFile(z).extractall(tmp)
                ncs = sorted(tmp.rglob("*.nc"))
                ds = xr.open_mfdataset([str(n) for n in ncs], combine="by_coords")
                o = monthly_stats(ds, a); o.attrs.update(version=a.version, year=y, month=m, source="C3S satellite-sea-level-global")
                o.to_netcdf(f.with_suffix(".tmp.nc")); f.with_suffix(".tmp.nc").rename(f); ds.close()
                print(y, m, "done", flush=True)
            except Exception as e:
                print(y, m, "FAILED:", e, flush=True)
            finally:
                shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
