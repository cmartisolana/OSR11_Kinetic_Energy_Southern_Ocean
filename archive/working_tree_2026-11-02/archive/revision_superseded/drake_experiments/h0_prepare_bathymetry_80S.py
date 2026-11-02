#!/usr/bin/env python3
"""Subset a global bedrock-elevation grid (ETOPO 2022 or GEBCO) to 80S-55S at ~1/8 deg for the Antarctic-slope
index (h0_antarctic_slope_index.py). Same procedure as revision/prepare_bathymetry.py, which stops at 70 S.
Usage: python h0_prepare_bathymetry_80S.py /path/to/global_file.nc [out.nc]"""
import sys
from pathlib import Path
import numpy as np, xarray as xr
src = Path(sys.argv[1]); out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parent / "data" / "bathymetry_80s.nc"
ds = xr.open_dataset(src)
var = next(v for v in ("z", "elevation", "Band1") if v in ds)
lat = "lat" if "lat" in ds.coords else "latitude"; lon = "lon" if "lon" in ds.coords else "longitude"
da = ds[var].sel({lat: slice(-80, -55)})
step = max(1, int(round(0.125 / float(abs(da[lat][1] - da[lat][0])))))
da = da.coarsen({lat: step, lon: step}, boundary="trim").mean().rename({lat: "latitude", lon: "longitude"}).rename("elevation")
out.parent.mkdir(parents=True, exist_ok=True); da.astype(np.float32).to_netcdf(out); print("saved", out, dict(da.sizes))
