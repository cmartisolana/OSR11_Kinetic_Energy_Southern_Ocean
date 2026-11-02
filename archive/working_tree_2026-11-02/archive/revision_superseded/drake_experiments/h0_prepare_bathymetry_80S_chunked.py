#!/usr/bin/env python3
"""Low-memory version of h0_prepare_bathymetry_80S.py for the GEBCO 2026 15-arcsec subset (80-45 S).
Block-averages 30x30 GEBCO cells (1/240 deg) to 1/8 deg, reading 30 rows at a time; cell centres then coincide
with the CMEMS 0.125 deg grid (-79.9375, ...). Output format identical to h0_prepare_bathymetry_80S.py.
Usage: python h0_prepare_bathymetry_80S_chunked.py <gebco nc> [out.nc]"""
import sys
from pathlib import Path
import numpy as np, netCDF4 as nc, xarray as xr
src = sys.argv[1]
out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parents[1] / "data" / "bathymetry_80s.nc"
d = nc.Dataset(src); lat = d["lat"][:]; lon = d["lon"][:]; z = d["elevation"]
step = int(round(0.125 / abs(lat[1] - lat[0])))
nr, nc_ = len(lat)//step, len(lon)//step
E = np.empty((nr, nc_), np.float32)
for i in range(nr):
    blk = z[i*step:(i+1)*step, :nc_*step].astype(np.float32)
    E[i] = blk.reshape(step, nc_, step).mean(axis=(0, 2))
la = lat[:nr*step].reshape(nr, step).mean(1); lo = lon[:nc_*step].reshape(nc_, step).mean(1)
da = xr.DataArray(E, coords=dict(latitude=la, longitude=lo), dims=("latitude", "longitude"), name="elevation",
                  attrs=dict(units="m", source=Path(src).name, note=f"{step}x{step} block mean of GEBCO 2026"))
out.parent.mkdir(parents=True, exist_ok=True); da.to_netcdf(out); print("saved", out, dict(da.sizes), la[[0, -1]], lo[[0, -1]])
