"""Subset a global bathymetry grid to the Southern Ocean for isobaths (A5d) and
meander-box selection (B1). Bathymetry is used only as context, not as ocean data.

1. Download a global bedrock-elevation NetCDF yourself, e.g. ETOPO 2022
   (60 arc-second, bedrock) from NOAA NCEI, or a GEBCO grid.
2. Run:
       python revision/prepare_bathymetry.py /path/to/global_file.nc
   -> revision/data/bathymetry_so.nc  (elevation in m, 70S-30S, 1/8 deg)

Purpose : subset and block-average a global bathymetry grid (GEBCO; ETOPO also works) to the Southern Ocean.
Inputs  : a global bedrock-elevation NetCDF downloaded by the user (GEBCO grid, https://www.gebco.net).
Outputs : outputs/base/bathymetry_so.nc (elevation, m, 70-30 S, ~1/8 deg). Distributed in outputs/base/.
Run     : python scripts/00_extract/00c_prepare_bathymetry.py /path/to/GEBCO_global.nc
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np
import xarray as xr

from osr11 import paths as PA  # noqa: E402
src = Path(sys.argv[1])
ds = xr.open_dataset(src)
var = next(v for v in ("z", "elevation", "Band1") if v in ds)
lat = "lat" if "lat" in ds.coords else "latitude"
lon = "lon" if "lon" in ds.coords else "longitude"
da = ds[var].sel({lat: slice(-70, -30)})
# block-average to ~1/8 deg
step = max(1, int(round(0.125 / float(abs(da[lat][1] - da[lat][0])))))
da = da.coarsen({lat: step, lon: step}, boundary="trim").mean()
da = da.rename({lat: "latitude", lon: "longitude"}).rename("elevation")
out = PA.BATHY
out.parent.mkdir(parents=True, exist_ok=True)
da.astype(np.float32).to_netcdf(out)
print("saved", out, dict(da.sizes))
