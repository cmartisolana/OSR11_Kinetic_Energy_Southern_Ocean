"""2001-2025 mean of the monthly fields (written for the publication repository; no generator existed before).

Purpose : outputs/base/clim_allsat.nc = unweighted mean over the 300 monthly files of adt_mean, u_mean and v_mean (each month
          counts once; NaN where a pixel is never valid). Used for the mean SAF/PF streamlines (dk8, dk12, dk34), the maps
          and the wind domains. Reproduces the working-repository file to float32 rounding (max |difference| 1.4e-6 m for
          ADT, 8.4e-7 m s-1 for u; identical NaN pattern), checked on 2 Nov 2026.
Inputs  : $OSR11_DATA/monthly_allsat/YYYY/MM.nc (00a)
Outputs : outputs/base/clim_allsat.nc (adt, u, v; float32)
Run     : python scripts/01_base/01d_climatology.py [--check]   (--check: compare with the existing file instead of writing)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import argparse

import numpy as np
import xarray as xr

from osr11 import paths as PA

ap = argparse.ArgumentParser(); ap.add_argument("--check", action="store_true"); a = ap.parse_args()
fs = PA.monthly_files("allsat")
S, N = {}, {}
for f in fs:
    with xr.open_dataset(f) as d:
        for k, v in (("adt", "adt_mean"), ("u", "u_mean"), ("v", "v_mean")):
            x = d[v].values.astype(np.float64); m = np.isfinite(x)
            if k not in S:
                S[k] = np.zeros(x.shape); N[k] = np.zeros(x.shape); lat, lon = d.latitude.values, d.longitude.values
            S[k][m] += x[m]; N[k][m] += 1
with np.errstate(invalid="ignore", divide="ignore"):
    M = {k: (S[k] / N[k]).astype(np.float32) for k in S}
ds = xr.Dataset({k: (("latitude", "longitude"), M[k]) for k in ("adt", "u", "v")}, coords={"latitude": lat, "longitude": lon})
if a.check:
    old = xr.open_dataset(PA.CLIM)
    for k in ("adt", "u", "v"):
        o, n = old[k].values.astype(float), ds[k].values.astype(float)
        print(k, "NaN pattern identical:", bool((np.isnan(o) == np.isnan(n)).all()), "max |diff|:", float(np.nanmax(np.abs(o - n))))
else:
    PA.BASE.mkdir(parents=True, exist_ok=True); ds.to_netcdf(PA.CLIM); print("saved", PA.CLIM)
