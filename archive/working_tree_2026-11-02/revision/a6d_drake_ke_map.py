"""Mean KE of the monthly-mean surface geostrophic flow in Drake Passage, 2001-2025 (Product 1; months with ice_frac > 0.10
excluded per cell) -> outputs/a6/a6d_drake_mean_ke.nc (background of Figure 3a).
Usage: python a6d_drake_ke_map.py <revision/data dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, xarray as xr
data, out = Path(sys.argv[1]), Path(sys.argv[2])
S = None; N = None
for yr in range(2001, 2026):
    for mo in range(1, 13):
        with xr.open_dataset(data / "monthly_allsat" / f"{yr}" / f"{mo:02d}.nc") as ds:
            d = ds.sel(longitude=slice(-76, -48), latitude=slice(-66, -52)).load()
        ke = 0.5 * (d.u_mean ** 2 + d.v_mean ** 2) * 1e4
        ke = ke.where(d.ice_frac.fillna(0) <= 0.10)
        S = ke.fillna(0) if S is None else S + ke.fillna(0)
        N = ke.notnull().astype(int) if N is None else N + ke.notnull().astype(int)
M = (S / N).where(N >= 60).rename("ke_mean"); M.attrs["units"] = "cm2 s-2"
M.to_dataset().to_netcdf(out / "a6d_drake_mean_ke.nc"); print("saved", float(M.max()))
