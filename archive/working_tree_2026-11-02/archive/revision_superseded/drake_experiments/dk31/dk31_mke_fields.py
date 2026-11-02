"""dk31a: annual kinetic energy of the mean flow (MKE) on the 1 deg grid of dk22, Product 1, 2001-2025.
MKE_y = 0.5 |U_y|^2 with U_y the annual mean of the monthly u_mean, v_mean at 1/8 deg, then block-averaged to 1 deg with
the same valid-cell rule as dk22 (all 64 pixels valid in every month and outside the sea-ice mask). Units cm2 s-2.
Usage: python dk31_mke_fields.py <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, xarray as xr
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
files = {y: [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for m in range(1, 13)] for y in range(2001, 2026)}
with xr.open_dataset(files[2001][0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
ny, nx = lat.size // 8, lon.size // 8
blk = lambda a: a.reshape(ny, 8, nx, 8)
valid = ~ice; UY, VY = [], []
for y, fl in files.items():
    su = np.zeros((lat.size, lon.size)); sv = np.zeros_like(su)
    for f in fl:
        with xr.open_dataset(f) as d: u, v = d.u_mean.values.astype(float), d.v_mean.values.astype(float)
        valid &= np.isfinite(u) & np.isfinite(v); su += np.nan_to_num(u); sv += np.nan_to_num(v)
    UY.append(su / 12); VY.append(sv / 12)
cellok = blk(valid).all(axis=(1, 3))
MKE = np.array([np.where(cellok, np.nanmean(blk(np.where(valid, 0.5 * (u ** 2 + v ** 2) * 1e4, np.nan)), axis=(1, 3)), np.nan) for u, v in zip(UY, VY)])
np.savez(OUT / "dk31_mke_fields.npz", MKE=MKE, lat=lat.reshape(ny, 8).mean(1), lon=lon.reshape(nx, 8).mean(1))
print(MKE.shape, np.nanmean(MKE))
