"""Climatological context maps used to fix the B1 sections for further meanders.

Only time-mean fields (2001-2025): mean speed |ubar|, mean EKE, mean ADT
contours and GEBCO depth. No time series are computed here.
Usage: python b1_regions.py <revision/data dir> <lon0> <lon1> <lat0> <lat1> <tag>
-> outputs/b1/regions_<tag>.png and outputs/b1/clim_<tag>.nc
"""
import sys
from pathlib import Path
import numpy as np, xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
data = Path(sys.argv[1]); lo0, lo1, la0, la1 = map(float, sys.argv[2:6]); tag = sys.argv[6]
su = sv = sa = s2 = None; U = []
for y in range(2001, 2026):
    for m in range(1, 13):
        with xr.open_dataset(data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc") as ds:
            d = ds.sel(longitude=slice(lo0, lo1), latitude=slice(la0, la1)).load()
        U.append((d.u_mean.values, d.v_mean.values)); sa = d.adt_mean if sa is None else sa + d.adt_mean
        s2 = d.var_intra if s2 is None else s2 + d.var_intra
uu = np.array([a for a, b in U]); vv = np.array([b for a, b in U])
ub, vb = np.nanmean(uu, 0), np.nanmean(vv, 0)
eke = 0.5 * np.nanmean((uu - ub) ** 2 + (vv - vb) ** 2, 0) + 0.5 * (s2 / 300).values
lat, lon = d.latitude.values, d.longitude.values
C = xr.Dataset({"spd": (("latitude", "longitude"), np.hypot(ub, vb)), "eke": (("latitude", "longitude"), eke * 1e4),
                "adt": (("latitude", "longitude"), (sa / 300).values), "u": (("latitude", "longitude"), ub),
                "v": (("latitude", "longitude"), vb)}, coords=dict(latitude=lat, longitude=lon))
C.to_netcdf(HERE / "outputs" / "b1" / f"clim_{tag}.nc")
bat = xr.open_dataset(data / "bathymetry_so.nc").elevation.sel(longitude=slice(lo0, lo1), latitude=slice(la0, la1))
fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
for ax, (fld, cm, vmax, lab) in zip(axs, [(C.spd, "viridis", 0.4, "mean speed (m s⁻¹)"),
                                          (C.eke, "magma", 600, "mean EKE (cm² s⁻²)"), (-bat, "Blues", 5500, "depth (m)")]):
    pc = ax.pcolormesh(fld.longitude, fld.latitude, fld, cmap=cm, vmin=0, vmax=vmax, shading="auto")
    ax.contour(lon, lat, C.adt, levels=np.arange(-1.2, 0.8, 0.1), colors="w" if ax is not axs[2] else "k", linewidths=0.4)
    ax.contour(bat.longitude, bat.latitude, bat, levels=[-3000], colors="r", linewidths=0.6)
    plt.colorbar(pc, ax=ax, label=lab, pad=0.01)
    ax.set_xticks(np.arange(np.ceil(lo0 / 5) * 5, lo1 + 1, 5)); ax.grid(alpha=0.3, lw=0.4)
fig.tight_layout(); fig.savefig(HERE / "outputs" / "b1" / f"regions_{tag}.png", dpi=130)
print("saved")
