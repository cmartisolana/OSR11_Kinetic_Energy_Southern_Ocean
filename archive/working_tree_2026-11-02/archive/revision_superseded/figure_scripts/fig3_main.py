"""Main Fig. 3 (v1 style): standing meanders from altimetry, 2001-2025.
Left: monthly path occurrence (fraction of months the core-level contour of ADT relative to its
regional mean crosses each 0.25° cell) with epoch-mean paths, over depth < 3000 m (grey) and land
(black); dashed lines: upstream / meander / downstream section limits.
(The path-length / downstream-EKE lag panels moved to Fig. 6, fig6_meander_eke.py, 2026-10-13.)
Usage: python fig3_main.py <revision/data dir>
"""
import sys
import numpy as np, pandas as pd, xarray as xr
import contourpy
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

data = sys.argv[1]
MEAN = [("PAR", "Pacific–Antarctic Ridge", (-156, -74, -64, -44), (-150, -135.13, -109.38, -80)),
        ("SEIR150", "Southeast Indian Ridge (150°E)", (122, 166, -64, -44), (125, 138, 140, 150, 158)),
        ("KERG", "Kerguelen Plateau (Fawn Trough)", (64, 100, -62, -42), (70, 77, 76, 84, 94))]
EP_LAB = ["2001–2008", "2009–2016", "2017–2025"]
pc_ = ccrs.PlateCarree()

fig = plt.figure(figsize=(8.4, 11.0))
gs = fig.add_gridspec(3, 1, hspace=0.36)
for r, (name, title, (lo0, lo1, la0, la1), secs) in enumerate(MEAN):
    ax = fig.add_subplot(gs[r, 0], projection=pc_)
    ax.set_extent([lo0, lo1, la0, la1], crs=pc_); ax.set_aspect("auto")
    O = xr.open_dataset(OUT / "b1" / f"occ_{name}.nc")
    bat = xr.open_dataset(f"{data}/bathymetry_so.nc").elevation.sel(longitude=slice(lo0, lo1), latitude=slice(la0, la1))
    ax.contourf(bat.longitude, bat.latitude, -bat, levels=[0, 3000, 9000], colors=["#D5DADF", "white"], transform=pc_, zorder=0)
    oc = np.where(O.occ > 0, O.occ, np.nan)
    im = ax.pcolormesh(O.glon, O.glat, oc, cmap="YlGnBu", vmin=0, vmax=0.6, shading="nearest", rasterized=True,
                       transform=pc_, zorder=1)
    ax.add_feature(cfeature.LAND, facecolor="#000000", edgecolor="#555555", linewidth=0.2, zorder=4)
    for e, col in enumerate(C_EPOCH):
        gen = contourpy.contour_generator(O.longitude.values, O.latitude.values,
                                          np.ma.masked_invalid(O.epoch_rel_adt.values[e]), name="serial")
        for k, L in enumerate([L for L in gen.lines(O.attrs["level"]) if np.ptp(L[:, 0]) > 2]):
            ax.plot(L[:, 0], L[:, 1], color=col, lw=1.4, transform=pc_, zorder=3, label=EP_LAB[e] if k == 0 else None)
    for xv in secs:
        ax.plot([xv, xv], [la0, la1], color="#222222", lw=0.8, ls=(0, (5, 3)), alpha=0.75, transform=pc_, zorder=3)
    gl = ax.gridlines(crs=pc_, draw_labels=True, linewidth=0.3, color="#7A8288", alpha=0.45, linestyle=":")
    gl.top_labels = False; gl.right_labels = False
    gl.xlabel_style = {"size": 8, "color": "#4A4A4A"}; gl.ylabel_style = {"size": 8, "color": "#4A4A4A"}
    ax.set_title(title, loc="left"); panel_label(ax, "abc"[r], x=-0.09, y=1.13)
    if r == 0:
        ax.legend(loc="lower left", fontsize=8, frameon=True, framealpha=0.9, edgecolor="#cccccc")
cax = fig.add_axes([0.2, 0.065, 0.6, 0.012])
cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label("Path occurrence (fraction of months)", fontsize=9)
cb.ax.tick_params(labelsize=8)
fig.savefig(FIG / "Fig3_standing_meanders.png"); fig.savefig(FIG / "Fig3_standing_meanders.pdf")
print("saved")
