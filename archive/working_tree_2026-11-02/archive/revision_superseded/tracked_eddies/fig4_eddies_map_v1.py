"""Candidate Figure 4 (v2.1, map version, draft 1): tracked-eddy distributions and trend maps.
a-d  Distributions (density, log scale) of amplitude, speed radius, rotational speed (all ACC-band observations) and
     lifetime (tracks born in the band), 2001-2011 (grey) and 2012-2023 (green); dashed: medians; below: late / early.
e-h  Trends 2001-2022 on 2 deg cells, % of the cell mean per decade: mean amplitude, mean rotational speed,
     sum(amplitude^2) per day, eddies present per day. Dots: FDR 10 % (draft AR(1)-effective-n test, see e7_maps.py).
     ACC band (mean ADT -0.6 and +0.2 m, black), 3000 m isobath (grey), sea ice (light blue).
Usage: python fig4_eddies_map_v1.py <e6 out dir> <e7 out dir> <repo outputs dir> <bathymetry nc> <revision dir> <figure dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

o6, o7, repo, bathy, rev, figdir = (Path(a) for a in sys.argv[1:7]); figdir.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(rev))
import style_nc as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import matplotlib.path as mpath  # noqa: E402
import cartopy.crs as ccrs  # noqa: E402

D = pd.read_csv(o6 / "e6_distributions.csv"); M = xr.open_dataset(o7 / "e7_cell_trends.nc")
adt = xr.open_dataset(repo / "mean_adt.nc").mean_adt; ice = xr.open_dataset(repo / "sea_ice_mask.nc").ice_mask
bat = xr.open_dataset(bathy).elevation
C_E, C_L = "#9A9A9A", st.C["te"]

fig = plt.figure(figsize=(st.W2, 185 * st.MM))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 2.3], hspace=0.35)
gtop = gs[0].subgridspec(2, 4, height_ratios=[3, 1], hspace=0.08, wspace=0.45)
panels = [("amplitude", "amplitude (cm)", 100, "a", (0, 35)), ("speed_radius", "speed radius (km)", 1e-3, "b", (20, 120)),
          ("speed_average", "rotational speed (cm s$^{-1}$)", 100, "c", (0, 70)), ("lifetime_d", "lifetime (d)", 1, "d", (0, 250))]
for k, (v, lab, sc, let, xl) in enumerate(panels):
    ax = fig.add_subplot(gtop[0, k]); axr = fig.add_subplot(gtop[1, k], sharex=ax)
    e = D[(D["var"] == v) & (D.period == "early") & D.left.notna()]; l = D[(D["var"] == v) & (D.period == "late") & D.left.notna()]
    x = (e.left.values + e.right.values) / 2 * sc
    ax.plot(x, e.density.values / sc, color=C_E, lw=0.9); ax.plot(x, l.density.values / sc, color=C_L, lw=0.9)
    for per, c in (("early", C_E), ("late", C_L)):
        ax.axvline(D[(D["var"] == v) & (D.period == per) & D.left.isna()].iloc[0]["median"] * sc, color=c, lw=0.6, ls="--")
    ax.set_yscale("log"); ax.tick_params(labelbottom=False); st.plabel(ax, let, lab)
    if k == 0: ax.set_ylabel("density")
    with np.errstate(divide="ignore", invalid="ignore"):
        r = l.density.values / e.density.values
    axr.plot(x, r, color=st.C["ink"], lw=0.8); axr.axhline(1, color="#9A9A9A", lw=0.5)
    axr.set_ylim(0.6, 1.6); axr.set_yticks([0.8, 1.0, 1.2, 1.4]); axr.set_xlabel(lab); ax.set_xlim(*xl)
    if k == 0: axr.set_ylabel("late / early")
fig.axes[0].legend(handles=[Line2D([], [], color=C_E, label="2001–2011"), Line2D([], [], color=C_L, label="2012–2023"),
                            Line2D([], [], color=st.C["ink"], ls="--", lw=0.6, label="median")], loc="upper right")

gmap = gs[1].subgridspec(2, 2, hspace=0.25, wspace=0.02)
maps = [("amp", "mean amplitude", "e", 15), ("spd", "mean rotational speed", "f", 15),
        ("amp2_per_day", "Σ amplitude² per day", "g", 30), ("present_per_day", "eddies present per day", "h", 15)]
proj = ccrs.SouthPolarStereo()
LON, LAT = np.meshgrid(M.lon.values, M.lat.values)
icec = ice.coarsen(latitude=int(round((M.lat.values[1]-M.lat.values[0])/float(ice.latitude[1]-ice.latitude[0]))), longitude=int(round((M.lon.values[1]-M.lon.values[0])/float(ice.longitude[1]-ice.longitude[0]))), boundary="trim").max()
icec = icec.interp(latitude=M.lat.values, longitude=M.lon.values, method="nearest").values > 0.5
for k, (v, lab, let, vmax) in enumerate(maps):
    ax = fig.add_subplot(gmap[k // 2, k % 2], projection=proj)
    ax.set_extent([-180, 180, -80, -35], crs=ccrs.PlateCarree())
    theta = np.linspace(0, 2 * np.pi, 200); verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    ax.set_boundary(mpath.Path(verts * 0.5 + 0.5), transform=ax.transAxes)
    pc = ax.pcolormesh(M.lon.values, M.lat.values, np.where(icec, np.nan, M[f"{v}_pct"].values), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       transform=ccrs.PlateCarree(), shading="nearest", zorder=1)
    f = M[f"{v}_fdr"].values.astype(bool) & ~icec
    ax.plot(LON[f], LAT[f], ".", color=st.C["ink"], ms=1.2, transform=ccrs.PlateCarree(), zorder=3, lw=0)
    ax.contour(adt.longitude, adt.latitude, adt, levels=[-0.6, 0.2], colors="#222222", linewidths=0.5, transform=ccrs.PlateCarree(), zorder=2)
    b = bat.coarsen(longitude=4, latitude=4, boundary="trim").mean()
    ax.contour(b.longitude, b.latitude, b, levels=[-3000], colors="#8A8A8A", linewidths=0.35, transform=ccrs.PlateCarree(), zorder=2)
    ax.contourf(ice.longitude, ice.latitude, np.where(ice.values > 0.5, 1.0, np.nan), levels=[0.5, 1.5], colors=[st.C["ice"]], transform=ccrs.PlateCarree(), zorder=2)
    ax.contourf(b.longitude, b.latitude, np.where(b.values >= 0, 1.0, np.nan), levels=[0.5, 1.5], colors=["#3A3A3A"], transform=ccrs.PlateCarree(), zorder=4)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.3, color="#9A9A9A", linestyle=":", alpha=0.6)
    gl.xlocator = plt.FixedLocator(range(-180, 181, 60)); gl.ylocator = plt.FixedLocator([-60, -40])
    ax.set_frame_on(False); st.plabel(ax, let, lab, dy=-2)
    cb = fig.colorbar(pc, ax=ax, orientation="horizontal", fraction=0.035, pad=0.01, aspect=35, shrink=0.7)
    cb.set_label("% of cell mean per decade"); cb.outline.set_linewidth(0.4)
fig.savefig(figdir / "Figure4_eddies_map_v1.png"); fig.savefig(figdir / "Figure4_eddies_map_v1.pdf")
print("saved")
