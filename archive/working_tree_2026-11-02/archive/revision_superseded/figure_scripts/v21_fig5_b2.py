"""Figure 5 (v2.1, box-free version): jet-path distortion and EKE around the ACC (b2_pathfollowing.py / b2_tests.py).
(a) Mean EKE 2001-2025 (log scale), the mean jet-core path (black; climatological core-level contour of relative ADT,
    central 10 deg of each 30 deg window), 3000 m isobath (grey); hatched: longitudes where the path is not tracked
    in >= 20 % of months.
(b) For each 30 deg window (centre on the x axis): correlation (best of lags 0-3, sign kept) between the monthly
    path length and EKE within 150 km of the path in the same window (black), and EKE within 150 km of the path
    15-35 deg further east (red). Filled: p < 0.05 (phase surrogates, max over lags).
Usage: OSR11_CACHE=<cache> python v21_fig5_b2.py <revision/data dir> <repo outputs dir>
"""
import os, sys
import numpy as np, pandas as pd, xarray as xr, contourpy
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.colors import LogNorm
from figstyle import *  # noqa
sys.path.insert(0, str(HERE)); import b1_meander as B  # noqa

data, repo = sys.argv[1], sys.argv[2]
CACHE = Path(os.environ.get("OSR11_CACHE", "/home/claude/w/cache"))
lat = np.load(CACHE / "lat.npy"); lon = np.load(CACHE / "lon.npy")
eke = np.load(CACHE / "eke.npy", mmap_mode="r")
E = np.zeros((lat.size, lon.size)); n = 0
for k in range(0, 300, 10):
    blk = np.asarray(eke[k:k + 10], np.float32); E += np.nansum(blk, 0); n += np.isfinite(blk).sum(0)
E = np.where(n > 0, E / np.maximum(n, 1), np.nan) * 1e4
ice = xr.open_dataset(f"{repo}/sea_ice_mask.nc").ice_mask.values > 0.5; E[ice] = np.nan
mean_adt = xr.open_dataset(f"{repo}/mean_adt.nc").mean_adt.values
T = pd.read_csv(OUT / "b2" / "b2_window_tests_v2.csv"); T0 = pd.read_csv(OUT / "b2" / "b2_window_tests.csv")
core = dict(zip(T0.centre, T0.core))
bat = xr.open_dataset(f"{data}/bathymetry_so.nc").elevation
pc = ccrs.PlateCarree()

fig = plt.figure(figsize=(11.0, 7.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1], hspace=0.12)
ax = fig.add_subplot(gs[0], projection=pc)
ax.set_extent([-180, 180, -65, -36], crs=pc); ax.set_aspect("auto")
im = ax.pcolormesh(lon, lat, E, cmap="YlGnBu", norm=LogNorm(vmin=20, vmax=2000), shading="nearest", rasterized=True, transform=pc, zorder=0)
ax.contour(bat.longitude, bat.latitude, -bat, levels=[3000], colors="#6b7280", linewidths=0.35, transform=pc, zorder=1)
ax.add_feature(cfeature.LAND, facecolor="#000000", edgecolor="#555555", linewidth=0.2, zorder=4)
# mean path: central 10 deg of each window
for c in np.arange(-180, 180, 10):
    rel = ((lon - c + 180) % 360) - 180; j = np.nonzero(np.abs(rel) <= 15)[0]; j = j[np.argsort(rel[j])]; lonw = rel[j][np.argsort(rel[j])] + c
    band = np.isfinite(mean_adt[:, j]) & (mean_adt[:, j] >= -0.6) & (mean_adt[:, j] <= 0.2)
    rows = np.nonzero(band.any(1))[0]; la0, la1 = max(lat[rows.min()] - 2, -64.9), min(lat[rows.max()] + 2, -36.0)
    i = np.nonzero((lat >= la0) & (lat <= la1))[0]
    rel_adt = mean_adt[np.ix_(i, j)]; rel_adt = rel_adt - np.nanmean(rel_adt)
    line = B.spanning_contour(rel_adt, lat[i], lonw, core[c]) if np.isfinite(core.get(c, np.nan)) else None
    cov = float(T.loc[T.centre == c, "coverage"].iloc[0])
    if cov < 0.8:
        ax.add_patch(plt.Rectangle((c - 5, -65), 10, 29, transform=pc, facecolor="none", edgecolor="#7A8288", hatch="///", linewidth=0, zorder=3))
        continue
    if line is None:
        continue
    keep = np.abs(line[:, 0] - c) <= 5
    L = line[keep]; L[:, 0] = ((L[:, 0] + 180) % 360) - 180
    brk = np.r_[True, np.abs(np.diff(L[:, 0])) < 180]
    ax.plot(np.where(brk, L[:, 0], np.nan), L[:, 1], color="#111111", lw=1.4, transform=pc, zorder=3)
gl = ax.gridlines(crs=pc, draw_labels=True, linewidth=0.3, color="#7A8288", alpha=0.45, linestyle=":", xlocs=np.arange(-180, 181, 30), ylocs=[-60, -50, -40])
gl.top_labels = False; gl.right_labels = False; gl.bottom_labels = False
gl.ylabel_style = {"size": 8, "color": "#4A4A4A"}
for name, x in (("Drake", -62), ("Kerguelen", 80), ("SEIR", 150), ("PAR", -120)):
    ax.text(x, -37.6, name, transform=pc, ha="center", va="center", fontsize=7.5, zorder=8, bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.8))
cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.03); cb.set_label(r"Mean EKE 2001–2025 (cm$^2$ s$^{-2}$)", fontsize=9); cb.ax.tick_params(labelsize=8)
ax.set_title("a) Mean EKE and jet-core path (hatched: path not tracked)", loc="left", fontsize=11, fontweight="bold")

b = fig.add_subplot(gs[1])
pos = ax.get_position(); pb = b.get_position(); b.set_position([pos.x0, pb.y0, pos.width, pb.height])
for link, lab, col in (("E", "EKE within 150 km of the path, same window", "#1F2933"), ("D", "EKE within 150 km of the path, 15–35° downstream", "#C0392B")):
    rr = T[[f"r_{link}_lag{k}" for k in range(4)]].values
    best = rr[np.arange(len(T)), np.nan_to_num(np.abs(rr)).argmax(1)]
    p = T[f"p_{link}"].values; ok = p < 0.05
    b.plot(T.centre, best, "-", color=col, lw=1.5, label=lab, zorder=2)
    b.plot(T.centre[ok], best[ok], "o", color=col, ms=5.5, zorder=3); b.plot(T.centre[~ok], best[~ok], "o", mfc="white", color=col, ms=5.5, zorder=3)
for c0 in T.centre[T.coverage < 0.8]:
    b.axvspan(c0 - 5, c0 + 5, color="#E9ECEF", lw=0, zorder=0)
b.axhline(0, color="#495057", lw=0.8); b.set_xlim(-180, 180); b.set_ylim(-0.3, 0.7)
b.set_xticks(np.arange(-180, 181, 30)); b.set_xticklabels([f"{abs(v)}°{'W' if v < 0 else 'E' if v > 0 else ''}" if abs(v) != 180 else "180°" for v in np.arange(-180, 181, 30)])
b.set_xlabel("Window centre"); b.set_ylabel("r (path length → EKE, best of lags 0–3)")
b.legend(fontsize=8, frameon=False, loc="upper left"); base_timeseries_style(b)
b.text(0.99, 0.96, "filled: p < 0.05", transform=b.transAxes, ha="right", va="top", fontsize=8, color="#4A4A4A")
b.set_title("b) Path length and EKE, 30° windows, 2001–2025", loc="left", fontsize=11, fontweight="bold")
fig.savefig(OUT / "figures_v21" / "Figure5.png"); fig.savefig(OUT / "figures_v21" / "Figure5.pdf")
print("saved")
