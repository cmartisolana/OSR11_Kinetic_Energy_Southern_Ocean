"""SI figure: mean jet-core path (+0.18 m relative-ADT contour, 150 km smoothing, as B3-B5) at the five meanders,
2001-2012 (dashed) and 2013-2025 (solid), over the 2001-2025 mean surface speed and the 3000 m isobath; upstream / ridge /
downstream regions of dk13 marked along the top. Separation of the two paths and random-split null from b5b_path_change.csv.
Usage: python figS_meander_paths_v326.py <revision dir> <b5 out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, B5, FIGD = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV))
import style_nc as S
sys.argv = [sys.argv[0], str(REV), "x", "x", str(B5)]
src = open(Path(__file__).parent / "b5_meander_shape.py").read()
exec(src[src.index("MEANDERS = {"):src.index("files = [")])                      # MEANDERS
exec(src[src.index("files = ["):src.index("af = OUT")].replace("OUT", "B5"))       # files, lat, lon
Zf = np.load(B5 / "b5_annual_fields.npz"); A = Zf["A"]
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
ACCM = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice
EXCL = (mean_adt > 0.45) | (lat[:, None] > -40) | ~np.isfinite(mean_adt)
import contourpy
from scipy.ndimage import gaussian_filter1d
exec(src[src.index("PAD = int(30"):src.index("def amp(")])                      # path()
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); spd = np.hypot(clim.u.values, clim.v.values) * 100
bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation
AM = pd.read_csv(B5 / "b5b_path_change.csv").set_index("meander")
Pe, Pl = path(np.nanmean(A[:12], 0)), path(np.nanmean(A[12:], 0))
plt = S.plt; C = S.C
fig, axs = plt.subplots(3, 2, figsize=(S.W2, 0.95 * S.W2)); axs = axs.ravel(); axs[-1].axis("off"); fig.subplots_adjust(hspace=0.62, wspace=0.18)
for ax, (nm, (u_, r_, d_)) in zip(axs[:5], MEANDERS.items()):
    a, b = u_[0] - 5, d_[1] + 5; mj = (lon >= a) & (lon <= b)
    lat_lo = -64; lat_hi = -40
    mi = (lat >= lat_lo) & (lat <= lat_hi)
    im = ax.pcolormesh(lon[mj], lat[mi], spd[np.ix_(mi, mj)], cmap="Greys", vmin=0, vmax=40, shading="auto", rasterized=True)
    bs = bat.sel(longitude=slice(a, b), latitude=slice(lat_lo, lat_hi))
    ax.contour(bs.longitude, bs.latitude, -bs.values, levels=[3000], colors="#B0763B", linewidths=0.4)
    for P, ls, lab in ((Pe, (0, (3, 1.5)), "2001–2012"), (Pl, "-", "2013–2025")):
        for lo_s, la_s in P:
            s = (lo_s >= a) & (lo_s <= b)
            if s.sum() > 2: ax.plot(lo_s[s], la_s[s], color=C["tr"], lw=0.9, ls=ls)
    for (x0, x1), col in ((u_, "#56B4E9"), (r_, "#555555"), (d_, "#E69F00")):
        ax.plot([x0, x1], [lat_hi - 0.6] * 2, color=col, lw=2.4, solid_capstyle="butt")
    q = AM.loc[nm]; ttl = f"{nm}\nmean separation {q.curve_dist_km:.0f} km (random splits: median {q.curve_null_median_km:.0f}, p = {q.curve_p:.2f})"
    S.plabel(ax, "abcde"[list(MEANDERS).index(nm)], ttl, dy=2); ax.tick_params(labelsize=5.5)
    ax.set_xlim(a, b); ax.set_ylim(lat_lo, lat_hi); ax.set_ylabel("Latitude (°)")
    ax.set_aspect(1 / np.cos(np.radians(52)), adjustable="box")
for ax in axs[3:5]: ax.set_xlabel("Longitude (°E)")
from matplotlib.lines import Line2D
axs[5].legend(handles=[Line2D([], [], color=C["tr"], ls=(0, (3, 1.5)), lw=0.9, label="mean path 2001–2012"), Line2D([], [], color=C["tr"], lw=0.9, label="mean path 2013–2025"),
                       Line2D([], [], color="#B0763B", lw=0.6, label="3000 m isobath")], loc="center", fontsize=6, ncol=1, title="Jet-core path (+0.18 m relative ADT); grey: mean speed;\nbars: upstream / ridge / downstream", title_fontsize=5.5)
S.save(fig, FIGD, "FigureS_meander_paths"); print("done")
