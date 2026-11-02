"""Fig. 2g-h (v3.42): isopycnal tilt in the Drake Passage region from Argo (Roemmich-Gilson, 1 deg, annual 2004-2025; dk32),
as maps with the extent and style of panels a-b.
(g) trend of the isopycnal tilt |grad z| of sigma0 = 27.5 (% of its 2004-2025 mean per decade; cells with mean tilt
    >= 20 m per 100 km); dots: p < 0.05 (n_eff from lag-1 autocorrelation); mean SAF and PF (altimetry); LM Gould line.
(h) mean depth of sigma0 = 27.5 with its 500 and 1000 dbar contours for 2004-2014 and 2015-2025; inset: trend of the
    cross-passage isopycnal slope along the LM Gould line (dk33).
Drawn as one row of the Drake-figure width; appended by fig2_append_v341.py.
Usage: python fig2gh_v342.py <revision dir> <drake out dir> <work out dir> <dk24 section npz> <out pdf>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
from scipy.ndimage import gaussian_filter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
REV, DO, WO, SEC, OUTF = (Path(a) for a in sys.argv[1:6])
plt.rcParams.update({   # style_nc.py, as fig2ef
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.5, "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 5.4,
    "legend.frameon": False, "legend.handlelength": 1.6, "legend.handletextpad": 0.4, "legend.columnspacing": 0.9,
    "legend.borderaxespad": 0.2, "axes.linewidth": 0.5, "axes.edgecolor": "#222222", "axes.labelcolor": "#222222",
    "text.color": "#222222", "xtick.color": "#222222", "ytick.color": "#222222", "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "pdf.fonttype": 42,
    "mathtext.default": "regular", "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans"})
PC = ccrs.PlateCarree()
I = np.load(WO / "dk32" / "dk32_iso.npz"); D33 = pd.read_csv(WO / "dk33" / "dk33_drake_isopycnals.csv").set_index("sigma0")
yrs, lat, lon = I["years"], I["lat"], I["lon"]; t = yrs - yrs.mean()
L = np.where(lon > 180, lon - 360, lon); sel = (L >= -90) & (L <= -35); o = np.argsort(L[sel]); LO = L[sel][o]
Z = I["p27.5"][:, :, sel][:, :, o]                       # (year, lat, lon) pressure of sigma0 = 27.5 (dbar ~ m)
dy = np.gradient(Z, lat, axis=1) / 111.195 * 100; dx = np.gradient(Z, LO, axis=2) / (111.195 * np.cos(np.radians(lat))[None, :, None]) * 100
G = np.hypot(dx, dy)                                     # tilt, m per 100 km
def trend(y):
    tt = t[:, None, None]; m = y.mean(0); b = (tt * (y - m)).sum(0) / (t ** 2).sum(); r = y - m - tt * b
    r1 = (r[1:] * r[:-1]).sum(0) / np.maximum((r * r).sum(0), 1e-12); ne = np.clip(len(t) * (1 - r1) / (1 + r1), 3, len(t))
    se = np.sqrt((r ** 2).sum(0) / (ne - 2) / (t ** 2).sum()); return 10 * b, 2 * stats.t.sf(np.abs(b / se), ne - 2)
bG, pG = trend(G); Gm = G.mean(0); rel = np.where(Gm >= 20, 100 * bG / Gm, np.nan); pG = np.where(Gm >= 20, pG, np.nan)
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc").sel(longitude=slice(-90, -35), latitude=slice(-66, -48))
Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation.sel(longitude=slice(-90, -35), latitude=slice(-66, -48)).coarsen(latitude=2, longitude=2, boundary="trim").mean()
TR = np.load(SEC); tlon, tlat = TR["lon"], TR["lat"]
PW, PH = 471.0592603718, 160.0
fig = plt.figure(figsize=(PW / 72, PH / 72))
EXT = [-86, -40, -65, -50]
def mapax(x0, x1, y0, y1, lab, ttl, ylab=True):
    ax = fig.add_axes([x0 / PW, y0 / PH, (x1 - x0) / PW, (y1 - y0) / PH], projection=PC)
    ax.set_extent(EXT, crs=PC); ax.set_aspect("auto")
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="#333333", edgecolor="none", zorder=5)
    ax.contour(bat.longitude, bat.latitude, -bat, levels=[2000], colors="#8A8A8A", linewidths=0.35, transform=PC, zorder=2)
    ax.set_xticks([-80, -70, -60, -50, -40], crs=PC); ax.set_xticklabels(["80°W", "70°W", "60°W", "50°W", "40°W"])
    ax.set_yticks([-65, -60, -55, -50], crs=PC); ax.set_yticklabels(["65°S", "60°S", "55°S", "50°S"] if ylab else [])
    ax.text(0.0, 1.03, f"$\\bf{{{lab}}}$  {ttl}", transform=ax.transAxes, fontsize=7.2, ha="left", va="bottom")
    return ax
def cbar(im, ax, label, ext):
    p = ax.get_position(); cax = fig.add_axes([p.x0 + 0.45 * p.width, p.y0 - 0.20, 0.5 * p.width, 0.035])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend=ext); cb.ax.tick_params(labelsize=6, width=0.4, length=2); cb.outline.set_linewidth(0.4)
    fig.text(p.x0 + 0.44 * p.width, p.y0 - 0.20 + 0.017, label, ha="right", va="center", fontsize=6.3)
Y0, Y1 = 50, 147
# g
ax = mapax(47, 214, Y0, Y1, "g", "Isopycnal tilt trend (σ$_0$ = 27.5), 2004–2025")
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
im = ax.pcolormesh(LO, lat, np.ma.masked_invalid(rel), transform=PC, cmap=cm_, vmin=-20, vmax=20, shading="nearest", rasterized=True, zorder=0)
jj, ii = np.nonzero(pG < 0.05); ax.scatter(LO[ii], lat[jj], s=2.0, c="k", lw=0, transform=PC, zorder=7)
ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="#1A1A1A", linewidths=1.0, linestyles="solid", transform=PC, zorder=6)
for x_, y_, s_ in ((-85, -55.6, "SAF"), (-85, -60.6, "PF")): ax.text(x_, y_, s_, transform=PC, fontsize=6.5, fontweight="bold", zorder=9, bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.75))
ax.plot(tlon, tlat, color="#E69F00", lw=1.4, transform=PC, zorder=8)
cbar(im, ax, "Tilt trend\n(% per decade)", "both")
# h
ax = mapax(251, 428, Y0, Y1, "h", "Depth of σ$_0$ = 27.5 kg m$^{-3}$, two periods", ylab=False)
pm = np.nanmean(Z, 0)
imd = ax.pcolormesh(LO, lat, np.ma.masked_invalid(pm), transform=PC, cmap="magma_r", vmin=0, vmax=2000, shading="nearest", rasterized=True, zorder=0)
CE, CL = "#00B4E6", "#111111"
for msk, col, ls, lw in ((yrs <= 2014, CE, "-", 1.1), (yrs >= 2015, CL, (0, (2.5, 1.2)), 0.6)):
    f = np.nanmean(Z[msk], 0); f = gaussian_filter(np.where(np.isfinite(f), f, 2500), 0.5)
    ax.contour(LO, lat, f, levels=[500, 1000], colors=col, linewidths=lw, linestyles=[ls], transform=PC, zorder=6)
ax.plot(tlon, tlat, color="#E69F00", lw=1.4, transform=PC, zorder=8)
ax.legend(handles=[Line2D([], [], color=CE, lw=1.1, label="2004–2014"), Line2D([], [], color=CL, lw=0.7, ls=(0, (2.5, 1.2)), label="2015–2025")],
          title="500 and 1000 dbar", title_fontsize=5.6, loc="upper left", fontsize=5.6, facecolor="white", framealpha=0.85, frameon=True, edgecolor="none").set_zorder(10)
r = D33.loc[27.5]
ax.text(0.985, 0.04, f"slope along the LM Gould line\n{r.slope_trend_pct:+.1f} % per decade ({r.slope_trend_pct - r.slope_ci_pct:+.1f} to {r.slope_trend_pct + r.slope_ci_pct:+.1f})".replace("-", "−"),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4, zorder=12, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
cbar(imd, ax, "Mean depth\n(dbar)", "max")
fig.savefig(OUTF); print("saved", OUTF, "| median tilt trend in band-ish cells:", np.nanmedian(rel))
