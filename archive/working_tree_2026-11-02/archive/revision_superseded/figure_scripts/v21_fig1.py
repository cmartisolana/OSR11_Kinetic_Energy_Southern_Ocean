"""v2.1 Figure 1: same layout as v2 Fig. 1 (a-f); inputs and statistics revised.
(a, b) mean KE and EKE 2001-2025 (all-sat, 1/8°); (c, d) deseasonalised Theil-Sen trends on 1° cells,
dots = significant after Benjamini-Hochberg FDR (alpha 0.10) with the revised test (a1_maps_bc.nc);
(e) ACC-band KE and EKE (raw monthly, 12-month running mean, trend); (f) trends by region (ACC band),
95 % CI and revised p (retest_all.py). Grey contours: mean ADT -0.6 and +0.2 m (ACC band).
Usage: python v21_fig1.py <cache dir>
"""
import sys
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

cache = sys.argv[1]
FIGV3 = OUT / "figures_v21"; FIGV3.mkdir(exist_ok=True)
ROOT_OUT = Path("/mnt/user-data/uploads/repository/outputs")
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
ice = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc").ice_mask.values > 0.5
lat = np.load(f"{cache}/lat.npy"); lon = np.load(f"{cache}/lon.npy")
ke = np.load(f"{cache}/ke.npy", mmap_mode="r"); eke = np.load(f"{cache}/eke.npy", mmap_mode="r")
kem = np.nanmean(np.asarray(ke, np.float32), 0) * 1e4; ekem = np.nanmean(np.asarray(eke, np.float32), 0) * 1e4
kem[ice] = np.nan; ekem[ice] = np.nan
M = xr.open_dataset(OUT / "a1_maps_bc.nc")
S = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv")
T = pd.read_csv(OUT / "a5" / "a5_trends_allsat.csv").query("period == '2001-2025'").set_index("series")
T = retest_override(T, lambda s, r: f"allsat|{s}|2001-2025")


def title(ax, t):
    ax.set_title(t, loc="left", fontsize=12, fontweight="bold")


proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(10.8, 13.0))
outer = fig.add_gridspec(2, 1, height_ratios=[2, 0.8], hspace=0.14)
gm = outer[0].subgridspec(2, 2, hspace=0.16, wspace=0.05)
gb = outer[1].subgridspec(1, 2, width_ratios=[1.35, 1], wspace=0.45)
axs = [fig.add_subplot(gm[i, j], projection=proj) for i in range(2) for j in range(2)]
LON, LAT = np.meshgrid(lon, lat)
vmax_mean = 1500
for ax, fld, ttl in ((axs[0], kem, "a) Mean KE (2001–2025)"), (axs[1], ekem, "b) Mean EKE (2001–2025)")):
    polar_ax_style(ax, ccrs, cfeature)
    im1 = ax.pcolormesh(LON, LAT, fld, transform=ccrs.PlateCarree(), cmap="YlGnBu", vmin=0, vmax=vmax_mean,
                        shading="auto", rasterized=True)
    cs = add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
    add_ice_mask(ax, ccrs, lon, lat, ice); add_sector_boundaries(ax, ccrs)
    ax.set_title(ttl, fontsize=12, fontweight="bold")
axs[1].clabel(cs, fmt="%.1f", fontsize=6, inline=True)
cb = fig.colorbar(im1, ax=axs[:2], shrink=0.42, pad=0.02); cb.set_label(r"Kinetic Energy (cm$^2$ s$^{-2}$)"); cb.ax.tick_params(labelsize=8)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F3F2EF")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
for ax, v, ttl in ((axs[2], "KE", "c) KE Trend"), (axs[3], "EKE", "d) EKE Trend")):
    polar_ax_style(ax, ccrs, cfeature)
    im3 = ax.pcolormesh(L1, A1, M[f"{v}_slope"].values, transform=ccrs.PlateCarree(), cmap=tc, vmin=-6, vmax=6,
                        shading="nearest", rasterized=True, zorder=0)
    add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
    add_ice_mask(ax, ccrs, lon, lat, ice); add_sector_boundaries(ax, ccrs)
    jj, ii = np.nonzero(M[f"{v}_fdr"].values == 1)
    ax.scatter(M.lon.values[ii], M.lat.values[jj], s=9, c="#111111", lw=0, transform=ccrs.PlateCarree(), zorder=6)
    ax.set_title(ttl, fontsize=12, fontweight="bold")
cb = fig.colorbar(im3, ax=axs[2:], shrink=0.42, pad=0.02, extend="both")
cb.set_label(r"Trend (cm$^2$ s$^{-2}$ yr$^{-1}$)"); cb.ax.tick_params(labelsize=8)

ax5 = fig.add_subplot(gb[0])
x = (S.year + (S.month - 0.5) / 12).values
for v, col in (("KE", STYLE["ke"]), ("EKE", STYLE["eke"])):
    y = S[f"ACC_{v}"].values; b = T.loc[f"ACC_{v}", "slope"]
    ts_panel(ax5, x, y, col, label=v, slope=b, intercept=theil_line(x, deseason(S, f"ACC_{v}"), b))
for i, v in enumerate(("KE", "EKE")):
    a = T.loc[f"ACC_{v}"]
    ax5.text(0.98, 0.13 - 0.07 * i, f"{v}: {a.slope:+.2f} cm$^2$ s$^{{-2}}$ yr$^{{-1}}$, p={a.p_cal:.3f} {stars(a.p_cal)}",
             ha="right", transform=ax5.transAxes, va="top", fontsize=9, color=STYLE["ke" if v == "KE" else "eke"])
ax5.set_xlabel("Year"); ax5.set_ylabel(r"Kinetic Energy (cm$^2$ s$^{-2}$)")
base_timeseries_style(ax5); ax5.legend(loc="upper left", ncol=2, fontsize=9, frameon=False)
title(ax5, "e) ACC KE and EKE Time Series (2001–2025)")
yl = ax5.get_ylim(); ax5.set_ylim(yl[0] - 0.12 * (yl[1] - yl[0]), yl[1] + 0.05 * (yl[1] - yl[0]))

ax7 = fig.add_subplot(gb[1])
rows = []
for reg, lab in (("ACC", "ACC"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian")):
    for v in ("KE", "EKE"):
        rows.append((f"{lab} {v}", {"all-satellite": retest_row(f"allsat|{reg}_{v}|2001-2025"),
                                    "two-satellite": retest_row(f"twosat|{reg}_{v}|2001-2025")}))
forest_panel(ax7, rows, series_colors={"all-satellite": "#1F2933", "two-satellite": "#7A8288"},
             series_markers={"all-satellite": "o", "two-satellite": "s"}, xlim=(-6, 16))
ax7.legend(loc="upper center", bbox_to_anchor=(0.35, -0.16), fontsize=7.5, frameon=False, ncol=3, handletextpad=0.3, columnspacing=1.0)
title(ax7, "f) Regional KE and EKE Trends")
fig.savefig(FIGV3 / "Figure1.png"); fig.savefig(FIGV3 / "Figure1.pdf")
print("saved")
