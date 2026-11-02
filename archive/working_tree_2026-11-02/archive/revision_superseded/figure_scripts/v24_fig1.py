"""v2.4 Figure 1 (six panels). Changes from v23_fig1.py (author comments, 2026-10-15):
(d) y-axis range fitted to the 12-month running means (was padded to ~2x the data range);
(e) split: (e) gridded KE and EKE trends by region, KE and EKE side by side in each region, coloured as in (d);
(f) tracked-eddy trends 2001-2023 with the gridded EKE and KE over the same months as reference.
Non-significant eddy parameters removed except the number of eddies present (null reference for
"more intense, not more numerous"); births and U/c remain in the text / SI S4. Set KEEP_NUMBER = False to drop it.
Maps (a-c) unchanged.
Usage: python v24_fig1.py <cache dir with lat.npy, lon.npy, ke.npy, eke.npy>
"""
import sys
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

KEEP_NUMBER = True
cache = sys.argv[1]
FIGV = OUT / "figures_v24"; FIGV.mkdir(exist_ok=True)
import contourpy
ROOT_OUT = HERE.parent / "outputs"
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
ice = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc").ice_mask.values > 0.5
bat = xr.open_dataset(HERE / "data" / "bathymetry_so.nc").elevation.coarsen(latitude=4, longitude=4, boundary="trim").mean()


def add_jet_path(ax):
    mean_adt = adt.values; ACCM = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice
    rel = mean_adt - np.nanmean(mean_adt[ACCM])
    rel = np.where((mean_adt > 0.45) | (adt.latitude.values[:, None] > -40), np.nan, rel)   # keep the ACC branch only
    gen = contourpy.contour_generator(adt.longitude.values, adt.latitude.values, np.ma.masked_invalid(rel), name="serial")
    for Lp in gen.lines(0.18):
        if np.ptp(Lp[:, 0]) < 5: continue
        ax.plot(Lp[:, 0], Lp[:, 1], color="#111111", lw=1.1, transform=ccrs.PlateCarree(), zorder=4)


def add_isobath(ax):
    ax.contour(bat.longitude, bat.latitude, -bat, levels=[3000], colors="#6b7280", linewidths=0.35, transform=ccrs.PlateCarree(), zorder=2)


lat = np.load(f"{cache}/lat.npy"); lon = np.load(f"{cache}/lon.npy")
eke = np.load(f"{cache}/eke.npy", mmap_mode="r")
ekem = np.nanmean(np.asarray(eke, np.float32), 0) * 1e4
ekem[ice] = np.nan
M = xr.open_dataset(OUT / "a1_maps_bc.nc")
S = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv")
T = pd.read_csv(OUT / "a5" / "a5_trends_allsat.csv").query("period == '2001-2025'").set_index("series")
T = retest_override(T, lambda s, r: f"allsat|{s}|2001-2025")


def title(ax, t):
    ax.set_title(t, loc="left", fontsize=12, fontweight="bold")


proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(11.5, 13.0))
outer = fig.add_gridspec(3, 1, height_ratios=[1.25, 0.62, 0.95], hspace=0.30)
gm = outer[0].subgridspec(1, 3, wspace=0.04)
gb = outer[1].subgridspec(1, 1)
gf = outer[2].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.62)
axs = [fig.add_subplot(gm[0, j], projection=proj) for j in range(3)]
LON, LAT = np.meshgrid(lon, lat)
ax = axs[0]; polar_ax_style(ax, ccrs, cfeature)
from matplotlib.colors import LogNorm
im1 = ax.pcolormesh(LON, LAT, ekem, transform=ccrs.PlateCarree(), cmap="YlGnBu", norm=LogNorm(vmin=20, vmax=2000), shading="auto", rasterized=True)
add_isobath(ax); add_jet_path(ax); add_ice_mask(ax, ccrs, lon, lat, ice); add_sector_boundaries(ax, ccrs)
ax.set_title("a) Mean EKE and jet-core path", fontsize=11, fontweight="bold")
cb = fig.colorbar(im1, ax=ax, orientation="horizontal", shrink=0.8, pad=0.04); cb.set_label(r"EKE (cm$^2$ s$^{-2}$)", fontsize=9); cb.ax.tick_params(labelsize=8)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F3F2EF")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
for ax, v, ttl in ((axs[1], "KE", "b) KE trend"), (axs[2], "EKE", "c) EKE trend")):
    polar_ax_style(ax, ccrs, cfeature)
    im3 = ax.pcolormesh(L1, A1, M[f"{v}_slope"].values, transform=ccrs.PlateCarree(), cmap=tc, vmin=-6, vmax=6,
                        shading="nearest", rasterized=True, zorder=0)
    add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values); add_isobath(ax)
    add_ice_mask(ax, ccrs, lon, lat, ice); add_sector_boundaries(ax, ccrs)
    jj, ii = np.nonzero(M[f"{v}_fdr"].values == 1)
    ax.scatter(M.lon.values[ii], M.lat.values[jj], s=9, c="#111111", lw=0, transform=ccrs.PlateCarree(), zorder=6)
    ax.set_title(ttl, fontsize=11, fontweight="bold")
cb = fig.colorbar(im3, ax=axs[1:], orientation="horizontal", shrink=0.4, pad=0.04, extend="both")
cb.set_label(r"Trend (cm$^2$ s$^{-2}$ yr$^{-1}$)", fontsize=9); cb.ax.tick_params(labelsize=8)

# (d) time series, y range fitted to the running means
ax5 = fig.add_subplot(gb[0])
x = (S.year + (S.month - 0.5) / 12).values
E = pd.read_csv(OUT / "eddies" / "e2_series.csv")
C_TE = "#2A9D8F"
lo_hi = []
for y, lab, col, ls, raw in ((S.ACC_KE.values, "KE", STYLE["ke"], "-", True), (S.ACC_EKE.values, "EKE", STYLE["eke"], "-", True),
                             (E.ACC_amp2_per_day.values, "Σ amplitude² of tracked eddies (2001–2023)", C_TE, (0, (3, 1.5)), False)):
    xx = x[:len(y)]; a = 100 * (y / np.nanmean(y) - 1); rm = running(a)
    if raw:
        ax5.plot(xx, a, color=col, lw=0.6, alpha=0.25, clip_on=True)
    ax5.plot(xx, rm, color=col, lw=2.0, ls=ls, label=lab)
    lo_hi += [np.nanmin(a if raw else rm), np.nanmax(a if raw else rm)]
ax5.axhline(0, color="#6b7280", lw=0.8, ls="--")
lo, hi = min(lo_hi), max(lo_hi); span = hi - lo
ax5.set_ylim(lo - 0.04 * span, hi + 0.22 * span)
for i, v in enumerate(("KE", "EKE")):
    a = T.loc[f"ACC_{v}"]
    ax5.text(0.99, 0.04 + 0.11 * (1 - i), f"{v}: {a.slope:+.2f} cm$^2$ s$^{{-2}}$ yr$^{{-1}}$, p={a.p_cal:.3f} {stars(a.p_cal)}",
             ha="right", transform=ax5.transAxes, va="bottom", fontsize=8.5, color=STYLE["ke" if v == "KE" else "eke"])
ax5.set_xlabel("Year", fontsize=10); ax5.set_ylabel("% of 2001–2025 mean", fontsize=10)
base_timeseries_style(ax5); ax5.legend(loc="upper left", ncol=3, fontsize=8.5, frameon=False); ax5.set_xlim(2000.8, 2025.5)
title(ax5, "d) ACC band: KE, EKE and tracked-eddy energy (12-month running means)")

# (e) gridded KE and EKE by region
ax7 = fig.add_subplot(gf[0, 0])
rows = []
for reg, lab in (("ACC", "ACC band"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian")):
    rows.append((lab, {"KE": retest_row(f"allsat|{reg}_KE|2001-2025"), "EKE": retest_row(f"allsat|{reg}_EKE|2001-2025")}))
forest_panel(ax7, rows, series_colors={"KE": STYLE["ke"], "EKE": STYLE["eke"]}, series_markers={"KE": "o", "EKE": "o"},
             xlim=(-2, 18), label_size=9, value_size=7.5, dodge=True)
ax7.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), fontsize=8, frameon=False, ncol=3, handletextpad=0.3, columnspacing=0.9)
title(ax7, "e) KE and EKE trends, 2001–2025")

# (f) tracked eddies 2001-2023 with gridded EKE / KE over the same months
ax8 = fig.add_subplot(gf[0, 1])
TE = pd.read_csv(OUT / "eddies" / "e2_e3_trends.csv").set_index("series")


def erow(k):
    r = TE.loc[k]; f = 1000 / r["mean"]
    return dict(pct=r.pct_per_decade, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p)


TD, GE, GK = "tracked eddies", "gridded EKE, same months", "gridded KE, same months"
erows = [("Σ amplitude²", {TD: erow("ACC_amp2_per_day")}), ("Mean amplitude", {TD: erow("ACC_amp")}),
         ("Rotational speed", {TD: erow("ACC_spd")})]
if KEEP_NUMBER:
    erows.append(("Number present", {TD: erow("ACC_present_per_day")}))
erows += [("Gridded EKE", {GE: erow("ACC_EKE_allsat")}), ("Gridded KE", {GK: erow("ACC_KE_allsat")})]
forest_panel(ax8, erows, series_colors={TD: C_TE, GE: STYLE["eke"], GK: STYLE["ke"]}, series_markers={TD: "D", GE: "o", GK: "o"},
             xlim=(-2, 18), label_size=9, value_size=7.5, dodge=False)
ax8.axhline(-(len(erows) - 2) + 0.5, color="#CED4DA", lw=0.8)
ax8.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), fontsize=8, frameon=False, ncol=2, handletextpad=0.3, columnspacing=0.9)
title(ax8, "f) Tracked eddies, 2001–2023")
fig.savefig(FIGV / "Figure1.png"); fig.savefig(FIGV / "Figure1.pdf")
print("saved")
for gl, d in rows + erows:
    for s, r in d.items(): print(f"{gl:20s} {s:22s} {r['pct']:+6.1f} [{r['pct_lo']:+6.1f}, {r['pct_hi']:+6.1f}] p={r['p']:.3f}")
