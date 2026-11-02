"""v2.2 Figure 1 (five panels): (a) mean EKE with mean jet-core path and 3000 m isobath; (b, c) KE and EKE trends on 1° cells;
(d) ACC-band KE, EKE and tracked-eddy sum of squared amplitudes (% of mean, 12-month running means); (e) forest of trends
incl. tracked-eddy properties 2001-2023. Based on v21_fig1.py.
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
FIGV3 = OUT / "figures_v23"; FIGV3.mkdir(exist_ok=True)
import contourpy
sys.path.insert(0, str(HERE)); import b1_meander as B
ROOT_OUT = Path("/mnt/user-data/uploads/repository/outputs")
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
ice = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc").ice_mask.values > 0.5
bat = xr.open_dataset("/mnt/user-data/uploads/repository/revision/data/bathymetry_so.nc").elevation.coarsen(latitude=4, longitude=4, boundary="trim").mean()
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
fig = plt.figure(figsize=(11.5, 13.5))
outer = fig.add_gridspec(3, 1, height_ratios=[1.25, 0.75, 1.0], hspace=0.28)
gm = outer[0].subgridspec(1, 3, wspace=0.04)
gb = outer[1].subgridspec(1, 1); gf = outer[2].subgridspec(1, 1)
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

ax5 = fig.add_subplot(gb[0])
x = (S.year + (S.month - 0.5) / 12).values
E = pd.read_csv(OUT / "eddies" / "e2_series.csv")
for y, lab, col, ls in ((S.ACC_KE.values, "KE", STYLE["ke"], "-"), (S.ACC_EKE.values, "EKE", STYLE["eke"], "-"),
                        (E.ACC_amp2_per_day.values, "Σ amplitude² of tracked eddies (2001–2023)", "#2A9D8F", (0, (3, 1.5)))):
    xx = x[:len(y)]; a = 100 * (y / np.nanmean(y) - 1)
    ax5.plot(xx, a, color=col, lw=0.6, alpha=0.25 if lab != "Σ amplitude² of tracked eddies (2001–2023)" else 0)
    ax5.plot(xx, running(a), color=col, lw=2.0, ls=ls, label=lab)
ax5.axhline(0, color="#6b7280", lw=0.8, ls="--")
for i, v in enumerate(("KE", "EKE")):
    a = T.loc[f"ACC_{v}"]
    ax5.text(0.98, 0.13 - 0.07 * i, f"{v}: {a.slope:+.2f} cm$^2$ s$^{{-2}}$ yr$^{{-1}}$, p={a.p_cal:.3f} {stars(a.p_cal)}",
             ha="right", transform=ax5.transAxes, va="top", fontsize=9, color=STYLE["ke" if v == "KE" else "eke"])
ax5.set_xlabel("Year"); ax5.set_ylabel("Anomaly (% of 2001–2025 mean)")
base_timeseries_style(ax5); ax5.legend(loc="upper left", ncol=3, fontsize=8.5, frameon=False); ax5.set_xlim(2000.8, 2025.5)
title(ax5, "d) ACC band: KE, EKE and tracked-eddy energy (12-month running means)")
yl = ax5.get_ylim(); ax5.set_ylim(yl[0] - 0.12 * (yl[1] - yl[0]), yl[1] + 0.25 * (yl[1] - yl[0]))

ax7 = fig.add_subplot(gf[0])
rows = []
for reg, lab in (("ACC", "ACC"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian")):
    for v in ("KE", "EKE"):
        rows.append((f"{lab} {v}", {"gridded KE / EKE, 2001–2025": retest_row(f"allsat|{reg}_{v}|2001-2025")}))
TE = pd.read_csv(OUT / "eddies" / "e2_e3_trends.csv").set_index("series")
def erow(k):
    r = TE.loc[k]; f = 1000 / r["mean"]
    return dict(pct=r.pct_per_decade, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p)
NL = pd.read_csv(OUT / "eddies" / "e5_nonlinearity_trends_lag7.csv").set_index("series")
for lab, k in (("Eddies: Σ amplitude²", "ACC_amp2_per_day"), ("Eddies: mean amplitude", "ACC_amp"), ("Eddies: rotational speed", "ACC_spd"),
               ("Eddies: number present", "ACC_present_per_day"), ("Eddies: births per month", "ACC_births")):
    rows.append((lab, {"tracked eddies 2001–2023": erow(k)}))
rows.append(("Eddies: nonlinearity U/c", {"tracked eddies 2001–2023": dict(pct=NL.loc["ratio_median", "pct_per_decade"], pct_lo=NL.loc["ratio_median", "pct_lo"], pct_hi=NL.loc["ratio_median", "pct_hi"], p=NL.loc["ratio_median", "p"])}))
forest_panel(ax7, rows, series_colors={"gridded KE / EKE, 2001–2025": "#1F2933", "tracked eddies 2001–2023": "#2A9D8F"},
             series_markers={"gridded KE / EKE, 2001–2025": "o", "tracked eddies 2001–2023": "D"}, xlim=(-6, 18), label_size=8.5, dodge=False)
ax7.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=8, frameon=False, ncol=3, handletextpad=0.3, columnspacing=1.0)
title(ax7, "e) Trends 2001–2025 (KE, EKE) and 2001–2023 (tracked eddies), % of the mean per decade")
fig.savefig(FIGV3 / "Figure1.png"); fig.savefig(FIGV3 / "Figure1.pdf")
print("saved")
