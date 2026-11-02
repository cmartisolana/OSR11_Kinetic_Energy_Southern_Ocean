"""Figure 1 (v2.6): Southern Ocean KE and EKE, 2001-2025, and the tracked-eddy decomposition (style_nc).
a mean EKE with the mean jet-core path and the 3000 m isobath; b, c KE and EKE trends on 1-degree cells (dots: FDR);
d ACC-band KE, EKE and tracked-eddy sum of squared amplitudes, 12-month running means in % of the mean;
e KE and EKE trends by region, 2001-2025; f tracked-eddy trends 2001-2023 with the gridded EKE over the same months.
Usage: OSR11_CACHE=<cache> python fig1.py
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, contourpy
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.colors import LogNorm
from style_nc import *  # noqa

FIGD = OUT / "figures_final"
cache = Path(os.environ.get("OSR11_CACHE", HERE / "cache"))
ROOT_OUT = HERE.parent / "outputs"
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
ice = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc").ice_mask.values > 0.5
bat = xr.open_dataset(HERE / "data" / "bathymetry_so.nc").elevation.coarsen(latitude=4, longitude=4, boundary="trim").mean()
lat = np.load(cache / "lat.npy"); lon = np.load(cache / "lon.npy")
ekem = np.nanmean(np.asarray(np.load(cache / "eke.npy", mmap_mode="r"), np.float32), 0) * 1e4; ekem[ice] = np.nan
M = xr.open_dataset(OUT / "a1_maps_bc.nc")
S = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv")
E = pd.read_csv(OUT / "eddies" / "e2_series.csv")
TE = pd.read_csv(OUT / "eddies" / "e2_e3_trends.csv").set_index("series")
PC = ccrs.PlateCarree()


def jet_path(ax):
    a = adt.values; acc = np.isfinite(a) & (a >= -0.6) & (a <= 0.2) & ~ice
    rel = a - np.nanmean(a[acc]); rel = np.where((a > 0.45) | (adt.latitude.values[:, None] > -40), np.nan, rel)
    for L in contourpy.contour_generator(adt.longitude.values, adt.latitude.values, np.ma.masked_invalid(rel), name="serial").lines(0.18):
        if np.ptp(L[:, 0]) >= 5:
            ax.plot(L[:, 0], L[:, 1], color="k", lw=0.7, transform=PC, zorder=5)


fig = plt.figure(figsize=(W2, 6.3))
g = fig.add_gridspec(3, 1, height_ratios=[1.25, 0.6, 0.8], hspace=0.42)
gm = g[0].subgridspec(1, 3, wspace=0.04)
gb = g[2].subgridspec(1, 2, wspace=0.9, width_ratios=[1, 1])
axm = [fig.add_subplot(gm[0, j], projection=ccrs.SouthPolarStereo()) for j in range(3)]
LON, LAT = np.meshgrid(lon, lat)
ax = axm[0]; polar(ax, ccrs, cfeature)
im = ax.pcolormesh(LON, LAT, ekem, transform=PC, cmap="YlGnBu", norm=LogNorm(20, 2000), shading="auto", rasterized=True)
ax.contour(bat.longitude, bat.latitude, -bat, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
jet_path(ax); add_ice_mask(ax, ccrs, lon, lat, ice); add_sector_boundaries(ax, ccrs)
for lo, s in ((-20, "Atlantic"), (75, "Indian"), (-140, "Pacific")):
    ax.text(lo, -41, s, transform=PC, ha="center", va="center", fontsize=5.3, zorder=8,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.75))
plabel(ax, "a", "Mean EKE")
hcbar(fig, im, ax, "EKE (cm$^2$ s$^{-2}$)")
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
for ax, v, l in ((axm[1], "KE", "b"), (axm[2], "EKE", "c")):
    polar(ax, ccrs, cfeature)
    im2 = ax.pcolormesh(L1, A1, M[f"{v}_slope"].values, transform=PC, cmap=tc, vmin=-6, vmax=6, shading="nearest", rasterized=True, zorder=0)
    add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
    add_ice_mask(ax, ccrs, lon, lat, ice); add_sector_boundaries(ax, ccrs)
    jj, ii = np.nonzero(M[f"{v}_fdr"].values == 1)
    ax.scatter(M.lon.values[ii], M.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
    plabel(ax, l, f"{v} trend")
hcbar(fig, im2, axm[1:], "Trend (cm$^2$ s$^{-2}$ yr$^{-1}$)", shrink=0.45, extend="both")

# d series
ax = fig.add_subplot(g[1])
x = (S.year + (S.month - 0.5) / 12).values
for y, lab, col, ls, raw in ((S.ACC_KE.values, "KE", C["ke"], "-", True), (S.ACC_EKE.values, "EKE", C["eke"], "-", True),
                             (E.ACC_amp2_per_day.values, "Σ squared amplitude of tracked eddies", C["te"], (0, (2.5, 1.2)), False)):
    xx = x[:len(y)]; a = 100 * (y / np.nanmean(y) - 1)
    if raw: ax.plot(xx, a, color=col, lw=0.35, alpha=0.3)
    ax.plot(xx, running(a), color=col, lw=1.0, ls=ls, label=lab)
zero(ax); ts_style(ax); ax.set_ylabel("Anomaly (% of mean)")
ax.legend(loc="upper left", ncol=3)
yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.18 * (yl[1] - yl[0]))
plabel(ax, "d", "ACC band", dx=-24)

# e regional trends
ax = fig.add_subplot(gb[0, 0])
rows = [(lab, {"KE": retest_row(f"allsat|{reg}_KE|2001-2025"), "EKE": retest_row(f"allsat|{reg}_EKE|2001-2025")})
        for reg, lab in (("ACC", "ACC band"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian"))]
forest(ax, rows, {"KE": C["ke"], "EKE": C["eke"]}, xlim=(-2, 18), xlabel="Trend 2001–2025 (% per decade)")
ax.legend(handles=sig_handles({"KE": C["ke"], "EKE": C["eke"]}), loc="upper center", bbox_to_anchor=(0.45, -0.3), ncol=3, fontsize=5.8)
plabel(ax, "e", "KE and EKE by region", dx=-40)


def erow(k):
    r = TE.loc[k]; f = 1000 / r["mean"]; return dict(pct=r.pct_per_decade, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p)


ax = fig.add_subplot(gb[0, 1])
TD, GE = "tracked eddies", "gridded EKE"
rows = [("Σ squared amplitude", {TD: erow("ACC_amp2_per_day")}), ("Mean amplitude", {TD: erow("ACC_amp")}),
        ("Rotational speed", {TD: erow("ACC_spd")}), ("Number present", {TD: erow("ACC_present_per_day")}),
        ("Gridded EKE", {GE: erow("ACC_EKE_allsat")})]
forest(ax, rows, {TD: C["te"], GE: C["eke"]}, {TD: "D", GE: "o"}, xlim=(-2, 18), dodge=False,
       xlabel="Trend 2001–2023 (% per decade)")
ax.axhline(-3.5, color=C["light"], lw=0.5)
ax.legend(handles=sig_handles({TD: C["te"], GE: C["eke"]}, {TD: "D"}), loc="upper center", bbox_to_anchor=(0.4, -0.3), ncol=3, fontsize=5.8)
plabel(ax, "f", "Tracked eddies, ACC band", dx=-58)
save(fig, FIGD, "Figure1")
