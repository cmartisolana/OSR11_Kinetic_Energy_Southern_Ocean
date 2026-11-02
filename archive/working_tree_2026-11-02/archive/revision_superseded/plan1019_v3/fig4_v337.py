"""Figure 4 (v3.37): EKE and its atmospheric drivers. Maps of the correlation of EKE with the local wind stress, the SAM and
Nino 3.4 (dk27), with the wind-stress trend and the driver series.
a  Trend of the eastward wind stress, 2001-2025 (% per decade; dk24), mean SAF and PF.
b  Annual series: ACC-band eastward wind stress (anomaly, %), SAM and Nino 3.4 (ONI), 2001-2025.
c-e Correlation of annual EKE (3 deg cells, detrended) with (c) the local eastward wind stress, (d) the SAM, (e) Nino 3.4,
   in the same year (upper) and one year earlier (lower); dots: p < 0.05 (effective sample size); rings: FDR 10 %.
Usage: python fig4_v337.py <revision dir> <drake out dir> <work out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, DO, WO, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
Z24 = np.load(WO / "dk24" / "dk24_patterns.npz", allow_pickle=True); M = np.load(WO / "dk27" / "dk27_maps.npz")
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
ice_ds = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc")
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); an = A.groupby("year").mean(numeric_only=True); yrs = an.index.values + 0.5
CSAM, CONI, CTAU = "#7B3294", "#2C7FB8", "#9A9A9A"
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")


def fronts(ax):
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="#1A1A1A", linewidths=0.45, transform=PC, zorder=6)


fig = plt.figure(figsize=(S.W2, 1.0 * S.W2))
gA = fig.add_gridspec(2, 1, height_ratios=[0.85, 1.25], hspace=0.20)
gT = gA[0].subgridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.22)
# a
ax = fig.add_subplot(gT[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
L2, A2 = np.meshgrid(Z24["lon"], Z24["lat"])
im = ax.pcolormesh(L2, A2, Z24["tau_trend_pct"], transform=PC, cmap=cm_, vmin=-20, vmax=20, shading="nearest", rasterized=True, zorder=0)
fronts(ax); add_ice_mask(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_ds.ice_mask.values > 0.5)
S.plabel(ax, "a", "Eastward wind-stress trend, 2001–2025")
S.hcbar(fig, im, ax, "wind-stress trend (% per decade)", shrink=0.7, extend="both")
# b drivers
ax = fig.add_subplot(gT[1]); S.zero(ax)
tau = 100 * (an.tau_accmask_ACC / an.tau_accmask_ACC.mean() - 1)
ax.bar(yrs, tau, 0.75, color="#D9D9D9", edgecolor="none", zorder=1, label="ACC-band wind stress (left)")
ax.set_ylim(-12, 12); ax.set_ylabel("wind-stress anomaly (%)"); S.ts_style(ax, (2000.6, 2025.6))
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True); ax2.spines["right"].set_linewidth(0.5)
ax2.plot(M["years"] + 0.5, M["sam_annual"], "-o", color=CSAM, lw=0.9, ms=2.0, label="SAM (right)")
ax2.plot(M["years"] + 0.5, M["oni_annual"], "-s", color=CONI, lw=0.9, ms=1.8, label="Niño 3.4, ONI (right)")
ax2.axhline(0, color="#BDBDBD", lw=0.4); ax2.set_ylim(-2.4, 2.4); ax2.set_ylabel("index", fontsize=5.5); ax2.tick_params(labelsize=5)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper left", ncol=3, fontsize=4.8)
r_ = np.corrcoef(an.tau_accmask_ACC - np.polyval(np.polyfit(yrs, an.tau_accmask_ACC, 1), yrs), M["sam_annual"] - np.polyval(np.polyfit(yrs, M["sam_annual"], 1), yrs))[0, 1]
ax.text(0.99, 0.03, f"detrended: r(wind, SAM) = {r_:.2f}", transform=ax.transAxes, ha="right", fontsize=4.8, color="#444444")
S.plabel(ax, "b", "Drivers: ACC-band wind stress, SAM and Niño 3.4, annual means", dy=4)
# c-e maps
gm = gA[1].subgridspec(2, 3, wspace=0.04, hspace=0.10)
L3, A3 = np.meshgrid(M["lon3"], M["lat3"])
COLS = (("wind", "c", "local wind stress"), ("sam", "d", "SAM"), ("oni", "e", "Niño 3.4"))
for j, (nm, lab, ttl) in enumerate(COLS):
    for i, L in enumerate((0, 1)):
        ax = fig.add_subplot(gm[i, j], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature, lons=(i == 0 and j == 0))
        imr = ax.pcolormesh(L3, A3, M[f"r_{nm}_{L}"], transform=PC, cmap=cm_, vmin=-0.7, vmax=0.7, shading="nearest", rasterized=True, zorder=0)
        jj, ii = np.nonzero(M[f"p_{nm}_{L}"] < 0.05); ax.scatter(M["lon3"][ii], M["lat3"][jj], s=0.6, c="k", lw=0, transform=PC, zorder=7)
        jj, ii = np.nonzero(M[f"fdr_{nm}_{L}"]); ax.scatter(M["lon3"][ii], M["lat3"][jj], s=7, facecolor="none", edgecolor="k", lw=0.5, transform=PC, zorder=8)
        fronts(ax)
        if i == 0: S.plabel(ax, lab, f"EKE vs {ttl}", dy=4)
        if j == 0: ax.text(-0.06, 0.5, "same year" if L == 0 else "EKE one year later", transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=5.6)
fig.canvas.draw(); pos = [a.get_position() for a in fig.axes[-6:]]
x0 = min(p.x0 for p in pos); x1 = max(p.x1 for p in pos); y0 = min(p.y0 for p in pos)
cax = fig.add_axes([x0 + 0.15 * (x1 - x0), y0 - 0.03, 0.7 * (x1 - x0), 0.009])
cb = fig.colorbar(imr, cax=cax, orientation="horizontal", extend="both")
cb.set_label("correlation (annual means, detrended; dots: p < 0.05; rings: significant after FDR control)", fontsize=5.4)
cb.ax.tick_params(labelsize=5, width=0.4, length=2); cb.outline.set_linewidth(0.4)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
