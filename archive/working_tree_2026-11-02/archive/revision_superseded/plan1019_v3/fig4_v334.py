"""Figure 4 (v3.34): wind stress and EKE (dk17/s1, dk24). Maps left (as Fig. 3), series and sector comparison right.
a  Trend of the eastward wind stress, 2001-2025 (least squares on annual means, % of the cell mean per decade), mean SAF/PF.
b  Eastward wind stress averaged over the ACC band and its three sectors: monthly ACC band (thin), 12-month running means,
   trends (calibrated Theil-Sen on monthly series, dk17/s1; solid: p < 0.05, dashed otherwise).
c  EKE trend against wind-stress trend in twelve 30 deg sectors of the ACC band (Zhang et al. 2021), coloured by ocean
   sector, 95 % intervals; least-squares line across sectors.
d  Correlation of the annual EKE of four regions (black arcs) with the wind stress one year earlier on 3 deg cells, both
   detrended (Hogg et al. 2022).
Usage: python fig4_v334.py <revision dir> <drake out dir> <dk24 out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, DO, D24, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
Z = np.load(D24 / "dk24_patterns.npz", allow_pickle=True); ST = pd.read_csv(D24 / "dk24_sector_trends.csv")
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
ice_ds = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc")
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); WT = pd.read_csv(DO / "dk17" / "s1_wind_trends.csv").set_index("series")
SECC = {"Atlantic": "#E69F00", "Indian": "#56B4E9", "Pacific": "#CC79A7"}


def sector_of(c):
    l = c % 360
    return "Indian" if 20 <= l < 147 else ("Pacific" if 147 <= l < 290 else "Atlantic")


def fronts(ax):
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="#1A1A1A", linewidths=0.5, transform=PC, zorder=6)


fig = plt.figure(figsize=(S.W2, 0.98 * S.W2))
gR = fig.add_gridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.24)
g1 = gR[0].subgridspec(2, 1, hspace=0.20, height_ratios=[1.0, 1.0])
gC = gR[1].subgridspec(2, 1, height_ratios=[0.85, 1.15], hspace=0.32)
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
# a wind trend map
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
L2, A2 = np.meshgrid(Z["lon"], Z["lat"])
im = ax.pcolormesh(L2, A2, Z["tau_trend_pct"], transform=PC, cmap=cm_, vmin=-20, vmax=20, shading="nearest", rasterized=True, zorder=0)
fronts(ax); add_ice_mask(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_ds.ice_mask.values > 0.5)
S.plabel(ax, "a", "Eastward wind-stress trend, 2001–2025")
S.hcbar(fig, im, ax, "wind-stress trend (% per decade)", shrink=0.7, extend="both")
# d regional lag maps (left, below a)
gd = g1[1].subgridspec(2, 2, wspace=0.04, hspace=0.16)
REG = {"SE Pacific (130-80W)": (-130, -80), "Kerguelen lee (70-100E)": (70, 100), "SE Indian-Macquarie (110-170E)": (110, 170), "Atlantic (50W-20E)": (-50, 20)}
TT = {"SE Pacific (130-80W)": "SE Pacific, 130–80°W", "Kerguelen lee (70-100E)": "Kerguelen lee, 70–100°E",
      "SE Indian-Macquarie (110-170E)": "SE Indian–Macquarie, 110–170°E", "Atlantic (50W-20E)": "Atlantic, 50°W–20°E"}
L4, A4 = np.meshgrid(Z["lon4"], Z["lat4"]); axs = []
for k, (nm, (a_, b_)) in enumerate(REG.items()):
    ax = fig.add_subplot(gd[k // 2, k % 2], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature, lons=False); axs.append(ax)
    imr = ax.pcolormesh(L4, A4, Z[f"reg_{nm}_lag1"], transform=PC, cmap=cm_, vmin=-0.7, vmax=0.7, shading="nearest", rasterized=True, zorder=0)
    fronts(ax); xs = np.linspace(a_, b_, 40)
    ax.plot(xs, np.full_like(xs, -36.5), color="#111111", lw=2.2, transform=PC, zorder=8, solid_capstyle="butt")
    ax.set_title(TT[nm], fontsize=4.8, pad=1.5)
S.plabel(axs[0], "d", "Regional EKE vs wind stress one year earlier", dy=12)
fig.canvas.draw()
p0, p3 = axs[2].get_position(), axs[3].get_position()
cax = fig.add_axes([p0.x0 + 0.02, p0.y0 - 0.035, p3.x1 - p0.x0 - 0.04, 0.008])
cb = fig.colorbar(imr, cax=cax, orientation="horizontal", extend="both"); cb.set_label("correlation (annual means, detrended)", fontsize=5.5)
cb.ax.tick_params(labelsize=5, width=0.4, length=2); cb.outline.set_linewidth(0.4)
# b wind-stress series
ax = fig.add_subplot(gC[0])
x = A.x.values; tau = A.tau_accmask_ACC.values
for sec, col in (("ACC", "#111111"), ("Pacific", SECC["Pacific"]), ("Indian", SECC["Indian"]), ("Atlantic", SECC["Atlantic"])):
    y = A[f"tau_accmask_{sec}"].values; r = WT.loc[f"tau_accmask_{sec}"]
    ax.plot(x, pd.Series(y).rolling(12, center=True).mean(), color=col, lw=1.1 if sec == "ACC" else 0.7, zorder=3 if sec == "ACC" else 2,
            label=f"{'ACC band' if sec == 'ACC' else sec}: {r.pct_dec:+.1f} % per decade ({r.pct_dec_lo:+.1f} to {r.pct_dec_hi:+.1f})")
    xx = np.array([2001, 2026]); ax.plot(xx, y.mean() + r.slope_per_yr * (xx - x.mean()), color=col, lw=0.7,
                                         ls="-" if round(r.p, 3) < 0.05 else (0, (3, 1.6)), zorder=4)
S.ts_style(ax, (2000.8, 2026.2)); ax.set_ylabel("eastward wind stress (N m$^{-2}$)")
ax.legend(loc="upper left", ncol=2, fontsize=4.5, columnspacing=0.8)
ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 0.03)
S.plabel(ax, "b", "Wind stress over the ACC band, 12-month running means and trends", dy=4)
# c sector scatter, coloured by ocean sector, with fit
ax = fig.add_subplot(gC[1]); S.zero(ax); S.zero(ax, "x")
E_ = ST[ST["var"] == "EKE"].set_index("lon"); T_ = ST[ST["var"] == "tau"].set_index("lon")
xx = np.linspace(T_.pct.min() - 1, T_.pct.max() + 1, 10); b1, b0 = np.polyfit(T_.pct, E_.pct, 1)
ax.plot(xx, b0 + b1 * xx, color="#777777", lw=0.7, ls=(0, (3, 1.6)), zorder=1)
OFF = {-15: (-15, 4), 45: (3, 4), 105: (-6, -8), 75: (-15, -2), -135: (-6, 5), 165: (4, -6), -105: (5, 2), -75: (5, -2), 135: (4, 3), -165: (4, -7), 15: (4, -6), -45: (3, 4)}
for c_ in E_.index:
    col = SECC[sector_of(c_)]; e, t = E_.loc[c_], T_.loc[c_]
    ax.plot([t.lo, t.hi], [e.pct, e.pct], color=col, lw=0.5, alpha=0.35, zorder=1); ax.plot([t.pct, t.pct], [e.lo, e.hi], color=col, lw=0.5, alpha=0.35, zorder=1)
    ax.plot(t.pct, e.pct, "o", ms=4.2, color=col, mfc=col if e.p < 0.05 else "white", mew=0.9, zorder=3)
    ax.annotate(f"{abs(c_):.0f}°{'W' if c_ < 0 else 'E'}", (t.pct, e.pct), xytext=OFF.get(int(c_), (4, 3)), textcoords="offset points",
                fontsize=4.6, color="#333333", zorder=4)
r = np.corrcoef(T_.pct, E_.pct)[0, 1]; sp = stats.spearmanr(T_.pct, E_.pct)
ax.text(0.02, 0.97, f"across the 12 sectors: r = {r:.2f}, Spearman ρ = {sp.correlation:.2f} (p = {sp.pvalue:.2f})\nslope {b1:.1f} % EKE per % wind stress",
        transform=ax.transAxes, fontsize=5.0, va="top")
ax.set_xlabel("wind-stress trend (% per decade)", labelpad=1); ax.set_ylabel("EKE trend (% per decade)")
ax.set_xlim(-8, 19); ax.set_ylim(-14, 27)
ax.legend(handles=[Line2D([], [], ls="", marker="o", ms=4.5, color=c, label=s) for s, c in SECC.items()] +
                  [Line2D([], [], ls="", marker="o", ms=4.5, color="#777777", mfc="white", label="open: EKE p ≥ 0.05")],
          loc="lower right", fontsize=4.8, ncol=2)
ax.text(0.98, 0.20, "labels: sector centre longitude (°)", transform=ax.transAxes, ha="right", fontsize=4.4, color="#777777")
S.plabel(ax, "c", "EKE and wind-stress trends in 30° sectors of the ACC band", dy=4)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
