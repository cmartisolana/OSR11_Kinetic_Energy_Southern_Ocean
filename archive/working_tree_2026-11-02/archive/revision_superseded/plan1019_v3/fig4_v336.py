"""Figure 4 (v3.36: d wind patterns of SE Pacific wind, SAM, Nino 3.4; e regional EKE vs four drivers at lag 1 yr; dk26)\nFigure 4 (v3.35: c regional lag maps for six regions, d wind coherence, e lag curves with local wind; dk25)\nFigure 4 (v3.34): wind stress and EKE (dk17/s1, dk24). Maps left (as Fig. 3), series and sector comparison right.
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
from matplotlib.patches import Patch
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


D25 = D24.parent / "dk25"; M25 = np.load(D25 / "dk25_maps.npz", allow_pickle=True); LC = pd.read_csv(D25 / "dk25_lagcurves.csv")
WC = np.load(D25 / "dk25_windcoherence.npz")
fig = plt.figure(figsize=(S.W2, 1.02 * S.W2))
gA = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.15], hspace=0.22)
gT = gA[0].subgridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.24); g1 = [gT[0], None]; gC = [gT[1]]
gB = gA[1].subgridspec(1, 2, width_ratios=[1.55, 1.0], wspace=0.22)
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
# a wind trend map
ax = fig.add_subplot(gT[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
L2, A2 = np.meshgrid(Z["lon"], Z["lat"])
im = ax.pcolormesh(L2, A2, Z["tau_trend_pct"], transform=PC, cmap=cm_, vmin=-20, vmax=20, shading="nearest", rasterized=True, zorder=0)
fronts(ax); add_ice_mask(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_ds.ice_mask.values > 0.5)
S.plabel(ax, "a", "Eastward wind-stress trend, 2001–2025")
S.hcbar(fig, im, ax, "wind-stress trend (% per decade)", shrink=0.7, extend="both")
# b wind-stress series
ax = fig.add_subplot(gT[1])
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
# c regional lag maps (Hogg et al. 2022): EKE of six regions vs the wind stress one year earlier, 3 deg cells
REG = {"Kerguelen Plateau": (65, 95), "Southeast Indian Ridge": (130, 150), "Macquarie Ridge": (155, 175),
       "Pacific–Antarctic Ridge": (-145, -100), "Southeast Pacific": (-100, -70), "Drake Passage–Scotia Sea": (-70, -35)}
RC = {"Kerguelen Plateau": "#56B4E9", "Southeast Indian Ridge": "#0072B2", "Macquarie Ridge": "#009E73",
      "Pacific–Antarctic Ridge": "#CC79A7", "Southeast Pacific": "#D55E00", "Drake Passage–Scotia Sea": "#E69F00"}
TTL = {"Kerguelen Plateau": "Kerguelen Plateau", "Southeast Indian Ridge": "SE Indian Ridge", "Macquarie Ridge": "Macquarie Ridge",
       "Pacific–Antarctic Ridge": "Pacific–Antarctic Ridge", "Southeast Pacific": "SE Pacific", "Drake Passage–Scotia Sea": "Drake Passage–Scotia Sea"}
gm = gB[0].subgridspec(2, 3, wspace=0.03, hspace=0.10); L3, A3 = np.meshgrid(M25["lon3"], M25["lat3"]); axs = []
for k, (nm, (a_, b_)) in enumerate(REG.items()):
    ax = fig.add_subplot(gm[k // 3, k % 3], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature, lons=False); axs.append(ax)
    imr = ax.pcolormesh(L3, A3, M25[f"r|{nm}|1"], transform=PC, cmap=cm_, vmin=-0.7, vmax=0.7, shading="nearest", rasterized=True, zorder=0)
    jj, ii = np.nonzero(M25[f"p|{nm}|1"] < 0.05); ax.scatter(M25["lon3"][ii], M25["lat3"][jj], s=0.4, c="k", lw=0, transform=PC, zorder=7)
    fronts(ax); xs = np.linspace(a_, b_, 40)
    ax.plot(xs, np.full_like(xs, -36.3), color=RC[nm], lw=2.6, transform=PC, zorder=8, solid_capstyle="butt")
    ax.set_title(TTL[nm], fontsize=4.9, pad=1.5, color="#222222")
S.plabel(axs[0], "c", "Correlation of regional EKE (coloured arc) with the wind stress one year earlier", dy=12)
fig.canvas.draw(); p0, p2 = axs[3].get_position(), axs[5].get_position()
cax = fig.add_axes([p0.x0 + 0.03, p0.y0 - 0.03, p2.x1 - p0.x0 - 0.06, 0.008])
cb = fig.colorbar(imr, cax=cax, orientation="horizontal", extend="both"); cb.set_label("correlation (annual means, detrended; dots: p < 0.05)", fontsize=5.3)
cb.ax.tick_params(labelsize=5, width=0.4, length=2); cb.outline.set_linewidth(0.4)
# d the wind pattern: one-point map of the SE Pacific wind, and the wind-stress patterns of the SAM and ENSO (dk26)
D26 = D24.parent / "dk26"; M26 = np.load(D26 / "dk26_maps.npz"); RE = pd.read_csv(D26 / "dk26_regional_eke.csv")
IC = pd.read_csv(D26 / "dk26_index_corr.csv")
gd = gB[1].subgridspec(2, 1, height_ratios=[0.62, 1.0], hspace=0.42); gdm = gd[0].subgridspec(1, 3, wspace=0.05)
ann_ = IC[IC.season == "annual"].set_index("pair")
for k, (fld, ttl) in enumerate(((WC["r"], "SE Pacific wind\n(black arc)"), (M26["r_sam"], f"SAM\npattern r = {float(M26['pattern_r_sam']):.2f}"),
                                 (M26["r_oni"], f"Niño 3.4\npattern r = {float(M26['pattern_r_oni']):.2f}".replace("-", "−")))):
    ax = fig.add_subplot(gdm[k], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature, lons=False)
    ax.pcolormesh(L3, A3, fld, transform=PC, cmap=cm_, vmin=-0.7, vmax=0.7, shading="nearest", rasterized=True, zorder=0)
    if k == 0:
        xs = np.linspace(-100, -70, 30); ax.plot(xs, np.full_like(xs, -36.3), color="#111111", lw=2.4, transform=PC, zorder=8, solid_capstyle="butt")
        S.plabel(ax, "d", "Wind stress correlated with:", dy=10)
    ax.set_title(ttl, fontsize=4.4, pad=1.5, linespacing=1.0)
# e regional EKE one year after: local wind, SE Pacific wind, SAM, Nino 3.4
ax = fig.add_subplot(gd[1]); S.zero(ax)
IX = [("local wind", "#111111"), ("SE Pacific wind", "#777777"), ("SAM", "#7B3294"), ("ONI", "#2C7FB8")]
w_ = 0.19
for gi, nm in enumerate(REG):
    for k, (ix, col) in enumerate(IX):
        if ix == "local wind":
            q = LC[(LC.region == nm) & (LC.wind == "local") & (LC.lag == 1)].iloc[0]
        else:
            q = RE[(RE.region == nm) & (RE["index"] == ix) & (RE.lag == 1)].iloc[0]
        x = gi + (k - 1.5) * w_
        ax.bar(x, q.r, w_ * 0.9, color=col if q.p < 0.05 else "white", edgecolor=col, lw=0.7, zorder=3)
ax.axhspan(-0.40, 0.40, color="#F1F1F1", lw=0, zorder=0)
ax.set_xticks(range(len(REG))); ax.set_xticklabels(["KP", "SEIR", "MR", "PAR", "SE Pac.", "DP–SS"], fontsize=5.2); ax.tick_params(axis="x", length=0)
ax.set_ylim(-0.7, 0.7); ax.set_ylabel("correlation")
ax.legend(handles=[Patch(facecolor=c, label=("Niño 3.4" if l == "ONI" else l)) for l, c in IX] + [Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")],
          loc="lower left", ncol=3, fontsize=4.3, handlelength=1.0, columnspacing=0.6)
S.plabel(ax, "e", "Regional EKE one year after the forcing", dy=4)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
