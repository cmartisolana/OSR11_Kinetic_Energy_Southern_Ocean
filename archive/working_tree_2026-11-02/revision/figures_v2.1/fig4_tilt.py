"""Figure 4 (v2.1, 28 Oct 2026, third pass): wind forcing and the isopycnal slope (eddy-saturation test; Thompson and Naveira
Garabato, 2014: under eddy saturation the mean isopycnal slope, and thus the baroclinic transport, is insensitive to the wind,
the extra momentum input being taken up by eddies and standing meanders). Maps on top, c and d full width. No EKE (Figure 3).
a  Eastward wind-stress trend 2001-2025 (Product 2, 10-3 N m-2 per decade; f3w_wind.py), dots: FDR 10 %.
b  Trend of the cross-stream slope of sigma0 = 27.3 (Roemmich-Gilson Argo, monthly, 2004-2025, % of the local mean per decade;
   calibrated monthly test, f4m_map.py); domain ACC band + 100 km, ice removed; dots: FDR 10 %.
c  Wind stress, EKE and slope on the common cells, anomalies in % of the 2004-2025 mean (12-month running means)
   with the 2004-2025 trends and their 95 % wedges (f4e_sat.py D1), values labelled at the right.
d  Sensitivity d ln X / d ln tau of KE, EKE and slope from interannual anomalies, wind leading by 1 yr (f4s_sens.py), against the
   reference values of the eddy-saturation studies (1: Meredith et al. 2012, Marshall et al. 2017; 0.7: Stewart et al. 2023; 0 for
   the transport: Boening et al. 2008); e the lagged correlations (f4e_sat.py D2) that fix the 1-yr lag. Trend ratios: SI.
Usage: python fig4_tilt.py <revision dir> <f4 out dir> <figure dir> <f4d out dir> <f4m out dir> <f3w out dir> <f4b out dir> <f4e out dir> <f4s out dir> <f4ae out dir>
       (F4_OFFLINE=1: land from the bathymetry, for machines without the Natural Earth files)
"""
import sys, os, types
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

REV, F4, FIGD, FD, FM, FW, FB, FE, FS, FA = (Path(a) for a in sys.argv[1:11]); FIGD.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV))
import style_nc as S  # noqa: E402
from style_nc import add_ice_mask  # noqa: E402
import cartopy.crs as ccrs, cartopy.feature as cfeature  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
OUT = REV / "outputs"; ROOT_OUT = REV.parent / "outputs"
M = xr.open_dataset(F4 / "f4_tilt_maps.nc"); MM_ = xr.open_dataset(FM / "f4m_map.nc")
MR = pd.read_csv(FM / "f4m_results.csv").set_index("series"); MS = np.load(FM / "f4m_series.npz")
L5m = pd.read_csv(FM / "f4m_lon5_trends.csv").sort_values("lon")
MW = xr.open_dataset(FW / "f3w_wind_map.nc"); WMo = pd.read_csv(FD / "f4d_wind_monthly.csv")
WB = pd.read_csv(FB / "f4b_bins.csv").sort_values("lon"); WR = pd.read_csv(FB / "f4b_results.csv").set_index("series")
RW4 = pd.read_csv(FM / "f4w_results.csv").iloc[0]
RT = pd.read_csv(F4 / "f4_regional_trends.csv").set_index("series")
SER = np.load(F4 / "f4_regional_series.npz"); yrs = SER["years"]
L5 = pd.read_csv(F4 / "f4_lon5_trends.csv"); L5 = L5[L5.surface.astype(str) == "27.3"].sort_values("lon")
DR = pd.read_csv(FD / "f4d_results.csv").set_index("series"); EB = pd.read_csv(FD / "f4d_eke_bins.csv").sort_values("lon")
RTR = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key"); KM = 111.195
DN = np.load(REV / "plan1019" / "drake" / "out" / "dk14" / "dk14_hp_density_circ.npz")
clim = xr.open_dataset(OUT / "clim_allsat.nc")
ice_ds = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc"); ice_g = ice_ds.ice_mask.values > 0.5
lat8, lon8 = ice_ds.latitude.values, ice_ds.longitude.values
bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation
bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
DK = REV / "plan1019" / "drake" / "out" / "dk12"
Z = np.load(DK / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z["lev_saf"]), float(Z["lev_pf"])
TB2 = pd.read_csv(DK / "dk12_band_trends_allsat.csv").sort_values("lon")

FEAT = [(75, -47.5, "Kerguelen\nPlateau"), (122, -46.0, "Southeast\nIndian Ridge"), (166, -49.0, "Campbell\nPlateau"),
        (-140, -60.5, "Pacific–Antarctic\nRidge"), (-66, -60.0, "Drake\nPassage")]
OBS = [(-142, "PAR"), (-65, "Drake"), (-40, "Scotia"), (37, "SWIR"), (78, "Kerguelen"), (130, "SEIR"), (160, "MR")]
LONX = (-180, 180); TICK = np.arange(-180, 181, 60)
CS, CW, CE = "#8E44AD", C["wind"], C["eke"]          # slope (Argo purple), wind grey, EKE green (as Fig. 3)
OFFLINE = os.environ.get("F4_OFFLINE") == "1"
_cf = types.SimpleNamespace(LAND=cfeature.ShapelyFeature([], PC)) if OFFLINE else cfeature


def sigp(p):
    return np.isfinite(p) and round(float(p), 3) < 0.05


def fp(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"


def lonfmt(ticks):
    return [("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in ticks]


def edges(c):
    d = np.diff(c).mean(); return np.r_[c - d / 2, c[-1] + d / 2]


def base_map(ax):
    S.polar(ax, ccrs, _cf)
    if OFFLINE:
        ax.contourf(bc.longitude, bc.latitude, (bc.values >= 0).astype(float), levels=[0.5, 1.5], colors=["#3A3A3A"], transform=PC, zorder=3)
    ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[LEV_PF, LEV_SAF], colors="#1A1A1A", linewidths=0.6, transform=PC, zorder=6)
    add_ice_mask(ax, ccrs, lon8, lat8, ice_g)
    ax.contour(M.lon.values, M.lat.values, M["buffer100"].values.astype(float), levels=[0.5], colors="#555555", linewidths=0.5,
               linestyles=[(0, (2, 1.5))], transform=PC, zorder=5)


def trend_map(ax, field, fdr, vmax, letter, title, cblabel, feat=False):
    base_map(ax)
    tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")
    im = ax.pcolormesh(le, ae, field, transform=PC, cmap=tc, vmin=-vmax, vmax=vmax, shading="flat", rasterized=True, zorder=0)
    jj, ii = np.nonzero(fdr)
    ax.scatter(M.lon.values[ii], M.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
    if feat:
        for lo_, la_, s_ in FEAT:
            ax.text(lo_, la_, s_, transform=PC, fontsize=4.6, ha="center", va="center", color="#222222", zorder=9, linespacing=0.9,
                    bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
    S.plabel(ax, letter, title)
    S.hcbar(fig, im, ax, cblabel, shrink=0.7, extend="both")


def depth_axis(ax):
    axd = ax.twinx(); axd.fill_between(TB2.lon, TB2.depth_mean / 1000, 7, color="#E3E3E3", step="mid", lw=0, zorder=0)
    axd.set_ylim(7, -3.5); axd.set_yticks([0, 1, 2, 3, 4, 5, 6]); axd.set_ylabel("Depth (km)", color="#777777", fontsize=5.5)
    axd.tick_params(axis="y", colors="#777777", labelsize=5); axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5)
    ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)


def profile(ax, x, y, lo, hi, p, col, mk, band=False, YL=(-12, 16)):
    if band: ax.fill_between(x, np.clip(lo, *YL), np.clip(hi, *YL), color=col, alpha=0.18, lw=0, zorder=1)
    ax.plot(x, y, color=col, lw=0.8 if band else 0.7, zorder=3)
    sg = np.round(p, 3) < 0.05
    ax.plot(x[sg], y[sg], mk, color=col, ms=2.0 if mk == "o" else 1.7, zorder=4)
    ax.plot(x[~sg], y[~sg], mk, mfc="white", color=col, ms=2.0 if mk == "o" else 1.7, mew=0.4, zorder=4)


def lon_axis(ax, ylabel, YL=(-12, 16), labels=True, xlabel=True):
    S.zero(ax); ax.set_xlim(*LONX); ax.set_xticks(TICK); ax.set_xticklabels(lonfmt(TICK) if xlabel else [])
    if xlabel: ax.set_xlabel("Longitude (flow from left to right)")
    ax.set_ylabel(ylabel); ax.set_ylim(*YL)
    if labels:
        for x_, s_ in OBS: ax.text(x_, 0.99, s_, ha="center", va="top", fontsize=4.8, color="#222222", transform=ax.get_xaxis_transform())
        ax.annotate("", xy=(175, 0.88), xytext=(155, 0.88), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"))
        ax.text(165, 0.90, "flow", ha="center", va="bottom", fontsize=5, color="#555555", transform=ax.get_xaxis_transform())


le, ae = edges(M.lon.values), edges(M.lat.values)
fig = plt.figure(figsize=(S.W2, 1.04 * S.W2))
gA = fig.add_gridspec(3, 1, height_ratios=[1.5, 0.62, 0.78], hspace=0.24)
g1 = gA[0].subgridspec(1, 2, wspace=0.10)
gC = [gA[1], gA[2]]

# a wind-stress trend map
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); base_map(ax)
tcw = plt.get_cmap("RdBu_r").copy(); tcw.set_bad("#F2F2F2")
imw = ax.pcolormesh(edges(MW.lon.values), edges(MW.lat.values), 1000 * MW["wind_trend_dec"].values, transform=PC, cmap=tcw, vmin=-12, vmax=12,
                    shading="flat", rasterized=True, zorder=0)
jj, ii = np.nonzero(MW["wind_fdr"].values > 0)
ax.scatter(MW.lon.values[ii], MW.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
for lo_, la_, s_ in FEAT:
    ax.text(lo_, la_, s_, transform=PC, fontsize=4.6, ha="center", va="center", color="#222222", zorder=9, linespacing=0.9,
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
S.plabel(ax, "a", "Eastward wind-stress trend, 2001–2025")
S.hcbar(fig, imw, ax, "Wind-stress trend (10$^{-3}$ N m$^{-2}$ per decade)", shrink=0.7, extend="both")

# b slope trend map
ax = fig.add_subplot(g1[1], projection=ccrs.SouthPolarStereo())
trend_map(ax, MM_["slope_pct"].values, MM_["slope_fdr"].values > 0, 8, "b", "Isopycnal-slope trend, σ$_0$ = 27.3, 2004–2025",
          "Slope trend (% of mean per decade)")
ax.legend(handles=[Line2D([], [], color="#1A1A1A", lw=0.6, label="mean SAF and PF"),
                   Line2D([], [], color="#555555", lw=0.5, ls=(0, (2, 1.5)), label="ACC band + 100 km")],
          loc="lower left", bbox_to_anchor=(-0.04, -0.04), fontsize=4.8)

# c wind stress, EKE and slope on the common cells (f4e_series.csv), anomalies in % of the 2004-2025 mean: 
#   12-month running means, 2004-2025; trend lines through the mean with their 95 % wedge (monthly test, f4e D1; solid:
#   p < 0.05); values labelled at the right end of each trend line.
ES = pd.read_csv(FE / "f4e_series.csv"); ER = pd.read_csv(FE / "f4e_results.csv"); D1 = ER[ER.test == "D1"].set_index("series")
ax = fig.add_subplot(gC[0]); xe = ES.x.values; c4 = ES.year.values >= 2004; xc = xe[c4].mean(); yrs4 = np.arange(2004, 2026)
xx = np.linspace(2004, 2026, 50); LABY = {"wind": 2.8, "EKE": 12.0, "slope": -6.4}
for k, col, nm, z in (("wind", CW, "wind stress", 2), ("EKE", CE, "EKE", 3), ("slope", CS, "isopycnal slope", 4)):
    v = ES[k].values; mu = np.nanmean(v[c4]); pa = 100 * (v / mu - 1); r = D1.loc[k]
    ax.fill_between(xx, r.pct_lo / 10 * (xx - xc), r.pct_hi / 10 * (xx - xc), color=col, alpha=0.10, lw=0, zorder=1)
    ax.plot(xe[c4], pd.Series(pa).rolling(12, center=True).mean().values[c4], color=col, lw=1.0 if k != "slope" else 1.3, alpha=0.9, zorder=z)
    ax.plot(xx, r.pct / 10 * (xx - xc), color=col, lw=0.9, ls="-" if sigp(r.p) else (0, (3, 1.6)), zorder=5)
    ax.text(2026.35, LABY[k], f"{nm}\n{r.pct:+.1f} % per decade\n({r.pct_lo:+.1f} to {r.pct_hi:+.1f}), {fp(r.p)}", color=col, fontsize=4.9,
            ha="left", va="center", linespacing=1.0)
ax.plot([2003.8, 2026.2], [0, 0], color="#9A9A9A", lw=0.5, zorder=0); ax.set_xlim(2003.8, 2030.6); ax.set_xticks([2005, 2010, 2015, 2020, 2025]); ax.spines["bottom"].set_bounds(2003.8, 2026.2)
ax.set_ylabel("Anomaly\n(% of the 2004–2025 mean)"); ax.set_ylim(-15, 19); ax.set_yticks([-10, -5, 0, 5, 10, 15])
for yg in (-10, -5, 5, 10, 15): ax.plot([2003.8, 2026.2], [yg, yg], color="#EDEDED", lw=0.4, zorder=0)
ax.text(0.01, 0.97, "12-month running means; straight lines and shading: trend 2004–2025 with its 95 % interval",
        transform=ax.transAxes, ha="left", va="top", fontsize=4.6, color="#555555")
S.plabel(ax, "c", "Wind stress, EKE and isopycnal slope over the same cells (ACC band + 100 km)", dy=4)

# d sensitivity to the wind stress, d ln X / d ln tau, the quantity of the eddy-saturation model studies (Meredith et al. 2012;
#   Munday et al. 2013; Marshall et al. 2017: EKE ~ tau, i.e. about 1; Stewart et al. 2023: about 0.7; transport/slope about 0,
#   Boening et al. 2008), from interannual anomalies with the wind leading by 1 yr (2005-2025; f4s_sens.py); e: the lagged
#   correlations that fix the lag (f4e D2, detrended annual, 2004-2025; shading: pointwise 95 % surrogate envelope).
#   The ratios of the 2004-2025 trends (f4e D1 ratio) are in the SI.
SR = pd.read_csv(FS / "f4s_results.csv").set_index("series"); D2 = ER[ER.test == "D2"]
gd = gC[1].subgridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.30); COLS = {"KE": C["ke"], "EKE": CE, "slope": CS}; NM = ["KE", "EKE", "isopycnal\nslope"]
REFH = [Line2D([], [], color="#777777", lw=0.6, ls=(0, (3, 1.6)), label="eddy energy linear in τ (Meredith et al., 2012; Marshall et al., 2017)"),
        Line2D([], [], color="#AAAAAA", lw=0.6, ls=(0, (1, 1.2)), label="eddy energy, 0.7 (Stewart et al., 2023)"),
        Line2D([], [], color="#9A9A9A", lw=0.5, label="transport and isopycnal slope, ≈ 0 (Böning et al., 2008)")]
ax = fig.add_subplot(gd[0])
for q, k in enumerate(("KE", "EKE", "slope")):
    r = SR.loc[k]; col = COLS[k]
    ax.bar(q, r.elasticity, 0.6, color=col if sigp(r.p) else "white", edgecolor=col, lw=0.9, zorder=3); ax.plot([q, q], [r.el_lo, r.el_hi], color="#333333", lw=0.7, zorder=4)
    ax.text(q, max(r.el_hi, r.elasticity) + 0.05, f"{r.elasticity:+.2f}", ha="center", va="bottom", fontsize=4.8)
ax.axhline(1.0, color="#777777", lw=0.6, ls=(0, (3, 1.6)), zorder=0); ax.axhline(0.7, color="#AAAAAA", lw=0.6, ls=(0, (1, 1.2)), zorder=0); S.zero(ax)
ax.set_xticks(range(3)); ax.set_xticklabels(NM, fontsize=5.5); ax.tick_params(axis="x", length=0); ax.set_xlim(-0.6, 2.6)
ax.set_ylabel("d ln X / d ln τ"); ax.set_ylim(-0.25, 2.7)
ax.legend(handles=REFH, loc="upper left", fontsize=4.4, title="reference values", title_fontsize=4.5, handlelength=1.6)
S.plabel(ax, "d", "Sensitivity to the wind, wind leading by 1 yr", dy=4)
ax = fig.add_subplot(gd[1])
for k, mk, nm in (("KE", "s", "KE"), ("EKE", "o", "EKE"), ("slope", "D", "isopycnal slope")):
    q = D2[D2.series == k].sort_values("lag"); col = COLS[k]
    if k == "EKE": ax.fill_between(q.lag, -q.env95, q.env95, color="#E6E6E6", lw=0, zorder=0, label="95 % surrogate envelope")
    sg = q.p.values < 0.05
    ax.plot(q.lag, q.r, color=col, lw=0.9, zorder=3, label=f"{nm} (family-wise {fp(q.p_fw.iloc[0])})")
    ax.plot(q.lag[sg], q.r[sg], mk, color=col, ms=3, zorder=4); ax.plot(q.lag[~sg], q.r[~sg], mk, mfc="white", color=col, ms=3, mew=0.6, zorder=4)
S.zero(ax); ax.set_xticks([0, 1, 2, 3]); ax.set_xlabel("Lag (yr, wind leading)"); ax.set_ylabel("r with wind stress\n(detrended annual, 2004–2025)")
ax.set_ylim(-0.6, 1.0); ax.legend(loc="upper right", fontsize=4.8, ncol=1)
S.plabel(ax, "e", "Lagged correlation with the wind", dy=4)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
