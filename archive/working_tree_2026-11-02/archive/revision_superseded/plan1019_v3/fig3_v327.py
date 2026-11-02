"""Main-text figures for manuscript v3.3, Product 1 (all-satellite) unless stated; style_nc house style.
Redesigned toward the physical displays of the reference studies (EKE-trend maps with fronts, annual anomaly
time series of forcing and response as in Meredith and Hogg 2006 / Hogg et al. 2015, trend profiles across the
passage and across streamlines as in Gutierrez-Villanueva et al. 2023, front-position maps and longitude profiles
with topography as in Chambers 2018 / Sokolov and Rintoul 2009). No regions or values of earlier studies are drawn.
Figure 1  Drake Passage: where the flow changed (map, forcing/response series, transect, streamlines).
Figure 2  Circumpolar: EKE trend map, forcing and response series, sector partition, tracked eddies.
Figure 3  Pathway and standing meanders: half-power positions (map, circumpolar series, longitude profile)
          and along-stream EKE, jet speed, sinuosity, with the depth under the path.
Trend lines in the time-series panels are the reported trends (Theil-Sen on deseasonalised monthly series),
drawn through the anomaly mean: solid if p < 0.05 (p rounded to three decimals), dashed otherwise.
Usage: python fig_v33.py <revision dir> <drake out dir> <plan1019 out dir> <figure dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, DO, PO, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import retest_override, add_ice_mask, add_acc_contours, add_sector_boundaries
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import LogNorm, TwoSlopeNorm
plt = S.plt; C = S.C
OUT = REV / "outputs"; ROOT_OUT = REV.parent / "outputs"
PC = ccrs.PlateCarree(); KM = 111.195
FEAT = [(75, -47.5, "Kerguelen\nPlateau"), (122, -46.0, "Southeast\nIndian Ridge"), (166, -49.0, "Campbell\nPlateau"),
        (-140, -60.5, "Pacific–Antarctic\nRidge"), (-66, -60.0, "Drake\nPassage")]


def sigp(p):
    return np.isfinite(p) and round(float(p), 3) < 0.05


def ann(df, col, mean=None, y0=2001, y1=2025):
    s = df[(df.year >= y0) & (df.year <= y1)].groupby("year")[col].mean()
    mu = df[col].mean() if mean is None else mean
    return 100 * (s / mu - 1)


def trend_line(ax, pct, p, col, x0=2001, x1=2025.99, xc=None, lw=0.9):
    xc = (x0 + x1) / 2 if xc is None else xc
    xx = np.array([x0, x1]); yy = pct / 10 * (xx - xc)
    ax.plot(xx, yy, color=col, lw=lw, ls="-" if sigp(p) else (0, (3, 1.6)), zorder=3)


def series(ax, s, col, lab, ms=1.8, lw=0.6, mk="o", ls="-"):
    ax.plot(s.index + 0.5, s.values, ls, marker=mk, color=col, lw=lw, ms=ms, label=lab, alpha=0.9, zorder=4)


def bar_trends(ax, groups, series_def, width=0.17, ylim=None, ylabel=None):
    """groups: [(label, {key: dict(pct, lo, hi, p)})]; series_def: [(key, label, colour)]. Filled bar if p < 0.05,
    white with coloured edge otherwise; whiskers are 95 % intervals."""
    for gi, (gl, d) in enumerate(groups):
        keys = [k for k, _, _ in series_def if k in d]
        for k_i, k in enumerate(keys):
            col = dict((a, c) for a, _, c in series_def)[k]; r = d[k]
            x = gi + (k_i - (len(keys) - 1) / 2) * width
            ax.bar(x, r["pct"], width * 0.9, color=col if sigp(r["p"]) else "white", edgecolor=col, lw=0.8, zorder=3)
            ax.plot([x, x], [r["lo"], r["hi"]], color="#333333", lw=0.6, zorder=4)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g for g, _ in groups]); ax.tick_params(axis="x", length=0)
    S.zero(ax)
    if ylim: ax.set_ylim(*ylim)
    if ylabel: ax.set_ylabel(ylabel)
    h = [Patch(facecolor=c, edgecolor=c, label=l) for _, l, c in series_def]
    h.append(Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05"))
    return h


def lonfmt(ticks):
    return [("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in ticks]


bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation
clim = xr.open_dataset(OUT / "clim_allsat.nc")
ice_ds = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc"); ice_g = ice_ds.ice_mask.values > 0.5
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
RT = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")


def rt(key):
    r = RT.loc[key]; f = 1000 / (r.slope * 1000 / r.pct_per_decade)
    return dict(pct=r.pct_per_decade, lo=r.ci_lo * f, hi=r.ci_hi * f, p=r.p_new)


trk = pd.read_csv(PO / "d1" / "d1b_track.csv")
# ================================================================== Figure 3 (v3.27): circumpolar energy and pathway
# a EKE trend map with the mean SAF and PF; b occurrence of half-power points (two periods) with the jet-core path;
# c moving-period trends of ACC-band EKE and KE (start years 2001-2016, end 2025; after Barcelo-Llull et al. 2025);
# d EKE trend between the SAF and PF by 2 deg longitude, with depth; e meridional displacement by 5 deg longitude
# (half-power latitude, jet-core path latitude, LM Gould half-power) with detection coverage. Hovmoller -> SI.
M = xr.open_dataset(OUT / "a1_maps_bc.nc")
DK = DO / "dk12"; TB2 = pd.read_csv(DK / "dk12_band_trends_allsat.csv")
Z = np.load(DK / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z["lev_saf"]), float(Z["lev_pf"])
MP = pd.read_csv(DO / "dk15" / "dk15_moving_periods_allsat.csv")
DN = np.load(DO / "dk14" / "dk14_hp_density_circ.npz")
LB = RT[RT.index.str.startswith("lonbin|")].copy(); LB["lon0"] = LB.index.str.split("|").str[1].astype(float); LB["c"] = LB.lon0 + 2.5
LB["tr"] = LB.slope * KM * 10; LB["lo"] = LB.ci_lo * KM * 10; LB["hi"] = LB.ci_hi * KM * 10
CL = pd.read_csv(OUT / "b3" / "b3_climatology.csv"); TLAT = pd.read_csv(OUT / "b3" / "b3_trends_lat.csv")
DF = pd.read_csv(OUT / "a2" / "detection_frequency_lon.csv"); DF["b"] = np.floor(DF.lon / 5) * 5 + 2.5; DFb = DF.groupby("b").freq.mean()
_dk2 = pd.read_csv(DO / "dk2" / "dk2_trends.csv"); _r = _dk2[(_dk2["product"] == "allsat") & (_dk2.period == "2001-2025")].iloc[0]
lat, lon = ice_ds.latitude.values, ice_ds.longitude.values
LONX = (-180, 180); TICK = np.arange(-180, 181, 30)
OBS = [(-142, "PAR"), (-65, "Drake"), (-40, "Scotia"), (37, "SWIR"), (78, "Kerguelen"), (138, "SEIR"), (160, "Macquarie")]
bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
fig = plt.figure(figsize=(S.W2, 1.38 * S.W2))
gR = fig.add_gridspec(4, 1, height_ratios=[1.3, 0.62, 0.95, 0.72], hspace=0.36)
g1 = gR[0].subgridspec(1, 2, wspace=0.12)
# a
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
im = ax.pcolormesh(L1, A1, 10 * M["EKE_slope"].values, transform=PC, cmap=tc, vmin=-60, vmax=60, shading="nearest", rasterized=True, zorder=0)
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[LEV_PF, LEV_SAF], colors="#1A1A1A", linewidths=0.6, linestyles="solid", transform=PC, zorder=6)
add_ice_mask(ax, ccrs, lon, lat, ice_g)
jj, ii = np.nonzero(M["EKE_fdr"].values == 1)
ax.scatter(M.lon.values[ii], M.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
for lo_, la_, s_ in FEAT:
    ax.text(lo_, la_, s_, transform=PC, fontsize=4.6, ha="center", va="center", color="#222222", zorder=9, linespacing=0.9,
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
S.plabel(ax, "a", "EKE trend, 2001–2025")
S.hcbar(fig, im, ax, "EKE trend (cm$^2$ s$^{-2}$ per decade)", shrink=0.7, extend="both")
# b
ax = fig.add_subplot(g1[1], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
dd = np.where(DN["all"] > 0, DN["all"], np.nan)
lc = 0.5 * (DN["lon_edges"][1:] + DN["lon_edges"][:-1]); ac = 0.5 * (DN["lat_edges"][1:] + DN["lat_edges"][:-1])
imd = ax.pcolormesh(DN["lon_edges"], DN["lat_edges"], dd, cmap="magma_r", vmin=0, vmax=15, shading="flat", transform=PC, rasterized=True, zorder=0)
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
from scipy.ndimage import gaussian_filter
for key, col, ls in (("early", "#56B4E9", (0, (3, 1.5))), ("late", "#0072B2", "-")):
    ax.contour(lc, ac, gaussian_filter(DN[key], 1.0), levels=[8], colors=col, linewidths=0.45, linestyles=[ls], transform=PC, zorder=4)
ax.plot(np.r_[CL.lon, CL.lon.iloc[0] + 360], np.r_[CL.lat_mean, CL.lat_mean.iloc[0]], color="#D7191C", lw=0.6, transform=PC, zorder=6)
ax.legend(handles=[Line2D([], [], color="#56B4E9", ls=(0, (3, 1.5)), lw=0.8, label="8 %, 2001–2015"), Line2D([], [], color="#0072B2", lw=0.8, label="8 %, 2016–2025"),
                   Line2D([], [], color="#D7191C", lw=0.8, label="jet-core path")], loc="upper left", bbox_to_anchor=(-0.1, 1.0), fontsize=4.6)
S.plabel(ax, "b", "Occurrence of half-power points, 2001–2025")
S.hcbar(fig, imd, ax, "Half-power points (% of months per cell)", shrink=0.7, extend="max")
# c moving-period trends
ax = fig.add_subplot(gR[1]); S.zero(ax); st = np.arange(2001, 2017); w_ = 0.36
for k, (v, lab, col) in enumerate((("ACC_EKE", "EKE", C["eke"]), ("ACC_KE", "KE", C["ke"]))):
    d = MP[MP["var"] == v].set_index("start").loc[st]; x = st + (k - 0.5) * w_
    sg = (d.p.round(3) < 0.05).values
    ax.bar(x[sg], d.pct.values[sg], w_ * 0.9, color=col, zorder=3, label=lab); ax.bar(x[~sg], d.pct.values[~sg], w_ * 0.9, color="white", edgecolor=col, lw=0.8, zorder=3)
    ax.vlines(x, d.pct_lo, d.pct_hi, color="#333333", lw=0.5, zorder=4)
ax.set_xticks(st); ax.set_xticklabels([f"{y}–\n2025" for y in st], fontsize=5.0); ax.set_xlim(2000.4, 2016.6)
ax.set_ylim(-25, 45); ax.set_ylabel("% per decade")
ax.legend(handles=[Patch(facecolor=C["eke"], label="EKE"), Patch(facecolor=C["ke"], label="KE"), Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")],
          loc="upper left", ncol=3, fontsize=5.0)
ax.text(0.995, 0.97, "intervals clipped", transform=ax.transAxes, ha="right", va="top", fontsize=4.6, color="#777777")
S.plabel(ax, "c", "ACC-band trends for periods ending in 2025", dy=4)
# d along-ACC EKE trend between the fronts + depth
g2 = gR[2].subgridspec(2, 1, height_ratios=[1, 0.38], hspace=0.08)
ax = fig.add_subplot(g2[0]); T2 = TB2.dropna(subset=["slope_dec"])
ax.fill_between(T2.lon, T2.lo, T2.hi, color=C["eke"], alpha=0.18, lw=0, zorder=1)
ax.plot(T2.lon, T2.slope_dec, color=C["eke"], lw=0.8, zorder=3)
sg = T2.p.round(3) < 0.05
ax.plot(T2.lon[sg], T2.slope_dec[sg], "o", color=C["eke"], ms=1.8, zorder=4); ax.plot(T2.lon[~sg], T2.slope_dec[~sg], "o", mfc="white", color=C["eke"], ms=1.8, mew=0.4, zorder=4)
S.zero(ax); ax.set_xlim(*LONX); ax.set_xticks(TICK); ax.set_xticklabels([]); ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ dec$^{-1}$)")
ax.annotate("", xy=(175, 0.93), xytext=(155, 0.93), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"))
ax.text(165, 0.97, "flow", ha="center", va="bottom", fontsize=5, color="#555555", transform=ax.get_xaxis_transform())
S.plabel(ax, "d", "EKE trend between the SAF and PF, and depth", dy=4)
ax = fig.add_subplot(g2[1]); TT = TB2.sort_values("lon")
ax.fill_between(TT.lon, -TT.depth_mean / 1000, -6, color="#8A8A8A", step="mid", lw=0)
ax.step(TT.lon, -TT.depth_min / 1000, where="mid", color="#222222", lw=0.6)
ax.set_ylim(-5.5, 1.6); ax.set_yticks([-5, -3, -1]); ax.set_yticklabels(["5", "3", "1"]); ax.set_ylabel("Depth\n(km)")
ax.set_xlim(*LONX); ax.set_xticks(TICK); ax.set_xticklabels([])
for x_, s_ in OBS: ax.text(x_, 1.55, s_, ha="center", va="top", fontsize=4.8, color="#222222")
# e meridional displacement by longitude
ax = fig.add_subplot(gR[3]); S.zero(ax)
axc = ax.twinx(); axc.bar(DFb.index, 100 * DFb.values, width=5, color=C["pos"], alpha=0.10, lw=0, zorder=0)
axc.set_ylim(0, 100); axc.set_ylabel("Detection (% of months)", color="#9C6A92", fontsize=5.5); axc.spines["right"].set_visible(True); axc.tick_params(colors="#9C6A92", labelsize=5)
ax.set_zorder(axc.get_zorder() + 1); ax.patch.set_visible(False)
def pts(ax, x, y, lo, hi, p, col, ms=2.4, mk="o", fdr=None, dx=0.0):
    sg = np.round(p, 3) < 0.05
    ax.vlines(x + dx, lo, hi, color=col, lw=0.6)
    ax.plot((x + dx)[sg], y[sg], mk, color=col, ms=ms, zorder=4); ax.plot((x + dx)[~sg], y[~sg], mk, mfc="white", color=col, ms=ms, mew=0.6, zorder=4)
    if fdr is not None: ax.plot((x + dx)[fdr], y[fdr], "o", mfc="none", mec="#111111", ms=5.5, mew=0.7, zorder=5)
pts(ax, LB.c.values, LB.tr.values, LB.lo.values, LB.hi.values, LB.p_new.values, C["pos"], dx=-0.8)
d = TLAT.dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]
pts(ax, d.lon.values, d.tr.values, d.lo.values, d.hi.values, d.p.values, "#555555", mk="D", ms=2.0, dx=0.8, fdr=d.fdr10.astype(bool).values if "fdr10" in d else None)
ax.errorbar([-64], [_r.trend_km_dec], yerr=[[_r.trend_km_dec - _r.ci_lo], [_r.ci_hi - _r.trend_km_dec]], fmt="*", color=C["obs"],
            mfc="white" if not sigp(_r.p) else C["obs"], ms=5.5, lw=0.7, zorder=6)
ax.set_ylim(-60, 60); ax.set_ylabel("km per decade"); ax.set_xlim(*LONX); ax.set_xticks(TICK); ax.set_xticklabels(lonfmt(TICK))
ax.set_xlabel("Longitude (flow from left to right)")
ax.legend(handles=[Line2D([], [], ls="-", marker="o", color=C["pos"], ms=2.4, lw=0.6, label="half-power latitude (5° bins)"),
                   Line2D([], [], ls="-", marker="D", color="#555555", ms=2.0, lw=0.6, label="jet-core path latitude"),
                   Line2D([], [], ls="-", marker="*", color=C["obs"], mfc="white", ms=5, lw=0.6, label="half-power, LM Gould section")],
          loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "e", "Meridional displacement (positive northward)", dy=12)
S.save(fig, FIGD, "Figure3"); plt.close(fig)
print("done")
