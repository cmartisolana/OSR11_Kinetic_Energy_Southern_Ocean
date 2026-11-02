"""Figure 3 (v2.1, 28 Oct 2026, third pass): EKE and the pathway, Product 1 (all-satellite), 2001-2025. No wind (author request:
the wind and the isopycnal slope are in Figure 4, the eddy-saturation figure). Maps on top, c and d full width.
 a  EKE trend (cm2 s-2 per decade; a1_maps_bc), dots: FDR 10 %; mean SAF and PF, 3000 m isobath, sea ice.
 b  Occurrence of half-power points 2001-2025 (dk14) with the 8 % contours of 2001-2015 and 2016-2025 and the circumpolar
    half-power latitude trend (retest baseline|ACC|lat_hp).
 c  ACC band (mean ADT -0.6 to +0.2 m): EKE and KE, monthly (thin) and 12-month running mean (thick), with the 2001-2025 trends
    (retest A5; solid: p < 0.05).
 d  EKE trend between the mean SAF and PF in 2 deg longitude bins (cm2 s-2 per decade; dk12, shading 95 %) with mean seafloor
    depth (right) and the main bathymetric features.
The earlier-estimates panel (former 3e) is SI Fig. S1d.
v2.1 (29 Oct 2026, author request, as Fig. 2): (b) 10 % contours (2016-2025 line thicker) with the area overlap and mean
  latitude shift of the two contours added to the trend box; (c) trends written at the right of the plot as in Fig. 4c
  (no legend); (d) mean seafloor depth between the fronts at the native 1/8 deg resolution, dark fill in a bottom strip.
Usage: python fig3_circumpolar.py <revision dir> <figure dir>     (F4_OFFLINE=1: land from the bathymetry)
"""
import sys, os, types
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

REV, FIGD = (Path(a) for a in sys.argv[1:3])
sys.path.insert(0, str(REV))
import style_nc as S  # noqa: E402
from style_nc import add_ice_mask  # noqa: E402
import cartopy.crs as ccrs, cartopy.feature as cfeature  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
OUT = REV / "outputs"; ROOT_OUT = REV.parent / "outputs"
OFFLINE = os.environ.get("F4_OFFLINE") == "1"
_cf = types.SimpleNamespace(LAND=cfeature.ShapelyFeature([], PC)) if OFFLINE else cfeature
FEAT = [(75, -47.5, "Kerguelen\nPlateau"), (122, -46.0, "Southeast\nIndian Ridge"), (166, -49.0, "Campbell\nPlateau"),
        (-140, -60.5, "Pacific–Antarctic\nRidge"), (-66, -60.0, "Drake\nPassage")]
OBS = [(-142, "PAR"), (-65, "Drake"), (-40, "Scotia"), (37, "SWIR"), (78, "Kerguelen"), (130, "SEIR"), (160, "MR")]
LONX = (-180, 180); TICK = np.arange(-180, 181, 60); CE, CK, CW = C["eke"], C["ke"], C["wind"]


def sigp(p):
    return np.isfinite(p) and round(float(p), 3) < 0.05


def lonfmt(ticks):
    return [("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in ticks]


def edges(c):
    d = np.diff(c).mean(); return np.r_[c - d / 2, c[-1] + d / 2]


bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation; bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
clim = xr.open_dataset(OUT / "clim_allsat.nc")
ice_ds = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc"); ice_g = ice_ds.ice_mask.values > 0.5
lat8, lon8 = ice_ds.latitude.values, ice_ds.longitude.values
RT = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")
M = xr.open_dataset(OUT / "a1_maps_bc.nc"); KM = 111.195
DN = np.load(REV / "plan1019" / "drake" / "out" / "dk14" / "dk14_hp_density_circ.npz")
DK = REV / "plan1019" / "drake" / "out" / "dk12"; TB2 = pd.read_csv(DK / "dk12_band_trends_allsat.csv").sort_values("lon")
Z = np.load(DK / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z["lev_saf"]), float(Z["lev_pf"])


def base_map(ax):
    S.polar(ax, ccrs, _cf)
    if OFFLINE:
        ax.contourf(bc.longitude, bc.latitude, (bc.values >= 0).astype(float), levels=[0.5, 1.5], colors=["#3A3A3A"], transform=PC, zorder=3)
    ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[LEV_PF, LEV_SAF], colors="#1A1A1A", linewidths=0.6, transform=PC, zorder=6)
    add_ice_mask(ax, ccrs, lon8, lat8, ice_g)


fig = plt.figure(figsize=(S.W2, 1.04 * S.W2))
gA = fig.add_gridspec(3, 1, height_ratios=[1.5, 0.72, 0.78], hspace=0.26)
g1 = gA[0].subgridspec(1, 2, wspace=0.10)
gC = [gA[1], gA[2]]
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")

# a EKE trend map
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); base_map(ax)
im = ax.pcolormesh(edges(M.lon.values), edges(M.lat.values), 10 * M["EKE_slope"].values, transform=PC, cmap=tc, vmin=-60, vmax=60,
                   shading="flat", rasterized=True, zorder=0)
jj, ii = np.nonzero(M["EKE_fdr"].values == 1)
ax.scatter(M.lon.values[ii], M.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
for lo_, la_, s_ in FEAT:
    ax.text(lo_, la_, s_, transform=PC, fontsize=4.6, ha="center", va="center", color="#222222", zorder=9, linespacing=0.9,
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
S.plabel(ax, "a", "EKE trend, 2001–2025")
S.hcbar(fig, im, ax, "EKE trend (cm$^2$ s$^{-2}$ per decade)", shrink=0.7, extend="both")

# b half-power occurrence
from scipy.ndimage import gaussian_filter
ax = fig.add_subplot(g1[1], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, _cf)
if OFFLINE:
    ax.contourf(bc.longitude, bc.latitude, (bc.values >= 0).astype(float), levels=[0.5, 1.5], colors=["#3A3A3A"], transform=PC, zorder=3)
dd = np.where(DN["all"] > 0, DN["all"], np.nan)
lc = 0.5 * (DN["lon_edges"][1:] + DN["lon_edges"][:-1]); ac = 0.5 * (DN["lat_edges"][1:] + DN["lat_edges"][:-1])
imd = ax.pcolormesh(DN["lon_edges"], DN["lat_edges"], dd, cmap="magma_r", vmin=0, vmax=15, shading="flat", transform=PC, rasterized=True, zorder=0)
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
CEa, CLa = "#00B4E6", "#111111"
LEVB = 10; SM = {k: gaussian_filter(DN[k], 1.0) for k in ("early", "late")}
for key, col, ls, lw in (("early", CEa, "-", 1.1), ("late", CLa, (0, (2.5, 1.2)), 0.8)):
    ax.contour(lc, ac, SM[key], levels=[LEVB], colors=col, linewidths=lw, linestyles=[ls], transform=PC, zorder=4 if key == "early" else 5)
# descriptive difference of the two contours (as Fig. 2b): overlap of the enclosed areas (intersection over union, cos-weighted)
#   and mean over longitudes of the change in occurrence-weighted latitude inside either area (late minus early; + = northward)
_me, _ml = SM["early"] >= LEVB, SM["late"] >= LEVB; _w = np.cos(np.radians(ac))[:, None] * np.ones_like(_me, float)
IOU = (_w * (_me & _ml)).sum() / (_w * (_me | _ml)).sum(); _U = _me | _ml; _sh = []
for j in range(lc.size):
    m_ = _U[:, j]; de_, dl_ = DN["early"][m_, j], DN["late"][m_, j]
    if m_.any() and de_.sum() > 0 and dl_.sum() > 0: _sh.append((ac[m_] * dl_).sum() / dl_.sum() - (ac[m_] * de_).sum() / de_.sum())
SHIFT = KM * np.mean(_sh); print(f"Fig 3b: 10 % areas IoU = {IOU:.2f}; mean latitude shift = {SHIFT:+.1f} km over {len(_sh)} columns")
add_ice_mask(ax, ccrs, lon8, lat8, ice_g)
r0 = RT.loc["baseline|ACC|lat_hp"]
ax.legend(handles=[Line2D([], [], color=CEa, lw=1.1, label="2001–2015"), Line2D([], [], color=CLa, ls=(0, (2.5, 1.2)), lw=0.9, label="2016–2025")],
          title=f"{LEVB} % contour", title_fontsize=4.8, loc="lower left", bbox_to_anchor=(-0.04, -0.04), fontsize=4.8)
ax.text(0.99, -0.03, (f"area overlap {100 * IOU:.0f} %\nmean shift {SHIFT:+.1f} km\nhalf-power latitude trend\n{10 * KM * r0.slope:+.1f} km per decade\n"
                     f"({10 * KM * r0.ci_lo:+.1f} to {10 * KM * r0.ci_hi:+.1f})").replace("-", "−"),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=4.8, color="#222222", zorder=12,
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
S.plabel(ax, "b", "Occurrence of half-power points, 2001–2025")
S.hcbar(fig, imd, ax, "Half-power points (% of months per cell)", shrink=0.7, extend="max")

# c ACC-band series
A5 = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv"); tA = A5.x.values; LABC = {"KE": 300, "EKE": 255}
def fp(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"
ax = fig.add_subplot(gC[0])
for v, lab, col in (("ACC_EKE", "EKE", CE), ("ACC_KE", "KE", CK)):
    y = A5[v].values; r = RT.loc[f"allsat|{v}|2001-2025"]
    ax.plot(tA, y, color=col, lw=0.35, alpha=0.45, zorder=2)
    ax.plot(tA, pd.Series(y).rolling(12, center=True).mean(), color=col, lw=1.1, zorder=3,
            )
    xx = np.array([2001, 2026]); ax.plot(xx, y.mean() + r.slope * (xx - tA.mean()), color=col, lw=0.8, ls="-" if sigp(r.p_new) else (0, (3, 1.6)), zorder=4)
    _m = 10 * r.slope / (r.pct_per_decade / 100)                     # period mean implied by the relative trend
    ax.text(2026.5, LABC[lab], f"{lab}\n{r.pct_per_decade:+.1f} % per decade\n({1000 * r.ci_lo / _m:+.1f} to {1000 * r.ci_hi / _m:+.1f}), {fp(r.p_new)}".replace("-", "−"),
            color=col, fontsize=4.9, ha="left", va="center", linespacing=1.0)
S.ts_style(ax, (2000.8, 2031.6)); ax.spines["bottom"].set_bounds(2000.8, 2026.2); ax.set_ylabel("cm$^2$ s$^{-2}$")
_all = np.r_[A5.ACC_EKE.values, A5.ACC_KE.values]; _lo, _hi = np.nanmin(_all), np.nanmax(_all); _pad = 0.06 * (_hi - _lo)
ax.set_ylim(_lo - _pad, _hi + _pad); ax.set_yticks(np.arange(np.ceil(_lo / 25) * 25, _hi + _pad, 25))
for yg in ax.get_yticks(): ax.plot([2000.8, 2026.2], [yg, yg], color="#EDEDED", lw=0.4, zorder=0)
S.plabel(ax, "c", "EKE and KE averaged over the ACC band", dy=4)

# d trends by longitude between the fronts, % per decade
ax = fig.add_subplot(gC[1]); axd = ax.twinx()
_ab = clim.adt.interp(latitude=bat.latitude, longitude=bat.longitude).values   # mean ADT on the 1/8 deg bathymetry grid
_band = (_ab >= LEV_PF) & (_ab <= LEV_SAF) & (bat.values < 0)
with np.errstate(all="ignore"): BDEP = np.nanmean(np.where(_band, -bat.values, np.nan), 0) / 1000
axd.fill_between(bat.longitude.values, BDEP, 7, color="#3A3A3A", lw=0, zorder=0)       # bottom strip, below the EKE data
axd.set_ylim(6.2, -6.8); axd.set_yticks([2, 3, 4, 5, 6]); axd.set_ylabel("Depth\n(km)", color="#333333", fontsize=5.2)
axd.yaxis.set_label_coords(1.035, 0.15); axd.tick_params(axis="y", colors="#333333", labelsize=5)
axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5); axd.spines["right"].set_bounds(6.2, 1.2)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
# EKE only (author request 28 Oct 2026: the wind is not drawn here; panel d shows where along the path EKE increases,
# relative to the bathymetric features, for the standing-meander / eddy-saturation argument)
T2 = TB2.dropna(subset=["slope_dec"])
ax.fill_between(T2.lon, T2.lo, T2.hi, color=CE, alpha=0.30, lw=0, zorder=1)
ax.plot(T2.lon, T2.slope_dec, color=CE, lw=0.8, zorder=3); sg = T2.p.round(3) < 0.05
ax.plot(T2.lon[sg], T2.slope_dec[sg], "o", color=CE, ms=1.8, zorder=4); ax.plot(T2.lon[~sg], T2.slope_dec[~sg], "o", mfc="white", color=CE, ms=1.8, mew=0.4, zorder=4)
S.zero(ax); ax.set_xlim(*LONX); ax.set_xticks(TICK); ax.set_xticklabels(lonfmt(TICK)); ax.set_xlabel("Longitude (flow from left to right)")
ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ per decade)"); ax.set_ylim(-240, 215); ax.set_yticks([-100, -50, 0, 50, 100, 150, 200])
ax.spines["left"].set_bounds(-100, 215)
for x_, s_ in OBS: ax.text(x_, 0.99, s_, ha="center", va="top", fontsize=4.8, color="#222222", transform=ax.get_xaxis_transform())
ax.annotate("", xy=(175, 0.80), xytext=(155, 0.80), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"))
ax.text(165, 0.82, "flow", ha="center", va="bottom", fontsize=5, color="#555555", transform=ax.get_xaxis_transform())
S.plabel(ax, "d", "EKE trend between the SAF and PF in 2° bins; mean seafloor depth (bottom)", dy=4)
S.save(fig, FIGD, "Figure3"); plt.close(fig)
print("done")
