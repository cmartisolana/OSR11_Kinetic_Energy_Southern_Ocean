"""Figure 4 (v3.41): isopycnal slope across the ACC from Argo, in the layout of Figure 3 (maps left, series right).
Roemmich-Gilson Argo climatology (1 deg, 0-2000 dbar), annual 2004-2025, TEOS-10 sigma0 (dk32).
a  Trend of the depth of sigma0 = 27.5 (dbar per decade); dots p < 0.05 in the ACC band; mean SAF, PF; 3000 m isobath.
b  Mean depth of sigma0 = 27.5 with its 500 and 1000 dbar contours for 2004-2014 (cyan) and 2015-2025 (black dashed);
   inset: ACC-band trend of the cross-ACC slope.
c  ACC-band cross-ACC isopycnal slope (m per 100 km; least squares of the isopycnal depth on latitude across the band per
   1 deg longitude, averaged), annual, sigma0 = 27.3, 27.5, 27.7, with 2004-2025 trends (95 %).
d  Trend of the cross-ACC slope of sigma0 = 27.5 by 10 deg longitude (95 %), with mean seafloor depth (right).
e  Trends 2004-2025 by sector (% per decade): wind stress, cross-ACC slope (27.5) and EKE; least squares, 95 % with n_eff
   from the lag-1 autocorrelation of the residuals; open: p >= 0.05.
Usage: python fig4_v341.py <revision dir> <drake out dir> <work out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
from scipy.ndimage import gaussian_filter
REV, DO, WO, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
I = np.load(WO / "dk32" / "dk32_iso.npz"); SLL = np.load(WO / "dk32" / "dk32_slope_lon.npz"); Z31 = np.load(WO / "dk31" / "dk31_series.npz")
DS = pd.read_csv(WO / "dk32" / "dk32_slopes.csv")
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
TB2 = pd.read_csv(DO / "dk12" / "dk12_band_trends_allsat.csv").sort_values("lon")
ice_ds = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc")
bc = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation.coarsen(latitude=4, longitude=4, boundary="trim").mean()
yrs = I["years"]; lat, lon = I["lat"], I["lon"]; band = I["band"]; t = yrs - yrs.mean()
CW, CE = C["wind"], C["eke"]; CS = {27.3: "#6A51A3", 27.5: "#2171B5", 27.7: "#6BAED6"}; M100 = 100 / 111.195   # dbar per deg -> m per 100 km
FEAT = [(75, -47.5, "Kerguelen\nPlateau"), (122, -46.0, "Southeast\nIndian Ridge"), (166, -49.0, "Campbell\nPlateau"),
        (-140, -60.5, "Pacific–Antarctic\nRidge"), (-66, -60.0, "Drake\nPassage")]
OBS = [(-142, "PAR"), (-65, "Drake"), (-40, "Scotia"), (37, "SWIR"), (78, "Kerguelen"), (130, "SEIR"), (160, "MR")]
def trend(y, axis=0):
    y = np.moveaxis(np.asarray(y, float), axis, 0); tt = t.reshape((-1,) + (1,) * (y.ndim - 1)); m = y.mean(0)
    b = (tt * (y - m)).sum(0) / (t ** 2).sum(); r = y - m - tt * b
    r1 = (r[1:] * r[:-1]).sum(0) / np.maximum((r * r).sum(0), 1e-12); ne = np.clip(len(t) * (1 - r1) / (1 + r1), 3, len(t))
    se = np.sqrt((r ** 2).sum(0) / (ne - 2) / (t ** 2).sum()); p = 2 * stats.t.sf(np.abs(b / se), ne - 2)
    return 10 * b, 10 * stats.t.ppf(0.975, ne - 2) * se, p
lon180 = np.where(lon > 180, lon - 360, lon); o = np.argsort(lon180); LO = lon180[o]
def mapbase(ax, fc="#1A1A1A"):
    S.polar(ax, ccrs, cfeature)
    ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors=fc, linewidths=0.6, linestyles="solid", transform=PC, zorder=6)
    add_ice_mask(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_ds.ice_mask.values > 0.5)
fig = plt.figure(figsize=(S.W2, 0.98 * S.W2))
gR = fig.add_gridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.24)
g1 = gR[0].subgridspec(2, 1, hspace=0.16)
gC = gR[1].subgridspec(3, 1, height_ratios=[1.0, 1.15, 0.9], hspace=0.45)
# a
Pz = I["p27.5"]; b, ci, p = trend(Pz)
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); mapbase(ax)
im = ax.pcolormesh(LO, lat, np.ma.masked_invalid(b[:, o]), transform=PC, cmap=cm_, vmin=-40, vmax=40, shading="nearest", rasterized=True, zorder=0)
jj, ii = np.nonzero((p[:, o] < 0.05) & band[:, o]); ax.scatter(LO[ii], lat[jj], s=1.0, c="k", lw=0, transform=PC, zorder=7)
for lo_, la_, s_ in FEAT:
    ax.text(lo_, la_, s_, transform=PC, fontsize=4.6, ha="center", va="center", color="#222222", zorder=9, linespacing=0.9,
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
S.plabel(ax, "a", "Depth trend of σ$_0$ = 27.5 kg m$^{-3}$, 2004–2025")
S.hcbar(fig, im, ax, "dbar per decade (positive: deeper)", shrink=0.7, extend="both")
# b
ax = fig.add_subplot(g1[1], projection=ccrs.SouthPolarStereo()); mapbase(ax, fc="#9A9A9A")
pm = np.nanmean(Pz, 0)
imd = ax.pcolormesh(LO, lat, np.ma.masked_invalid(pm[:, o]), transform=PC, cmap="magma_r", vmin=0, vmax=2000, shading="nearest", rasterized=True, zorder=0)
CEo, CLt = "#00B4E6", "#111111"; e_, l_ = yrs <= 2014, yrs >= 2015
for msk, col, ls, lw in ((e_, CEo, "-", 1.0), (l_, CLt, (0, (2.5, 1.2)), 0.55)):
    f = gaussian_filter(np.where(np.isfinite(np.nanmean(Pz[msk], 0)), np.nanmean(Pz[msk], 0), 2500), 0.6)[:, o]
    ax.contour(LO, lat, f, levels=[500, 1000], colors=col, linewidths=lw, linestyles=[ls], transform=PC, zorder=5)
r0 = DS[(DS.region == "ACC band") & (DS.sigma0 == 27.5)].iloc[0]
ax.legend(handles=[Line2D([], [], color=CEo, lw=1.0, label="2004–2014"), Line2D([], [], color=CLt, ls=(0, (2.5, 1.2)), lw=0.9, label="2015–2025")],
          title="500 and 1000 dbar", title_fontsize=4.8, loc="lower left", bbox_to_anchor=(-0.04, -0.04), fontsize=4.8)
ax.text(0.99, -0.03, f"cross-ACC slope trend\n{r0.slope_trend_pct:+.1f} % per decade\n({r0.slope_trend_pct - r0.slope_ci_pct:+.1f} to {r0.slope_trend_pct + r0.slope_ci_pct:+.1f})".replace("-", "−"),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=4.8, color="#222222", zorder=12, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
S.plabel(ax, "b", "Depth of σ$_0$ = 27.5 kg m$^{-3}$, two periods")
S.hcbar(fig, imd, ax, "Mean depth, 2004–2025 (dbar)", shrink=0.7, extend="max")
# c ACC-band slope series in m per 100 km
ax = fig.add_subplot(gC[0])
for s in (27.3, 27.5, 27.7):
    sl = SLL[f"sl{s}"]; cols = np.isfinite(sl).all(0); y = M100 * sl[:, cols].mean(1); bb, cc, pp = trend(y)
    ax.plot(yrs + 0.5, y, "-o", color=CS[s], lw=1.1 if s == 27.5 else 0.7, ms=1.8, zorder=3,
            label=f"σ$_0$ = {s}: {bb:+.2f} ({bb - cc:+.2f} to {bb + cc:+.2f}) m per 100 km per decade".replace("-", "−"))
    xx = np.array([2004, 2026]); ax.plot(xx, y.mean() + bb / 10 * (xx - (yrs.mean() + 0.5)), color=CS[s], lw=0.7, ls="-" if pp < 0.05 else (0, (3, 1.6)), zorder=2)
S.ts_style(ax, (2003.8, 2026.2)); ax.set_ylabel("Slope (m per 100 km)")
ax.set_ylim(80, 150); ax.legend(loc="upper left", fontsize=4.8, ncol=1)
S.plabel(ax, "c", "Cross-ACC isopycnal slope averaged over the ACC band", dy=4)
# d slope trend by 10 deg longitude, with depth
ax = fig.add_subplot(gC[1])
axd = ax.twinx(); axd.fill_between(TB2.lon, TB2.depth_mean / 1000, 7, color="#E3E3E3", step="mid", lw=0, zorder=0)
axd.set_ylim(7, -3.5); axd.set_yticks([0, 1, 2, 3, 4, 5, 6]); axd.set_ylabel("Depth (km)", color="#777777", fontsize=5.5)
axd.tick_params(axis="y", colors="#777777", labelsize=5); axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
l = np.where(SLL["lon"] > 180, SLL["lon"] - 360, SLL["lon"]); edges = np.arange(-180, 181, 10)
sl = SLL["sl27.5"]; xs, ys, lo_, hi_, ps = [], [], [], [], []
for a_, b_ in zip(edges[:-1], edges[1:]):
    c = (l >= a_) & (l < b_) & np.isfinite(sl).all(0)
    if c.sum() < 4: continue
    y = sl[:, c].mean(1); bb, cc, pp = trend(100 * (y / y.mean() - 1)); xs.append((a_ + b_) / 2); ys.append(bb); lo_.append(bb - cc); hi_.append(bb + cc); ps.append(pp)
xs, ys, ps = map(np.array, (xs, ys, ps))
ax.fill_between(xs, lo_, hi_, color=CS[27.5], alpha=0.18, lw=0, zorder=1); ax.plot(xs, ys, color=CS[27.5], lw=0.8, zorder=3)
sg = ps < 0.05; ax.plot(xs[sg], ys[sg], "o", color=CS[27.5], ms=1.8, zorder=4); ax.plot(xs[~sg], ys[~sg], "o", mfc="white", color=CS[27.5], ms=1.8, mew=0.4, zorder=4)
S.zero(ax); ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 60)); ax.set_xticklabels(["180°", "120°W", "60°W", "0°", "60°E", "120°E", "180°"])
ax.set_xlabel("Longitude (flow from left to right)"); ax.set_ylabel("Slope trend\n(% per decade)"); ax.set_ylim(-20, 20)
for x_, s_ in OBS: ax.text(x_, 0.99, s_, ha="center", va="top", fontsize=4.8, color="#222222", transform=ax.get_xaxis_transform())
S.plabel(ax, "d", "Cross-ACC slope trend of σ$_0$ = 27.5 (left) and mean seafloor depth (right)", dy=4)
# e sector bars
ax = fig.add_subplot(gC[2]); S.zero(ax); w_ = 0.24; m4 = np.arange(2001, 2026) >= 2004
REGS = ["ACC band", "Pacific", "Indian", "Atlantic"]
for gi, rg in enumerate(REGS):
    q = DS[(DS.region == rg) & (DS.sigma0 == 27.5)].iloc[0]
    items = [trend(Z31[f"{rg}|tau_raw"][m4]) + (CW,), (q.slope_trend_pct, q.slope_ci_pct, q.slope_p, CS[27.5]), trend(Z31[f"{rg}|eke_raw"][m4]) + (CE,)]
    for k, (v, c_, p_, col) in enumerate(items):
        x = gi + (k - 1) * w_
        ax.bar(x, v, w_ * 0.9, color=col if p_ < 0.05 else "white", edgecolor=col, lw=0.8, zorder=3)
        ax.vlines(x, v - c_, v + c_, color="#333333", lw=0.5, zorder=4)
ax.set_xticks(range(4)); ax.set_xticklabels(REGS); ax.tick_params(axis="x", length=0)
ax.set_ylabel("Trend 2004–2025\n(% per decade)"); ax.set_ylim(-4, 20)
ax.legend(handles=[Patch(facecolor=CW, label="wind stress"), Patch(facecolor=CS[27.5], label="isopycnal slope (σ$_0$ = 27.5)"), Patch(facecolor=CE, label="EKE"),
                   Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper left", ncol=4, fontsize=4.7)
S.plabel(ax, "e", "Wind stress, isopycnal slope and EKE trends by sector", dy=4)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
