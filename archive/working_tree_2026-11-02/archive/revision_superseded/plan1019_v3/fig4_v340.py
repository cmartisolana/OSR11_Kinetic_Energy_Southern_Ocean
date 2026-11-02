"""Figure 4 (draft v3.40): eddy saturation from Argo + altimetry + wind, 2004-2025 (dk32, dk31).
a  Trend of the depth (pressure) of the sigma0 = 27.5 surface (dbar per decade), RG Argo climatology, ACC band cells
   outlined by the mean SAF and PF; dots p < 0.05 (n_eff from lag-1 autocorrelation).
b  Trend of the local meridional isopycnal slope of sigma0 = 27.5 (centred difference in latitude), % per decade of its
   2004-2025 mean, ACC band cells only.
c  ACC band annual anomalies (% of the 2004-2025 mean): cross-ACC isopycnal slope (27.3, 27.5, 27.7), eastward wind
   stress and EKE.
d  Trend of the cross-ACC isopycnal slope (27.5) by 20 deg longitude, with 95 % intervals.
e  Trends 2004-2025 by sector (% per decade): wind stress, isopycnal slope (27.3/27.5/27.7) and EKE; open: p >= 0.05.
Usage: python fig4_v340.py <revision dir> <drake out dir> <work out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, DO, WO, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
I = np.load(WO / "dk32" / "dk32_iso.npz"); SLL = np.load(WO / "dk32" / "dk32_slope_lon.npz"); Z31 = np.load(WO / "dk31" / "dk31_series.npz")
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
yrs = I["years"]; lat, lon = I["lat"], I["lon"]; band = I["band"]; t = yrs - yrs.mean()
CW, CE = C["wind"], C["eke"]; CS = {27.3: "#6A51A3", 27.5: "#2171B5", 27.7: "#6BAED6"}
def trend(y, axis=0):
    y = np.moveaxis(np.asarray(y, float), axis, 0); tt = t.reshape((-1,) + (1,) * (y.ndim - 1)); m = y.mean(0)
    b = (tt * (y - m)).sum(0) / (t ** 2).sum(); r = y - m - tt * b
    r1 = (r[1:] * r[:-1]).sum(0) / np.maximum((r * r).sum(0), 1e-12); ne = np.clip(len(t) * (1 - r1) / (1 + r1), 3, len(t))
    se = np.sqrt((r ** 2).sum(0) / (ne - 2) / (t ** 2).sum()); p = 2 * stats.t.sf(np.abs(b / se), ne - 2)
    return 10 * b, 10 * stats.t.ppf(0.975, ne - 2) * se, p
lon180 = np.where(lon > 180, lon - 360, lon); o = np.argsort(lon180); LO = lon180[o]
def fronts(ax): ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="#1A1A1A", linewidths=0.6, transform=PC, zorder=6)
fig = plt.figure(figsize=(S.W2, 0.98 * S.W2))
gR = fig.add_gridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.24)
g1 = gR[0].subgridspec(2, 1, hspace=0.16)
gC = gR[1].subgridspec(3, 1, height_ratios=[1.0, 0.9, 1.0], hspace=0.5)
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
# a depth trend sigma0 27.5
Pz = I["p27.5"]; b, ci, p = trend(Pz)
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
im = ax.pcolormesh(LO, lat, np.ma.masked_invalid(b[:, o]), transform=PC, cmap=cm_, vmin=-40, vmax=40, shading="nearest", rasterized=True, zorder=0)
fronts(ax); jj, ii = np.nonzero((p[:, o] < 0.05) & band[:, o]); ax.scatter(LO[ii], lat[jj], s=1.2, c="k", lw=0, transform=PC, zorder=7)
S.plabel(ax, "a", "Depth trend of σ$_0$ = 27.5 kg m$^{-3}$, 2004–2025")
S.hcbar(fig, im, ax, "dbar per decade (dots: p < 0.05, ACC band)", shrink=0.8, extend="both")
# b local meridional slope trend
G = np.gradient(Pz, lat, axis=1)    # dbar per degree; positive where the surface deepens northward
bG, ciG, pG = trend(G); rel = np.where(band & (np.abs(G.mean(0)) > 20), 100 * bG / G.mean(0), np.nan)
ax = fig.add_subplot(g1[1], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
im2 = ax.pcolormesh(LO, lat, np.ma.masked_invalid(rel[:, o]), transform=PC, cmap=cm_, vmin=-15, vmax=15, shading="nearest", rasterized=True, zorder=0)
fronts(ax); jj, ii = np.nonzero((pG[:, o] < 0.05) & np.isfinite(rel[:, o])); ax.scatter(LO[ii], lat[jj], s=1.2, c="k", lw=0, transform=PC, zorder=7)
S.plabel(ax, "b", "Isopycnal slope trend (σ$_0$ = 27.5)")
S.hcbar(fig, im2, ax, "% per decade (+: steeper; dots: p < 0.05)", shrink=0.8, extend="both")
# c ACC band series
ax = fig.add_subplot(gC[0]); S.zero(ax)
m4 = np.arange(2001, 2026) >= 2004
for key, col, lab, mk in (("tau_raw", CW, "wind stress", "s"), ("eke_raw", CE, "EKE", "o")):
    y = Z31[f"ACC band|{key}"][m4]; y = 100 * ((1 + y / 100) / (1 + y / 100).mean() - 1)
    ax.plot(yrs + 0.5, y, "-" + mk, color=col, lw=0.9, ms=1.8, zorder=3, label=lab)
for s in (27.3, 27.5, 27.7):
    sl = SLL[f"sl{s}"]; cols = np.isfinite(sl).all(0); y = sl[:, cols].mean(1); y = 100 * (y / y.mean() - 1)
    ax.plot(yrs + 0.5, y, "-", color=CS[s], lw=1.3 if s == 27.5 else 0.8, zorder=4, label=f"isopycnal slope, σ$_0$ = {s}")
S.ts_style(ax, (2003.8, 2026.2)); ax.set_ylabel("anomaly (% of mean)")
ax.legend(loc="upper left", ncol=3, fontsize=4.6)
S.plabel(ax, "c", "ACC band: wind stress, EKE and cross-ACC isopycnal slope", dy=4)
# d slope trend vs longitude (27.5)
ax = fig.add_subplot(gC[1]); S.zero(ax)
sl = SLL["sl27.5"]; l = np.where(SLL["lon"] > 180, SLL["lon"] - 360, SLL["lon"]); edges = np.arange(-180, 181, 20)
for s in (27.3, 27.5, 27.7):
    sl = SLL[f"sl{s}"]; xs, ys, ls_, hs_ = [], [], [], []
    for a_, b_ in zip(edges[:-1], edges[1:]):
        c = (l >= a_) & (l < b_) & np.isfinite(sl).all(0)
        if c.sum() < 5: continue
        y = sl[:, c].mean(1); bb, cc, pp = trend(100 * (y / y.mean() - 1)); xs.append((a_ + b_) / 2); ys.append(bb); ls_.append(bb - cc); hs_.append(bb + cc)
    xs = np.array(xs) + {27.3: -2.5, 27.5: 0, 27.7: 2.5}[s]
    if s == 27.5: ax.vlines(xs, ls_, hs_, color=CS[s], lw=0.6)
    ax.plot(xs, ys, "-o", color=CS[s], lw=1.1 if s == 27.5 else 0.6, ms=2.2 if s == 27.5 else 1.4)
ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 60)); ax.set_xticklabels(["180°", "120°W", "60°W", "0°", "60°E", "120°E", "180°"])
ax.set_ylabel("% per decade"); ax.set_xlabel("Longitude")
S.plabel(ax, "d", "Cross-ACC isopycnal slope trend by 20° longitude", dy=4)
# e sector bars
D = pd.read_csv(WO / "dk32" / "dk32_slopes.csv")
ax = fig.add_subplot(gC[2]); S.zero(ax)
REGS = ["ACC band", "Pacific", "Indian", "Atlantic"]; w_ = 0.15
for gi, rg in enumerate(REGS):
    items = []
    y = Z31[f"{rg}|tau_raw"][m4]; items.append(("wind", *trend(y), CW))
    for s in (27.3, 27.5, 27.7):
        q = D[(D.region == rg) & (D.sigma0 == s)].iloc[0]; items.append((f"s{s}", q.slope_trend_pct, q.slope_ci_pct, q.slope_p, CS[s]))
    y = Z31[f"{rg}|eke_raw"][m4]; items.append(("eke", *trend(y), CE))
    for k, (nm, v, ci_, p_, col) in enumerate(items):
        x = gi + (k - 2) * w_
        ax.bar(x, v, w_ * 0.9, color=col if p_ < 0.05 else "white", edgecolor=col, lw=0.8, zorder=3)
        ax.vlines(x, v - ci_, v + ci_, color="#333333", lw=0.5, zorder=4)
ax.set_xticks(range(4)); ax.set_xticklabels(REGS, fontsize=5.6); ax.tick_params(axis="x", length=0)
ax.set_ylabel("% per decade"); ax.set_ylim(-4, 20)
ax.legend(handles=[Patch(facecolor=CW, label="wind stress")] + [Patch(facecolor=CS[s], label=f"isopycnal slope {s}") for s in (27.3, 27.5, 27.7)] +
          [Patch(facecolor=CE, label="EKE"), Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper right", ncol=3, fontsize=4.5)
S.plabel(ax, "e", "Trends 2004–2025: wind, isopycnal slope and EKE", dy=4)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("band median depth trend 27.5:", np.nanmedian(np.where(band, b, np.nan)), "| median slope trend %:", np.nanmedian(rel))
