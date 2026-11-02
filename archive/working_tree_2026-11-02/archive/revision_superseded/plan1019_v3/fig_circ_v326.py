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


# ================================================================== Figure 3 (v3.26): circumpolar
# a EKE trend map (1 deg, FDR dots) with the mean SAF and PF of Fig. 2a; b EKE trend between the SAF and PF in 2 deg bins
# with depth (Fig. 2c around the ACC); c annual EKE between the fronts, % of each bin's mean (Fig. 2d around the ACC);
# d EKE, KE, along-stream speed (and two-satellite EKE) trends for the ACC band and sectors. No tracked eddies.
M = xr.open_dataset(OUT / "a1_maps_bc.nc")
DK = DO / "dk12"
TB2 = pd.read_csv(DK / "dk12_band_trends_allsat.csv"); HV = pd.read_csv(DK / "dk12_band_annual_allsat.csv", index_col=0)
Z = np.load(DK / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z["lev_saf"]), float(Z["lev_pf"])
lat, lon = ice_ds.latitude.values, ice_ds.longitude.values
LONX = (-180, 180); TICK = np.arange(-180, 181, 30)
OBS = [(-142, "PAR"), (-65, "Drake"), (-40, "Scotia"), (37, "SWIR"), (78, "Kerguelen"), (138, "SEIR"), (160, "Macquarie")]
fig = plt.figure(figsize=(S.W2, 1.12 * S.W2))
gR = fig.add_gridspec(3, 1, height_ratios=[1.25, 0.95, 0.8], hspace=0.28)
g1 = gR[0].subgridspec(1, 2, width_ratios=[1, 1.05], wspace=0.28)
# a
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo())
S.polar(ax, ccrs, cfeature)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
im = ax.pcolormesh(L1, A1, 10 * M["EKE_slope"].values, transform=PC, cmap=tc, vmin=-60, vmax=60, shading="nearest", rasterized=True, zorder=0)
bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[LEV_PF, LEV_SAF], colors="#1A1A1A", linewidths=0.6, linestyles="solid", transform=PC, zorder=6)
add_ice_mask(ax, ccrs, lon, lat, ice_g); add_sector_boundaries(ax, ccrs)
jj, ii = np.nonzero(M["EKE_fdr"].values == 1)
ax.scatter(M.lon.values[ii], M.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
for lo_, la_, s_ in FEAT:
    ax.text(lo_, la_, s_, transform=PC, fontsize=4.6, ha="center", va="center", color="#222222", zorder=9, linespacing=0.9,
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
for lo_, s_ in ((-110, "Pacific"), (65, "Indian"), (-15, "Atlantic")):
    ax.text(lo_, -38.5, s_, transform=PC, fontsize=5.2, ha="center", va="center", color="#555555", style="italic", zorder=9)
S.plabel(ax, "a", "EKE trend, 2001–2025")
S.hcbar(fig, im, ax, "EKE trend (cm$^2$ s$^{-2}$ per decade)", shrink=0.7, extend="both")
# d (top right): sector bars
ax = fig.add_subplot(g1[1])
groups = []
for reg, gl in (("ACC", "ACC band"), ("Pacific", "Pacific"), ("Indian", "Indian"), ("Atlantic", "Atlantic")):
    groups.append((gl, {"EKE": rt(f"allsat|{reg}_EKE|2001-2025"), "EKE2": rt(f"twosat|{reg}_EKE|2001-2025"),
                        "KE": rt(f"allsat|{reg}_KE|2001-2025"), "SPD": rt(f"allsat|{reg}_SPD|2001-2025")}))
h = bar_trends(ax, groups, [("EKE", "EKE", C["eke"]), ("EKE2", "EKE, two-satellite", C["two"]), ("KE", "KE", C["ke"]),
                            ("SPD", "along-stream speed", C["tr"])], width=0.19, ylim=(-4, 17), ylabel="Trend 2001–2025 (% per decade)")
ax.legend(handles=h, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "d", "Eddies and mean flow by sector", dy=22)
# b: along-ACC trend between the fronts + depth
g2 = gR[1].subgridspec(2, 1, height_ratios=[1, 0.38], hspace=0.08)
ax = fig.add_subplot(g2[0]); T2 = TB2.dropna(subset=["slope_dec"])
ax.fill_between(T2.lon, T2.lo, T2.hi, color=C["eke"], alpha=0.18, lw=0, zorder=1)
ax.plot(T2.lon, T2.slope_dec, color=C["eke"], lw=0.8, zorder=3)
sg = T2.p.round(3) < 0.05
ax.plot(T2.lon[sg], T2.slope_dec[sg], "o", color=C["eke"], ms=1.8, zorder=4); ax.plot(T2.lon[~sg], T2.slope_dec[~sg], "o", mfc="white", color=C["eke"], ms=1.8, mew=0.4, zorder=4)
S.zero(ax); ax.set_xlim(*LONX); ax.set_xticks(TICK); ax.set_xticklabels([]); ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ dec$^{-1}$)")
ax.annotate("", xy=(175, 0.93), xytext=(155, 0.93), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"))
ax.text(165, 0.97, "flow", ha="center", va="bottom", fontsize=5, color="#555555", transform=ax.get_xaxis_transform())
S.plabel(ax, "b", "EKE trend between the SAF and PF, and depth", dy=4)
ax = fig.add_subplot(g2[1]); TT = TB2.sort_values("lon")
ax.fill_between(TT.lon, -TT.depth_mean / 1000, -6, color="#8A8A8A", step="mid", lw=0)
ax.step(TT.lon, -TT.depth_min / 1000, where="mid", color="#222222", lw=0.6)
ax.set_ylim(-5.5, 1.6); ax.set_yticks([-5, -3, -1]); ax.set_yticklabels(["5", "3", "1"]); ax.set_ylabel("Depth\n(km)")
ax.set_xlim(*LONX); ax.set_xticks(TICK); ax.set_xticklabels(lonfmt(TICK))
for x_, s_ in OBS: ax.text(x_, 1.55, s_, ha="center", va="top", fontsize=4.8, color="#222222")
# c: Hovmoller
ax = fig.add_subplot(gR[2]); lc = HV.columns.astype(float).values; yy = HV.index.values
imh = ax.pcolormesh(np.append(lc - 1, lc[-1] + 1), np.append(yy, yy[-1] + 1), HV.values, cmap="RdBu_r", vmin=-40, vmax=40, shading="flat", rasterized=True)
ax.set_xlim(*LONX); ax.set_xticks(TICK); ax.set_xticklabels(lonfmt(TICK)); ax.set_xlabel("Longitude (flow from left to right)")
ax.set_ylim(2001, 2026); ax.set_yticks(np.array([2005, 2010, 2015, 2020, 2025]) + 0.5); ax.set_yticklabels(["2005", "2010", "2015", "2020", "2025"])
for x_, s_ in OBS: ax.axvline(x_, color="#555555", lw=0.3, ls=":")
caxh = ax.inset_axes([1.01, 0.08, 0.012, 0.84])
cbh = fig.colorbar(imh, cax=caxh, extend="both"); cbh.ax.tick_params(labelsize=5, length=2, width=0.4)
cbh.set_label("EKE anomaly (% of mean)", fontsize=5.2); cbh.outline.set_linewidth(0.4)
S.plabel(ax, "c", "Annual EKE between the SAF and PF", dy=4)
S.save(fig, FIGD, "Figure3"); plt.close(fig)
print("done")
