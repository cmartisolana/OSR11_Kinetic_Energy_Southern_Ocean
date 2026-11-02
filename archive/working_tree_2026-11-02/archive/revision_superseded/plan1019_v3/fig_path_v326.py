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

# ================================================================== Figure 4 (v3.26): pathway and standing meanders
# Redesigned after the referee comments: a where half-power points occur (two periods), b position trends by longitude
# with detection coverage, c EKE and jet speed along the path, d sinuosity, e depth and meander regions,
# f the standing-meander chain at five meanders (dk13). The circumpolar-mean half-power series is in the text and S3.
DN = np.load(DO / "dk14" / "dk14_hp_density_circ.npz")
CH = pd.read_csv(DO / "dk13" / "dk13_meander_chain_allsat.csv")
LB = RT[RT.index.str.startswith("lonbin|")].copy(); LB["lon0"] = LB.index.str.split("|").str[1].astype(float); LB["c"] = LB.lon0 + 2.5
LB["tr"] = LB.slope * KM * 10; LB["lo"] = LB.ci_lo * KM * 10; LB["hi"] = LB.ci_hi * KM * 10
CL = pd.read_csv(OUT / "b3" / "b3_climatology.csv")
T = {v: pd.read_csv(OUT / "b3" / f"b3_trends_{v}.csv") for v in ("eke", "speed", "sinu", "lat")}
T["kappa"] = pd.read_csv(DO / "b4" / "b4_trends_rms_kappa.csv")
SHALLOW = CL.lon[CL.depth_min < 3000].values
DF = pd.read_csv(OUT / "a2" / "detection_frequency_lon.csv"); DF["b"] = np.floor(DF.lon / 5) * 5 + 2.5; DFb = DF.groupby("b").freq.mean()
_dk2 = pd.read_csv(DO / "dk2" / "dk2_trends.csv"); _r = _dk2[(_dk2["product"] == "allsat") & (_dk2.period == "2001-2025")].iloc[0]
bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
MEANDERS = {"Kerguelen Plateau": ((55, 65), (65, 80), (80, 95)), "Southeast Indian Ridge": ((120, 130), (130, 140), (140, 150)),
            "Macquarie Ridge": ((150, 155), (155, 160), (160, 175)), "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
            "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35))}
SHORT = {"Kerguelen Plateau": "KP", "Southeast Indian Ridge": "SEIR", "Macquarie Ridge": "MR", "Pacific–Antarctic Ridge": "PAR",
         "Drake Passage–Scotia Sea": "DP"}
UCOL, RCOL, DCOL = "#56B4E9", "#555555", "#E69F00"

fig = plt.figure(figsize=(S.W2, 1.3 * S.W2))
gs = fig.add_gridspec(5, 1, height_ratios=[2.3, 1, 1, 0.9, 0.6], hspace=0.5)
gtop = gs[0].subgridspec(1, 2, width_ratios=[1.0, 1.1], wspace=0.25)
# a
ax = fig.add_subplot(gtop[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
dd = np.where(DN["all"] > 0, DN["all"], np.nan)
lc = 0.5 * (DN["lon_edges"][1:] + DN["lon_edges"][:-1]); ac = 0.5 * (DN["lat_edges"][1:] + DN["lat_edges"][:-1])
im = ax.pcolormesh(DN["lon_edges"], DN["lat_edges"], dd, cmap="magma_r", vmin=0, vmax=15, shading="flat", transform=PC, rasterized=True, zorder=0)
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
from scipy.ndimage import gaussian_filter
for key, col, ls in (("early", "#56B4E9", (0, (3, 1.5))), ("late", "#0072B2", "-")):
    ax.contour(lc, ac, gaussian_filter(DN[key], 1.0), levels=[8], colors=col, linewidths=0.45, linestyles=[ls], transform=PC, zorder=4)
ax.plot(np.r_[CL.lon, CL.lon.iloc[0] + 360], np.r_[CL.lat_mean, CL.lat_mean.iloc[0]], color="#D7191C", lw=0.6, transform=PC, zorder=6)
for nm, (u, r_, d_) in MEANDERS.items():
    ax.text(np.mean(r_), -37.5, SHORT[nm], transform=PC, fontsize=5.0, ha="center", va="center", fontweight="bold", zorder=9,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
ax.legend(handles=[Line2D([], [], color="#56B4E9", ls=(0, (3, 1.5)), lw=0.8, label="8 %, 2001–2015"),
                   Line2D([], [], color="#0072B2", lw=0.8, label="8 %, 2016–2025"), Line2D([], [], color="#D7191C", lw=0.8, label="jet-core path")],
          loc="upper left", bbox_to_anchor=(-0.12, 1.0), fontsize=4.6)
S.hcbar(fig, im, ax, "Half-power points (% of months per cell)", shrink=0.7, extend="max", pad=0.05)
S.plabel(ax, "a", "Occurrence of half-power points, 2001–2025")
# f: meander chain
ax = fig.add_subplot(gtop[1])
KCOL = "#CC79A7"
LINKS = [("upstream speed", UCOL, "o"), ("sinuosity, ridge and downstream", RCOL, "s"), ("curvature, ridge and downstream", KCOL, "v"),
         ("downstream EKE", C["eke"], "D"), ("downstream speed", C["tr"], "^")]
names = list(MEANDERS)
for gi, nm in enumerate(names):
    for li, (lk, col, mk) in enumerate(LINKS):
        r_ = CH[(CH.meander == nm) & (CH.link == lk)].iloc[0]; y = gi + (li - 2) * 0.16
        ax.hlines(y, r_.lo, r_.hi, color=col, lw=0.7)
        ax.plot(r_.pct, y, mk, color=col, mfc=col if sigp(r_.p) else "white", ms=3.0, mew=0.6, zorder=4)
        if r_.fdr10: ax.plot(r_.pct, y, "o", mfc="none", mec="#111111", ms=6, mew=0.6, zorder=5)
S.zero(ax, "x"); ax.set_yticks(range(len(names))); ax.set_yticklabels([n.replace("–", "–\n") if len(n) > 20 else n for n in names], fontsize=5.5)
ax.invert_yaxis(); ax.set_xlabel("Trend 2001–2025 (% per decade)"); ax.tick_params(axis="y", length=0)
ax.legend(handles=[Line2D([], [], ls="-", marker=mk, color=col, ms=3, lw=0.7, label=lk) for lk, col, mk in LINKS],
          loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "f", "Standing-meander chain at five meanders", dy=24)

def bands(ax):
    for l in SHALLOW: ax.axvspan(l - 2.5, l + 2.5, color="#EDEDED", lw=0, zorder=0)
    ax.axvspan(-60, -50, facecolor="none", edgecolor="#BDBDBD", hatch="////", lw=0, zorder=0)
def pts(ax, x, y, lo, hi, p, col, ms=2.6, mk="o", fdr=None, dx=0.0):
    sg = np.round(p, 3) < 0.05
    ax.vlines(x + dx, lo, hi, color=col, lw=0.6)
    ax.plot((x + dx)[sg], y[sg], mk, color=col, ms=ms, zorder=4)
    ax.plot((x + dx)[~sg], y[~sg], mk, mfc="white", color=col, ms=ms, mew=0.6, zorder=4)
    if fdr is not None: ax.plot((x + dx)[fdr], y[fdr], "o", mfc="none", mec="#111111", ms=5.5, mew=0.7, zorder=5)
TICK = np.arange(-180, 181, 30)
axs = [fig.add_subplot(gs[k]) for k in range(1, 5)]
# b: position trends + detection coverage
ax = axs[0]; bands(ax); S.zero(ax)
axc = ax.twinx(); axc.bar(DFb.index, 100 * DFb.values, width=5, color=C["pos"], alpha=0.10, lw=0, zorder=0)
axc.set_ylim(0, 100); axc.set_ylabel("Detection (% of months)", color="#9C6A92", fontsize=5.5); axc.spines["right"].set_visible(True); axc.tick_params(colors="#9C6A92", labelsize=5)
ax.set_zorder(axc.get_zorder() + 1); ax.patch.set_visible(False)
pts(ax, LB.c.values, LB.tr.values, LB.lo.values, LB.hi.values, LB.p_new.values, C["pos"], dx=-0.8)
d = T["lat"].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]
pts(ax, d.lon.values, d.tr.values, d.lo.values, d.hi.values, d.p.values, "#555555", mk="D", ms=2.2, dx=0.8, fdr=d.fdr10.astype(bool).values if "fdr10" in d else None)
ax.errorbar([-64], [_r.trend_km_dec], yerr=[[_r.trend_km_dec - _r.ci_lo], [_r.ci_hi - _r.trend_km_dec]], fmt="*", color=C["obs"],
            mfc="white" if not sigp(_r.p) else C["obs"], ms=5.5, lw=0.7, zorder=6)
ax.set_ylim(-60, 60); ax.set_ylabel("km per decade")
ax.legend(handles=[Line2D([], [], ls="-", marker="o", color=C["pos"], ms=2.6, lw=0.6, label="half-power latitude (5° bins)"),
                   Line2D([], [], ls="-", marker="D", color="#555555", ms=2.2, lw=0.6, label="jet-core path latitude"),
                   Line2D([], [], ls="-", marker="*", color=C["obs"], mfc="white", ms=5, lw=0.6, label="half-power, LM Gould section")],
          loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "b", "Meridional displacement (positive northward)", dx=-28, dy=12)
# c: EKE and speed along the path
ax = axs[1]; bands(ax); S.zero(ax)
for v, col, dx, mk in (("eke", C["eke"], -0.8, "o"), ("speed", C["tr"], 0.8, "s")):
    d = T[v].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]
    pts(ax, d.lon.values, d.tr.values, d.lo.values, d.hi.values, d.p.values, col, mk=mk, ms=2.3, dx=dx, fdr=d.fdr10.astype(bool).values if "fdr10" in d else None)
ax.set_ylabel("% per decade")
ax.legend(handles=[Line2D([], [], ls="-", marker="o", color=C["eke"], ms=2.4, lw=0.6, label="EKE within 150 km of the path"),
                   Line2D([], [], ls="-", marker="s", color=C["tr"], ms=2.4, lw=0.6, label="jet speed along the path")],
          loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "c", "EKE and jet speed along the jet-core path", dx=-28, dy=12)
# d: sinuosity
ax = axs[2]; bands(ax); S.zero(ax)
for v, col, dx, mk in (("sinu", C["sin"], -0.8, "o"), ("kappa", "#CC79A7", 0.8, "v")):
    d = T[v].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]
    pts(ax, d.lon.values, d.tr.values, d.lo.values, d.hi.values, d.p.values, col, mk=mk, ms=2.3, dx=dx, fdr=d.fdr10.astype(bool).values if "fdr10" in d else None)
ax.set_ylabel("% per decade")
ax.legend(handles=[Line2D([], [], ls="-", marker="o", color=C["sin"], ms=2.4, lw=0.6, label="sinuosity"),
                   Line2D([], [], ls="-", marker="v", color="#CC79A7", ms=2.4, lw=0.6, label="RMS curvature")],
          loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "d", "Shape of the jet-core path", dx=-28, dy=12)
for ax in axs[:3]:
    ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels([])
# e: depth and meander regions
ax = axs[3]; dep = CL.sort_values("lon")
ax.fill_between(dep.lon, -dep.depth_min / 1000, -6.5, color="#C8C8C8", step="mid", lw=0, zorder=0)
ax.fill_between(dep.lon, -dep.depth / 1000, -6.5, color="#6A6A6A", step="mid", lw=0, zorder=1)
bands(ax); ax.set_ylim(-6.0, 0.6); ax.set_yticks([-5, -3, -1]); ax.set_yticklabels(["5", "3", "1"]); ax.set_ylabel("Depth (km)")
for nm, (u, r_, d_) in MEANDERS.items():
    for (a, b), col in ((u, UCOL), (r_, RCOL), (d_, DCOL)): ax.plot([a, b], [0.25, 0.25], color=col, lw=2.2, solid_capstyle="butt", zorder=3)
    ax.text(np.mean(r_), -0.35, SHORT[nm], ha="center", va="top", fontsize=5.0, fontweight="bold", zorder=4)
ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels(lonfmt(TICK)); ax.set_xlabel("Longitude (flow from left to right)")
ax.legend(handles=[Line2D([], [], color=c_, lw=2.2, label=l_) for c_, l_ in ((UCOL, "upstream"), (RCOL, "ridge"), (DCOL, "downstream"))],
          loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=3, fontsize=5.0)
S.plabel(ax, "e", "Depth under the jet-core path and meander regions (bars)", dx=-28, dy=12)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
