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

# ================================================================== Figure 4 (v3.27): transport and meander flexing
# a annual anomalies of the circumpolar surface transport proxy (a8) and along-stream speed (a5), with trends;
# b standing-meander chain at five meanders (dk13 + transport from a8b); c transport-proxy trend by 5 deg longitude (a8b);
# d shape of the jet-core path by 5 deg longitude (sinuosity b3, RMS curvature b4); e depth under the path and meander
# regions. Half-power results are in Fig. 3.
CH = pd.read_csv(DO / "dk13" / "dk13_meander_chain_allsat.csv")
TB8 = pd.read_csv(DO / "a8b" / "a8b_trends_allsat.csv")
A8 = pd.read_csv(OUT / "a8" / "a8_series_allsat.csv"); A8t = pd.read_csv(OUT / "a8" / "a8_trends_allsat.csv").set_index("series")
A5 = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv"); A5["year"] = 2001 + A5.i // 12
CL = pd.read_csv(OUT / "b3" / "b3_climatology.csv")
T = {v: pd.read_csv(OUT / "b3" / f"b3_trends_{v}.csv") for v in ("sinu",)}; T["kappa"] = pd.read_csv(DO / "b4" / "b4_trends_rms_kappa.csv")
SHALLOW = CL.lon[CL.depth_min < 3000].values
MEANDERS = {"Kerguelen Plateau": ((55, 65), (65, 80), (80, 95)), "Southeast Indian Ridge": ((120, 130), (130, 140), (140, 150)),
            "Macquarie Ridge": ((150, 155), (155, 160), (160, 175)), "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
            "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35))}
SHORT = {"Kerguelen Plateau": "KP", "Southeast Indian Ridge": "SEIR", "Macquarie Ridge": "MR", "Pacific–Antarctic Ridge": "PAR", "Drake Passage–Scotia Sea": "DP"}
UCOL, RCOL, DCOL, KCOL = "#56B4E9", "#555555", "#E69F00", "#CC79A7"
TICK = np.arange(-180, 181, 30)
def bands(ax):
    for l in SHALLOW: ax.axvspan(l - 2.5, l + 2.5, color="#EDEDED", lw=0, zorder=0)
    ax.axvspan(-60, -50, facecolor="none", edgecolor="#BDBDBD", hatch="////", lw=0, zorder=0)
def pts(ax, x, y, lo, hi, p, col, ms=2.4, mk="o", fdr=None, dx=0.0):
    sg = np.round(p, 3) < 0.05
    ax.vlines(x + dx, lo, hi, color=col, lw=0.6)
    ax.plot((x + dx)[sg], y[sg], mk, color=col, ms=ms, zorder=4); ax.plot((x + dx)[~sg], y[~sg], mk, mfc="white", color=col, ms=ms, mew=0.6, zorder=4)
    if fdr is not None: ax.plot((x + dx)[fdr], y[fdr], "o", mfc="none", mec="#111111", ms=5.5, mew=0.7, zorder=5)
def mbars(ax, y=1.0):
    for nm, (u_, r_, d_) in MEANDERS.items():
        for (a, b), col in ((u_, UCOL), (r_, RCOL), (d_, DCOL)): ax.plot([a, b], [y, y], color=col, lw=2.2, solid_capstyle="butt", transform=ax.get_xaxis_transform(), clip_on=False)
        ax.text(np.mean(r_), y + 0.04, SHORT[nm], transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.0, fontweight="bold")

fig = plt.figure(figsize=(S.W2, 1.22 * S.W2))
gs = fig.add_gridspec(4, 1, height_ratios=[1.45, 0.9, 0.8, 0.5], hspace=0.42)
gtop = gs[0].subgridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.42)
# a: circumpolar transport proxy and along-stream speed
ax = fig.add_subplot(gtop[0])
series(ax, ann(A8, "ACC_T"), C["tr"], "surface transport proxy"); r = A8t.loc["ACC_T"]; trend_line(ax, r.pct_per_decade, r.p, C["tr"])
series(ax, ann(A5, "ACC_SPD"), "#56B4E9", "along-stream speed", mk="s"); r2 = rt("allsat|ACC_SPD|2001-2025"); trend_line(ax, r2["pct"], r2["p"], "#56B4E9")
S.zero(ax); S.ts_style(ax); ax.set_ylim(-6, 6); ax.set_ylabel("Anomaly (% of mean)")
ax.legend(loc="upper left", fontsize=5.0)
S.plabel(ax, "a", "Mean flow, ACC band", dy=4)
# b: meander chain
ax = fig.add_subplot(gtop[1])
LINKS = [("upstream speed", UCOL, "o"), ("transport, ridge and downstream", C["tr"], "P"), ("sinuosity, ridge and downstream", RCOL, "s"),
         ("curvature, ridge and downstream", KCOL, "v"), ("downstream EKE", C["eke"], "D"), ("downstream speed", UCOL, "^")]
names = list(MEANDERS)
for gi, nm in enumerate(names):
    for li, (lk, col, mk) in enumerate(LINKS):
        if lk.startswith("transport"):
            r_ = TB8[TB8.key == nm].iloc[0]; pct, lo, hi, p, fd = r_.pct, r_.pct_lo, r_.pct_hi, r_.p, False
        else:
            r_ = CH[(CH.meander == nm) & (CH.link == lk)].iloc[0]; pct, lo, hi, p, fd = r_.pct, r_.lo, r_.hi, r_.p, bool(r_.fdr10)
        y = gi + (li - 2.5) * 0.13
        ax.hlines(y, lo, hi, color=col, lw=0.7); ax.plot(pct, y, mk, color=col, mfc=col if sigp(p) else "white", ms=2.8, mew=0.6, zorder=4)
        if fd: ax.plot(pct, y, "o", mfc="none", mec="#111111", ms=6, mew=0.6, zorder=5)
S.zero(ax, "x"); ax.set_yticks(range(len(names))); ax.set_yticklabels([n.replace("–", "–\n") if len(n) > 20 else n for n in names], fontsize=5.5)
ax.invert_yaxis(); ax.set_xlabel("Trend 2001–2025 (% per decade)"); ax.tick_params(axis="y", length=0)
ax.legend(handles=[Line2D([], [], ls="-", marker=mk, color=col, ms=3, lw=0.7, label=lk) for lk, col, mk in LINKS],
          loc="lower left", bbox_to_anchor=(-0.05, 1.0), ncol=2, fontsize=4.8, columnspacing=0.6)
S.plabel(ax, "b", "Standing-meander chain at five meanders", dy=34)
# c: transport by longitude
ax = fig.add_subplot(gs[1]); bands(ax); S.zero(ax)
d = TB8[TB8.lon.notna()].dropna(subset=["p"]); d = d[~d.lon.between(-60, -50)]
pts(ax, d.lon.values, d.pct.values, d.pct_lo.values, d.pct_hi.values, d.p.values, C["tr"], mk="P", fdr=d.fdr10.astype(bool).values)
ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels([]); ax.set_ylabel("% per decade"); mbars(ax)
S.plabel(ax, "c", "Surface transport proxy across the ACC band, by longitude", dx=-28, dy=14)
# d: path shape
ax = fig.add_subplot(gs[2]); bands(ax); S.zero(ax)
for v, col, dx, mk in (("sinu", C["sin"], -0.8, "o"), ("kappa", KCOL, 0.8, "v")):
    d = T[v].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]
    pts(ax, d.lon.values, d.tr.values, d.lo.values, d.hi.values, d.p.values, col, mk=mk, ms=2.3, dx=dx, fdr=d.fdr10.astype(bool).values if "fdr10" in d else None)
ax.set_ylabel("% per decade"); ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels([])
ax.legend(handles=[Line2D([], [], ls="-", marker="o", color=C["sin"], ms=2.4, lw=0.6, label="sinuosity"),
                   Line2D([], [], ls="-", marker="v", color=KCOL, ms=2.4, lw=0.6, label="RMS curvature")], loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=5.0)
S.plabel(ax, "d", "Shape of the jet-core path (flexing)", dx=-28, dy=12)
# e: depth
ax = fig.add_subplot(gs[3]); dep = CL.sort_values("lon")
ax.fill_between(dep.lon, -dep.depth_min / 1000, -6.5, color="#C8C8C8", step="mid", lw=0, zorder=0)
ax.fill_between(dep.lon, -dep.depth / 1000, -6.5, color="#6A6A6A", step="mid", lw=0, zorder=1)
bands(ax); ax.set_ylim(-6.0, 0); ax.set_yticks([-5, -3, -1]); ax.set_yticklabels(["5", "3", "1"]); ax.set_ylabel("Depth (km)")
ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels(lonfmt(TICK)); ax.set_xlabel("Longitude (flow from left to right)")
ax.legend(handles=[Line2D([], [], color=c_, lw=2.2, label=l_) for c_, l_ in ((UCOL, "upstream"), (RCOL, "ridge"), (DCOL, "downstream"))],
          loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=5.0, ncol=3)
S.plabel(ax, "e", "Depth under the jet-core path (bars in c: meander regions)", dx=-28)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
