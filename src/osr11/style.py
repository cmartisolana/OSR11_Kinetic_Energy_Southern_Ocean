"""Shared figure style for the OSR11 manuscript and Supporting Information (final revision, 2 Nov 2026).
Replaces revision/style_nc.py + revision/figstyle.py (no dependence on figstyle).
Publication repository: copied from figures_final/scripts/style.py; retest_row and bathy_path read from osr11.paths.

Canvas: every figure is drawn on a FIXED canvas 170 mm wide (Copernicus text width; W2) or 83 mm (single column; W1) and
saved without cropping (bbox_inches=None), so all figures print at the same scale and the same font sizes. Axes are placed in
millimetres (axmm), panel letters in figure coordinates (letter), so letters align in columns and rows.
Type: Arial/Helvetica (Liberation Sans fallback), ticks 6 pt, labels and titles 7 pt, panel letters 8 pt bold; nothing below 6 pt.
Colours, one per quantity (Okabe-Ito based): KE vermilion, EKE bluish green, cross-track velocity blue, half-power (frontal)
position reddish purple, isopycnal slope indigo, wind grey, in situ orange, two-satellite product grey.
Significance: filled symbols p < 0.05, open p >= 0.05; trend lines solid p < 0.05, dashed otherwise; shading or bars 95 %
intervals; map dots (white, black edge) cells passing the 10 % false-discovery-rate control.
Maps: land #3A3A3A, sea ice light blue (#8EC5FF, alpha 0.35), 3000 m isobath grey 0.3 pt, mean SAF and PF solid black 0.7 pt
(both solid: explicit linestyles, so negative ADT levels are not dashed by matplotlib), ACC band + 100 km grey dashed.
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

MM = 1 / 25.4
W2 = 170.0          # mm, full text width
W1 = 83.0           # mm, single column
FS_TICK, FS_LAB, FS_TTL, FS_LET, FS_ANN = 6.0, 7.0, 7.0, 8.0, 6.0

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS_LAB, "axes.labelsize": FS_LAB, "axes.titlesize": FS_TTL, "axes.titleweight": "normal",
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK, "legend.fontsize": FS_ANN, "legend.title_fontsize": FS_ANN,
    "legend.frameon": False, "legend.handlelength": 1.8, "legend.handletextpad": 0.5, "legend.columnspacing": 1.2,
    "legend.borderaxespad": 0.2, "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": False, "axes.edgecolor": "#222222", "axes.labelcolor": "#222222", "text.color": "#222222",
    "xtick.color": "#222222", "ytick.color": "#222222", "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.pad": 2.0, "ytick.major.pad": 2.0, "axes.labelpad": 2.5,
    "lines.linewidth": 1.0, "lines.markersize": 3.5, "patch.linewidth": 0.5, "axes.unicode_minus": True,
    "figure.dpi": 150, "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
    # mathtext in the same face as the text (stixsans fallback is avoided: 'custom' with the text font)
    "mathtext.default": "regular", "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})
from matplotlib import font_manager as _fm
if not any("Arial" in f.name for f in _fm.fontManager.ttflist):          # container / Linux: Liberation Sans = Arial metrics
    plt.rcParams.update({"mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic", "mathtext.bf": "Liberation Sans:bold"})

C = dict(ke="#D55E00", eke="#009E73", vel="#0072B2", pos="#A2508A", slope="#6A51A3", wind="#7F7F7F", obs="#E69F00",
         two="#9A9A9A", ink="#222222", grey="#7F7F7F", light="#E6E6E6", gridline="#EDEDED",
         land="#3A3A3A", ice="#8EC5FF", ice_flat="#D8EBFF", isobath="#8A8A8A", front="#1A1A1A", buffer="#555555",
         early="#00B4E6", late="#111111", masked="#F2F2F2")
LW = dict(smooth=1.0, raw=0.35, trend=0.8, profile=0.8, front=0.7, isobath=0.3, buffer=0.5, early=1.2, late=0.9)
MS = 2.2            # marker size of profiles (pt)
TREND_CMAP = "RdBu_r"; OCC_CMAP = "magma_r"
EKE_TREND_LIM = 80; OCC_LIM = 20
EKE_TREND_TICKS = [-80, -40, 0, 40, 80]; OCC_TICKS = [0, 5, 10, 15, 20]
LEVB = 10           # % occurrence contour of Figs. 2b and 3b
# bathymetric/topographic keys shared by the polar maps (rim ticks) and the longitude panels (Figs. 3d, 4): (longitude, key)
KEYS = [(-142, "PAR"), (-65, "DP"), (-40, "NSR"), (37, "SWIR"), (78, "KP"), (130, "SEIR"), (160, "MR")]
KEY_NAMES = {"PAR": "Pacific–Antarctic Ridge", "DP": "Drake Passage", "NSR": "North Scotia Ridge", "SWIR": "Southwest Indian Ridge",
             "KP": "Kerguelen Plateau", "SEIR": "Southeast Indian Ridge", "MR": "Macquarie Ridge"}


# ---------------------------------------------------------------- text helpers
def um(s):
    """true minus signs: hyphen-minus before a digit (not inside words or year ranges) -> U+2212"""
    return re.sub(r"(?<![\w\)])-(?=[\d.])", "−", s)


def fnum(x, nd=1, sign=True):
    return um(f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}")


def fp(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"


def sigp(p):
    return bool(np.isfinite(p) and round(float(p), 3) < 0.05)


def latfmt(v):
    return "0°" if v == 0 else f"{abs(v):g}°{'S' if v < 0 else 'N'}"


def lonfmt(v):
    return "180°" if abs(v) == 180 else "0°" if v == 0 else f"{abs(v):g}°{'W' if v < 0 else 'E'}"


def running(y, w=12):
    return pd.Series(np.asarray(y, float)).rolling(w, center=True, min_periods=w // 2).mean().values


def runfull(y, w=12):
    return pd.Series(np.asarray(y, float)).rolling(w, center=True).mean().values


# ---------------------------------------------------------------- canvas, axes and labels in mm
def figure(width_mm, height_mm):
    fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
    fig._mm = (width_mm, height_mm)
    return fig


def axmm(fig, x, y, w, h, **kw):
    """axes at (x, y) mm from the lower-left corner, w x h mm"""
    W, H = fig._mm
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def letter(fig, x, y, let, title=None, gap=3.6):
    """panel letter (8 pt bold) with its lower-left corner at (x, y) mm; optional title (7 pt) after it"""
    W, H = fig._mm
    fig.text(x / W, y / H, let, fontsize=FS_LET, fontweight="bold", ha="left", va="baseline")
    if title:
        fig.text((x + gap) / W, y / H, title, fontsize=FS_TTL, ha="left", va="baseline")


def cbar(fig, mappable, x, y, w, label, extend="neither", ticks=None, h=2.2, ticklabels=None):
    """one colourbar geometry for all maps: horizontal, h mm high, label below"""
    cax = axmm(fig, x, y, w, h)
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", extend=extend, extendfrac=0.04)
    cb.set_label(label, fontsize=FS_LAB, labelpad=1.5)
    cb.ax.tick_params(labelsize=FS_TICK, width=0.4, length=2, pad=1.5); cb.outline.set_linewidth(0.4)
    if ticks is not None:
        cb.set_ticks(ticks); cb.set_ticklabels(ticklabels or [um(f"{t:g}") for t in ticks])
    return cb


def zero(ax, axis="y", lim=None):
    if lim is None:
        (ax.axhline if axis == "y" else ax.axvline)(0, color="#9A9A9A", lw=0.5, zorder=0.5)
    else:
        ax.plot(lim, [0, 0], color="#9A9A9A", lw=0.5, zorder=0.5)


def hgrid(ax, ys, xlim):
    for yg in ys:
        ax.plot(xlim, [yg, yg], color=C["gridline"], lw=0.4, zorder=0)


def yticks_minus(ax):
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: um(f"{v:g}")))


def xticks_minus(ax):
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: um(f"{v:g}")))


def ts_axis(ax, x0=2000.8, x1=2026.2, xmax=None, ticks=(2005, 2010, 2015, 2020, 2025)):
    """time axis; xmax > x1 leaves room for labels at the right, the bottom spine stops at x1"""
    ax.set_xlim(x0, xmax or x1); ax.set_xticks(list(ticks)); ax.spines["bottom"].set_bounds(x0, x1)


def profile(ax, x, y, p, col, mk="o", ms=MS, lw=LW["profile"], z=3, ls="-"):
    """line with filled (p < 0.05) / open markers"""
    x, y, p = np.asarray(x, float), np.asarray(y, float), np.asarray(p, float)
    if lw: ax.plot(x, y, color=col, lw=lw, ls=ls, zorder=z)
    s = np.array([sigp(v) for v in p], bool)
    ax.plot(x[s], y[s], mk, color=col, mec=col, ms=ms, mew=0.6, ls="none", zorder=z + 1)
    ax.plot(x[~s], y[~s], mk, color=col, mfc="white", mec=col, ms=ms, mew=0.6, ls="none", zorder=z + 1)


def trend_line(ax, x0, x1, xc, yc, slope, lo, hi, p, col, wedge=True, z=4):
    """Trend through (xc, yc) with its 95 % wedge (slope, lo, hi per unit x); solid if p < 0.05, dashed otherwise"""
    xx = np.linspace(x0, x1, 60)
    if wedge:
        ax.fill_between(xx, yc + lo * (xx - xc), yc + hi * (xx - xc), color=col, alpha=0.12, lw=0, zorder=1)
    ax.plot(xx, yc + slope * (xx - xc), color=col, lw=LW["trend"], ls="-" if sigp(p) else (0, (3, 1.6)), zorder=z)


def forest(ax, rows, colors, markers=None, xlim=None, dodge=True, xlabel=None, gap=0.6, ms=3.4):
    """rows: [(label, {series: dict(pct, pct_lo, pct_hi, p)})], top first; markers filled if p < 0.05"""
    markers = markers or {}
    series = list(dict.fromkeys(s for _, d in rows for s in d))
    ns = len(series); off = {s: (i - (ns - 1) / 2) * (0.26 if dodge and ns > 1 else 0) for i, s in enumerate(series)}
    ys = -np.arange(len(rows), dtype=float)
    for (lab, d), y in zip(rows, ys):
        for s, r in d.items():
            yy = y - off[s]; c = colors[s]; m = markers.get(s, "o")
            ax.plot([r["pct_lo"], r["pct_hi"]], [yy, yy], color=c, lw=0.8, solid_capstyle="butt", zorder=2)
            sg = sigp(r["p"])
            ax.plot(r["pct"], yy, m, color=c, mfc=c if sg else "white", mec=c, mew=0.8, ms=ms, zorder=3)
    ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0)
    ax.set_yticks(ys); ax.set_yticklabels([l for l, _ in rows])
    ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_ylim(ys[-1] - gap, ys[0] + gap)
    if xlim: ax.set_xlim(*xlim)
    if xlabel: ax.set_xlabel(xlabel)
    xticks_minus(ax)
    return ys


def sig_handles(series_colors, markers=None, ns_label="open: p ≥ 0.05", ms=3.4):
    markers = markers or {}
    h = [Line2D([], [], ls="-", lw=0.8, marker=markers.get(s, "o"), color=c, mfc=c, mec=c, ms=ms, label=s)
         for s, c in series_colors.items()]
    h.append(Line2D([], [], ls="none", marker="o", color="#555555", mfc="white", mew=0.8, ms=ms, label=ns_label))
    return h


# ---------------------------------------------------------------- maps
def polar(ax, ccrs, cfeature, keys=True, lonlabels=True, extent_lat=-80, land=True):
    """South-polar map 80-35 S in a circle, land, dotted graticule; bathymetric keys and meridian labels OUTSIDE the rim"""
    import matplotlib.path as mpath
    PC = ccrs.PlateCarree()
    ax.set_extent([-180, 180, extent_lat, -35], crs=PC)
    th = np.linspace(0, 2 * np.pi, 400); v = np.vstack([np.sin(th), np.cos(th)]).T
    ax.set_boundary(mpath.Path(v * 0.5 + 0.5), transform=ax.transAxes)
    if land:
        ax.add_feature(cfeature.LAND, facecolor=C["land"], edgecolor="none", zorder=3)
    gl = ax.gridlines(crs=PC, draw_labels=False, linewidth=0.25, color="#9A9A9A", alpha=0.6, linestyle=":")
    gl.xlocator = matplotlib.ticker.FixedLocator(range(-180, 181, 60)); gl.ylocator = matplotlib.ticker.FixedLocator([-70, -60, -50, -40])
    if "geo" in ax.spines: ax.spines["geo"].set_linewidth(0.5)
    rim = _rim(ax, ccrs)
    if lonlabels:
        for lo, s in ((0, "0°"), (90, "90°E"), (180, "180°"), (-90, "90°W")):
            rim(lo, s, color="#777777", tick=False, r=1.035)
    if keys:
        for lo, s in KEYS:
            rim(lo, s, color="#222222", tick=True, r=1.065)


def _rim(ax, ccrs):
    """returns f(lon, text) that writes text radially outside the circular boundary at that longitude"""
    PC = ccrs.PlateCarree()
    def f(lon, s, color="#222222", tick=True, r=1.06):
        x, y = ax.projection.transform_point(lon, -35.0, PC)
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); cx, cy = (x0 + x1) / 2, (y0 + y1) / 2; R = (x1 - x0) / 2
        ang = np.arctan2(y - cy, x - cx); ca, sa = np.cos(ang), np.sin(ang)
        if tick:
            ax.plot([cx + R * ca, cx + 1.025 * R * ca], [cy + R * sa, cy + 1.025 * R * sa], color=color, lw=0.6,
                    transform=ax.transData, clip_on=False, zorder=12)
        ha = "center" if abs(ca) < 0.3 else ("left" if ca > 0 else "right")
        va = "center" if abs(sa) < 0.3 else ("bottom" if sa > 0 else "top")
        ax.text(cx + r * R * ca, cy + r * R * sa, s, transform=ax.transData, ha=ha, va=va, fontsize=FS_ANN, color=color,
                clip_on=False, zorder=12)
    return f


def add_ice(ax, ccrs, lon, lat, ice):
    ov = np.where(ice, 1.0, np.nan)
    ax.contourf(lon, lat, ov, levels=[0.5, 1.5], colors=[C["ice"]], alpha=0.35, transform=ccrs.PlateCarree(), zorder=2)


def add_isobath(ax, lon, lat, depth_pos, transform=None, z=2):
    kw = dict(levels=[3000], colors=C["isobath"], linewidths=LW["isobath"], linestyles="solid", zorder=z)
    if transform is not None: kw["transform"] = transform
    ax.contour(lon, lat, depth_pos, **kw)


def add_fronts(ax, lon, lat, adt, lev_pf, lev_saf, transform=None, lw=None, z=6):
    kw = dict(levels=sorted([lev_pf, lev_saf]), colors=C["front"], linewidths=lw or LW["front"], linestyles="solid", zorder=z)
    if transform is not None: kw["transform"] = transform
    ax.contour(lon, lat, adt, **kw)


def sig_dots(ax, x, y, transform=None, s=5):
    kw = dict(s=s, c="white", edgecolors="k", linewidths=0.35, zorder=7)
    if transform is not None: kw["transform"] = transform
    ax.scatter(x, y, **kw)


def map_legend_handles(fdr=True, buffer=False, occ=False, fronts=True):
    h = []
    if fronts: h.append(Line2D([], [], color=C["front"], lw=LW["front"], label="mean SAF and PF"))
    if buffer: h.append(Line2D([], [], color=C["buffer"], lw=LW["buffer"], ls=(0, (2, 1.5)), label="ACC band + 100 km"))
    h.append(Line2D([], [], color=C["isobath"], lw=0.6, label="3000 m isobath"))
    h.append(Patch(facecolor=C["ice_flat"], edgecolor="none", label="sea ice"))
    if fdr: h.append(Line2D([], [], ls="none", marker="o", mfc="white", mec="k", mew=0.35, ms=2.6, label="10 % FDR"))
    if occ:
        h += [Line2D([], [], color=C["early"], lw=LW["early"], label=f"{LEVB} % contour 2001–2015"),
              Line2D([], [], color=C["late"], lw=LW["late"], ls=(0, (2.5, 1.2)), label=f"{LEVB} % contour 2016–2025")]
    return h


def lon_axis(ax, ticks=np.arange(-180, 181, 60), xlabel="Longitude (flow from left to right)"):
    ax.set_xlim(-180, 180); ax.set_xticks(ticks); ax.set_xticklabels([lonfmt(t) for t in ticks])
    if xlabel: ax.set_xlabel(xlabel)


def keys_on_axis(ax, y=0.99):
    for x_, s_ in KEYS:
        ax.text(x_, y, s_, ha="center", va="top", fontsize=FS_ANN, color="#222222", transform=ax.get_xaxis_transform())


def flow_arrow(ax, x0, x1, y=0.84):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"))
    ax.text((x0 + x1) / 2, y + 0.02, "flow", ha="center", va="bottom", fontsize=FS_ANN, color="#555555", transform=ax.get_xaxis_transform())


def edges(c):
    d = np.diff(c).mean(); return np.r_[c - d / 2, c[-1] + d / 2]


# ---------------------------------------------------------------- statistics tables
def retest_row(rev, key):
    """slope, ci (native units), pct per decade with CI, p from the revised test (revision/outputs/retest/retest_summary.csv)"""
    from osr11 import paths as _paths   # publication repository: always outputs/retest/ ('rev' kept for compatibility)
    D = pd.read_csv(_paths.OUT / "retest" / "retest_summary.csv").set_index("key")
    r = D.loc[key]; mean = r.slope * 1000 / r.pct_per_decade if r.pct_per_decade else np.nan
    f = 1000 / mean if np.isfinite(mean) and mean else np.nan
    return dict(slope=r.slope, lo=r.ci_lo, hi=r.ci_hi, pct=r.pct_per_decade, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p_new)


# ---------------------------------------------------------------- output
def save(fig, outdir, name):
    outdir = Path(outdir); (outdir / "png").mkdir(parents=True, exist_ok=True); (outdir / "pdf").mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "png" / f"{name}.png", dpi=600, bbox_inches=None)
    fig.savefig(outdir / "pdf" / f"{name}.pdf", bbox_inches=None)
    print("saved", outdir / "png" / f"{name}.png", f"{fig._mm[0]:.0f} x {fig._mm[1]:.0f} mm")


def bathy_path(rev=None):
    """outputs/base/bathymetry_so.nc (publication repository; distributed)."""
    from osr11 import paths as _paths
    return _paths.BATHY
