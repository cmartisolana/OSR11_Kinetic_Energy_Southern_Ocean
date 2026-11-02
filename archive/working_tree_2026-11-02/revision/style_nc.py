"""Common figure style for the manuscript and Supporting Information (v2.6), after Nature-family figure guidelines:
180 mm (double column) or 88 mm (single column) wide, sans-serif type at 5-7 pt, lowercase bold panel letters,
thin axes with outward ticks and no grid, frameless legends, no statistics printed inside panels (they are in the
captions and text), one colour per quantity across all figures (colour-blind-safe, Okabe-Ito based):
  KE vermilion, EKE bluish green, tracked eddies bluish green dashed / diamonds, transport blue, frontal position
  reddish purple, wind grey, in situ orange, two-satellite product grey.
Significance: filled symbols p < 0.05, open symbols p >= 0.05; bars are 95 % confidence intervals.
"""
from figstyle import (OUT, HERE, running, theil_line, add_ice_mask, add_acc_contours, add_sector_boundaries,  # noqa: F401
                      retest_override, retest_row)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

MM = 1 / 25.4
W2 = 180 * MM      # double column
W1 = 88 * MM       # single column

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.5, "axes.labelsize": 6.5, "axes.titlesize": 6.5, "axes.titleweight": "normal",
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "legend.frameon": False,
    "legend.handlelength": 1.6, "legend.handletextpad": 0.4, "legend.columnspacing": 1.0, "legend.borderaxespad": 0.2,
    "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "axes.edgecolor": "#222222", "axes.labelcolor": "#222222", "text.color": "#222222",
    "xtick.color": "#222222", "ytick.color": "#222222", "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.size": 1.5, "ytick.minor.size": 1.5, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "lines.linewidth": 1.0, "lines.markersize": 3.5, "patch.linewidth": 0.5,
    "figure.dpi": 150, "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.default": "regular", "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans",
    "mathtext.it": "Liberation Sans:italic",
})

C = dict(ke="#D55E00", eke="#009E73", te="#009E73", tr="#0072B2", pos="#A2508A", wind="#7F7F7F", obs="#E69F00",
         sin="#56B4E9", two="#9A9A9A", ink="#222222", grey="#7F7F7F", light="#E6E6E6", ice="#DCEBF7")


def plabel(ax, letter, title=None, dx=-0.0, dy=4, x=0.0):
    """Lowercase bold panel letter at the upper-left corner of the axes, optional short plain title after it."""
    ax.annotate(letter, xy=(x, 1), xycoords="axes fraction", xytext=(dx, dy), textcoords="offset points",
                fontsize=8, fontweight="bold", ha="left", va="bottom")
    if title:
        ax.annotate(title, xy=(x, 1), xycoords="axes fraction", xytext=(dx + 9, dy + 0.5), textcoords="offset points",
                    fontsize=6.5, ha="left", va="bottom")


def zero(ax, axis="y"):
    (ax.axhline if axis == "y" else ax.axvline)(0, color="#9A9A9A", lw=0.5, zorder=0)


def ts_style(ax, xlim=(2000.8, 2025.6)):
    ax.set_xlim(*xlim); ax.set_xticks([2005, 2010, 2015, 2020, 2025])


def forest(ax, rows, colors, markers=None, xlim=None, dodge=True, xlabel=None, gap=0.6, ms=3.6):
    """rows: [(label, {series: dict(pct, pct_lo, pct_hi, p)})], top first. Markers filled if p < 0.05, open otherwise.
    Returns the y positions of the row centres."""
    markers = markers or {}
    series = list(dict.fromkeys(s for _, d in rows for s in d))
    ns = len(series); off = {s: (i - (ns - 1) / 2) * (0.26 if dodge and ns > 1 else 0) for i, s in enumerate(series)}
    ys = -np.arange(len(rows), dtype=float)
    for (lab, d), y in zip(rows, ys):
        for s, r in d.items():
            yy = y - off[s]; c = colors[s]; m = markers.get(s, "o")
            ax.plot([r["pct_lo"], r["pct_hi"]], [yy, yy], color=c, lw=0.8, solid_capstyle="butt", zorder=2)
            sig = np.isfinite(r["p"]) and r["p"] < 0.05
            ax.plot(r["pct"], yy, m, color=c, mfc=c if sig else "white", mec=c, mew=0.8, ms=ms, zorder=3)
    ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0)
    ax.set_yticks(ys); ax.set_yticklabels([l for l, _ in rows])
    ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_ylim(ys[-1] - gap, ys[0] + gap)
    if xlim: ax.set_xlim(*xlim)
    if xlabel: ax.set_xlabel(xlabel)
    return ys


def sig_handles(series_colors, markers=None, ns_label="open: p ≥ 0.05"):
    markers = markers or {}
    h = [Line2D([], [], ls="-", lw=0.8, marker=markers.get(s, "o"), color=c, mfc=c, mec=c, ms=3.6, label=s)
         for s, c in series_colors.items()]
    h.append(Line2D([], [], ls="none", marker="o", color="#555555", mfc="white", mew=0.8, ms=3.6, label=ns_label))
    return h


def polar(ax, ccrs, cfeature, extent_lat=-80, lons=True):
    import matplotlib.path as mpath
    ax.set_extent([-180, 180, extent_lat, -35], crs=ccrs.PlateCarree())
    th = np.linspace(0, 2 * np.pi, 200); v = np.vstack([np.sin(th), np.cos(th)]).T
    ax.set_boundary(mpath.Path(v * 0.5 + 0.5), transform=ax.transAxes)
    ax.add_feature(cfeature.LAND, facecolor="#3A3A3A", edgecolor="none", zorder=3)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.25, color="#9A9A9A", alpha=0.6, linestyle=":")
    gl.xlocator = plt.FixedLocator(range(-180, 181, 60)); gl.ylocator = plt.FixedLocator([-70, -60, -50, -40])
    if lons:
        for lo, s in ((0, "0°"), (180, "180°")):
            ax.text(lo, -36.5, s, transform=ccrs.PlateCarree(), ha="center", va="center", fontsize=5.5, color="#555555",
                    zorder=6, bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.8))
    ax.spines["geo"].set_linewidth(0.5) if "geo" in ax.spines else None


def hcbar(fig, mappable, ax, label, **kw):
    kw.setdefault("orientation", "horizontal"); kw.setdefault("shrink", 0.75); kw.setdefault("pad", 0.03); kw.setdefault("aspect", 25)
    cb = fig.colorbar(mappable, ax=ax, **kw); cb.set_label(label, fontsize=6); cb.ax.tick_params(labelsize=5.5, width=0.4, length=2)
    cb.outline.set_linewidth(0.4); return cb


def save(fig, folder, name):
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{name}.png", dpi=600, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(folder / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    print("saved", folder / name)
