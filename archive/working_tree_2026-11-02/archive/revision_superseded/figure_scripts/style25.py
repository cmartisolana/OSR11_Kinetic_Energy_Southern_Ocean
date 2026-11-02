"""Common style for the final (v2.5) main-text and Supporting Information figures.

Figures are drawn at print size for State of the Planet / Copernicus (two-column width 17.4 cm ~ 7.2 in), so that the
type is 7-9 pt in the printed page and no rescaling is needed. Panel titles are "(a) Title", left-aligned, bold.
Colours: KE red, EKE green, transport blue, wind grey, position slate blue, tracked eddies teal (as before).
Significance: one level, p < 0.05, filled symbols; open symbols p >= 0.05; 95 % confidence intervals everywhere.
"""
from figstyle import *  # noqa: F401,F403  (OUT, STYLE, helpers, retest_row, ...)
import figstyle as _fs
import matplotlib.pyplot as plt
import numpy as np

W_FULL = 7.2      # in, two-column width
W_HALF = 3.5      # in, one-column width

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.2, "grid.linewidth": 0.4,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "savefig.dpi": 400, "figure.dpi": 150, "pdf.fonttype": 42, "ps.fonttype": 42,
})
STYLE.update({"smooth_lw": 1.4, "raw_lw": 0.5, "trend_lw": 1.0})
C = dict(ke=STYLE["ke"], eke=STYLE["eke"], te="#2A9D8F", pos=STYLE["lat"], tr="#1D6FA3", obs="#D9822B", wind="#7A8288",
         ink="#1F2933", grey="#6b7280", light="#E3E6EA")


def ptitle(ax, letter, text=None, **kw):
    """(a) Title, left-aligned, bold."""
    s = f"({letter})" + (f" {text}" if text else "")
    ax.set_title(s, loc="left", fontsize=kw.pop("fontsize", 8.5), fontweight="bold", pad=kw.pop("pad", 4), **kw)


def forest(ax, rows, **kw):
    """forest_panel at print size."""
    kw.setdefault("label_size", 7.5); kw.setdefault("value_size", 6.5)
    _fs.forest_panel(ax, rows, **kw)
    for ln in ax.lines:
        if ln.get_marker() not in (None, "None", ""):
            ln.set_markersize(4.5); ln.set_markeredgewidth(1.0)
        else:
            ln.set_linewidth(min(ln.get_linewidth(), 1.3))
    ax.xaxis.label.set_size(7.5)
    return ax


def legend_below(ax, ncol=3, y=-0.2, **kw):
    kw.setdefault("fontsize", 7); kw.setdefault("frameon", False)
    return ax.legend(loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol, handletextpad=0.3, columnspacing=1.0,
                     borderaxespad=0, **kw)


def save(fig, folder, name):
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{name}.png", dpi=400); fig.savefig(folder / f"{name}.pdf")
    print("saved", folder / name)


def polar_style(ax, ccrs, cfeature, extent_lat=-80, label_lons=True):
    """Polar map, circular boundary, light graticule labelled at 0/90/180/270 and 40/60 S."""
    import matplotlib.path as mpath
    ax.set_extent([-180, 180, extent_lat, -35], crs=ccrs.PlateCarree())
    th = np.linspace(0, 2 * np.pi, 200); v = np.vstack([np.sin(th), np.cos(th)]).T
    ax.set_boundary(mpath.Path(v * 0.5 + 0.5), transform=ax.transAxes)
    ax.add_feature(cfeature.LAND, facecolor="#2B2B2B", edgecolor="none", zorder=3)
    ax.coastlines("50m", linewidth=0.3, color="#2B2B2B", zorder=3)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.3, color="#7A8288", alpha=0.5, linestyle=":")
    gl.xlocator = plt.FixedLocator(range(-180, 181, 60)); gl.ylocator = plt.FixedLocator([-70, -60, -50, -40])
    if label_lons:
        for lo, s in ((0, "0°"), (180, "180°")):
            ax.text(lo, -36.3, s, transform=ccrs.PlateCarree(), ha="center", va="center", fontsize=6, color="#4A4A4A",
                    zorder=6, bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
    ax.set_frame_on(False)
