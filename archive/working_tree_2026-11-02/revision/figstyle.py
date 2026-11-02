"""Shared style for the main-text figures (v3), reproducing the v1/v2 figure style
(notebooks 01, 03a, 03b, 07: STYLE dict, rcParams, polar map styling, colours)."""
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.path as mpath

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

STYLE = {
    "raw_alpha": 0.18, "raw_lw": 0.7, "smooth_lw": 1.9, "trend_lw": 1.4,
    "grid_alpha": 0.18, "grid_lw": 0.5, "axis_label_size": 11, "title_size": 11,
    "tick_size": 9, "annot_size": 9, "panel_size": 11,
    "neutral": "#7A8288", "ke": "#e96c6a", "eke": "#2A9D8F", "trend_sig": "#1F2933",
    "lat": "#5285A0", "width": "#9D2A90", "ns": "#D0D4D8",
    "transport": "#1f77b4", "transport_obs": "#264653", "wind": "#d36b20", "box": "#2a9d8f",
}
REGION_ORDER = ["Indian", "Pacific", "Atlantic", "ACC"]
REGION_COLORS = {"Indian": "#2A9D8F", "Pacific": "#457B9D", "Atlantic": "#E07A5F", "ACC": "#2F3E46"}
C_EPOCH = ["#457B9D", "#7A8288", "#E07A5F"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10, "axes.labelsize": STYLE["axis_label_size"], "axes.titlesize": STYLE["title_size"],
    "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
    "axes.facecolor": "white", "axes.grid": True, "axes.linewidth": 0.6, "figure.facecolor": "white",
    "grid.alpha": 0.18, "grid.linewidth": 0.5, "xtick.labelsize": STYLE["tick_size"],
    "ytick.labelsize": STYLE["tick_size"], "legend.fontsize": 8, "legend.framealpha": 0.9,
    "legend.edgecolor": "#cccccc", "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "savefig.facecolor": "white",
})


def panel_label(ax, label, x=0.015, y=0.985):
    ax.text(x, y, f"({label})", transform=ax.transAxes, ha="left", va="top",
            fontsize=STYLE["panel_size"], fontweight="semibold", color="#1F2933", zorder=20)


def base_timeseries_style(ax):
    ax.grid(alpha=STYLE["grid_alpha"], lw=STYLE["grid_lw"])
    ax.tick_params(length=3.0, width=0.7, color="#67727A")


def stars(p):
    if not np.isfinite(p):
        return ""
    return "*" if p < 0.05 else ""      # single level (v2.2): p < 0.05


def ptxt(p):
    return "p<0.01" if p < 0.01 else f"p={p:.2f}"


def running(y, w=12):
    return pd.Series(y).rolling(w, center=True, min_periods=w // 2).mean().values


def deseason(df, col):
    return (df[col] - df.groupby("month")[col].transform("mean") + df[col].mean()).values


def ts_panel(ax, x, y, color, label=None, slope=None, intercept=None):
    """v1 time-series style: faint raw monthly, 12-month running mean, dashed trend."""
    ax.plot(x, y, color=color, alpha=STYLE["raw_alpha"], marker=".", ms=2, lw=STYLE["raw_lw"])
    ax.plot(x, running(y), color=color, lw=STYLE["smooth_lw"], label=label)
    if slope is not None:
        xf = np.array([np.nanmin(x), np.nanmax(x)])
        ax.plot(xf, intercept + slope * xf, color=color, ls="--", lw=STYLE["trend_lw"])


def theil_line(x, y, slope):
    """intercept so that a Theil-Sen line with the given slope passes through the median point."""
    ok = np.isfinite(y)
    return np.nanmedian(y[ok] - slope * x[ok])


def clim_allsat(data):
    f = OUT / "clim_allsat.nc"
    if f.exists():
        return xr.open_dataset(f)
    raise FileNotFoundError(f)


# ---------------- polar maps (notebook 03a) ----------------
def polar_ax_style(ax, ccrs, cfeature, extent_lat=-82):
    ax.set_extent([-180, 180, extent_lat, -35], crs=ccrs.PlateCarree())
    theta = np.linspace(0, 2 * np.pi, 200)
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    ax.set_boundary(mpath.Path(verts * 0.5 + 0.5), transform=ax.transAxes)
    ax.add_feature(cfeature.LAND, facecolor="#000000", edgecolor="#555555", linewidth=0.2, zorder=2)
    ax.coastlines("50m", linewidth=0.45, color="#333333")
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.30, color="#7A8288",
                      alpha=0.45, linestyle=":")
    gl.xlocator = plt.FixedLocator(range(-180, 181, 60)); gl.ylocator = plt.FixedLocator([-60, -40])
    gl.top_labels = False; gl.right_labels = False
    gl.xlabel_style = {"size": 7, "color": "#4A4A4A"}; gl.ylabel_style = {"size": 7, "color": "#4A4A4A"}
    ax.set_frame_on(False)


def add_acc_contours(ax, ccrs, lon, lat, adt):
    return ax.contour(lon, lat, adt, levels=[-0.6, 0.2], colors=["#5f5f5f", "#2f2f2f"], linewidths=[0.6, 0.9],
                      alpha=0.65, transform=ccrs.PlateCarree(), zorder=2)


def add_ice_mask(ax, ccrs, lon, lat, ice):
    ov = np.where(ice, 1.0, np.nan)
    ax.contourf(lon, lat, ov, levels=[0.5, 1.5], colors=["#8EC5FF"], alpha=0.35, transform=ccrs.PlateCarree(), zorder=2)
    ax.contour(lon, lat, np.where(np.isfinite(ov), 1.0, np.nan), levels=[0.5], colors=["#595959"],
               linewidths=[0.35], alpha=0.5, transform=ccrs.PlateCarree(), zorder=4)


def add_sector_boundaries(ax, ccrs):
    lats = np.linspace(-82, -35, 200)
    for lon in (20.0, 147.0, -70.0):
        ax.plot(np.full_like(lats, lon), lats, color="#222222", lw=0.85, ls=(0, (5, 3)),
                transform=ccrs.PlateCarree(), zorder=5, alpha=0.75)


# ---------------- revised statistics (retest_all.py, 2026-10-13) ----------------
def retest_override(df, keyfun, cols=None):
    """Replace p-values and CIs in a results table by those of the bias-corrected, calibrated test
    (outputs/retest/retest_summary.csv). keyfun(index, row) -> retest key or None. Returns a copy."""
    cols = cols or {"p_cal": "p_new", "ci_cal_lo": "ci_lo", "ci_cal_hi": "ci_hi"}
    D = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")
    df = df.copy(); n = 0
    for idx, row in df.iterrows():
        k = keyfun(idx, row)
        if k is not None and k in D.index:
            for c_old, c_new in cols.items():
                df.loc[idx, c_old] = D.loc[k, c_new]
            n += 1
    print(f"retest_override: {n} of {len(df)} rows updated")
    return df


# ---------------- summary ("forest") panels for trend comparisons (v2.1) ----------------
def retest_row(key):
    """slope, ci_lo, ci_hi (native units), pct per decade with CI, p — from the revised test."""
    D = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")
    r = D.loc[key]; mean = r.slope * 1000 / r.pct_per_decade if r.pct_per_decade else np.nan
    f = 1000 / mean if np.isfinite(mean) and mean else np.nan
    return dict(slope=r.slope, lo=r.ci_lo, hi=r.ci_hi, pct=r.pct_per_decade, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p_new)


def forest_panel(ax, rows, unit=r"% dec$^{-1}$", xlabel=None, group_gap=0.6, series_colors=None, series_markers=None,
                 label_size=9, value_size=8, xlim=None, dodge=True):
    """rows: list of (group_label, {series_name: dict(pct, pct_lo, pct_hi, p)}) in display order (top first).
    Draws dot + 95 % range per series, filled when p < 0.05, open otherwise; value text at the right."""
    series = list(dict.fromkeys(s for _, d in rows for s in d))
    ns = len(series); off = np.linspace(0.22, -0.22, ns) if (ns > 1 and dodge) else [0.0] * ns
    series_colors = series_colors or {s: list(REGION_COLORS.values())[i] for i, s in enumerate(series)}
    series_markers = series_markers or {s: m for s, m in zip(series, ["o", "s", "D", "^"])}
    y0 = 0.0; ys = []; labels = []
    for gl, d in rows:
        ys.append(y0); labels.append(gl)
        for s, o in zip(series, off):
            if s not in d:
                continue
            r = d[s]; sig = r["p"] < 0.05; col = series_colors[s]
            ax.plot([r["pct_lo"], r["pct_hi"]], [y0 + o, y0 + o], color=col, lw=1.6, solid_capstyle="round", zorder=3)
            ax.plot(r["pct"], y0 + o, series_markers[s], ms=6.5, color=col, mfc=col if sig else "white", mew=1.4, zorder=4)
        y0 -= 1.0 + group_gap * 0
    ax.axvline(0, color="#495057", lw=0.8, zorder=2)
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=label_size)
    ax.set_ylim(min(ys) - 0.7, max(ys) + 0.7)
    ax.set_xlabel(xlabel or f"Trend ({unit})", fontsize=9); ax.tick_params(axis="x", labelsize=8)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", alpha=0); ax.grid(axis="x", alpha=STYLE["grid_alpha"], lw=STYLE["grid_lw"])
    if xlim:
        ax.set_xlim(*xlim)
    xr = ax.get_xlim()
    for (gl, d), y in zip(rows, ys):
        txt = "  ".join(f"{d[s]['pct']:+.1f}{stars(d[s]['p']) or ' (n.s.)'}" for s in series if s in d)
        ax.text(xr[1] + 0.02 * (xr[1] - xr[0]), y, txt, va="center", ha="left", fontsize=value_size, color="#1F2933", clip_on=False)
    for s in series:
        ax.plot([], [], series_markers[s], color=series_colors[s], mfc=series_colors[s], mew=1.4, ms=6.5, label=s)
    ax.plot([], [], "o", color="#7A8288", mfc="white", mew=1.4, ms=6.5, label="not significant (p ≥ 0.05)")
