"""Shared conventions of the Supporting Information figures S1-S6 (SI v3, 3 Nov 2026).
Publication repository: copied from figures_final/scripts_SI_v3/si_common.py; imports osr11.style and resolves the
output and figure directories through osr11.paths (command-line arguments optional).
Imports the main-figure style module (osr11.style) and adds
  PROD        the one convention that distinguishes the altimetry products in every SI panel:
              Product 1 (all-satellite, DUACS DT2024 1/8 deg) = circle marker, solid line;
              Product 4 (two-satellite climate product, 1/4 deg) = square marker, dashed line.
              Colours stay with the quantity (KE orange, EKE green, half-power position purple, velocity blue, in situ orange).
  forest()    dot plot of trends with 95 % intervals: one colour per quantity, one marker/linestyle per product,
              filled marker p < 0.05, open p >= 0.05.
  retest_row  trend, 95 % interval and p of one key of outputs/retest/retest_summary.csv (calibrated test, Sect. S1)
Scripts take two optional arguments: <outputs dir> (default $OSR11_OUTPUTS = outputs/) and <figure dir> (default
$OSR11_FIGURES = rendered/); derived CSVs written by the scripts go to <outputs dir> (see each script).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from osr11 import style as S, paths as PA
from matplotlib.lines import Line2D

KM = 111.195                     # km per degree of latitude
PROD = {"P1": dict(marker="o", ls="-", name="Product 1 (all-satellite)", key="allsat", short="Product 1"),
        "P4": dict(marker="s", ls=(0, (2.2, 1.3)), name="Product 4 (two-satellite)", key="twosat", short="Product 4")}
SIG_NS = "open: p ≥ 0.05"


def args():
    """(outputs dir, figure dir, outputs dir): from the command line when given, else osr11.paths (OUT, FIG)."""
    OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else PA.OUT
    FIGD = Path(sys.argv[2]) if len(sys.argv) > 2 else PA.FIG
    return OUT, FIGD, OUT


def retest_row(OUT, key, scale=None):
    """trend (pct per decade unless scale given: slope*scale per decade), 95 % interval, p of one retest_summary key"""
    D = pd.read_csv(Path(OUT) / "retest" / "retest_summary.csv").set_index("key"); r = D.loc[key]
    if scale is None:
        mean = r.slope * 1000 / r.pct_per_decade if r.pct_per_decade else np.nan
        f = 1000 / mean if np.isfinite(mean) and mean else np.nan
        return dict(val=r.pct_per_decade, lo=r.ci_lo * f, hi=r.ci_hi * f, p=r.p_new)
    return dict(val=r.slope * scale, lo=r.ci_lo * scale, hi=r.ci_hi * scale, p=r.p_new)


def forest(ax, rows, xlim=None, xlabel=None, dodge=0.26, gap=0.6, ms=3.4, sep_after=()):
    """rows: [(label, [(colour, prod key or None, dict(val, lo, hi, p)), ...])], top first; entries of a row are dodged
    vertically in order. sep_after: row indices after which a light separator is drawn."""
    ys = -np.arange(len(rows), dtype=float)
    for (lab, ent), y in zip(rows, ys):
        n = len(ent)
        for i, (col, pk, r) in enumerate(ent):
            yy = y - (i - (n - 1) / 2) * (dodge if n > 1 else 0)
            st = PROD[pk] if pk else dict(marker="o", ls="-")
            ax.plot([r["lo"], r["hi"]], [yy, yy], color=col, lw=0.8, ls=st["ls"], solid_capstyle="butt", zorder=2)
            sg = S.sigp(r["p"])
            ax.plot(r["val"], yy, st["marker"], color=col, mfc=col if sg else "white", mec=col, mew=0.8, ms=ms, zorder=3)
    ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0)
    for i in sep_after:
        ax.axhline(ys[i] - 0.5, color=S.C["light"], lw=0.5, zorder=0)
    ax.set_yticks(ys); ax.set_yticklabels([l for l, _ in rows]); ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False); ax.set_ylim(ys[-1] - gap, ys[0] + gap)
    if xlim: ax.set_xlim(*xlim)
    if xlabel: ax.set_xlabel(xlabel)
    S.xticks_minus(ax)
    return ys


def prod_handles(color="#222222", ms=3.4, ns=True, keys=("P1", "P4")):
    h = [Line2D([], [], ls=PROD[k]["ls"], lw=0.8, marker=PROD[k]["marker"], color=color, mfc=color, mec=color, ms=ms,
                label=PROD[k]["name"]) for k in keys]
    if ns:
        h.append(Line2D([], [], ls="none", marker="o", color="#555555", mfc="white", mew=0.8, ms=ms, label=SIG_NS))
    return h


def note(fig, x, y, text, color="#222222", ha="left", va="top"):
    """6 pt annotation at (x, y) mm"""
    W, H = fig._mm
    fig.text(x / W, y / H, S.um(text), fontsize=S.FS_ANN, color=color, ha=ha, va=va, linespacing=1.15)
