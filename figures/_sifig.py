"""Shared set-up of the Supplement figure scripts S1, S2, S4-S7 (final layout, 10 Nov 2026).

Imported first by every figS*.py in this folder. It
  - puts the repository src/ on sys.path and imports osr11.style (S) and osr11.si_common (Q);
  - raises every font of the Supplement figures to >= 7 pt at print size (ticks, legends, annotations 7 pt; labels and
    titles 7 pt; panel letters 8 pt bold), as required by the Supplement figure standard (the main figures keep 6 pt ticks);
  - points cartopy to the local Natural Earth shapefiles (natural_earth/physical/, no download);
  - resolves <outputs dir> and <figure dir> from the command line ([outputs dir] [figure dir]) or from $OSR11_OUTPUTS /
    $OSR11_FIGURES (osr11.paths defaults: outputs/, rendered/).
Nothing is written outside <figure dir>.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
if len(sys.argv) > 1:                                   # outputs dir from the command line: also seen by osr11.paths
    os.environ["OSR11_OUTPUTS"] = str(Path(sys.argv[1]).resolve())
if len(sys.argv) > 2:
    os.environ["OSR11_FIGURES"] = str(Path(sys.argv[2]).resolve())

from osr11 import style as S, si_common as Q, paths as P   # noqa: E402,F401  (Q is used by the figure scripts as F.Q)

FS = 7.0
S.FS_TICK = S.FS_ANN = S.FS_LAB = S.FS_TTL = FS
S.plt.rcParams.update({"font.size": FS, "axes.labelsize": FS, "axes.titlesize": FS, "xtick.labelsize": FS,
                       "ytick.labelsize": FS, "legend.fontsize": FS, "legend.title_fontsize": FS})
OUT, FIGD = P.OUT, P.FIG
C = S.C
plt = S.plt
GREY_TXT = "#5A5A5A"     # secondary text (periods, notes): darker than C['grey'] so that 7 pt stays legible in print


def cartopy_local():
    """cartopy reading the Natural Earth shapefiles shipped in natural_earth/physical (osr11.style.local_natural_earth)"""
    return S.local_natural_earth()


def note(fig, x, y, text, color="#222222", ha="left", va="top", **kw):
    W, H = fig._mm
    return fig.text(x / W, y / H, S.um(text), fontsize=FS, color=color, ha=ha, va=va, linespacing=1.15, **kw)


def check(fig):
    """report every text and legend whose extent leaves the canvas (mm) - a guard against clipped labels"""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); W, H = fig._mm; dpi = fig.dpi; bad = []
    arts = [t for t in fig.findobj(match=lambda a: hasattr(a, "get_window_extent") and a.get_visible()
                                   and a.__class__.__name__ in ("Text", "Legend", "Annotation"))]
    for a in arts:
        if a.__class__.__name__ in ("Text", "Annotation") and not a.get_text().strip():
            continue
        try:
            bb = a.get_window_extent(r)
        except Exception:
            continue
        x0, x1, y0, y1 = bb.x0 / dpi * 25.4, bb.x1 / dpi * 25.4, bb.y0 / dpi * 25.4, bb.y1 / dpi * 25.4
        if x0 < -0.05 or y0 < -0.05 or x1 > W + 0.05 or y1 > H + 0.05:
            bad.append((a.__class__.__name__, getattr(a, "get_text", lambda: "legend")()[:40], round(x0, 1), round(x1, 1), round(y0, 1), round(y1, 1)))
    print("check: outside canvas:", bad if bad else "none")
    return bad
