"""SI Figure S3: cross-track velocity trend along the LM Gould transect in passage coordinates (25 km bins of distance from the
northern end), as in Gutierrez-Villanueva et al. (2023, passage coordinates): in situ ADCP velocity at 70-118 m, Oct 2005-Apr 2019
(orange, 95 % bars), altimetric surface geostrophic velocity (Product 1) over the same period (green dashed) and over 2001-2025
(green solid, 95 % shading). Grey fill: 2001-2025 mean altimetric cross-track velocity (right axis); light blue: seasonally
ice-covered part of the section (408-850 km). Labels: mean positions of the SAF, PF and SACCF (fixed distances 40, 400, 790 km, as
in the earlier version of this figure). Filled symbols: p < 0.05 (calibrated Theil-Sen, dk5_trends.csv, coord = passage).
Replaces revision/figures_v2.1/figS3_velocity_transect.py (which wrote figS/new_S12.png at 300 dpi).

Publication repository: ported from figures_final/scripts/figS3_velocity_transect.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/drake/dk5/dk5_trends.csv
Outputs : rendered/{png,pdf}/FigureS3
Run     : python figures/figS3_velocity_transect.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd
from osr11 import style as S, paths as P  # noqa: E402
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else P.FIG
plt = S.plt; C = S.C
K5 = pd.read_csv(P.OUT / "drake" / "dk5" / "dk5_trends.csv")
def prof(src, per): return K5[(K5.coord == "passage") & (K5.src == src) & (K5["var"] == "u") & (K5.period == per)].sort_values("bin")

H = 70.0
fig = S.figure(S.W1, H)
ax = S.axmm(fig, 13, 21, 58, 43); ax2 = ax.twinx()
mu = prof("allsat", "2001-2025")
ax2.axvspan(408, 850, color=C["ice_flat"], lw=0, zorder=-1)
ax2.fill_between(mu.bin, 0, mu["mean"], color="#BDBDBD", alpha=0.35, lw=0, zorder=0); ax2.set_ylim(0, 45)
ax2.set_ylabel("Mean velocity (cm s$^{-1}$)", color="#8A8A8A"); ax2.spines["right"].set_visible(True); ax2.tick_params(colors="#8A8A8A")
ax2.spines["right"].set_color("#8A8A8A")
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
qi = prof("insitu", "2005-2019"); qa = prof("allsat", "2005-2019"); ql = prof("allsat", "2001-2025")
ax.fill_between(ql.bin, 10 * ql.ci_lo, 10 * ql.ci_hi, color=C["eke"], alpha=0.18, lw=0, zorder=1)
S.profile(ax, ql.bin, 10 * ql.slope, ql.p, C["eke"], ms=2.2, z=3)
S.profile(ax, qa.bin, 10 * qa.slope, qa.p, C["eke"], ms=2.2, z=3, ls=(0, (3, 1.6)))
ax.vlines(qi.bin, 10 * qi.ci_lo, 10 * qi.ci_hi, color=C["obs"], lw=0.5, alpha=0.8, zorder=2)
S.profile(ax, qi.bin, 10 * qi.slope, qi.p, C["obs"], mk="D", ms=2.4, z=4)
S.zero(ax); ax.set_xlim(0, 850); ax.set_xlabel("Distance from the northern end (km)")
ax.set_ylabel("Velocity trend (cm s$^{-1}$ per decade)"); S.yticks_minus(ax)
for x_, lab in ((40, "SAF"), (400, "PF"), (790, "SACCF")):
    ax.text(x_, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=S.FS_ANN, color="#555555", fontweight="bold")
fig.legend(handles=[Line2D([], [], color=C["obs"], marker="D", ms=2.4, lw=S.LW["profile"], label="in situ ADCP, 70–118 m, Oct 2005–Apr 2019"),
                    Line2D([], [], color=C["eke"], marker="o", ms=2.2, lw=S.LW["profile"], ls=(0, (3, 1.6)), label="altimetry, same period"),
                    Line2D([], [], color=C["eke"], marker="o", ms=2.2, lw=S.LW["profile"], label="altimetry, 2001–2025"),
                    Patch(facecolor=C["ice_flat"], label="seasonal sea ice")],
           loc="lower left", bbox_to_anchor=(2 / S.W1, 0.0), ncol=2, fontsize=S.FS_ANN, columnspacing=0.8, handlelength=2.0)
S.save(fig, FIGD, "FigureS3"); plt.close(fig)
print("S3 in situ bins", len(qi), "sig", int((qi.p.round(3) < 0.05).sum()), "; altimetry 2001-2025 sig", int((ql.p.round(3) < 0.05).sum()))
