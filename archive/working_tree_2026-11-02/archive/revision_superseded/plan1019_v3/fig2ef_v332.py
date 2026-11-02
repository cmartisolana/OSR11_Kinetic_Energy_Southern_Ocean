"""Fig. 2e-f (v3.32): eddy-saturation diagnostics in Drake Passage, drawn on a transparent page of the same size as
the v2.0 Drake figure, at the positions of the old panels (e, f), to be overlaid by fig2_compose_v332.py.
(e) 2001-2025 trends (% per decade, 95 % intervals; monthly calibrated test, dk1_trends.csv) of eastward wind stress,
    along-stream speed, KE and EKE (Product 1) upstream, in the passage box and downstream.
(f) d ln X / d ln tau as the ratio of the 2001-2025 trends of annual ln X and ln tau (dk23_trend_ratios.csv), Products 1
    and 4; reference lines at 0 and 1 (EKE proportional to tau, Marshall et al. 2017).
Usage: python fig2ef_v332.py <dk1 dir> <dk23 dir> <out pdf>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
DK1, DK23, OUTF = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
plt.rcParams.update({   # style_nc.py
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.5, "axes.labelsize": 6.5, "axes.titlesize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "legend.frameon": False, "legend.handlelength": 1.6, "legend.handletextpad": 0.4,
    "legend.columnspacing": 1.0, "legend.borderaxespad": 0.2, "axes.linewidth": 0.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.edgecolor": "#222222", "axes.labelcolor": "#222222", "text.color": "#222222",
    "xtick.color": "#222222", "ytick.color": "#222222", "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "lines.linewidth": 1.0, "patch.linewidth": 0.5, "pdf.fonttype": 42, "mathtext.default": "regular",
    "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic"})
C = dict(ke="#D55E00", eke="#009E73", tr="#0072B2", wind="#7F7F7F", two="#9A9A9A")
sigp = lambda p: np.isfinite(p) and round(float(p), 3) < 0.05
PW, PH = 471.0592603718, 476.1464860667                     # v2.0 page (pt)
AXE = (37.676328, 23.170957, 224.201362, 156.512689)       # old panel e axes (pt)
AXF = (246.584366, 23.170957, 433.109399, 156.512689)      # old panel f axes (pt)
REG = [("up", "Upstream\n100–80°W"), ("box", "Drake Passage\n80–52°W"), ("down", "Downstream\n55–35°W")]

fig = plt.figure(figsize=(PW / 72, PH / 72)); fig.patch.set_alpha(0)
def axes(b):
    return fig.add_axes([b[0] / PW, b[1] / PH, (b[2] - b[0]) / PW, (b[3] - b[1]) / PH])
def label(x, letter, title):
    fig.text(x / PW, 168.2 / PH, letter, fontsize=8, fontweight="bold", va="baseline")
    fig.text((x + 9) / PW, 168.7 / PH, title, fontsize=6.5, va="baseline")

# e: trends by region
TR = pd.read_csv(DK1 / "dk1_trends.csv"); TR = TR[TR.period == "2001-2025"].set_index("name")
SER = [("wind_{r}_TAUX", "eastward wind stress", C["wind"]), ("allsat_{r}_SPD", "along-stream speed", C["tr"]),
       ("allsat_{r}_KE", "KE", C["ke"]), ("allsat_{r}_EKE", "EKE", C["eke"])]
ax = axes(AXE); ax.set_facecolor("white"); w = 0.18
for gi, (r, _) in enumerate(REG):
    for k, (key, _, col) in enumerate(SER):
        q = TR.loc[key.format(r=r)]; x = gi + (k - 1.5) * w
        ax.bar(x, q.pct, w * 0.88, color=col if sigp(q.p) else "white", edgecolor=col, lw=0.8, zorder=3)
        ax.plot([x, x], [q.pct_lo, q.pct_hi], color="#333333", lw=0.6, zorder=4)
ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=0)
ax.set_xlim(-0.5, 2.5); ax.set_ylim(-8, 32)
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in REG]); ax.tick_params(axis="x", length=0)
ax.set_ylabel("Trend 2001–2025 (% per decade)")
h = [Patch(facecolor=c, edgecolor=c, label=l) for _, l, c in SER] + [Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")]
ax.legend(handles=h, loc="upper left", ncol=3, fontsize=5.4, handlelength=1.1, columnspacing=0.8)
label(13.676328, "e", "Wind stress, mean flow and eddy energy, 2001–2025")

# f: sensitivity to the wind stress (trend ratios); two-satellite product shown for EKE, where the products differ
RT = pd.read_csv(DK23 / "dk23_trend_ratios.csv")
VAR = [("SPD", "along-stream speed", C["tr"]), ("KE", "KE", C["ke"]), ("EKE", "EKE", C["eke"])]
AXF2 = (AXF[0], AXF[1], AXF[2], 141.5)
ax = axes(AXF2); ax.set_facecolor("white"); YL = (-1.0, 3.0)
ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=0)
ax.axhline(1, color="#555555", lw=0.6, ls=(0, (1, 1.5)), zorder=0)
pts = [(v, "allsat", col, "o", 3.4) for v, _, col in VAR] + [("EKE", "twosat", C["two"], "D", 2.6)]
off = [-0.27, -0.09, 0.09, 0.27]
for gi, (r, _) in enumerate(REG):
    for k, (v, prod, c, mk, ms) in enumerate(pts):
        q = RT[(RT.region == r) & (RT["var"] == v) & (RT["product"] == prod)].iloc[0]; x = gi + off[k]
        lo, hi = max(q.lo, YL[0]), min(q.hi, YL[1])
        ax.plot([x, x], [lo, hi], color=c, lw=0.7, zorder=2)
        if q.lo < YL[0]: ax.plot(x, YL[0] + 0.04, "v", color=c, ms=2.2, zorder=2, clip_on=False)
        if q.hi > YL[1]: ax.plot(x, YL[1] - 0.04, "^", color=c, ms=2.2, zorder=2, clip_on=False)
        ax.plot(x, q.ratio, mk, color=c, mfc=c if q.lo > 0 else "white", ms=ms, mew=0.7, zorder=3)
ax.set_xlim(-0.5, 2.5); ax.set_ylim(*YL)
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in REG]); ax.tick_params(axis="x", length=0)
ax.spines["right"].set_visible(True); ax.spines["left"].set_visible(False)
ax.yaxis.tick_right(); ax.yaxis.set_label_position("right"); ax.set_ylabel("d ln X / d ln τ")
h = [Line2D([], [], marker="o", ls="-", lw=0.7, color=c, ms=3.4, label=l) for _, l, c in VAR] + \
    [Line2D([], [], marker="D", ls="-", lw=0.7, color=C["two"], ms=2.6, label="EKE, two-satellite"),
     Line2D([], [], color="#555555", lw=0.6, ls=(0, (1, 1.5)), label="EKE linear in τ"),
     Line2D([], [], marker="o", ls="", color="#555555", mfc="white", ms=3.4, mew=0.7, label="open: interval includes 0")]
ax.legend(handles=h, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5.2, handlelength=1.3, columnspacing=0.7, borderaxespad=0.1)
label(246.584366, "f", "Sensitivity to the wind stress (ratio of 2001–2025 trends)")
fig.savefig(OUTF, transparent=True); print("ok")
