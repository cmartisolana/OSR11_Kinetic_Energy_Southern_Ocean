"""Figure (draft, v1): response of the eddies and of the mean flow to the wind stress (eddy saturation), 2001-2025 (dk31).
a  ACC band, annual anomalies (% of mean, linearly detrended): eastward wind stress of the PREVIOUS year, EKE and the
   kinetic energy of the annual-mean flow (MKE).
b  Sensitivity (% per % wind stress, detrended annual means) of EKE and MKE to the wind of single earlier years (lag 0-4)
   and to the wind averaged over the last N years (N = 1-6): ACC band with 95 % intervals (filled: p < 0.05), sectors thin.
   Reference: 1 = EKE proportional to the wind stress (eddy-saturated models); 0 = no response.
c  Trends 2001-2025 (% per decade) of the wind stress, MKE and EKE, ACC band and sectors; 95 % intervals with an
   effective sample size from the lag-1 autocorrelation of the residuals; open bars: p >= 0.05.
Usage: python fig_saturation_v1.py <revision dir> <dk31 out> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
REV, D31, FIGD = (Path(a) for a in sys.argv[1:4])
sys.path.insert(0, str(REV))
import style_nc as S
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C
D = pd.read_csv(D31 / "dk31_response.csv"); Z = np.load(D31 / "dk31_series.npz")
CW, CE, CM = C["wind"], C["eke"], C["ke"]
REGS = ["ACC band", "Pacific", "Indian", "Atlantic"]
yrs = np.arange(2001, 2026); t = yrs - yrs.mean()
def trend_ci(y):
    b, a = np.polyfit(t, y, 1); res = y - a - b * t; r1 = np.corrcoef(res[1:], res[:-1])[0, 1]; ne = max(3, 25 * (1 - r1) / (1 + r1))
    se = np.sqrt((res ** 2).sum() / (ne - 2) / (t ** 2).sum()); p = 2 * stats.t.sf(abs(b / se), ne - 2)
    return 10 * b, 10 * stats.t.ppf(0.975, ne - 2) * se, p
fig = plt.figure(figsize=(S.W2, 0.78 * S.W2))
g = fig.add_gridspec(2, 1, height_ratios=[0.85, 1.0], hspace=0.42)
gB = g[1].subgridspec(1, 3, width_ratios=[0.75, 1.0, 1.05], wspace=0.12)
# a
ax = fig.add_subplot(g[0]); S.zero(ax)
s = "ACC band"; tau, e, k = Z[f"{s}|tau"], Z[f"{s}|eke"], Z[f"{s}|mke"]
ax.bar(yrs[1:] + 0.5, tau[:-1], 0.75, color="#D9D9D9", lw=0, zorder=1, label="wind stress, previous year")
ax.plot(yrs + 0.5, e, "-o", color=CE, lw=1.1, ms=2.2, zorder=4, label="EKE")
ax.plot(yrs + 0.5, k, "-s", color=CM, lw=0.9, ms=1.9, zorder=3, label="mean-flow KE")
q = D[(D.region == s) & (D.kind == "lag") & (D.k == 1)].set_index("var")
ax.text(0.99, 0.04, f"r with previous-year wind: EKE {q.loc['EKE', 'r']:.2f} (p = {q.loc['EKE', 'p']:.2f}), mean-flow KE {q.loc['MKE', 'r']:.2f} (p = {q.loc['MKE', 'p']:.2f})",
        transform=ax.transAxes, ha="right", fontsize=5.0, color="#333333")
S.ts_style(ax, (2001.0, 2026.0)); ax.set_ylabel("anomaly, detrended (%)")
ax.legend(loc="upper left", ncol=3, fontsize=5.0)
S.plabel(ax, "a", "ACC band: EKE, mean-flow KE and the previous year's wind stress (detrended)", dy=4)
# b (two axes: lag, window)
axL = fig.add_subplot(gB[0]); axW = fig.add_subplot(gB[1], sharey=axL)
for ax_, kind, ks in ((axL, "lag", range(0, 5)), (axW, "window", range(1, 7))):
    S.zero(ax_); ax_.axhline(1, color="#BDBDBD", lw=0.6, ls=(0, (3, 2)), zorder=0)
    for var, col, dx in (("EKE", CE, -0.12), ("MKE", CM, 0.12)):
        for reg, ls in zip(REGS[1:], ("-", (0, (3, 1.4)), (0, (1, 1)))):
            qq = D[(D.region == reg) & (D.kind == kind) & (D["var"] == var)]
            ax_.plot(qq.k + dx, qq.sens, color=col, lw=0.5, ls=ls, alpha=0.6, zorder=2)
        qq = D[(D.region == "ACC band") & (D.kind == kind) & (D["var"] == var)]
        ax_.vlines(qq.k + dx, qq.lo, qq.hi, color=col, lw=0.6, zorder=3)
        sg = qq.p < 0.05
        ax_.plot(qq.k + dx, qq.sens, color=col, lw=1.0, zorder=4)
        ax_.plot((qq.k + dx)[sg], qq.sens[sg], "o", color=col, ms=2.6, zorder=5)
        ax_.plot((qq.k + dx)[~sg], qq.sens[~sg], "o", color=col, mfc="white", mew=0.7, ms=2.6, zorder=5)
    ax_.set_xticks(list(ks))
axL.set_xlabel("wind of year y − L  (lag L, yr)"); axW.set_xlabel("wind averaged over the last N years")
axL.set_ylabel("sensitivity (% per % wind stress)"); axL.set_ylim(-1.5, 4.5)
plt.setp(axW.get_yticklabels(), visible=False)
axL.text(4.3, 1.05, "EKE ~ τ", fontsize=4.8, color="#777777", ha="right", va="bottom")
axW.legend(handles=[Line2D([], [], color=CE, marker="o", ms=2.5, lw=1.0, label="EKE, ACC band"), Line2D([], [], color=CM, marker="o", ms=2.5, lw=1.0, label="mean-flow KE, ACC band"),
                    Line2D([], [], color="#888888", lw=0.5, label="Pacific / Indian / Atlantic"), Line2D([], [], color="#555555", marker="o", mfc="white", ls="", ms=2.5, label="open: p ≥ 0.05")],
           loc="upper left", fontsize=4.6)
S.plabel(axL, "b", "Response to the past wind stress (interannual)", dy=4)
# c trends
ax = fig.add_subplot(gB[2]); S.zero(ax); w_ = 0.26
for gi, reg in enumerate(REGS):
    for j, (v, col) in enumerate((("tau", CW), ("mke", CM), ("eke", CE))):
        b, ci, p = trend_ci(Z[f"{reg}|{v}_raw"]); x = gi + (j - 1) * w_
        ax.bar(x, b, w_ * 0.9, color=col if p < 0.05 else "white", edgecolor=col, lw=0.8, zorder=3)
        ax.vlines(x, b - ci, b + ci, color="#333333", lw=0.5, zorder=4)
ax.set_xticks(range(4)); ax.set_xticklabels(REGS, fontsize=5.4); ax.tick_params(axis="x", length=0)
ax.set_ylabel("trend 2001–2025 (% per decade)"); ax.yaxis.set_label_position("right"); ax.yaxis.tick_right()
ax.legend(handles=[Patch(facecolor=CW, label="wind stress"), Patch(facecolor=CM, label="mean-flow KE"), Patch(facecolor=CE, label="EKE"),
                   Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper left", fontsize=4.6, ncol=2)
S.plabel(ax, "c", "Trends of the wind, the mean flow and the eddies", dy=4)
ax.set_ylim(-4, 22)
S.save(fig, FIGD, "FigureX_eddy_saturation"); plt.close(fig)
for reg in REGS: print(reg, [tuple(np.round(trend_ci(Z[f"{reg}|{v}_raw"]), 3)) for v in ("tau", "mke", "eke")])
