"""SI Figure S2 (new script, 2 Nov 2026; no earlier script produced this figure): Drake Passage EKE year by year and by period,
from existing outputs only (no new analysis).
a  annual EKE between the mean SAF and PF in 2 deg longitude bins, 86-40 W, anomaly in % of each bin's 2001-2025 mean, Product 1
   (dk12_band_annual_allsat.csv; the band and bins of main Figs. 2c and 3d); dashed: the narrow passage 68-58 W; horizontal line:
   the 2001-2015 / 2016-2025 split.
b  the same averaged over the bins upstream (86-68 W), in the passage (68-58 W) and downstream (58-40 W) (unweighted bin means);
   horizontal segments: means of 2001-2015 and 2016-2025; labels: increase of the 2016-2025 regional mean EKE over the 2001-2015
   regional mean (ratio of the regional means of dk8_along_acc_fronts_periods_allsat.csv, the values quoted in Sect. 3.1).
c  EKE trends 2001-2025 (% per decade, 95 %, calibrated test) in the three fixed regional domains of dk1 (upstream 100-80 W in the
   ACC band; Drake box 80-52 W, 63-54 S; downstream 55-35 W in the ACC band; ice-masked), Products 1 and 4 (dk1_trends.csv).
d  2016-2025 minus 2001-2015 mean EKE in the same domains, % of the 2001-2025 mean, Products 1 and 4 (from dk1_series.csv).
A Product 4 version of a-b needs dk12_circumpolar_band.py run on data/monthly_twosat (author's Mac).

Publication repository: ported from figures_final/scripts/figS2_drake_periods.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/drake/dk12/dk12_band_annual_allsat.csv, outputs/drake/dk1/*, outputs/drake/dk8/dk8_along_acc_fronts_periods_allsat.csv
Outputs : rendered/{png,pdf}/FigureS2
Run     : python figures/figS2_drake_periods.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd
from osr11 import style as S, paths as P  # noqa: E402
from matplotlib.colors import TwoSlopeNorm
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else P.FIG
plt = S.plt; C = S.C
DO = P.OUT / "drake"
A = pd.read_csv(DO / "dk12" / "dk12_band_annual_allsat.csv", index_col=0); lon = A.columns.astype(float).values
sel = (lon > -86) & (lon < -40); A = A.loc[:, sel]; lon = lon[sel]; yrs = A.index.values.astype(int)
DOM = [("upstream, 86–68°W", (-86, -68), "#7FC9B0"), ("passage, 68–58°W", (-68, -58), C["eke"]), ("downstream, 58–40°W", (-58, -40), "#00513A")]
T1 = pd.read_csv(DO / "dk1" / "dk1_trends.csv").set_index("series"); S1 = pd.read_csv(DO / "dk1" / "dk1_series.csv")
A_, B_ = "Product 1 (all-satellite)", "Product 4 (two-satellite)"; COLS = {A_: C["ink"], B_: C["two"]}; MK = {A_: "o", B_: "s"}
REG = [("up", "Upstream, 100–80°W"), ("box", "Drake box, 80–52°W"), ("down", "Downstream, 55–35°W")]

H = 140.0
fig = S.figure(S.W2, H)
# a ---------------------------------------------------------------- Hovmoller
ax = S.axmm(fig, 16, 66, 62, 66)
im = ax.pcolormesh(S.edges(lon), np.r_[yrs - 0.5, yrs[-1] + 0.5], A.values, cmap=S.TREND_CMAP, norm=TwoSlopeNorm(0, -40, 40),
                   shading="flat", rasterized=True)
for x_ in (-68, -58): ax.axvline(x_, color="#555555", lw=0.6, ls=(0, (2, 2)))
ax.axhline(2015.5, color="#222222", lw=0.6)
ax.set_xlim(-86, -40); ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels([S.lonfmt(v) for v in (-80, -70, -60, -50, -40)])
ax.set_ylim(2025.5, 2000.5); ax.set_yticks([2001, 2005, 2010, 2015, 2020, 2025]); ax.set_ylabel("Year")
ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
S.letter(fig, 2, H - 5, "a", "Annual EKE between the SAF and PF, Product 1")
S.cbar(fig, im, 16 + 0.1 * 62, 55, 0.8 * 62, "EKE anomaly (% of bin mean)", extend="both", ticks=[-40, -20, 0, 20, 40])
# b ---------------------------------------------------------------- domain means and period means
ax = S.axmm(fig, 100, 66, 64, 66); diffs = {}; labs = []; ratios = {}
P8 = pd.read_csv(P.OUT / "drake" / "dk8" / "dk8_along_acc_fronts_periods_allsat.csv").pivot(index="lon", columns="period", values="mean")
for i, (nm, (a, b), col) in enumerate(DOM):
    m = (lon > a) & (lon < b); s = A.loc[:, m].mean(axis=1).values
    e, l = s[yrs <= 2015].mean(), s[yrs >= 2016].mean(); diffs[nm] = l - e
    ax.plot(yrs, s, "-o", color=col, lw=0.9, ms=2.2, zorder=3, label=nm)
    ax.plot([2000.6, 2015.4], [e, e], color=col, lw=2.0, alpha=0.6, solid_capstyle="butt", zorder=2)
    ax.plot([2015.6, 2025.4], [l, l], color=col, lw=2.0, alpha=0.6, solid_capstyle="butt", zorder=2)
    q8 = P8[(P8.index > a) & (P8.index < b)]; rel = 100 * (q8["2016-2025"].mean() / q8["2001-2015"].mean() - 1); ratios[nm] = rel
    labs.append([l, S.um(f"{rel:+.0f} %"), col])
labs.sort(key=lambda q: -q[0])
for i in range(1, len(labs)): labs[i][0] = min(labs[i][0], labs[i - 1][0] - 2.6)
for yv, t_, col in labs: ax.text(2025.9, yv, t_, color=col, fontsize=S.FS_ANN, ha="left", va="center")
ax.axvline(2015.5, color="#BDBDBD", lw=0.5, ls=(0, (2, 2)))
S.zero(ax); S.ts_axis(ax, 2000.5, 2025.5, 2028.0, ticks=(2001, 2005, 2010, 2015, 2020, 2025))
ax.set_ylabel("EKE anomaly (% of bin mean)"); S.yticks_minus(ax); ax.legend(loc="upper left", fontsize=S.FS_ANN)
yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.15 * (yl[1] - yl[0]))
S.letter(fig, 86, H - 5, "b", "Upstream, passage and downstream means, Product 1")
print("S2b late-minus-early (% points of 2001-2025 bin mean):", {k: round(v, 1) for k, v in diffs.items()}, "; labels, 2016-2025 vs 2001-2015 regional mean (dk8 periods, %):", {k: round(v, 1) for k, v in ratios.items()})
# c ---------------------------------------------------------------- dk1 trends, both products
def trow(p, r):
    q = T1.loc[f"{p}_{r}_EKE|2001-2025"]; return dict(pct=q.pct, pct_lo=q.pct_lo, pct_hi=q.pct_hi, p=q.p)
rows = [(lab, {A_: trow("allsat", r), B_: trow("twosat", r)}) for r, lab in REG]
ax = S.axmm(fig, 40, 17, 40, 16)
S.forest(ax, rows, COLS, MK, xlim=(-5, 25), xlabel="EKE trend 2001–2025 (% per decade)")
S.letter(fig, 2, 38, "c", "Regional domains, trend")
# d ---------------------------------------------------------------- dk1 period differences, both products
ax = S.axmm(fig, 126, 17, 40, 16); ys = -np.arange(3, dtype=float); out = {}
for j, (p, nm) in enumerate((("allsat", A_), ("twosat", B_))):
    for i, (r, lab) in enumerate(REG):
        v = S1[f"{p}_{r}_EKE"].values; mu = v.mean(); d = 100 * (v[S1.year >= 2016].mean() - v[S1.year <= 2015].mean()) / mu
        out[(p, r)] = d; yy = ys[i] - (j - 0.5) * 0.26
        ax.plot([0, d], [yy, yy], color=COLS[nm], lw=0.8); ax.plot(d, yy, MK[nm], color=COLS[nm], ms=3.4)
ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0); ax.set_yticks(ys); ax.set_yticklabels([lab for _, lab in REG])
ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False); ax.set_ylim(-2.6, 0.6); ax.set_xlim(-5, 30)
ax.set_xlabel("Difference (% of mean)"); S.xticks_minus(ax)
S.letter(fig, 86, 38, "d", "Regional domains, period difference")
fig.legend(handles=S.sig_handles(COLS, MK), loc="center", ncol=3, bbox_to_anchor=(0.5, 3.0 / H))
S.save(fig, FIGD, "FigureS2"); plt.close(fig)
print("S2c trends:", {k: {p: round(v["pct"], 1) for p, v in d.items()} for k, d in rows})
print("S2d period differences (%):", {k: round(v, 1) for k, v in out.items()})
