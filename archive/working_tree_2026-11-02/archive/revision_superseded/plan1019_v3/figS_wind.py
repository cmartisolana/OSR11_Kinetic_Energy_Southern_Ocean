"""SI Figure S11: eastward wind stress (Product 2), 2001-2025: monthly deseasonalised anomalies, annual means and trends."""
import sys, numpy as np, pandas as pd
from pathlib import Path
REV = Path("/mnt/user-data/uploads/repository/revision"); sys.path.insert(0, str(REV))
import style_nc as S
plt = S.plt; C = S.C
S1 = pd.read_csv("/home/claude/work/drake/out/dk1/dk1_series.csv"); T1 = pd.read_csv("/home/claude/work/drake/out/dk1/dk1_trends.csv")
T1 = T1[T1.period == "2001-2025"].set_index("name")
W6 = pd.read_csv(REV / "outputs/a6/a6_wind_stress_series.csv"); W6["year"] = W6.year.astype(int); W6["month"] = W6.month.astype(int)
RT = pd.read_csv(REV / "outputs/retest/retest_summary.csv").set_index("key")
r = RT.loc["tau_Circumpolar_45_60S|2001-2025"]; f = 1000 / (r.slope * 1000 / r.pct_per_decade)
rows = [(S1, "wind_up_TAUX", "Upstream of Drake Passage (100–80°W, 63–54°S)", T1.loc["wind_up_TAUX", ["pct", "pct_lo", "pct_hi", "p"]].values),
        (S1, "wind_box_TAUX", "Drake box (80–52°W, 63–54°S)", T1.loc["wind_box_TAUX", ["pct", "pct_lo", "pct_hi", "p"]].values),
        (S1, "wind_down_TAUX", "Downstream (55–35°W, 63–54°S)", T1.loc["wind_down_TAUX", ["pct", "pct_lo", "pct_hi", "p"]].values),
        (W6, "Circumpolar_45_60S", "Circumpolar, 45–60°S", np.array([r.pct_per_decade, r.ci_lo * f, r.ci_hi * f, r.p_new]))]
fig = plt.figure(figsize=(S.W2, 0.6 * S.W2)); gs = fig.add_gridspec(4, 2, width_ratios=[1, 0.1], hspace=0.35, wspace=0.05)
for k, (df, col, ttl, (v, lo, hi, p)) in enumerate(rows):
    ax = fig.add_subplot(gs[k, 0]); axt = fig.add_subplot(gs[k, 1])
    y = df[col].astype(float); a = 100 * (y - y.groupby(df.month.values).transform("mean")) / y.mean()
    x = df.year + (df.month - 0.5) / 12
    ax.plot(x, a, color=C["wind"], lw=0.35, alpha=0.5); an = a.groupby(df.year.values).mean()
    ax.plot(an.index + 0.5, an.values, "-o", color="#444444", lw=0.9, ms=2.2, mec="white", mew=0.3)
    S.zero(ax); S.ts_style(ax); ax.set_ylim(-60, 60); ax.set_ylabel("%")
    ax.text(0.01, 0.95, ttl, transform=ax.transAxes, fontsize=5.5, va="top", bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.5))
    if k < 3: ax.set_xticklabels([])
    sg = round(float(p), 3) < 0.05
    axt.plot([0, 0], [lo, hi], color="#444444", lw=0.9); axt.plot(0, v, "o", color="#444444", mfc="#444444" if sg else "white", ms=3.2)
    axt.axhline(0, color="#9A9A9A", lw=0.5); axt.set_xlim(-1, 1); axt.set_xticks([]); axt.set_ylim(-5, 20)
    axt.yaxis.tick_right(); axt.spines["left"].set_visible(False); axt.spines["right"].set_visible(True); axt.spines["bottom"].set_visible(False)
    axt.tick_params(axis="y", labelsize=5); axt.yaxis.set_label_position("right"); axt.set_ylabel("% dec$^{-1}$", fontsize=5)
    if k == 0: axt.set_title("trend", fontsize=5.5, pad=2)
fig.savefig("figS/new_S11.png", dpi=300, bbox_inches="tight", pad_inches=0.02); print("ok")
