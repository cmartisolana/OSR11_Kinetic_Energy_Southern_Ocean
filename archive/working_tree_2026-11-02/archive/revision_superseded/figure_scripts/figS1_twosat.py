"""SI Figure S1: all-satellite (Product 1) vs two-satellite (Product 4) trends, % of the mean per decade, 2001-2025:
KE, EKE and along-stream speed for the ACC band and its sectors, and the circumpolar surface transport proxy."""
import pandas as pd
from figstyle import *  # noqa
FIG = OUT / "figures_v23"; FIG.mkdir(exist_ok=True)
A, B = "all-satellite (Product 1)", "two-satellite (Product 4)"
rows = []
for reg in ("ACC", "Pacific", "Atlantic", "Indian"):
    for v in ("KE", "EKE", "SPD"):
        try:
            rows.append((f"{reg} {'along-stream speed' if v == 'SPD' else v}", {A: retest_row(f"allsat|{reg}_{v}|2001-2025"), B: retest_row(f"twosat|{reg}_{v}|2001-2025")}))
        except KeyError:
            pass
for prod, lab in (("allsat", A), ("twosat", B)):
    T = pd.read_csv(OUT / "a8" / f"a8_trends_{prod}.csv").set_index("series")
    for s, nm in (("ACC_T", "ACC surface transport proxy"),):
        r = T.loc[s]; d = dict(pct=r.pct_per_decade, pct_lo=r.pct_lo, pct_hi=r.pct_hi, p=r.p)
        if prod == "allsat": rows.append((nm, {lab: d}))
        else: rows[-1][1][lab] = d
fig, ax = plt.subplots(figsize=(8.5, 8.0))
forest_panel(ax, rows, series_colors={A: "#1F2933", B: "#7A8288"}, series_markers={A: "o", B: "s"}, xlim=(-6, 18))
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=8, frameon=False, ncol=3)
ax.set_title("Trends 2001–2025 in the two altimetry products (% of the mean per decade)", loc="left", fontsize=11, fontweight="bold")
fig.savefig(FIG / "FigureS1_products.png"); fig.savefig(FIG / "FigureS1_products.pdf"); print("saved", len(rows))
