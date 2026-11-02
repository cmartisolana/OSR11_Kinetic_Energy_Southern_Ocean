"""Figure S1 (lean SI, 26 Oct 2026): the two-satellite product. Replaces si_figures.py S1: the withdrawn surface transport
proxy is removed, and the mean-KE map is dropped (it needs the 1.6 GB monthly cache; the ACC band is shown in Fig. 3a).
a  ACC-band KE and EKE, 12-month running means in % of the 2001-2025 mean, Products 1 and 4.
b  ACC-band trends 2001-2025: KE, EKE, its low-frequency and sub-monthly parts, along-stream speed; both products.
c  Sector trends of KE, EKE and along-stream speed; both products.
d  (added 28 Oct 2026; former main-text Fig. 3e) EKE trends of earlier studies over their own periods and of this study (Product 1,
   2001-2025) over the same regions (plan1019 m1): Hogg et al. (2015) sectors, 1993-2012, +/- 90 % [AUTHOR CHECK values];
   Martinez-Moreno et al. (2021) 35-65 S, 1993-2020, no interval.
All trends and intervals from the revised test (outputs/retest/retest_summary.csv). Filled: p < 0.05.
Usage: python figS1_twosat_v2.py <figure dir>   (run from revision/)
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from style_nc import *  # noqa

FIGD = Path(sys.argv[1])
A_, B_ = "All-satellite (Product 1)", "Two-satellite (Product 4)"
COLS = {A_: C["ink"], B_: C["two"]}; MK = {A_: "o", B_: "s"}
pa = lambda y: 100 * (y / np.nanmean(y) - 1)
S = {p: pd.read_csv(OUT / "a5" / f"a5_series_{p}.csv") for p in ("allsat", "twosat")}
x = 2001 + (np.arange(300) + 0.5) / 12

fig = plt.figure(figsize=(W2, 6.3))
g = fig.add_gridspec(3, 1, height_ratios=[0.9, 1.25, 0.85], hspace=0.75)
gb = g[1].subgridspec(1, 2, wspace=0.85)
ax = fig.add_subplot(g[0])
for v, col in (("KE", C["ke"]), ("EKE", C["eke"])):
    for p, ls in (("allsat", "-"), ("twosat", (0, (3, 1.5)))):
        ax.plot(x, running(pa(S[p][f"ACC_{v}"].values)), color=col, ls=ls, lw=1.0,
                label=f"{v}, {'Product 1' if p == 'allsat' else 'Product 4'}")
zero(ax); ts_style(ax); ax.set_ylabel("% of mean")
ax.legend(loc="upper left", ncol=2)
yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.25 * (yl[1] - yl[0]))
plabel(ax, "a", "ACC band, 12-month running means", dx=-24)
terms = [("KE", "KE"), ("EKE", "EKE"), ("LF", "EKE, low-frequency part"), ("SUB", "EKE, sub-monthly part"),
         ("SPD", "Along-stream speed")]
rows = [(lab, {A_: retest_row(f"allsat|ACC_{t}|2001-2025"), B_: retest_row(f"twosat|ACC_{t}|2001-2025")}) for t, lab in terms]
ax = fig.add_subplot(gb[0])
forest(ax, rows, COLS, MK, xlim=(-2, 25), xlabel="Trend 2001–2025 (% per decade)")
plabel(ax, "b", "ACC band", dx=-80)
rows = []
for reg in ("Pacific", "Atlantic", "Indian"):
    for v in ("KE", "EKE", "SPD"):
        rows.append((f"{reg}: {'speed' if v == 'SPD' else v}",
                     {A_: retest_row(f"allsat|{reg}_{v}|2001-2025"), B_: retest_row(f"twosat|{reg}_{v}|2001-2025")}))
ax = fig.add_subplot(gb[1])
forest(ax, rows, COLS, MK, xlim=(-6, 20), xlabel="Trend 2001–2025 (% per decade)", ms=3.2)
for yy in (-2.5, -5.5): ax.axhline(yy, color=C["light"], lw=0.5)
plabel(ax, "c", "Sectors", dx=-60)
fig.legend(handles=sig_handles(COLS, MK), loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.325))
from matplotlib.patches import Patch
PO_ = HERE / "plan1019" / "outputs"
M1 = pd.read_csv(PO_ / "m1" / "m1_trends.csv").set_index("series")
GR = [("Indian sector", "Hogg_Indian_EKE", (6.0, 1.6)), ("Pacific sector", "Hogg_Pacific_EKE", (6.4, 1.8)),
      ("Atlantic sector", "Hogg_Atlantic_EKE", None), ("35–65°S", "box35_65_EKE", (5.1, None)), ("ACC band", "ACC_EKE", None)]
CITE = {"Hogg_Indian_EKE": "Hogg et al.\n(2015)", "Hogg_Pacific_EKE": "Hogg et al.\n(2015)", "Hogg_Atlantic_EKE": "box of Hogg\net al. (2015)", "box35_65_EKE": "Martínez-Moreno\net al. (2021)"}
CPUB = "#222222"; ax = fig.add_subplot(g[2]); zero(ax); w_ = 0.36
for gi, (lab, key, pub) in enumerate(GR):
    r = M1.loc[f"allsat|{key}|2001-2025"]; sg = np.isfinite(r.p) and round(float(r.p), 3) < 0.05
    xs = gi + (w_ / 2 if pub is not None else 0)
    ax.bar(xs, r.pct_per_decade, w_ * 0.9, color=C["eke"] if sg else "white", edgecolor=C["eke"], lw=0.8, zorder=3)
    ax.vlines(xs, r.pct_ci_lo, r.pct_ci_hi, color="#333333", lw=0.5, zorder=4)
    if pub is not None:
        v, e = pub; xp = gi - w_ / 2
        ax.bar(xp, v, w_ * 0.9, color=CPUB, zorder=3)
        if e is not None: ax.vlines(xp, v - e, v + e, color="#777777", lw=0.6, zorder=4)
ax.axvline(3.5, color="#BDBDBD", lw=0.5, ls=(0, (2, 2)))
ax.set_xticks(range(len(GR))); ax.set_xticklabels([q[0] + ("\n" + CITE[q[1]] if q[1] in CITE else "") for q in GR], fontsize=5.5, linespacing=1.05)
ax.tick_params(axis="x", length=0); ax.set_ylim(-2, 22); ax.set_ylabel("EKE trend (% per decade)")
ax.legend(handles=[Patch(facecolor=CPUB, label="earlier study (its own period)"), Patch(facecolor=C["eke"], label="this study, Product 1, 2001–2025"),
                   Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper left", ncol=3, fontsize=5.5)
plabel(ax, "d", "EKE trends: earlier estimates and this study", dx=-24)
save(fig, FIGD, "FigureS1")
