"""SI Figure S1: all-satellite (Product 1) and two-satellite (Product 4) products, 2001-2025.
a  ACC-band KE and EKE, 12-month running means in % of the 2001-2025 mean (a5_series_{allsat,twosat}.csv).
b  ACC-band trends (% per decade, 95 %): KE, EKE, its low-frequency and sub-monthly parts, along-stream speed; both products.
c  sector trends of KE, EKE and along-stream speed; both products.
d  EKE trends of earlier studies over their own periods (black; values typed from the papers: Hogg et al., 2015, Indian 6.0 +- 1.6
   and Pacific 6.4 +- 1.8 % per decade, 1993-2012, 90 %; Martinez-Moreno et al., 2021, 35-65 S, 5.1 % per decade, 1993-2020,
   no interval) and of Product 1, 2001-2025, over the same regions (m1_trends.csv), plus the Atlantic box of Hogg et al. and the
   ACC band. [AUTHOR CHECK the typed literature values.]
All trends and intervals from the revised calibrated test (outputs/retest/retest_summary.csv; m1 for d). Filled: p < 0.05.

Publication repository: ported from figures_final/scripts/figS1_twosat.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/a5/a5_series_{allsat,twosat}.csv, outputs/retest/retest_summary.csv, outputs/m1/m1_trends.csv
Outputs : rendered/{png,pdf}/FigureS1
Run     : python figures/figS1_products.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd
from osr11 import style as S, paths as P  # noqa: E402
from matplotlib.patches import Patch
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else P.FIG; OUT = P.OUT
plt = S.plt; C = S.C
A_, B_ = "Product 1 (all-satellite)", "Product 4 (two-satellite)"
COLS = {A_: C["ink"], B_: C["two"]}; MK = {A_: "o", B_: "s"}
pa = lambda y: 100 * (y / np.nanmean(y) - 1)
SER = {p: pd.read_csv(OUT / "a5" / f"a5_series_{p}.csv") for p in ("allsat", "twosat")}
x = 2001 + (np.arange(300) + 0.5) / 12
rr = lambda k: S.retest_row(P.OUT, k)

H = 168.0
fig = S.figure(S.W2, H)
# a
ax = S.axmm(fig, 18, 128, 148, 32)
for v, col in (("KE", C["ke"]), ("EKE", C["eke"])):
    for p, ls in (("allsat", "-"), ("twosat", (0, (3, 1.5)))):
        ax.plot(x, S.running(pa(SER[p][f"ACC_{v}"].values)), color=col, ls=ls, lw=S.LW["smooth"],
                label=f"{v}, {'Product 1' if p == 'allsat' else 'Product 4'}")
S.zero(ax); S.ts_axis(ax); ax.set_ylabel("% of mean"); S.yticks_minus(ax)
ax.legend(loc="upper left", ncol=4)
yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.25 * (yl[1] - yl[0]))
S.letter(fig, 2, H - 5, "a", "ACC band, 12-month running means")
# b
terms = [("KE", "KE"), ("EKE", "EKE"), ("LF", "EKE, low-frequency part"), ("SUB", "EKE, sub-monthly part"), ("SPD", "Along-stream speed")]
rows = [(lab, {A_: rr(f"allsat|ACC_{t}|2001-2025"), B_: rr(f"twosat|ACC_{t}|2001-2025")}) for t, lab in terms]
ax = S.axmm(fig, 36, 66, 46, 44)
S.forest(ax, rows, COLS, MK, xlim=(-2, 25), xlabel="Trend 2001–2025 (% per decade)")
S.letter(fig, 2, 115, "b", "ACC band")
# c
rows = []
for reg in ("Pacific", "Atlantic", "Indian"):
    for v in ("KE", "EKE", "SPD"):
        rows.append((f"{reg}: {'speed' if v == 'SPD' else v}", {A_: rr(f"allsat|{reg}_{v}|2001-2025"), B_: rr(f"twosat|{reg}_{v}|2001-2025")}))
ax = S.axmm(fig, 120, 66, 46, 44)
S.forest(ax, rows, COLS, MK, xlim=(-6, 20), xlabel="Trend 2001–2025 (% per decade)", ms=3.0)
for yy in (-2.5, -5.5): ax.axhline(yy, color=C["light"], lw=0.5)
S.letter(fig, 88, 115, "c", "Sectors")
fig.legend(handles=S.sig_handles(COLS, MK), loc="center", ncol=3, bbox_to_anchor=(0.5, 54.5 / H))
# d
M1 = pd.read_csv(P.OUT / "m1" / "m1_trends.csv").set_index("series")
GR = [("Indian sector", "Hogg_Indian_EKE", (6.0, 1.6)), ("Pacific sector", "Hogg_Pacific_EKE", (6.4, 1.8)),
      ("Atlantic sector", "Hogg_Atlantic_EKE", None), ("35–65°S", "box35_65_EKE", (5.1, None)), ("ACC band", "ACC_EKE", None)]
CITE = {"Hogg_Indian_EKE": "Hogg et al. (2015)", "Hogg_Pacific_EKE": "Hogg et al. (2015)", "Hogg_Atlantic_EKE": "box of Hogg et al. (2015)",
        "box35_65_EKE": "Martínez-Moreno et al. (2021)"}
CPUB = "#222222"; ax = S.axmm(fig, 18, 12, 148, 30); S.zero(ax); w_ = 0.36
for gi, (lab, key, pub) in enumerate(GR):
    r = M1.loc[f"allsat|{key}|2001-2025"]; sg = S.sigp(r.p)
    xs = gi + (w_ / 2 if pub is not None else 0)
    ax.bar(xs, r.pct_per_decade, w_ * 0.9, color=C["eke"] if sg else "white", edgecolor=C["eke"], lw=0.8, zorder=3)
    ax.vlines(xs, r.pct_ci_lo, r.pct_ci_hi, color="#333333", lw=0.5, zorder=4)
    if pub is not None:
        v, e = pub; xp = gi - w_ / 2
        ax.bar(xp, v, w_ * 0.9, color=CPUB, zorder=3)
        if e is not None: ax.vlines(xp, v - e, v + e, color="#777777", lw=0.6, zorder=4)
ax.axvline(3.5, color="#BDBDBD", lw=0.5, ls=(0, (2, 2)))
ax.set_xticks(range(len(GR))); ax.set_xticklabels([q[0] + ("\n" + CITE[q[1]] if q[1] in CITE else "") for q in GR], linespacing=1.1)
ax.tick_params(axis="x", length=0); ax.set_ylim(-2, 22); ax.set_ylabel("EKE trend (% per decade)"); S.yticks_minus(ax)
ax.legend(handles=[Patch(facecolor=CPUB, label="earlier study (its own period)"), Patch(facecolor=C["eke"], label="this study, Product 1, 2001–2025"),
                   Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper left", ncol=3)
S.letter(fig, 2, 47, "d", "EKE trends: earlier estimates and this study")
S.save(fig, FIGD, "FigureS1"); plt.close(fig)
for k in ("ACC_EKE", "ACC_KE"):
    print(k, {p: round(rr(f"{p}|{k}|2001-2025")["pct"], 2) for p in ("allsat", "twosat")})
