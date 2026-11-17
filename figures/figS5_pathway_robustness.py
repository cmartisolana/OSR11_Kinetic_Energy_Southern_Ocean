"""Supplement Fig. S5 (final layout, 10 Nov 2026): sensitivity of the half-power method, Product 1, 2001-2025.

Same rows and numbers as figures_old/figS5_pathway_robustness.py; layout and presentation only: one panel, rows grouped
(Threshold / Envelope / Minimum width / Retention / Topography) with the group name and the baseline value of the varied
parameter in a left column; the adopted relative threshold (main text; 1826 tracks) is the first row, in bold on a light
highlight; the fixed 200 cm2 s-2 threshold of Chambers (2018) (900 tracks), of which the other rows are one-parameter
variants, is in bold and its 95 % interval is the shaded vertical band (legend); the trend and interval of each row are
written at the right. No left margin. Markers stay open because no trend is significant (filled would mean p < 0.05).
Values (km per decade, positive northward; calibrated Theil-Sen, Sect. S1.7):
  fixed-threshold baseline and parameter variants: outputs/retest/retest_summary.csv, keys <set>|ACC|lat_hp (slope in deg per
  year x 111.195 km x 10); relative threshold: outputs/x7/x7_trends.csv (circumpolar, lat_hp); retention 75 % / 50 %:
  outputs/a2/a2_retention_trends.csv; rough / smooth halves: outputs/x1/x1_bathy_split.csv.
Usage: python figures/figS5_pathway_robustness.py [outputs dir] [figure dir]   (defaults: $OSR11_OUTPUTS, rendered/)
"""
import _sifig as F
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.transforms as mtrans

S, Q, C, plt, OUT, FIGD = F.S, F.Q, F.C, F.plt, F.OUT, F.FIGD
KM = Q.KM
R = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")


def kmrow(key):
    r = R.loc[key]; return dict(val=r.slope * KM * 10, lo=r.ci_lo * KM * 10, hi=r.ci_hi * KM * 10, p=r.p_new)


RT = pd.read_csv(OUT / "a2" / "a2_retention_trends.csv").set_index("min_cov")
B = pd.read_csv(OUT / "x1" / "x1_bathy_split.csv").set_index("subset")
X7 = pd.read_csv(OUT / "x7" / "x7_trends.csv"); X7 = X7[(X7.region == "circumpolar") & (X7.series == "lat_hp")].iloc[0]
rrow = lambda r: dict(val=r.slope, lo=r.lo, hi=r.hi, p=r.p)   # noqa: E731
U = "cm$^2$ s$^{-2}$"
base = kmrow("baseline|ACC|lat_hp")
# (group, group note, [(label, row, style)]) ; style: 'adopted', 'baseline' or ''
GROUPS = [
    ("Threshold", "", [
        (f"relative: 20 % of maximum KE (adopted; {int(X7.n_tracks)} tracks)", dict(val=X7.slope, lo=X7.lo, hi=X7.hi, p=X7.p), "adopted"),
        (f"fixed 200 {U} (Chambers, 2018; 900 tracks)", base, "baseline"),
        (f"fixed 150 {U}", kmrow("thr150|ACC|lat_hp"), ""), (f"fixed 250 {U}", kmrow("thr250|ACC|lat_hp"), ""),
        (f"fixed 300 {U}", kmrow("thr300|ACC|lat_hp"), ""), (f"fixed 400 {U}", kmrow("thr400|ACC|lat_hp"), "")]),
    ("Envelope", "baseline:\nto 25 % of peak", [
        ("to 20 % of peak", kmrow("frac20|ACC|lat_hp"), ""), ("to 33 % of peak", kmrow("frac33|ACC|lat_hp"), ""),
        ("to 50 % of peak", kmrow("frac50|ACC|lat_hp"), "")]),
    ("Minimum width", "baseline:\n100 km", [
        ("50 km", kmrow("minw50|ACC|lat_hp"), ""), ("200 km", kmrow("minw200|ACC|lat_hp"), ""), ("300 km", kmrow("minw300|ACC|lat_hp"), "")]),
    ("Retention", "baseline:\n≥ 90 % of months", [
        (f"≥ 75 % of months ({int(RT.loc[0.75].n_tracks)} tracks)", rrow(RT.loc[0.75]), ""),
        (f"≥ 50 % of months ({int(RT.loc[0.50].n_tracks)} tracks)", rrow(RT.loc[0.50]), "")]),
    ("Topography", "baseline tracks,\nsplit at median\nroughness", [
        (f"rough half ({int(B.loc['rough_half'].n_tracks)} tracks)", rrow(B.loc["rough_half"]), ""),
        (f"smooth half ({int(B.loc['smooth_half'].n_tracks)} tracks)", rrow(B.loc["smooth_half"]), "")]),
]
assert abs(RT.loc[0.90].slope - base["val"]) < 0.01      # retention 90 % = the fixed-threshold baseline

ROW, GAP = 1.0, 0.7
ys, y = [], 0.0
for g, _, rows in GROUPS:
    ys.append([y - i * ROW for i in range(len(rows))]); y = ys[-1][-1] - ROW - GAP
ymin, ymax = ys[-1][-1] - 0.6, 0.6
RMM = 4.3                                   # mm per row
AH = (ymax - ymin) * RMM
H = AH + 25.0
X_GRP, X_ROW, AX0, AXW = 2.0, 33.0, 102.0, 45.0
fig = S.figure(S.W2, H)
ax = S.axmm(fig, AX0, 13.0, AXW, AH)
XL = (-10, 12)
ax.set_ylim(ymin, ymax); ax.set_xlim(*XL)
ax.axvspan(base["lo"], base["hi"], color=C["pos"], alpha=0.10, lw=0, zorder=0)
ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0.5)
tr = mtrans.blended_transform_factory(fig.transFigure, ax.transData)
trv = mtrans.blended_transform_factory(ax.transAxes, ax.transData)
out = []
for (g, gnote, rows), yy in zip(GROUPS, ys):
    ax.text(X_GRP / S.W2, yy[0], g, transform=tr, ha="left", va="center", fontsize=F.FS, fontweight="bold")
    if gnote:
        ax.text(X_GRP / S.W2, yy[0] - 0.45, S.um(gnote), transform=tr, ha="left", va="top", fontsize=F.FS, color=F.GREY_TXT, linespacing=1.1)
    for (lab, r, st), yv in zip(rows, yy):
        bold = st in ("adopted", "baseline")
        ax.text(X_ROW / S.W2, yv, lab, transform=tr, ha="left", va="center", fontsize=F.FS, fontweight="bold" if bold else "normal")
        c = C["pos"]; lw = 1.1 if bold else 0.8; ms = 3.8 if bold else 3.2
        ax.plot([r["lo"], r["hi"]], [yv, yv], color=c, lw=lw, solid_capstyle="butt", zorder=2)
        sg = S.sigp(r["p"])
        ax.plot(r["val"], yv, "o", color=c, mfc=c if sg else "white", mec=c, mew=0.9 if bold else 0.8, ms=ms, zorder=3)
        ax.text(1.03, yv, S.um(f"{r['val']:+.1f} [{r['lo']:+.1f}, {r['hi']:+.1f}]"), transform=trv, ha="left", va="center",
                fontsize=F.FS, color="#222222" if bold else F.GREY_TXT, fontweight="bold" if bold else "normal")
        out.append((g, lab, r))
    if g != GROUPS[-1][0]:
        ax.plot([0, 1], [yy[-1] - (ROW + GAP) / 2] * 2, transform=trv, color=C["light"], lw=0.5, zorder=0, clip_on=False)
# highlight of the adopted row across the label column, the axes and the value column (figure coordinates)
y_ad = ys[0][0]
y0f = ax.transData.transform((0, y_ad - 0.5))[1]; y1f = ax.transData.transform((0, y_ad + 0.5))[1]
inv = fig.transFigure.inverted()
fy0, fy1 = inv.transform((0, y0f))[1], inv.transform((0, y1f))[1]
fig.patches.append(plt.Rectangle(((X_ROW - 1.5) / S.W2, fy0), (S.W2 - X_ROW + 0.5) / S.W2, fy1 - fy0, transform=fig.transFigure,
                                 facecolor="#F1F1F1", edgecolor="none", zorder=-10))
ax.set_ylim(ymin, ymax); ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.patch.set_visible(False)
ax.set_xlim(*XL); ax.set_xticks([-10, -5, 0, 5, 10]); S.xticks_minus(ax)
ax.set_xlabel("Half-power latitude trend, 2001–2025\n(km per decade, positive northward)")
ax.text(1.03, ymax + 0.15, "trend [95 % int.]", transform=trv, ha="left", va="bottom", fontsize=F.FS, color=F.GREY_TXT)
fig.legend(handles=[Line2D([], [], color=C["pos"], lw=0.8, marker="o", mfc="white", mew=0.8, ms=3.2, label="trend and 95 % interval (open: p ≥ 0.05; none significant)"),
                    Patch(facecolor=C["pos"], alpha=0.10, edgecolor="none", label="95 % interval of the fixed 200 cm$^2$ s$^{-2}$ threshold (Chambers, 2018)"),
                    Patch(facecolor="#F1F1F1", edgecolor="none", label="adopted in the main text")],
           loc="upper left", bbox_to_anchor=(X_GRP / S.W2, (H - 1.5) / H), ncol=2, handlelength=2.0, columnspacing=1.6, labelspacing=0.4,
           borderaxespad=0)
F.check(fig); S.save(fig, FIGD, "FigureS5"); plt.close(fig)
for g, lab, r in out:
    print(f"S5 {g:13s} {lab.replace('$^2$', '2').replace('$^{-2}$', '-2'):55s} {r['val']:+6.2f} [{r['lo']:+6.2f}, {r['hi']:+6.2f}] p={r['p']:.3f}")
vals = [r["val"] for g, lab, r in out if g in ("Threshold", "Envelope", "Minimum width") and "relative" not in lab and "fixed 200" not in lab]
print("S5 alternative fixed thresholds, envelopes, minimum widths: %+.1f to %+.1f km per decade" % (min(vals), max(vals)))
