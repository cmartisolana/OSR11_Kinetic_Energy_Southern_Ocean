"""Supplement Fig. S4 (final layout, 10 Nov 2026): comparison with previous studies.

Same data and numbers as figures_old/figS4_previous_studies.py (panels a-d, same rows); layout and presentation only:
panels a-c stacked with one label column and aligned rows (study + region on the first line, published period and our
period on a second, grey line); the marker semantics are stated once in a legend at the top (black diamond = published
value with its interval where given; colour = this study over the overlap with the published period; light colour = this
study over 2001-2025; circle/solid = Product 1, square/dashed = Product 4; open = p >= 0.05); the purple/grey key of c sits
next to c and the line key of d above d. x ranges widened so that no interval is clipped (a: -30..32, b: -25..35,
c: -50..55); clipped intervals would still end in an arrow head.
Values: outputs/m1/m1_previous_studies.csv (05k_previous_studies.py). d: figures_old/x11/x11_lag_corr.csv
(x11_lag_previous_studies.py; variant "as Hogg et al. (2015): not detrended"), the 3-year values of Hogg et al. (2015) as
diamonds (r = 0.50 Pacific, 0.32 Indian).
  a  EKE trends in % per decade: Hogg et al. (2015) Indian and Pacific boxes; Martinez-Moreno et al. (2021) ACC region
     (ours: 35-65 S area integral, not their mask); Barcelo-Llull et al. (2025) high-EKE regions, all-sat and two-sat
     (ours: ACC band).
  b  EKE trends in cm2 s-2 per decade: Zhang et al. (2021), Hogg boxes (Pacific box 62-40 S as published by them) and the
     180-150 W sector (ours: between the mean SAF and PF).
  c  Front-position trends in km per decade (positive northward): Chambers (2018) half-power point, 90-170 E and 160-90 W;
     SSH-contour displacement of Sokolov and Rintoul (2009) as summarised by Gille (2014), against our fixed ADT contour
     and half-power latitude (2001-2025).
  d  Correlation of 12-month running-mean EKE in the boxes of Hogg et al. (2015) with the 12-month running-mean SAM index
     (Marshall, 2003), SAM leading, 2001-2025; dots: maxima (filled: pointwise p < 0.05).
Usage: python figures/figS4_previous_studies.py [outputs dir] [figure dir]   (defaults: $OSR11_OUTPUTS, rendered/)
"""
import _sifig as F
import numpy as np, pandas as pd
import matplotlib.transforms as mtrans
from matplotlib.lines import Line2D

S, Q, C, plt, OUT, FIGD = F.S, F.Q, F.C, F.plt, F.OUT, F.FIGD
T = pd.read_csv(OUT / "m1" / "m1_previous_studies.csv")
LIGHT = {C["eke"]: "#8FD6C2", C["pos"]: "#D3A9C8", C["grey"]: "#BDBDBD"}
PUB = "#111111"
MS_PUB, MS_OUR = 3.2, 3.2


def get(study, region, series, period):
    r = T[(T.study == study) & (T.region == region) & (T.our_series == series) & (T.our_period == period)]
    assert len(r) == 1, (study, region, series, period, len(r))
    return r.iloc[0]


NCLIP = [0]


def seg(ax, lo, hi, y, col, ls, xl):
    """interval clipped to the axis, with an arrow head where clipped"""
    a, b = max(lo, xl[0]), min(hi, xl[1])
    ax.plot([a, b], [y, y], color=col, lw=0.8, ls=ls, solid_capstyle="butt", zorder=2)
    for end, cl, d in ((lo, a, -1), (hi, b, 1)):
        if end != cl:
            NCLIP[0] += 1
            ax.annotate("", xy=(cl, y), xytext=(cl - d * 0.04 * (xl[1] - xl[0]), y),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color=col, mutation_scale=4, shrinkA=0, shrinkB=0), zorder=2)


def panel(fig, ax, groups, xl, xlabel, xticks, dodge=0.17, lab_x=3.0):
    """groups: [(line 1, line 2, published row, [(colour, prod, row), ...])], top first; labels in the figure label column"""
    ys = -np.arange(len(groups), dtype=float)
    tr = mtrans.blended_transform_factory(fig.transFigure, ax.transData)
    for (l1, l2, pub, ent), y in zip(groups, ys):
        n = len(ent) + 1; off = lambda i: y + ((n - 1) / 2 - i) * dodge   # noqa: E731
        yp = off(0)
        if np.isfinite(pub.published_lo):
            seg(ax, pub.published_lo, pub.published_hi, yp, PUB, "-", xl)
        ax.plot(pub.published_value, yp, "D", color=PUB, ms=MS_PUB, zorder=3)
        for i, (col, pk, r) in enumerate(ent, start=1):
            yy = off(i); st = Q.PROD[pk]
            seg(ax, r.our_lo, r.our_hi, yy, col, st["ls"], xl)
            sg = S.sigp(r.p)
            ax.plot(r.our_value, yy, st["marker"], color=col, mfc=col if sg else "white", mec=col, mew=0.8, ms=MS_OUR, zorder=3)
        ax.text(lab_x / S.W2, y + 0.04, l1, transform=tr, ha="left", va="bottom", fontsize=F.FS, color="#222222")
        ax.text(lab_x / S.W2, y - 0.04, S.um(l2), transform=tr, ha="left", va="top", fontsize=F.FS, color=F.GREY_TXT)
    for y in ys[:-1]:
        ax.axhline(y - 0.5, color=C["light"], lw=0.5, zorder=0)
    ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0)
    ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_ylim(ys[-1] - 0.5, ys[0] + 0.5)
    ax.set_xlim(*xl); ax.set_xticks(xticks); ax.set_xlabel(xlabel); S.xticks_minus(ax)
    return ys


def eke_entries(study, region, matched, prods=("P1", "P4"), col=C["eke"]):
    e = []
    for per, c in ((matched, col), ("2001-2025", LIGHT[col])):
        for pk in prods:
            e.append((c, pk, get(study, region, Q.PROD[pk]["name"].split(" (")[0], per)))
    return e


HO, MM, BL, ZH = "Hogg et al. (2015)", "Martinez-Moreno et al. (2021)", "Barcelo-Llull et al. (2025)", "Zhang et al. (2021)"
RI, RP = "Indian box 40-150E, 57-44S", "Pacific box 150E-72W, 62-48S"
RM = "ACC region (their mask: 35-60S without boundary currents and gyres); ours: 35-65S area integral, not their ACC-region mask"
RB = "high-EKE regions (global, 5 % of ocean area); ours: ACC band (not the same region)"
ZI, ZP, ZS = "Hogg Indian box 40-150E, 57-44S", "Hogg Pacific box 150E-72W, 62-40S (as published by Zhang et al. 2021, p. 11)", "30 deg sector 180-150W (ours: between mean SAF and PF)"
CH, GI = "Chambers (2018)", "Gille (2014), after Sokolov and Rintoul (2009)"
HP, CT = "Product 1, half-power latitude", "Product 1, fixed ADT contour"

ga = [("Hogg et al. (2015), Indian box", "1993–2012 · ours 2001–2012", get(HO, RI, "Product 1", "2001-2012"), eke_entries(HO, RI, "2001-2012")),
      ("Hogg et al. (2015), Pacific box", "1993–2012 · ours 2001–2012", get(HO, RP, "Product 1", "2001-2012"), eke_entries(HO, RP, "2001-2012")),
      ("Martínez-Moreno et al. (2021), ACC region", "1993–2020 · ours 35–65°S, 2001–2020", get(MM, RM, "Product 1", "2001-2020"), eke_entries(MM, RM, "2001-2020")),
      ("Barceló-Llull et al. (2025), high EKE, all-sat.", "1993–2022 · ours ACC band, 2001–2022", get(BL, RB, "Product 1", "2001-2022"), eke_entries(BL, RB, "2001-2022", prods=("P1",))),
      ("Barceló-Llull et al. (2025), high EKE, two-sat.", "1993–2022 · ours ACC band, 2001–2022", get(BL, RB, "Product 4", "2001-2022"), eke_entries(BL, RB, "2001-2022", prods=("P4",)))]
gb = [("Zhang et al. (2021), Indian box", "1993–2019 · ours 2001–2019", get(ZH, ZI, "Product 1", "2001-2019"), eke_entries(ZH, ZI, "2001-2019")),
      ("Zhang et al. (2021), Pacific box 62–40°S", "1993–2019 · ours 2001–2019", get(ZH, ZP, "Product 1", "2001-2019"), eke_entries(ZH, ZP, "2001-2019", prods=("P1",))),
      ("Zhang et al. (2021), 180–150°W", "1993–2019 · ours 2001–2019, between the fronts", get(ZH, ZS, "Product 1", "2001-2019"), eke_entries(ZH, ZS, "2001-2019", prods=("P1",)))]


def pe(study, region, matched, series=HP, col=C["pos"]):
    e = []
    if matched:
        e.append((col, "P1", get(study, region, series, matched)))
    e.append((LIGHT[col], "P1", get(study, region, series, "2001-2025")))
    return e


gc = [("Chambers (2018), 90–170°E", "1993–2016 · ours 2001–2016", get(CH, "90-170E", HP, "2001-2016"), pe(CH, "90-170E", "2001-2016")),
      ("Chambers (2018), 160–90°W", "1993–2016 · ours 2001–2016", get(CH, "160-90W (200-270E)", HP, "2001-2016"), pe(CH, "160-90W (200-270E)", "2001-2016")),
      ("Sokolov and Rintoul (2009), SSH contours", "1993–2007, in Gille (2014) · ours 2001–2025", get(GI, "circumpolar, SSH contours", CT, "2001-2025"),
       pe(GI, "circumpolar, SSH contours", None, CT, C["grey"]) + pe(GI, "circumpolar, SSH contours", None))]

# ---------------------------------------------------------------- geometry (mm)
H = 228.0; AX0, AXW = 66.0, 100.0; GA, GB = 10.0, 9.0          # axes left/width; group heights
fig = S.figure(S.W2, H)
# top legend: marker semantics stated once
h = [Line2D([], [], ls="-", lw=0.8, marker="D", color=PUB, ms=MS_PUB, label="published value (its interval where given)"),
     Line2D([], [], ls="none", marker="o", color="#444444", mfc="white", mew=0.8, ms=MS_OUR, label="open: p ≥ 0.05"),
     Line2D([], [], ls="-", lw=0.8, marker="o", color=C["eke"], ms=MS_OUR, label="this study, overlap with the published period"),
     Line2D([], [], ls="-", lw=0.8, marker="o", color=LIGHT[C["eke"]], ms=MS_OUR, label="this study, 2001–2025 (light colour)"),
     Line2D([], [], ls=Q.PROD["P1"]["ls"], lw=0.8, marker="o", color="#444444", ms=MS_OUR, label="Product 1 (all-satellite)"),
     Line2D([], [], ls=Q.PROD["P4"]["ls"], lw=0.8, marker="s", color="#444444", ms=MS_OUR, label="Product 4 (two-satellite)")]
fig.legend(handles=h, loc="upper left", bbox_to_anchor=(2 / S.W2, (H - 1.5) / H), ncol=3, handlelength=2.4, columnspacing=1.6,
           labelspacing=0.45, borderaxespad=0)
# a ---------------------------------------------------------------- relative EKE trends
ya_top = H - 19.0; ha = GA * len(ga)
ax = S.axmm(fig, AX0, ya_top - ha, AXW, ha)
panel(fig, ax, ga, (-30, 32), "EKE trend (% per decade)", [-30, -20, -10, 0, 10, 20, 30])
S.letter(fig, 2, ya_top + 2.5, "a", "EKE, relative trends")
# b ---------------------------------------------------------------- absolute EKE trends (Zhang et al. 2021)
yb_top = ya_top - ha - 17.0; hb = GB * len(gb)
ax = S.axmm(fig, AX0, yb_top - hb, AXW, hb)
panel(fig, ax, gb, (-25, 35), "EKE trend (cm$^2$ s$^{−2}$ per decade)", [-20, -10, 0, 10, 20, 30])
S.letter(fig, 2, yb_top + 2.5, "b", "EKE, absolute trends")
# c ---------------------------------------------------------------- front position
yc_top = yb_top - hb - 17.0; hc = GB * len(gc)
ax = S.axmm(fig, AX0, yc_top - hc, AXW, hc)
panel(fig, ax, gc, (-50, 55), "Trend (km per decade, positive northward)", [-40, -20, 0, 20, 40], dodge=0.2)
S.letter(fig, 2, yc_top + 2.5, "c", "Front position")
ax.legend(handles=[Line2D([], [], ls="-", lw=0.8, marker="o", color=C["pos"], ms=MS_OUR, label="half-power latitude"),
                   Line2D([], [], ls="-", lw=0.8, marker="o", color=C["grey"], ms=MS_OUR, label="fixed ADT contour")],
          loc="lower right", bbox_to_anchor=(1.0, 1.0 + 1.2 / hc), ncol=2, handlelength=2.0, columnspacing=1.4, borderaxespad=0)
# d ---------------------------------------------------------------- lag behind the SAM index, method of Hogg et al. (2015)
X11 = pd.read_csv(OUT / "x11" / "x11_lag_corr.csv")
X11 = X11[X11.variant.str.startswith("as Hogg")]
yd_let = yc_top - hc - 16.0; yd_top = yd_let - 11.5; hd = yd_top - 12.0
ax = S.axmm(fig, 16.0, 12.0, AX0 + AXW - 16.0, hd); S.zero(ax)
LS = {"Pacific": "-", "Indian": (0, (4, 1.6)), "Atlantic": (0, (1, 1.2))}
mx = {}
for reg, ls in LS.items():
    q = X11[X11.region == reg]; ax.plot(q.lag_months / 12, q.r, color=C["eke"], lw=1.0, ls=ls, zorder=3, label=f"{reg} box")
    k = q.r.idxmax(); mx[reg] = (int(q.lag_months[k]), round(q.r[k], 2), round(q.p_pointwise[k], 3))
    ax.plot(q.lag_months[k] / 12, q.r[k], "o", color=C["eke"], mfc=C["eke"] if q.p_pointwise[k] < 0.05 else "white", mec=C["eke"],
            mew=0.8, ms=MS_OUR + 0.4, zorder=4)
for reg, rv in (("Pacific", 0.50), ("Indian", 0.32)):
    ax.plot(3.0, rv, "D", color=PUB, ms=MS_PUB + 0.4, zorder=5); ax.text(3.07, rv, f"{reg}", fontsize=F.FS, va="center", ha="left", color=PUB)
ax.set_xlim(0, 4); ax.set_xticks([0, 1, 2, 3, 4]); ax.set_ylim(-0.3, 0.65); ax.set_yticks([-0.25, 0, 0.25, 0.5]); S.yticks_minus(ax)
ax.set_xlabel("Lag (years, SAM leading)"); ax.set_ylabel("r, EKE with SAM")
ax.legend(handles=ax.get_legend_handles_labels()[0] + [Line2D([], [], ls="none", marker="D", color=PUB, ms=MS_PUB + 0.4, label="Hogg et al. (2015), 1993–2012, lag 3 years"),
                   Line2D([], [], ls="none", marker="o", color=C["eke"], mfc=C["eke"], mew=0.8, ms=MS_OUR + 0.4, label="maximum, pointwise p < 0.05"),
                   Line2D([], [], ls="none", marker="o", color=C["eke"], mfc="white", mew=0.8, ms=MS_OUR + 0.4, label="maximum, p ≥ 0.05")],
          loc="lower left", bbox_to_anchor=(0.0, 1.0 + 2.0 / hd), ncol=3, handlelength=2.4, columnspacing=1.6, labelspacing=0.35, borderaxespad=0)
S.letter(fig, 2, yd_let, "d", "Lag of EKE behind the SAM index, 12-month running means, method of Hogg et al. (2015)")
F.check(fig); S.save(fig, FIGD, "FigureS4"); plt.close(fig)
print("S4 clipped interval ends:", NCLIP[0], "| d axes height (mm):", round(hd, 1))
for nm, gg in (("a", ga), ("b", gb), ("c", gc)):
    for l1, l2, pub, ent in gg:
        print(f"S4{nm} {l1} | published {pub.published_value:+.1f}", "| ours:", ", ".join(f"{r.our_series.split(',')[0]} {r.our_period} {r.our_value:+.1f} [{r.our_lo:+.1f}, {r.our_hi:+.1f}] p={r.p:.3f}" for _, _, r in ent))
print("S4d maxima (lag months, r, pointwise p):", mx)
