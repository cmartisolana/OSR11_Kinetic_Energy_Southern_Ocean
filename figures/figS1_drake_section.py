"""Supplement Fig. S1 (final layout, 10 Nov 2026): Drake Passage section along the LM Gould transect.

Same data and numbers as figures_old/figS1_drake_section.py (5 Nov 2026); layout and presentation only:
170 x 128 mm, a and b side by side, c full width; legends above the panels or in free areas; fronts drawn the same way in
every panel (dashed vertical at the mean SAF / PF level, bold label above the axis); c oriented south (left) to north
(right) as main-text Fig. 2d, with the SAF and PF zones and the five 0.05 m bins shaded as in Fig. 2d and the four bins
between the fronts marked.

  a  Mean cross-track velocity per 0.25 deg latitude bin, 248 transects Feb 2005-Dec 2019: in situ ADCP 70-118 m and
     Products 1 (solid) and 4 (dashed) sampled on them (outputs/drake/dk3/dk3_velocity_profile_{allsat,twosat}.csv).
  b  Cross-track velocity trend per 25 km bin (passage coordinates; outputs/drake/dk5/dk5_trends.csv, var u): in situ
     2005-2019 (diamonds, light bars = 95 % intervals, clipped at +-10 with a caret), Product 1 Oct 2005-Apr 2019 (dashed)
     and 2001-2025 (solid, shading = 95 % interval). Bin centres to latitude with outputs/drake/dk39/dk39_transport_bins.csv.
  c  Cross-track velocity trend in the 17 ADT bins of main-text Fig. 2d with fixed edges (circles) and with edges following
     the 1.8 cm per decade rise of the section-mean ADT (squares) (review check c04: c04_results.csv in $OSR11_REVIEW,
     written by scripts/03_drake/03w_adt_bin_drift.py). Zones (psi edges of Gutierrez-Villanueva et al. 2023 mapped to
     ADT with outputs/drake/dk5/dk5_stream_edges.csv): SAF zone psi -0.4 to -0.6, PF zone psi -1.0 to -1.25, between the
     fronts psi -0.6 to -1.0 (the four bins of Fig. 2d between the SAF and PF zones).
Fronts: latitude where the 2001-2025 mean ADT on the line (outputs/drake/dk24/dk24_section_fields_allsat.npz) crosses the
SAF and PF levels of Fig. 2 (outputs/drake/dk34/dk34_drake_fields.npz). Light blue: seasonally ice-covered part of the line
(bins with fewer than 300 valid months, outputs/drake/dk25/dk25_profile_trends.csv). Filled symbols: p < 0.05.
The Pearson r of the trend profiles is printed (the old script also wrote it to outputs/drake/dk5/dk5_trend_profile_r.csv;
this script writes nothing outside the figure directory).
Usage: python figures/figS1_drake_section.py [outputs dir] [figure dir]   (defaults: $OSR11_OUTPUTS, rendered/)
"""
import _sifig as F
import numpy as np, pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple
import matplotlib.ticker as mticker

S, Q, C, plt, OUT, FIGD = F.S, F.Q, F.C, F.plt, F.OUT, F.FIGD
DO = OUT / "drake"
REV = OUT / "c04"   # ADT-bin drift check (scripts/03_drake/03w_adt_bin_drift.py)
P3 = {k: pd.read_csv(DO / "dk3" / f"dk3_velocity_profile_{k}.csv") for k in ("allsat", "twosat")}
SMY = pd.read_csv(DO / "dk3" / "dk3_summary.csv", header=None, index_col=0)[1]
K5 = pd.read_csv(DO / "dk5" / "dk5_trends.csv"); T39 = pd.read_csv(DO / "dk39" / "dk39_transport_bins.csv"); lat_of = dict(zip(T39.bin, T39.lat))
V = pd.read_csv(DO / "dk25" / "dk25_profile_trends.csv"); ice_lat = V.latb[V.months < 300].max() + 0.125
F34 = np.load(DO / "dk34" / "dk34_drake_fields.npz"); Z = np.load(DO / "dk24" / "dk24_section_fields_allsat.npz")
LEV = {"SAF": float(F34["lev_saf"]), "PF": float(F34["lev_pf"])}
_a = np.nanmean(Z["adt"], 0); _la = Z["lat"]
EDG = pd.read_csv(DO / "dk5" / "dk5_stream_edges.csv"); psi2adt = dict(zip(EDG.psi.round(3), EDG.adt_edge))
ZONE = {"SAF": (psi2adt[-0.6], psi2adt[-0.4]), "PF": (psi2adt[-1.25], psi2adt[-1.0])}
BETWEEN = (psi2adt[-1.0], psi2adt[-0.6])


def _cross(lev):
    j = np.where(np.diff(np.sign(_a - lev)) != 0)[0][0]
    return _la[j] + (lev - _a[j]) * (_la[j + 1] - _la[j]) / (_a[j + 1] - _a[j])


FRONTS = {k: _cross(v) for k, v in LEV.items()}
XL = (-63.6, -54.9); XT = [-63, -61, -59, -57, -55]
CI, CA = C["obs"], C["vel"]
CI_LIGHT = "#F3CF8A"                  # light in situ orange for the bars
FLS = (0, (4, 2))                     # front line style (all panels)
SHADE_Z, SHADE_N = "#EFEFEF", "#DCDCDC"   # front zones (as Fig. 2d), 0.05 m bins (as Fig. 2d)


def fronts(ax, xs, z=1.5):
    for fn, x in xs.items():
        ax.axvline(x, color=C["front"], lw=0.6, ls=FLS, zorder=z)
        ax.text(x, 1.012, fn, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=F.FS, fontweight="bold")


def latax(ax):
    ax.axvspan(XL[0], ice_lat, color=C["ice_flat"], lw=0, zorder=0); ax.set_xlim(*XL); ax.set_xticks(XT)
    ax.set_xticklabels([S.latfmt(v) for v in XT]); S.zero(ax); ax.set_xlabel("Latitude")
    fronts(ax, FRONTS); S.yticks_minus(ax)


def k5(src, per):
    q = K5[(K5.coord == "passage") & (K5.src == src) & (K5["var"] == "u") & (K5.period == per)].sort_values("bin").copy()
    q["lat"] = q.bin.map(lat_of); return q


H = 128.0
fig = S.figure(S.W2, H)
TOP_Y, PH = 74.0, 38.0                 # row 1 axes bottom and height (mm)
# a ---------------------------------------------------------------- mean cross-track velocity
ax = S.axmm(fig, 14, TOP_Y, 63, PH); latax(ax)
p1 = P3["allsat"]; p4 = P3["twosat"]
ax.fill_between(p1.latb, 0, p1.ualt_mean, color=CA, alpha=0.13, lw=0, zorder=2)
ax.plot(p1.latb, p1.ualt_mean, color=CA, lw=0.9, ls=Q.PROD["P1"]["ls"], zorder=3)
ax.plot(p4.latb, p4.ualt_mean, color=CA, lw=0.9, ls=Q.PROD["P4"]["ls"], zorder=3)
ax.plot(p1.latb, p1.uin_mean, color=CI, lw=0.9, marker="D", ms=2.0, mew=0, zorder=4)
ax.set_ylim(-5, 48); ax.set_yticks([0, 10, 20, 30, 40]); ax.set_ylabel("Cross-track velocity (cm s$^{-1}$)")
S.letter(fig, 2, H - 5, "a", "Mean velocity, 248 transects, Feb 2005–Dec 2019")
ax.legend(handles=[Line2D([], [], color=CI, marker="D", ms=2.0, mew=0, lw=0.9, label="in situ ADCP, 70–118 m"),
                   Line2D([], [], color=CA, lw=0.9, label="altimetry, Product 1"),
                   Patch(facecolor=C["ice_flat"], edgecolor="none", label="seasonal ice"),
                   Line2D([], [], color=CA, lw=0.9, ls=Q.PROD["P4"]["ls"], label="altimetry, Product 4")],
          loc="lower left", bbox_to_anchor=(-0.02, 1.07), ncol=2, handlelength=2.0, columnspacing=1.0, labelspacing=0.3,
          borderaxespad=0)
# b ---------------------------------------------------------------- velocity trends per 25 km bin
YB0, YB = -10.0, 12.0
ax = S.axmm(fig, 103, TOP_Y, 63, PH); latax(ax)
qi, qa, ql = k5("insitu", "2005-2019"), k5("allsat", "2005-2019"), k5("allsat", "2001-2025")
# in situ: light, wide bars behind (clipped at the axis, caret at clipped ends), diamonds on top
lo_i, hi_i = 10 * qi.ci_lo.values, 10 * qi.ci_hi.values
ax.vlines(qi.lat, np.maximum(lo_i, YB0), np.minimum(hi_i, YB), color=CI_LIGHT, lw=1.6, zorder=1.2)
for x_, l_, h_ in zip(qi.lat, lo_i, hi_i):
    if h_ > YB: ax.plot(x_, YB - 0.25, "^", ms=2.2, color=CI_LIGHT, mew=0, zorder=1.3, clip_on=False)
    if l_ < YB0: ax.plot(x_, YB0 + 0.25, "v", ms=2.2, color=CI_LIGHT, mew=0, zorder=1.3, clip_on=False)
ax.fill_between(ql.lat, 10 * ql.ci_lo, 10 * ql.ci_hi, color=CA, alpha=0.16, lw=0, zorder=1.4)
S.profile(ax, ql.lat, 10 * ql.slope, ql.p, CA, ms=2.0, z=3)
S.profile(ax, qa.lat, 10 * qa.slope, qa.p, CA, ms=2.0, z=3, ls=(0, (3, 1.6)), lw=0.7)
S.profile(ax, qi.lat, 10 * qi.slope, qi.p, CI, mk="D", ms=2.6, z=5, lw=0)
ax.set_ylim(YB0, YB); ax.set_yticks([-10, -5, 0, 5, 10]); S.yticks_minus(ax)
ax.set_ylabel("Velocity trend (cm s$^{-1}$ per decade)")
S.letter(fig, 90, H - 5, "b", "Velocity trend per 25 km bin")
ax.legend(handles=[(Line2D([], [], color=CI_LIGHT, marker="|", ms=7, mew=1.6, ls=""), Line2D([], [], color=CI, marker="D", ms=2.6, ls="")),
                   (Patch(facecolor=CA, alpha=0.16, edgecolor="none"), Line2D([], [], color=CA, marker="o", ms=2.0, lw=0.8)),
                   Line2D([], [], color=CA, marker="o", ms=2.0, lw=0.7, ls=(0, (3, 1.6))),
                   Line2D([], [], color="#555555", marker="o", mfc="white", ls="", ms=2.4)],
          labels=["in situ, 2005–2019", "Product 1, 2001–2025", "Product 1, Oct 2005–Apr 2019", Q.SIG_NS],
          handler_map={tuple: HandlerTuple(ndivide=1)},
          loc="lower left", bbox_to_anchor=(-0.19, 1.07), ncol=2, handlelength=1.8, columnspacing=0.9, labelspacing=0.3,
          borderaxespad=0)
# c ---------------------------------------------------------------- ADT bins: fixed edges vs edges following the ADT rise
ax = S.axmm(fig, 14, 13, 152, 38); S.zero(ax)
D = pd.read_csv(REV / "c04_results.csv"); v0 = D[D.variant.str.startswith("fixed")].copy(); v1 = D[D.variant.str.contains("secular")]
# fixed edges: the values plotted in main-text Fig. 2d (dk5_trends.csv, stream coordinates, Product 1, 2001-2025), so that the
# circles are identical to Fig. 2d; the drifting-edge squares are from the review check c04 (10 Nov 2026)
_k5 = K5[(K5.coord == "stream") & (K5.src == "allsat") & (K5["var"] == "u") & (K5.period == "2001-2025")].copy(); _k5.index = _k5["bin"].round(3).values
_b = v0.bin_psic.round(3).values
v0["slope_dec"] = 10 * _k5.loc[_b, "slope"].values; v0["lo"] = 10 * _k5.loc[_b, "ci_lo"].values
v0["hi"] = 10 * _k5.loc[_b, "ci_hi"].values; v0["p"] = _k5.loc[_b, "p"].values
xc = 0.5 * (v0.adt_lo.values + v0.adt_hi.values)
narrow = (v0.adt_hi.values - v0.adt_lo.values) < 0.07
XLc = (-1.2, 0.45)
ax.axvspan(v0.adt_lo.values[narrow].min(), v0.adt_hi.values[narrow].max(), color=SHADE_N, lw=0, zorder=0)
for zn, (z0, z1) in ZONE.items():
    ax.axvspan(z0, z1, color=SHADE_Z, lw=0, zorder=0)
ax.text(v0.adt_lo.values[narrow].mean() * 0 + 0.5 * (v0.adt_lo.values[narrow].min() + v0.adt_hi.values[narrow].max()), -7.6,
        "0.05 m bins", ha="center", va="bottom", fontsize=F.FS, color=F.GREY_TXT)
fronts(ax, {k: v for k, v in LEV.items()})
# bracket over the bins between the fronts
yb = 8.6
ax.plot([BETWEEN[0] + 0.006, BETWEEN[0] + 0.006, BETWEEN[1] - 0.006, BETWEEN[1] - 0.006], [yb - 0.5, yb, yb, yb - 0.5], color="#222222", lw=0.6, zorder=6)
ax.text(0.5 * sum(BETWEEN), yb + 0.25, "between the fronts", ha="center", va="bottom", fontsize=F.FS, zorder=6)
for v, mk, c, dx in ((v0, "o", C["vel"], -0.011), (v1, "s", C["grey"], 0.011)):
    for xi, r in zip(xc + dx, v.itertuples()):
        ax.plot([xi, xi], [r.lo, r.hi], color=c, lw=0.7, zorder=2)
        ax.plot(xi, r.slope_dec, mk, color=c, mfc=c if S.sigp(r.p) else "white", mec=c, mew=0.7, ms=3.0, zorder=3)
ax.set_xlim(*XLc); ax.set_xticks(np.arange(-1.2, 0.41, 0.2)); ax.set_ylim(-8, 10); ax.set_yticks([-5, 0, 5, 10])
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: S.um(f"{round(v, 2) + 0:g}"))); S.yticks_minus(ax)
ax.set_xlabel("Mean ADT at the bin centre (m; south ← → north)"); ax.set_ylabel("Velocity trend (cm s$^{-1}$ per decade)")
ax.legend(handles=[Line2D([], [], marker="o", color=C["vel"], lw=0.7, ms=3.0, label="fixed edges (as Fig. 2d)"),
                   Line2D([], [], marker="s", color=C["grey"], lw=0.7, ms=3.0, label="edges follow the ADT rise"),
                   Line2D([], [], color="#555555", marker="o", mfc="white", ls="", ms=3.0, label=Q.SIG_NS)],
          loc="center right", bbox_to_anchor=(1.0, (38 + 7.6) / 38), ncol=3, handlelength=1.6, columnspacing=1.2, borderaxespad=0)
S.letter(fig, 2, 13 + 38 + 6.6, "c", "Velocity trend in ADT bins, fixed and drifting edges")
F.check(fig); S.save(fig, FIGD, "FigureS1"); plt.close(fig)

print("S1a r(mean profiles):", {k: round(float(SMY[f"{k}|u_mean_profile_r"]), 3) for k in ("allsat", "twosat")},
      "| r bins:", round(float(SMY["allsat|u_r_bins_all"]), 2), "n", int(SMY["allsat|u_n_bins_all"]), "rms", round(float(SMY["allsat|u_rms_diff_cm_s"]), 1))
print("S1b in situ bins", len(qi), "sig", int((qi.p.round(3) < 0.05).sum()), "; altimetry 2001-2025 sig", int((ql.p.round(3) < 0.05).sum()),
      "; r(trend profiles, in situ vs Product 1 Oct 2005-Apr 2019) =", round(np.corrcoef(qi.slope, qa.slope)[0, 1], 2),
      "; vs 2001-2025 =", round(np.corrcoef(qi.slope, ql.slope)[0, 1], 2),
      "; in situ bars clipped at -10/+12:", int((hi_i > YB).sum() + (lo_i < YB0).sum()), "ends")
print("S1 fronts (lat):", {k: round(v, 2) for k, v in FRONTS.items()}, "levels (m):", {k: round(v, 3) for k, v in LEV.items()}, "ice_lat", ice_lat)
print("S1c zones (ADT m): SAF", tuple(round(z, 3) for z in ZONE["SAF"]), "PF", tuple(round(z, 3) for z in ZONE["PF"]),
      "between the fronts", tuple(round(z, 3) for z in BETWEEN))
inb = (xc > BETWEEN[0]) & (xc < BETWEEN[1])
for nm, v in (("fixed", v0), ("drifting", v1)):
    print(f"S1c {nm} edges, bins between the fronts (centre ADT: slope [lo, hi], p):",
          "; ".join(f"{x:+.2f}: {r.slope_dec:+.1f} [{r.lo:+.1f}, {r.hi:+.1f}], p={r.p:.3f}" for x, r in zip(xc[inb], v[inb].itertuples())))
