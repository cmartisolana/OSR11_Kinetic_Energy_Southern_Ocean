"""SI Figure S4: LM Gould in situ record (Gutierrez-Villanueva et al., 2023) and altimetry (Product 1) on the same transect, in fixed
latitude (25 km) bins; in situ Oct 2005-Apr 2019, altimetry on the same months unless stated.
a  mean cross-track velocity per 0.25 deg bin: in situ ADCP 70-118 m and altimetric surface velocity on the same transects (dk3).
b  mean eastward transport per 25 km bin: in situ 0-760 m (Sv, left) and altimetric surface transport per unit depth
   (10^3 m2 s-1, right; Oct 2005-Apr 2019 mean, computed here from dk24_section_fields_allsat.npz with the dk39 bin definition;
   the right axis is scaled so that the two totals coincide).
c  annual EKE north of 58.5 S, 2005-2019: in situ 0-780 m depth-integrated and altimetric surface EKE on the same transects, each
   in % of its own 2005-2019 mean (dk3_eke_annual_north_allsat.csv).
d  velocity trend per bin, Oct 2005-Apr 2019, in situ and altimetry (dk5 passage).
e  transport trend per bin, same period: in situ (Sv per decade, left) and altimetric surface transport per unit depth
   (10^3 m2 s-1 per decade, right; dk5, dk39).
f  EKE trend per SSH streamline pair, Oct 2005-Apr 2019, % per decade: in situ 0-780 m (dk5 stream) and altimetry (dk24).
Trends: calibrated Theil-Sen, 95 % (bars in situ, shading altimetry); filled p < 0.05. Correlation values are printed to stdout and
given in the caption (no statistics inside the panels).

Publication repository: ported from figures_final/scripts/figS4_drake_validation.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/drake/{dk3,dk5,dk24,dk25,dk34,dk39}
Outputs : rendered/{png,pdf}/FigureS4
Run     : python figures/figS4_drake_validation.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd
from osr11 import style as S, paths as P  # noqa: E402
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else P.FIG
plt = S.plt; C = S.C
DR = DO = P.OUT / "drake"
P3 = pd.read_csv(DO / "dk3" / "dk3_velocity_profile_allsat.csv"); SMY = pd.read_csv(DO / "dk3" / "dk3_summary.csv", header=None, index_col=0)[1]
EA = pd.read_csv(DO / "dk3" / "dk3_eke_annual_north_allsat.csv")
K5 = pd.read_csv(DO / "dk5" / "dk5_trends.csv"); T39 = pd.read_csv(DO / "dk39" / "dk39_transport_bins.csv"); lat_of = dict(zip(T39.bin, T39.lat))
A24 = pd.read_csv(DO / "dk24" / "dk24_stream_trends.csv"); adt_of = dict(zip(A24.psic.round(3), A24.adt_c))
V = pd.read_csv(DO / "dk25" / "dk25_profile_trends.csv"); ice_lat = V.latb[V.months < 300].max() + 0.125
F34 = np.load(DO / "dk34" / "dk34_drake_fields.npz"); Z = np.load(DO / "dk24" / "dk24_section_fields_allsat.npz")
_a = np.nanmean(Z["adt"], 0); _la = Z["lat"]

# --- Oct 2005-Apr 2019 mean surface transport per unit depth per 25 km bin (dk39 definitions; the missing 'mean_0519' column)
lo, la, dsx = Z["lon"], Z["lat"], Z["dsx"]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
un = Z["u"] * nx + Z["v"] * ny; ice = Z["ice"] > 0.10
o = np.argsort(-la); dist = np.zeros(la.size); dist[o] = np.r_[0, np.cumsum(dsx[o][:-1] + dsx[o][1:]) / 2] / 1000
xm = 2001 + (np.arange(300) + 0.5) / 12; ov = (xm >= 2005 + 9 / 12) & (xm <= 2019 + 4 / 12)
m0519, mall = {}, {}
for b0 in np.arange(0, 850, 25):
    s = (dist >= b0) & (dist < b0 + 25)
    if not s.any(): continue
    q = (un[:, s] * dsx[s]).sum(1); q[ice[:, s].any(1)] = np.nan
    m0519[b0 + 12.5] = np.nanmean(q[ov]); mall[b0 + 12.5] = np.nanmean(q)
T39["mean_0519"] = T39.bin.map(m0519)
chk = np.nanmax(np.abs(T39.bin.map(mall) - T39["mean"]) / np.abs(T39["mean"]))
print(f"check: 2001-2025 bin means reproduce dk39 'mean' to {100 * chk:.2e} %")


def _cross(lev):
    j = np.where(np.diff(np.sign(_a - lev)) != 0)[0][0]; return _la[j] + (lev - _a[j]) * (_la[j + 1] - _la[j]) / (_a[j + 1] - _a[j])
FRONTS = [("SAF", _cross(float(F34["lev_saf"]))), ("PF", _cross(float(F34["lev_pf"])))]
XL = (-63.6, -54.9); XT = [-63, -61, -59, -57, -55]
CI, CA = C["obs"], C["eke"]


def latax(ax):
    ax.axvspan(XL[0], ice_lat, color=C["ice_flat"], lw=0, zorder=0); ax.set_xlim(*XL); ax.set_xticks(XT)
    ax.set_xticklabels([S.latfmt(v) for v in XT]); ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=1)
    for fn, fl in FRONTS:
        ax.axvline(fl, color=C["front"], lw=0.6, ls=(0, (4, 2)), zorder=1.5)
        ax.text(fl, 1.0, fn, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=S.FS_ANN, fontweight="bold")
    S.yticks_minus(ax)


def k5(coord, src, var, per):
    return K5[(K5.coord == coord) & (K5.src == src) & (K5["var"] == var) & (K5.period == per)].sort_values("bin").copy()


def right(ax2, col):
    ax2.spines["right"].set_visible(True); ax2.spines["right"].set_linewidth(0.5); ax2.spines["right"].set_color(col)
    ax2.tick_params(axis="y", colors=col)


H = 126.0; PW, PH = 36.0, 38.0; XS = (14, 70, 126); YS = (73, 15)
fig = S.figure(S.W2, H)
axs = {l: S.axmm(fig, XS[i % 3], YS[i // 3], PW, PH) for i, l in enumerate("abcdef")}
for i, (l, t) in enumerate(zip("abcdef", ("Mean cross-track velocity", "Mean eastward transport", "Annual EKE north of 58.5°S",
                                            "Velocity trend", "Transport trend", "EKE trend per ADT bin"))):
    S.letter(fig, XS[i % 3] - 11, YS[i // 3] + PH + 6.5, l, t)
# a
ax = axs["a"]; latax(ax)
ax.fill_between(P3.latb, 0, P3.ualt_mean, color=CA, alpha=0.2, lw=0, zorder=2); ax.plot(P3.latb, P3.ualt_mean, color=CA, lw=0.9, zorder=3)
ax.plot(P3.latb, P3.uin_mean, color=CI, lw=0.9, marker="D", ms=1.8, zorder=4)
ax.set_ylim(-5, 48); ax.set_ylabel("cm s$^{-1}$")
# b
qi = k5("passage", "insitu", "transport", "2005-2019"); qi["lat"] = qi.bin.map(lat_of); FAC = (T39["mean_0519"].sum() / 1e3) / qi["mean"].sum()
ax = axs["b"]; latax(ax)
ax.fill_between(qi.lat, 0, qi["mean"], color=CI, alpha=0.25, lw=0, zorder=2); ax.plot(qi.lat, qi["mean"], color=CI, lw=0.9, zorder=3)
ax.set_ylim(0, 8); ax.set_ylabel("In situ (Sv)", color=CI); ax.tick_params(axis="y", colors=CI)
ax2 = ax.twinx(); ax2.plot(T39.lat, T39["mean_0519"] / 1e3, color=CA, lw=0.9, zorder=3); ax2.set_ylim(0, 8 * FAC)
ax2.set_ylabel("Altimetry (10$^3$ m$^2$ s$^{-1}$)", color=CA); right(ax2, CA)
_m = qi.merge(T39[["bin", "mean_0519"]], on="bin"); r_tr = np.corrcoef(_m["mean"], _m["mean_0519"])[0, 1]
# c  (one axis: annual EKE of both records in % of their own 2005-2019 mean)
ax = axs["c"]; yy = EA.year + 0.5; pc_ = lambda v: 100 * (v / v.mean() - 1)
ax.plot(yy, pc_(EA.eke_in), "-D", color=CI, lw=0.9, ms=1.8); ax.plot(yy, pc_(EA.eke_alt), "-o", color=CA, lw=0.9, ms=1.8)
ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=0); ax.set_ylabel("% of 2005–2019 mean"); S.yticks_minus(ax)
ax.set_xlim(2004.8, 2020.2); ax.set_xticks([2005, 2010, 2015, 2020])
# d
ax = axs["d"]; latax(ax)
q = k5("passage", "insitu", "u", "2005-2019"); q["lat"] = q.bin.map(lat_of)
ax.vlines(q.lat, 10 * q.ci_lo, 10 * q.ci_hi, color=CI, lw=0.5, alpha=0.7, zorder=2); S.profile(ax, q.lat, 10 * q.slope, q.p, CI, mk="D", ms=2.0, lw=0)
q = k5("passage", "allsat", "u", "2005-2019"); q["lat"] = q.bin.map(lat_of)
ax.fill_between(q.lat, 10 * q.ci_lo, 10 * q.ci_hi, color=CA, alpha=0.18, lw=0, zorder=1); S.profile(ax, q.lat, 10 * q.slope, q.p, CA, ms=1.8)
ax.set_ylim(-9, 9); ax.set_ylabel("cm s$^{-1}$ per decade")
# e
ax = axs["e"]; latax(ax)
ax.vlines(qi.lat, 10 * qi.ci_lo, 10 * qi.ci_hi, color=CI, lw=0.5, alpha=0.7, zorder=2); S.profile(ax, qi.lat, 10 * qi.slope, qi.p, CI, mk="D", ms=2.0, lw=0)
ax.set_ylim(-2.2, 2.2); ax.set_ylabel("In situ (Sv per decade)", color=CI); ax.tick_params(axis="y", colors=CI)
ax2 = ax.twinx(); ax2.set_ylim(-2.2 * FAC, 2.2 * FAC); right(ax2, CA); S.yticks_minus(ax2)
ax2.fill_between(T39.lat, T39.ci_lo_0519 / 1e3, T39.ci_hi_0519 / 1e3, color=CA, alpha=0.18, lw=0, zorder=1)
S.profile(ax2, T39.lat, T39.slope_0519 / 1e3, T39.p_0519, CA, ms=1.8)
ax2.set_ylabel("Altimetry (10$^3$ m$^2$ s$^{-1}$ per decade)", color=CA)
# f
ax = axs["f"]; ED = pd.read_csv(DO / "dk5" / "dk5_stream_edges.csv"); ed = dict(zip(ED.psi.round(3), ED.adt_edge))
for zn, lo_, hi_ in (("SAF", -0.6, -0.4), ("PF", -1.25, -1.0)):
    ax.axvspan(ed[lo_], ed[hi_], color="#EDEDED", lw=0, zorder=0)
    ax.text((ed[lo_] + ed[hi_]) / 2, 1.0, zn, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=S.FS_ANN, fontweight="bold")
ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=1)
q = k5("stream", "insitu", "eke", "2005-2019"); q["adt"] = q.bin.round(3).map(adt_of); q = q.dropna(subset=["adt"]).sort_values("adt")
ax.vlines(q.adt, 1000 * q.ci_lo / q["mean"], 1000 * q.ci_hi / q["mean"], color=CI, lw=0.5, alpha=0.7, zorder=2)
S.profile(ax, q.adt, 1000 * q.slope / q["mean"], q.p, CI, mk="D", ms=2.0, lw=0)
q = A24[(A24["var"] == "eke") & (A24.period == "2005-2019")].sort_values("adt_c")
ax.fill_between(q.adt_c, q.pct_lo, q.pct_hi, color=CA, alpha=0.18, lw=0, zorder=1); S.profile(ax, q.adt_c, q.pct, q.p, CA, ms=1.8)
ax.set_xlim(-1.16, 0.44); ax.set_ylim(-80, 120); ax.set_ylabel("% per decade"); ax.set_xlabel("Mean ADT at the bin centre (m)")
S.yticks_minus(ax); S.xticks_minus(ax)
for l in "abde": axs[l].set_xlabel("Latitude")
axs["c"].set_xlabel("Year")
fig.legend(handles=[Line2D([], [], color=CI, marker="D", ms=2.0, lw=0.9, label="in situ (Gutierrez-Villanueva et al., 2023)"),
                    Line2D([], [], color=CA, marker="o", ms=1.8, lw=0.9, label="altimetry, Product 1"),
                    Line2D([], [], color="#555555", marker="o", mfc="white", ls="", ms=2.2, label="open: p ≥ 0.05"),
                    Patch(facecolor=C["ice_flat"], label="seasonal sea ice")], loc="center", ncol=4, bbox_to_anchor=(0.5, 3.0 / H))
S.save(fig, FIGD, "FigureS4"); plt.close(fig)
print("S4 correlations: mean velocity profile r = %.2f; transect bins r = %.2f; rms diff = %.1f cm/s; annual EKE r = %.2f; transects EKE r = %.2f; mean transport bins r = %.2f; FAC = %.3f"
      % (float(SMY["allsat|u_mean_profile_r"]), float(SMY["allsat|u_r_bins_all"]), float(SMY["allsat|u_rms_diff_cm_s"]),
         float(SMY["allsat|eke_r_annual_north"]), float(SMY["allsat|eke_r_transects_north"]), r_tr, FAC))
