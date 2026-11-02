"""Figure S (v3.57): agreement between the LM Gould in situ record and altimetry at the surface across Drake Passage, to
support the altimetric conclusions of Fig. 2. In situ: ADCP cross-track velocity 70-118 m, 0-760 m transport, 0-780 m
depth-integrated EKE (Gutierrez-Villanueva et al. 2023), Oct 2005-Apr 2019; altimetry on the same transects / months.
a  mean cross-track velocity per 0.25 deg bin (dk3); b mean transport per 25-km bin: in situ 0-760 m vs altimetric surface
   transport per unit depth, right axis scaled so the totals coincide (dk5, dk39); c annual EKE north of 58.5 S, in situ vs
   altimetry (dk3); d velocity trend per 25-km bin, Oct 2005-Apr 2019 (dk5); e transport trend per 25-km bin (dk5, dk39);
f  EKE trend per SSH streamline pair, Oct 2005-Apr 2019, % per decade (dk5, dk24). Trends: calibrated Theil-Sen, 95 %.
Usage: python figS_drake_validation_v357.py <revision dir> <drake dir> <work out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, DR, WO, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from matplotlib.lines import Line2D
plt = S.plt; C = S.C
sigp = lambda p: np.isfinite(p) and round(float(p), 3) < 0.05
P3 = pd.read_csv(DR / "out" / "dk3" / "dk3_velocity_profile_allsat.csv"); SM = pd.read_csv(DR / "out" / "dk3" / "dk3_summary.csv", header=None, index_col=0)[1]
EA = pd.read_csv(DR / "out" / "dk3" / "dk3_eke_annual_north_allsat.csv")
K5 = pd.read_csv(DR / "out" / "dk5" / "dk5_trends.csv"); T39 = pd.read_csv(WO / "dk39" / "dk39_transport_bins.csv"); lat_of = dict(zip(T39.bin, T39.lat))
A24 = pd.read_csv(DR / "out" / "dk24" / "dk24_stream_trends.csv"); adt_of = dict(zip(A24.psic.round(3), A24.adt_c))
V = pd.read_csv(DR / "out" / "dk25" / "dk25_profile_trends.csv"); ice_lat = V.latb[V.months < 300].max() + 0.125
F34 = np.load(WO / "dk34" / "dk34_drake_fields.npz"); _sec = np.load(DR / "out" / "dk24" / "dk24_section_fields_allsat.npz"); _a = np.nanmean(_sec["adt"], 0); _la = _sec["lat"]
def _cross(lev):
    j = np.where(np.diff(np.sign(_a - lev)) != 0)[0][0]; return _la[j] + (lev - _a[j]) * (_la[j + 1] - _la[j]) / (_a[j + 1] - _a[j])
FRONTS = [("SAF", _cross(float(F34["lev_saf"]))), ("PF", _cross(float(F34["lev_pf"])))]
XL = (-63.6, -54.9); XT = [-63, -61, -59, -57, -55]; XTL = ["63°S", "61°S", "59°S", "57°S", "55°S"]
CI, CA = C["obs"], C["eke"]
def latax(ax, lab, ttl, fronts=True):
    ax.axvspan(XL[0], ice_lat, color="#DCEBF7", lw=0, zorder=0); ax.set_xlim(*XL); ax.set_xticks(XT); ax.set_xticklabels(XTL); ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=1)
    if fronts:
        for fn, fl in FRONTS: ax.axvline(fl, color="#1A1A1A", lw=0.7, ls=(0, (4, 2)), zorder=1.5); ax.text(fl, 1.0, fn, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.2, fontweight="bold")
    S.plabel(ax, lab, ttl, dy=9)
def prof(ax, x, y, p, col, ls="-", ms=2.2, lw=0.9, z=3, mk="o"):
    x, y, p = np.asarray(x), np.asarray(y), np.asarray(p); ax.plot(x, y, color=col, ls=ls, lw=lw, zorder=z); s = np.array([sigp(v) for v in p])
    ax.plot(x[s], y[s], mk, color=col, ms=ms, zorder=z + 1); ax.plot(x[~s], y[~s], mk, color=col, mfc="white", ms=ms, mew=0.6, zorder=z + 1)
def k5(coord, src, var, per): q = K5[(K5.coord == coord) & (K5.src == src) & (K5["var"] == var) & (K5.period == per)].sort_values("bin").copy(); return q
fig = plt.figure(figsize=(S.W2, 0.62 * S.W2)); g = fig.add_gridspec(2, 3, wspace=0.55, hspace=0.5)
# a mean velocity
ax = fig.add_subplot(g[0, 0]); latax(ax, "a", "Mean cross-track velocity, same transects")
ax.fill_between(P3.latb, 0, P3.ualt_mean, color=CA, alpha=0.2, lw=0, zorder=2); ax.plot(P3.latb, P3.ualt_mean, color=CA, lw=0.9, zorder=3, label="altimetry, surface")
ax.plot(P3.latb, P3.uin_mean, color=CI, lw=0.9, marker="D", ms=1.6, zorder=4, label="in situ ADCP, 70–118 m")
ax.set_ylim(-5, 48); ax.set_ylabel("cm s$^{-1}$"); ax.legend(loc="upper left", fontsize=4.6)
ax.text(0.98, 0.42, f"mean profiles: r = {float(SM['allsat|u_mean_profile_r']):.2f}\ntransect bins: r = {float(SM['allsat|u_r_bins_all']):.2f}, rms diff {float(SM['allsat|u_rms_diff_cm_s']):.0f} cm s$^{{-1}}$",
        transform=ax.transAxes, ha="right", va="top", fontsize=4.6, color="#333333", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
# b mean transport
qi = k5("passage", "insitu", "transport", "2005-2019"); qi["lat"] = qi.bin.map(lat_of); FAC = (T39["mean_0519"].sum() / 1e3) / qi["mean"].sum()
ax = fig.add_subplot(g[0, 1]); latax(ax, "b", "Mean eastward transport per 25-km bin")
ax.fill_between(qi.lat, 0, qi["mean"], color=CI, alpha=0.25, lw=0, zorder=2); ax.plot(qi.lat, qi["mean"], color=CI, lw=0.9, zorder=3, label="in situ, 0–760 m")
ax.set_ylim(0, 8); ax.set_ylabel("Sv", color=CI); ax.tick_params(axis="y", colors=CI)
ax2 = ax.twinx(); ax2.plot(T39.lat, T39["mean_0519"] / 1e3, color=CA, lw=0.9, zorder=3, label="altimetry, surface"); ax2.set_ylim(0, 8 * FAC)
ax2.set_ylabel("10$^3$ m$^2$ s$^{-1}$", color=CA, fontsize=5.5); ax2.tick_params(labelsize=5, colors=CA); ax2.spines["right"].set_visible(True); ax2.spines["right"].set_linewidth(0.5)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=4.6)
_m = qi.merge(T39[["bin", "mean_0519"]], on="bin"); ax.text(0.55, 0.96, f"r = {np.corrcoef(_m['mean'], _m['mean_0519'])[0, 1]:.2f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=4.6, color="#333333")
# c annual EKE series
ax = fig.add_subplot(g[0, 2]); yy = EA.year + 0.5
ax.plot(yy, EA.eke_in, "-D", color=CI, lw=0.9, ms=2.0, label="in situ, 0–780 m (left)"); ax.set_ylabel("in situ EKE (J m$^{-1}$)", color=CI); ax.tick_params(axis="y", colors=CI)
ax2 = ax.twinx(); ax2.plot(yy, EA.eke_alt, "-o", color=CA, lw=0.9, ms=2.0, label="altimetry, surface (right)"); ax2.set_ylabel("altimetric EKE (cm$^2$ s$^{-2}$)", color=CA, fontsize=5.5)
ax2.tick_params(labelsize=5, colors=CA); ax2.spines["right"].set_visible(True); ax2.spines["right"].set_linewidth(0.5)
ax.set_xlim(2004.8, 2020.2); ax.set_xticks([2005, 2010, 2015, 2020])
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=4.6)
ax.text(0.98, 0.04, f"annual: r = {float(SM['allsat|eke_r_annual_north']):.2f}; transects: r = {float(SM['allsat|eke_r_transects_north']):.2f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=4.6, color="#333333")
S.plabel(ax, "c", "Annual EKE north of 58.5°S, 2005–2019", dy=4)
# d velocity trend 2005-2019
ax = fig.add_subplot(g[1, 0]); latax(ax, "d", "Velocity trend, Oct 2005–Apr 2019")
q = k5("passage", "insitu", "u", "2005-2019"); q["lat"] = q.bin.map(lat_of); ax.vlines(q.lat, 10 * q.ci_lo, 10 * q.ci_hi, color=CI, lw=0.5, alpha=0.7, zorder=2); prof(ax, q.lat, 10 * q.slope, q.p, CI, mk="D", ms=2.0, lw=0)
q = k5("passage", "allsat", "u", "2005-2019"); q["lat"] = q.bin.map(lat_of); ax.fill_between(q.lat, 10 * q.ci_lo, 10 * q.ci_hi, color=CA, alpha=0.18, lw=0, zorder=1); prof(ax, q.lat, 10 * q.slope, q.p, CA, ms=1.8)
ax.set_ylim(-9, 9); ax.set_ylabel("cm s$^{-1}$ per decade")
ax.legend(handles=[Line2D([], [], color=CI, marker="D", ls="", ms=2.0, label="in situ (bars: 95 %)"), Line2D([], [], color=CA, marker="o", ms=1.8, label="altimetry (shading: 95 %)"), Line2D([], [], color="#555555", marker="o", mfc="white", ls="", ms=2, label="open: p ≥ 0.05")], loc="upper left", fontsize=4.4)
# e transport trend 2005-2019
ax = fig.add_subplot(g[1, 1]); latax(ax, "e", "Transport trend per bin, Oct 2005–Apr 2019")
ax.vlines(qi.lat, 10 * qi.ci_lo, 10 * qi.ci_hi, color=CI, lw=0.5, alpha=0.7, zorder=2); prof(ax, qi.lat, 10 * qi.slope, qi.p, CI, mk="D", ms=2.0, lw=0)
ax.set_ylim(-0.9, 0.9); ax.set_ylabel("Sv per decade", color=CI); ax.tick_params(axis="y", colors=CI)
ax2 = ax.twinx(); ax2.set_ylim(-0.9 * FAC, 0.9 * FAC); ax2.spines["right"].set_visible(True); ax2.spines["right"].set_linewidth(0.5)
ax2.fill_between(T39.lat, T39.ci_lo_0519 / 1e3, T39.ci_hi_0519 / 1e3, color=CA, alpha=0.18, lw=0, zorder=1); prof(ax2, T39.lat, T39.slope_0519 / 1e3, T39.p_0519, CA, ms=1.8)
ax2.set_ylabel("10$^3$ m$^2$ s$^{-1}$ dec$^{-1}$", color=CA, fontsize=5.5); ax2.tick_params(labelsize=5, colors=CA)
# f EKE trend per streamline pair 2005-2019
ax = fig.add_subplot(g[1, 2])
for zn, lo_, hi_ in (("SAF", -0.6, -0.4), ("PF", -1.25, -1.0)):
    ED = pd.read_csv(DR / "out" / "dk5" / "dk5_stream_edges.csv"); ed = dict(zip(ED.psi.round(3), ED.adt_edge)); ax.axvspan(ed[lo_], ed[hi_], color="#EDEDED", lw=0, zorder=0)
    ax.text((ed[lo_] + ed[hi_]) / 2, 1.0, zn, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.2, fontweight="bold")
ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=1)
q = k5("stream", "insitu", "eke", "2005-2019"); q["adt"] = q.bin.round(3).map(adt_of); q = q.dropna(subset=["adt"]).sort_values("adt")
pc = lambda d: (1000 * d.slope / d["mean"], 1000 * d.ci_lo / d["mean"], 1000 * d.ci_hi / d["mean"])
y, lo_, hi_ = pc(q); ax.vlines(q.adt, lo_, hi_, color=CI, lw=0.5, alpha=0.7, zorder=2); prof(ax, q.adt, y, q.p, CI, mk="D", ms=2.0, lw=0)
q = A24[(A24["var"] == "eke") & (A24.period == "2005-2019")].sort_values("adt_c"); ax.fill_between(q.adt_c, q.pct_lo, q.pct_hi, color=CA, alpha=0.18, lw=0, zorder=1); prof(ax, q.adt_c, q.pct, q.p, CA, ms=1.8)
ax.set_xlim(-1.16, 0.44); ax.set_ylim(-80, 120); ax.set_ylabel("% per decade"); ax.set_xlabel("ADT of the streamline pair (m)   south ←  → north")
ax.legend(handles=[Line2D([], [], color=CI, marker="D", ls="", ms=2.0, label="in situ, 0–780 m"), Line2D([], [], color=CA, marker="o", ms=1.8, label="altimetry, surface")], loc="upper left", fontsize=4.4)
S.plabel(ax, "f", "EKE trend per streamline pair, Oct 2005–Apr 2019", dy=9)
S.save(fig, FIGD, "FigureS_drake_insitu_validation"); plt.close(fig); print("done")
