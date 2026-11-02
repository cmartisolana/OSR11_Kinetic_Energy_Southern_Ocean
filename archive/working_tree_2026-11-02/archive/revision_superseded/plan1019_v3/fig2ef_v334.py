"""Fig. 2e-f (v3.34): trends across the streamlines of the LM Gould transect (method of Sect. 2.5), altimetry (Product 1)
with the in situ values of Gutierrez-Villanueva et al. (2023) recomputed with the same test (dk5_trends.csv).
x: streamline pairs of GV23 mapped onto the ADT (centre of each pair, m), south to the left.
(e) upper: mean EKE and KE per pair; lower: trends (% of the mean per decade) of EKE and KE 2001-2025 (EKE 95 % interval
    shaded), altimetric EKE over Oct 2005-Apr 2019 (dashed) and in situ EKE 0-780 m over Oct 2005-Apr 2019 (orange).
(f) upper: mean surface transport per unit depth per pair; lower: trends of the surface transport 2001-2025 (shaded 95 %),
    over Oct 2005-Apr 2019 (dashed) and in situ 0-760 m transport over Oct 2005-Apr 2019 (orange).
Drawn on a transparent page of the v2.0 figure size at the old panel positions; composed by fig2_compose_v332.py.
Usage: python fig2ef_v333.py <dk24 dir> <dk5 dir> <out pdf>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
DK24, DK5, OUTF = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
plt.rcParams.update({   # style_nc.py
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.5, "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 5.4,
    "legend.frameon": False, "legend.handlelength": 1.6, "legend.handletextpad": 0.4, "legend.columnspacing": 0.9,
    "legend.borderaxespad": 0.2, "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#222222", "axes.labelcolor": "#222222", "text.color": "#222222", "xtick.color": "#222222",
    "ytick.color": "#222222", "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "lines.linewidth": 1.0,
    "pdf.fonttype": 42, "mathtext.default": "regular", "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans"})
C = dict(ke="#D55E00", eke="#009E73", tr="#0072B2", obs="#E69F00")
sigp = lambda p: np.isfinite(p) and round(float(p), 3) < 0.05
PW, PH = 471.0592603718, 476.1464860667
X0E, X1E, X0F, X1F, Y0, YM, Y1 = 37.676328, 224.201362, 246.584366, 433.109399, 23.170957, 112.0, 156.512689
A = pd.read_csv(DK24 / "dk24_stream_trends.csv")
G = pd.read_csv(DK5 / "dk5_trends.csv"); G = G[(G.coord == "stream") & (G.src == "insitu")].copy()
adt_of = dict(zip(A.psic.round(3), A.adt_c)); G["adt_c"] = G.bin.round(3).map(adt_of); G["pct"] = 1000 * G.slope / G["mean"]
G["pct_lo"] = 1000 * G.ci_lo / G["mean"]; G["pct_hi"] = 1000 * G.ci_hi / G["mean"]
def sel(var, per):
    return A[(A["var"] == var) & (A.period == per)].sort_values("adt_c")
E = pd.read_csv(DK5 / "dk5_stream_edges.csv"); ed = dict(zip(E.psi.round(3), E.adt_edge))
ZONES = [("SAF", ed[-0.6], ed[-0.4]), ("PF", ed[-1.25], ed[-1.0])]
XL = (-1.16, 0.44)
fig = plt.figure(figsize=(PW / 72, PH / 72)); fig.patch.set_alpha(0)
def ax_at(x0, x1, y0, y1):
    a = fig.add_axes([x0 / PW, y0 / PH, (x1 - x0) / PW, (y1 - y0) / PH]); a.set_facecolor("white"); a.set_xlim(*XL); return a
def zones(ax, labels=False):
    for nm, a, b in ZONES:
        ax.axvspan(a, b, color="#EDEDED", lw=0, zorder=0)
        if labels: ax.text(0.5 * (a + b), 1.02, nm, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.5, color="#555555", fontweight="bold")
def prof(ax, x, y, lo, hi, p, col, ls="-", ms=2.6, mk="o", shade=False, lw=0.9, z=3):
    if shade: ax.fill_between(x, lo, hi, color=col, alpha=0.15, lw=0, zorder=1)
    ax.plot(x, y, color=col, ls=ls, lw=lw, zorder=z)
    s = np.array([sigp(v) for v in p])
    ax.plot(np.asarray(x)[s], np.asarray(y)[s], mk, color=col, ms=ms, zorder=z + 1)
    ax.plot(np.asarray(x)[~s], np.asarray(y)[~s], mk, color=col, mfc="white", ms=ms, mew=0.6, zorder=z + 1)
def insitu(ax, var, yl):
    g = G[G["var"] == var].sort_values("adt_c"); x = g.adt_c.values + 0.012
    lo, hi = np.clip(g.pct_lo, *yl), np.clip(g.pct_hi, *yl)
    ax.vlines(x, lo, hi, color=C["obs"], lw=0.5, alpha=0.5, zorder=2)
    for xx, l0, h0 in zip(x, g.pct_lo, g.pct_hi):
        if l0 < yl[0]: ax.plot(xx, yl[0] + 0.01 * (yl[1] - yl[0]), "v", color=C["obs"], alpha=0.6, ms=1.6, zorder=2, clip_on=False)
        if h0 > yl[1]: ax.plot(xx, yl[1] - 0.01 * (yl[1] - yl[0]), "^", color=C["obs"], alpha=0.6, ms=1.6, zorder=2, clip_on=False)
    prof(ax, x, g.pct, None, None, g.p, C["obs"], ls="", ms=2.6, mk="D", z=4)
def xax(ax):
    ax.set_xlabel("ADT of the streamline pair (m)   south ← → north", labelpad=1.5)
def label(x, letter, title):
    fig.text(x / PW, 168.2 / PH, letter, fontsize=8, fontweight="bold", va="baseline")
    fig.text((x + 9) / PW, 168.7 / PH, title, fontsize=6.5, va="baseline")

# ---------------- e: EKE and KE
top = ax_at(X0E, X1E, YM + 6, Y1); zones(top, True)
e, k = sel("eke", "2001-2025"), sel("ke", "2001-2025")
top.fill_between(e.adt_c, 0, e["mean"], color=C["eke"], alpha=0.25, lw=0); top.plot(e.adt_c, e["mean"], color=C["eke"], lw=0.8)
top.plot(k.adt_c, k["mean"], color=C["ke"], lw=0.8)
top.set_ylim(0, 900); top.set_yticks([0, 400, 800]); top.set_xticklabels([]); top.set_ylabel("Mean\n(cm$^2$ s$^{-2}$)", fontsize=5.8)
ax = ax_at(X0E, X1E, Y0, YM); zones(ax); YLE = (-45, 75); ax.set_ylim(*YLE); ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=0)
prof(ax, e.adt_c, e.pct, e.pct_lo, e.pct_hi, e.p, C["eke"], shade=True)
prof(ax, k.adt_c, k.pct, None, None, k.p, C["ke"])
e2 = sel("eke", "2005-2019"); prof(ax, e2.adt_c, e2.pct, None, None, e2.p, C["eke"], ls=(0, (3, 1.5)), ms=2.0, lw=0.7, z=2)
insitu(ax, "eke", YLE)
ax.set_ylabel("Trend (% per decade)"); xax(ax)
ax.legend(handles=[Line2D([], [], color=C["eke"], marker="o", ms=2.6, label="EKE, 2001–2025"),
                   Line2D([], [], color=C["ke"], marker="o", ms=2.6, label="KE, 2001–2025"),
                   Line2D([], [], color=C["eke"], ls=(0, (3, 1.5)), lw=0.7, marker="o", ms=2.0, label="EKE, 2005–2019"),
                   Line2D([], [], color=C["obs"], ls="", marker="D", ms=2.6, label="in situ EKE 0–780 m, 2005–2019"),
                   Line2D([], [], color="#555555", ls="", marker="o", mfc="white", ms=2.6, mew=0.6, label="open: p ≥ 0.05")],
          loc="upper left", ncol=2, fontsize=5.0, handlelength=1.5, columnspacing=0.7)
label(13.676328, "e", "EKE and KE across streamlines, LM Gould transect")

# ---------------- f: surface agreement with the in situ record (dk3, Product 1) and the extended altimetric record
P = pd.read_csv(DK5 / "dk3_velocity_profile_allsat.csv"); Q = pd.read_csv(DK5 / "dk3_eke_annual_north_allsat.csv")
top = fig.add_axes([X0F / PW, 104 / PH, (X1F - X0F) / PW, (Y1 - 104) / PH]); top.set_facecolor("white")
top.axhline(0, color="#9A9A9A", lw=0.5, zorder=0)
top.plot(P.latb, P.ualt_mean, color=C["eke"], lw=0.9, marker="o", ms=1.8, label="altimetry (surface)")
top.plot(P.latb, P.uin_mean, color=C["obs"], lw=0.9, marker="D", ms=1.8, label="in situ ADCP, 70–118 m")
top.set_xlim(-63.6, -54.9); top.set_ylim(-5, 45); top.set_yticks([0, 20, 40])
top.set_xticks([-63, -61, -59, -57, -55]); top.set_xticklabels(["63°S", "61°S", "59°S", "57°S", "55°S"])
top.spines["right"].set_visible(True); top.spines["left"].set_visible(False); top.yaxis.tick_right(); top.yaxis.set_label_position("right")
top.set_ylabel("Mean cross-track\nvelocity (cm s$^{-1}$)", fontsize=5.8)
top.legend(loc="upper left", fontsize=5.0, handlelength=1.5)
Z = np.load(DK24 / "dk24_section_fields_allsat.npz"); u, v, vi, ice, lat = Z["u"], Z["v"], Z["vi"], Z["ice"] > 0.10, Z["lat"]
ek = (0.5 * ((u - np.nanmean(u, 0)) ** 2 + (v - np.nanmean(v, 0)) ** 2) + 0.5 * vi) * 1e4; ek[ice] = np.nan
yrs = np.repeat(np.arange(2001, 2026), 12); m = np.nanmean(ek[:, lat > -58.5], 1)
full = pd.Series(m).groupby(yrs).mean(); ref = full.loc[2005:2019].mean()
ax = fig.add_axes([X0F / PW, Y0 / PH, (X1F - X0F) / PW, (84 - Y0) / PH]); ax.set_facecolor("white")
ax.axhline(0, color="#9A9A9A", lw=0.5, zorder=0); ax.axvspan(2004.5, 2019.5, color="#F2F2F2", lw=0, zorder=0)
ax.plot(full.index, 100 * (full / ref - 1), color=C["eke"], lw=0.6, alpha=0.5, zorder=2, label="altimetry, all months 2001–2025")
ax.plot(Q.year, 100 * (Q.eke_alt / Q.eke_alt.mean() - 1), color=C["eke"], lw=0.9, marker="o", ms=1.8, zorder=3, label="altimetry, months of the transects")
ax.plot(Q.year, 100 * (Q.eke_in / Q.eke_in.mean() - 1), color=C["obs"], lw=0.9, marker="D", ms=1.8, zorder=4, label="in situ EKE 0–780 m")
ax.set_xlim(2000.5, 2025.5); ax.set_xticks([2005, 2010, 2015, 2020, 2025]); ax.set_ylim(-40, 60)
ax.spines["right"].set_visible(True); ax.spines["left"].set_visible(False); ax.yaxis.tick_right(); ax.yaxis.set_label_position("right")
ax.set_ylabel("Annual EKE north\nof 58.5°S (% of mean)", fontsize=5.8)
ax.legend(loc="upper left", fontsize=5.0, handlelength=1.5, ncol=1)
label(246.584366, "f", "Surface agreement with the in situ record")
fig.savefig(OUTF, transparent=True); print("ok")
