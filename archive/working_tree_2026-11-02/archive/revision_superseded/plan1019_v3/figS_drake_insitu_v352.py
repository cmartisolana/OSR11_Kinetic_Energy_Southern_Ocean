"""Figure S (v3.52): in situ validation along the LM Gould line (panel e of Fig. 2 v3.52 as a standalone SI figure). Original header:
Figure 2 (v3.52: maps, along-band EKE trend, transect in situ vs altimetry; option B layout; two periods 2001-2017 / 2018-2025 in d and g, dk37): Drake Passage, rebuilt from the data (replaces the composite of the v2.0 PDF).
Row 1 (maps, 86-40 W): (a) EKE trend 2001-2025 (Product 1, 2x2 pixel average; dk34); (b) occurrence of half-power points
  with the 6 % contours of 2001-2015 and 2016-2025 (dk10); (c) trend of the isopycnal tilt |grad z| of sigma0 = 27.5
  (Roemmich-Gilson Argo, 1 deg, 2004-2025, % of the mean per decade, cells with mean tilt >= 20 m per 100 km; dots p < 0.05).
Row 2: (d) as Fig. 3d: EKE trend between the SAF and PF (cm2 s-2 per decade, dk8) with mean seafloor depth (right);
  (e) annual EKE between the SAF and PF (% of the mean of each bin; dk34).
Row 3, profiles along the LM Gould line (south to the left; upper: mean, lower: trend): (f) cross-track velocity (dk3, dk25),
  (g) surface transport per unit depth accumulated from the north (dk35), (h) EKE (dk35) with the mean depth of sigma0 = 27.3-27.7
  (Argo 2004-2025, dk33) in the upper panel.
Usage: python fig2_drake_v352.py <revision dir> <drake dir> <work out dir> <dk10 npz> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
from scipy.ndimage import gaussian_filter
REV, DR, WO, D10, FIGD = (Path(a) for a in sys.argv[1:6])
sys.path.insert(0, str(REV))
import style_nc as S
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import TwoSlopeNorm
plt = S.plt; C = S.C
sigp = lambda p: np.isfinite(p) and round(float(p), 3) < 0.05
F34 = np.load(WO / "dk34" / "dk34_drake_fields.npz"); DN = np.load(D10); I = np.load(WO / "dk32" / "dk32_iso.npz")
SEC = np.load(WO / "dk33" / "dk33_drake_section.npz"); T33 = pd.read_csv(WO / "dk33" / "dk33_drake_isopycnals.csv").set_index("sigma0")
AL = pd.read_csv(DR / "dk8" / "dk8_along_acc_fronts_allsat.csv")
P3 = pd.read_csv(DR / "out" / "dk3" / "dk3_velocity_profile_allsat.csv"); V = pd.read_csv(DR / "out" / "dk25" / "dk25_profile_trends.csv")
trk = pd.read_csv(REV / "plan1019" / "outputs" / "d1" / "d1b_track.csv")
box = dict(longitude=slice(-87, -33), latitude=slice(-66, -46))
bs = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation.sel(**box)
cl = xr.open_dataset(REV / "outputs" / "clim_allsat.nc").sel(**box)
LEVF = [float(F34["lev_pf"]), float(F34["lev_saf"])]
EXT = (-86, -40, -65, -50); ASP = 1 / np.cos(np.radians(57)); PASS = (-68, -58)
CE, CL_ = "#00B4E6", "#111111"; CS = {27.3: "#6A51A3", 27.5: "#2171B5", 27.7: "#6BAED6"}
def base(ax, ylab=True, fronts=True, track=True):
    ax.contour(bs.longitude, bs.latitude, -bs.values, levels=[3000], colors="#8A8A8A", linewidths=0.35, zorder=2)
    ax.contourf(bs.longitude, bs.latitude, bs.values, levels=[0, 1e5], colors=["#3A3A3A"], zorder=6)
    if fronts:
        ax.contour(cl.longitude, cl.latitude, cl.adt.values, levels=LEVF, colors="#1A1A1A", linewidths=0.9, linestyles="solid", zorder=5)
        ax.text(-85.5, -55.5, "SAF", fontsize=5.5, fontweight="bold", zorder=8); ax.text(-85.5, -60.6, "PF", fontsize=5.5, fontweight="bold", zorder=8)
    if track:
        ax.plot(trk.lon[trk.perm], trk.lat[trk.perm], color=C["obs"], lw=1.4, zorder=8, solid_capstyle="butt")
        ax.plot(trk.lon[~trk.perm], trk.lat[~trk.perm], color=C["obs"], lw=1.0, ls=(0, (2, 1.2)), zorder=8)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3]); ax.set_aspect(ASP)
    ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels(["80°W", "70°W", "60°W", "50°W", "40°W"])
    ax.set_yticks([-65, -60, -55, -50]); ax.set_yticklabels(["65°S", "60°S", "55°S", "50°S"] if ylab else [])
    ax.tick_params(labelsize=5.5)
def hbar(im, ax, label, ext):
    cax = ax.inset_axes([0.12, -0.27, 0.76, 0.045]); cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend=ext)
    cb.ax.tick_params(labelsize=5, length=2, width=0.4); cb.outline.set_linewidth(0.4); cb.set_label(label, fontsize=5.4)
fig = plt.figure(figsize=(0.55 * S.W2, 0.9 * S.W2)); gS = fig.add_gridspec(1, 1)
# row 3 (v3.45): profiles along the LM Gould line, south to the left; upper: mean, lower: trend; seasonally ice-covered shaded
P35 = pd.read_csv(WO / "dk35" / "dk35_transect_profiles.csv")
ice_lat = V.latb[V.months < 300].max() + 0.125; XL = (-63.6, -54.9); XT = [-63, -61, -59, -57, -55]; XTL = ["63°S", "61°S", "59°S", "57°S", "55°S"]
def pair(cell, lab, ttl, right=False):
    gg = cell.subgridspec(2, 1, height_ratios=[0.75, 1], hspace=0.1); top = fig.add_subplot(gg[0]); bot = fig.add_subplot(gg[1])
    for a_ in (top, bot):
        a_.axvspan(XL[0], ice_lat, color="#DCEBF7", lw=0, zorder=0); a_.set_xlim(*XL); a_.set_xticks(XT); a_.axhline(0, color="#9A9A9A", lw=0.5, zorder=1)
        if right: a_.spines["right"].set_visible(True); a_.spines["left"].set_visible(False); a_.yaxis.tick_right(); a_.yaxis.set_label_position("right")
    top.set_xticklabels([]); bot.set_xticklabels(XTL)
    S.plabel(top, lab, ttl, dy=4); return top, bot
def prof(ax, x, y, p, col, ls="-", ms=2.2, lw=0.9, z=3, mk="o"):
    x, y, p = np.asarray(x), np.asarray(y), np.asarray(p); ax.plot(x, y, color=col, ls=ls, lw=lw, zorder=z); s = np.array([sigp(v) for v in p])
    ax.plot(x[s], y[s], mk, color=col, ms=ms, zorder=z + 1); ax.plot(x[~s], y[~s], mk, color=col, mfc="white", ms=ms, mew=0.6, zorder=z + 1)
# row 3, option B (v3.46): two wide panels along the LM Gould line (south to the left). (f) mean structure: cross-track
#   velocity (altimetry vs in situ ADCP) and EKE with the mean isopycnals (Argo); (g) trends: velocity, transport accumulated
#   from the north, EKE (stacked, shared latitude axis).
_sec = np.load(DR / "out" / "dk24" / "dk24_section_fields_allsat.npz"); _a = np.nanmean(_sec["adt"], 0); _la = _sec["lat"]
def _cross(lev):
    j = np.where(np.diff(np.sign(_a - lev)) != 0)[0][0]; return _la[j] + (lev - _a[j]) * (_la[j + 1] - _la[j]) / (_a[j + 1] - _a[j])
FRONTS = [("SAF", _cross(float(F34["lev_saf"]))), ("PF", _cross(float(F34["lev_pf"])))]   # mean (2001-2025) crossings of the map front levels
def strips(cell, n, ratios, lab, ttl):
    gg = cell.subgridspec(n, 1, height_ratios=ratios, hspace=0.12); axs = [fig.add_subplot(gg[i]) for i in range(n)]
    for i, a_ in enumerate(axs):
        a_.axvspan(XL[0], ice_lat, color="#DCEBF7", lw=0, zorder=0); a_.set_xlim(*XL); a_.set_xticks(XT); a_.axhline(0, color="#9A9A9A", lw=0.5, zorder=1)
        for fn, fl in FRONTS:
            a_.axvline(fl, color="#1A1A1A", lw=0.7, ls=(0, (4, 2)), zorder=1.5)
            if i == 0: a_.text(fl, 1.0, fn, transform=a_.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.2, fontweight="bold", color="#1A1A1A")
        a_.set_xticklabels(XTL if i == n - 1 else [])
    S.plabel(axs[0], lab, ttl, dy=9); return axs
# e (v3.52): along the LM Gould line in passage coordinates (25-km bins from the north, GV23; drawn at the latitude of the bin
#   centre on the mean track): in situ (0-760 m transport, 70-118 m velocity; Oct 2005-Apr 2019) vs altimetry (surface
#   transport per unit depth, velocity, EKE) over the same period and over 2001-2025. Trends: calibrated Theil-Sen (dk5, dk39).
K5 = pd.read_csv(DR / "out" / "dk5" / "dk5_trends.csv"); K5 = K5[K5.coord == "passage"]
T39 = pd.read_csv(WO / "dk39" / "dk39_transport_bins.csv"); lat_of = dict(zip(T39.bin, T39.lat))
def k5(src, var, per): q = K5[(K5.src == src) & (K5["var"] == var) & (K5.period == per)].sort_values("bin").copy(); q["lat"] = q.bin.map(lat_of); return q
e1, e2, e3, e4 = strips(gS[0], 4, [0.8, 1, 1, 1], "a", "Along the LM Gould line: in situ (Oct 2005–Apr 2019) and altimetry")
CI, CA = C["obs"], C["eke"]; L25, L0519 = "-", (0, (3, 1.5))
# 1 mean transport per bin
qi = k5("insitu", "transport", "2005-2019")
e1.fill_between(qi.lat, 0, qi["mean"], color=CI, alpha=0.25, lw=0, zorder=2); e1.plot(qi.lat, qi["mean"], color=CI, lw=0.9, zorder=3)
e1.set_ylabel("Transport\n(Sv per bin)", color=CI); e1.set_ylim(0, 8); e1.tick_params(axis="y", colors=CI)
ea = e1.twinx(); ea.plot(T39.lat, T39["mean"] / 1e3, color=CA, lw=0.9, zorder=3); ea.set_ylim(0, 8 * (T39["mean"].sum() / 1e3) / qi["mean"].sum())
ea.set_ylabel("Surface transport\n(10$^3$ m$^2$ s$^{-1}$)", color=CA, fontsize=5.5); ea.tick_params(labelsize=5, colors=CA); ea.spines["right"].set_visible(True); ea.spines["right"].set_linewidth(0.5)
e1.legend(handles=[Line2D([], [], color=CI, lw=0.9, label="in situ, 0–760 m"), Line2D([], [], color=CA, lw=0.9, label="altimetry, surface")], loc="upper left", fontsize=4.5, ncol=1)
# 2 transport trends: in situ in Sv per decade (left); altimetry in 10^3 m2 s-1 per decade (right), right axis scaled as in strip 1
FAC = (T39["mean"].sum() / 1e3) / qi["mean"].sum()        # 10^3 m2 s-1 of surface transport per Sv of 0-760 m transport
e2.vlines(qi.lat, 10 * qi.ci_lo, 10 * qi.ci_hi, color=CI, lw=0.5, alpha=0.7, zorder=2); prof(e2, qi.lat, 10 * qi.slope, qi.p, CI, mk="D", ms=2.0, lw=0)
e2.set_ylabel("Transport trend\n(Sv dec$^{-1}$)", color=CI); e2.set_ylim(-0.9, 1.3); e2.tick_params(axis="y", colors=CI)
eb = e2.twinx(); eb.set_ylim(-0.9 * FAC, 1.3 * FAC); eb.spines["right"].set_visible(True); eb.spines["right"].set_linewidth(0.5)
eb.fill_between(T39.lat, T39.ci_lo_0125 / 1e3, T39.ci_hi_0125 / 1e3, color=CA, alpha=0.18, lw=0, zorder=1)
prof(eb, T39.lat, T39.slope_0519 / 1e3, T39.p_0519, CA, ls=L0519, ms=1.5, lw=0.7, z=2); prof(eb, T39.lat, T39.slope_0125 / 1e3, T39.p_0125, CA, ms=1.8)
eb.set_ylabel("Surface transport\n(10$^3$ m$^2$ s$^{-1}$ dec$^{-1}$)", color=CA, fontsize=5.5); eb.tick_params(labelsize=5, colors=CA)
e2.legend(handles=[Line2D([], [], color=CI, marker="D", ls="", ms=2.0, label="in situ, Oct 2005–Apr 2019"), Line2D([], [], color=CA, ls=L0519, lw=0.7, marker="o", ms=1.5, label="altimetry, Oct 2005–Apr 2019"),
                   Line2D([], [], color=CA, marker="o", ms=1.8, label="altimetry, 2001–2025 (shading: 95 %)")], loc="upper left", fontsize=4.3, ncol=3, frameon=True, facecolor="white", framealpha=0.85, edgecolor="none")
# 3 velocity trends
q = k5("insitu", "u", "2005-2019"); e3.vlines(q.lat, 10 * q.ci_lo, 10 * q.ci_hi, color=CI, lw=0.5, alpha=0.7, zorder=2); prof(e3, q.lat, 10 * q.slope, q.p, CI, mk="D", ms=2.0, lw=0)
q = k5("allsat", "u", "2005-2019"); prof(e3, q.lat, 10 * q.slope, q.p, CA, ls=L0519, ms=1.5, lw=0.7, z=2)
q = k5("allsat", "u", "2001-2025"); e3.fill_between(q.lat, 10 * q.ci_lo, 10 * q.ci_hi, color=CA, alpha=0.18, lw=0, zorder=1); prof(e3, q.lat, 10 * q.slope, q.p, CA, ms=1.8)
e3.set_ylabel("Velocity trend\n(cm s$^{-1}$ dec$^{-1}$)"); e3.set_ylim(-9, 9)
# 4 EKE trends
q = k5("allsat", "eke", "2005-2019"); prof(e4, q.lat, 1000 * q.slope / q["mean"], q.p, CA, ls=L0519, ms=1.5, lw=0.7, z=2)
q = k5("allsat", "eke", "2001-2025"); e4.fill_between(q.lat, 1000 * q.ci_lo / q["mean"], 1000 * q.ci_hi / q["mean"], color=CA, alpha=0.18, lw=0, zorder=1); prof(e4, q.lat, 1000 * q.slope / q["mean"], q.p, CA, ms=1.8)
e4.set_ylabel("EKE trend\n(% per decade)"); e4.set_ylim(-40, 60)
for a_ in (e3, e4): a_.spines["right"].set_visible(True); a_.spines["left"].set_visible(False); a_.yaxis.tick_right(); a_.yaxis.set_label_position("right")
e4.text(XL[0] + 0.1, 0.04, "seasonally ice-covered", transform=e4.get_xaxis_transform(), fontsize=4.6, color="#4A7BA6", va="bottom")
S.save(fig, FIGD, "FigureS_drake_insitu_validation"); plt.close(fig)
print("done")
