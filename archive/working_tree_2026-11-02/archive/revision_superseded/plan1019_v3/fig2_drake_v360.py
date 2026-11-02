"""Figure 2 (v3.58: maps, along-band EKE trend, velocity trend per streamline pair; option B layout; two periods 2001-2017 / 2018-2025 in d and g, dk37): Drake Passage, rebuilt from the data (replaces the composite of the v2.0 PDF).
Row 1 (maps, 86-40 W): (a) EKE trend 2001-2025 (Product 1, 2x2 pixel average; dk34); (b) occurrence of half-power points
  with the 6 % contours of 2001-2015 and 2016-2025 (dk10); (c) trend of the isopycnal tilt |grad z| of sigma0 = 27.5
  (Roemmich-Gilson Argo, 1 deg, 2004-2025, % of the mean per decade, cells with mean tilt >= 20 m per 100 km; dots p < 0.05).
Row 2: (d) as Fig. 3d: EKE trend between the SAF and PF (cm2 s-2 per decade, dk8) with mean seafloor depth (right);
  (e) annual EKE between the SAF and PF (% of the mean of each bin; dk34).
Row 3, profiles along the LM Gould line (south to the left; upper: mean, lower: trend): (f) cross-track velocity (dk3, dk25),
  (g) surface transport per unit depth accumulated from the north (dk35), (h) EKE (dk35) with the mean depth of sigma0 = 27.3-27.7
  (Argo 2004-2025, dk33) in the upper panel.
Usage: python fig2_drake_v360.py <revision dir> <drake dir> <work out dir> <dk10 npz> <figure dir>"""
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
        ax.plot(trk.lon[trk.perm], trk.lat[trk.perm], color="#00A651", lw=1.4, zorder=8, solid_capstyle="butt")
        ax.plot(trk.lon[~trk.perm], trk.lat[~trk.perm], color="#00A651", lw=1.0, ls=(0, (2, 1.2)), zorder=8)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3]); ax.set_aspect(ASP)
    ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels(["80°W", "70°W", "60°W", "50°W", "40°W"])
    ax.set_yticks([-65, -60, -55, -50]); ax.set_yticklabels(["65°S", "60°S", "55°S", "50°S"] if ylab else [])
    ax.tick_params(labelsize=5.5)
def hbar(im, ax, label, ext):
    cax = ax.inset_axes([0.12, -0.27, 0.76, 0.045]); cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend=ext)
    cb.ax.tick_params(labelsize=5, length=2, width=0.4); cb.outline.set_linewidth(0.4); cb.set_label(label, fontsize=5.4)
fig = plt.figure(figsize=(S.W2, 0.66 * S.W2))
gA = fig.add_gridspec(2, 1, height_ratios=[0.62, 0.5], hspace=0.42)
g1 = gA[0].subgridspec(1, 3, wspace=0.08); g2 = gA[1].subgridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.45)
gL = g2[0].subgridspec(2, 1, height_ratios=[1.0, 0.75], hspace=0.5)
# a
ax = fig.add_subplot(g1[0])
im = ax.pcolormesh(F34["slon"], F34["slat"], F34["slope_map"], cmap="RdBu_r", norm=TwoSlopeNorm(0, -80, 80), shading="auto", rasterized=True, zorder=1)
base(ax)
for x_, y_, t_ in ((-71, -63.6, "Drake Passage"), (-45, -51.3, "Scotia Sea")):
    ax.text(x_, y_, t_, fontsize=5.0, ha="center", style="italic", color="#333333", zorder=9, bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.7))
hbar(im, ax, "EKE trend (cm$^2$ s$^{-2}$ per decade)", "both"); S.plabel(ax, "a", "EKE trend, 2001–2025", dy=3)
# b
ax = fig.add_subplot(g1[1])
dd = np.where(DN["all"] > 0, DN["all"], np.nan)
imd = ax.pcolormesh(DN["lon_edges"], DN["lat_edges"], dd, cmap="magma_r", vmin=0, vmax=20, shading="flat", rasterized=True, zorder=1)
lc = 0.5 * (DN["lon_edges"][1:] + DN["lon_edges"][:-1]); ac = 0.5 * (DN["lat_edges"][1:] + DN["lat_edges"][:-1])
for key, col, ls in (("early", "#56B4E9", (0, (3, 1.5))), ("late", "#0072B2", "-")):
    ax.contour(lc, ac, gaussian_filter(DN[key], 1.0), levels=[6], colors=col, linewidths=0.8, linestyles=[ls], zorder=4)
base(ax, ylab=False, fronts=False)
ax.legend(handles=[Line2D([], [], color="#56B4E9", ls=(0, (3, 1.5)), lw=0.8, label="2001–2015"), Line2D([], [], color="#0072B2", lw=0.8, label="2016–2025")],
          title="6 % contour", title_fontsize=4.6, fontsize=4.6, loc="upper left", framealpha=0.85, frameon=True, edgecolor="none").set_zorder(10)
hbar(imd, ax, "Half-power points (% of months per cell)", "max"); S.plabel(ax, "b", "Half-power points, 2001–2025", dy=3)
# c isopycnal tilt trend (RG 1 deg)
yrs, lat, lon = I["years"], I["lat"], I["lon"]; t = yrs - yrs.mean()
L = np.where(lon > 180, lon - 360, lon); sel = (L >= -90) & (L <= -30); o = np.argsort(L[sel]); LO = L[sel][o]
Z = I["p27.5"][:, :, sel][:, :, o]
dy = np.gradient(Z, lat, axis=1) / 111.195 * 100; dx = np.gradient(Z, LO, axis=2) / (111.195 * np.cos(np.radians(lat))[None, :, None]) * 100
G = np.hypot(dx, dy); tt = t[:, None, None]; Gm = G.mean(0); b = (tt * (G - Gm)).sum(0) / (t ** 2).sum(); r = G - Gm - tt * b
r1 = (r[1:] * r[:-1]).sum(0) / np.maximum((r * r).sum(0), 1e-12); ne = np.clip(len(t) * (1 - r1) / (1 + r1), 3, len(t))
pG = 2 * stats.t.sf(np.abs(b / np.sqrt((r ** 2).sum(0) / (ne - 2) / (t ** 2).sum())), ne - 2)
rel = np.where(Gm >= 20, 1000 * b / Gm, np.nan)
ax = fig.add_subplot(g1[2])
cmt = plt.get_cmap("PuOr_r").copy(); cmt.set_bad("#F2F2F2")
imt = ax.pcolormesh(LO, lat, np.ma.masked_invalid(rel), cmap=cmt, vmin=-20, vmax=20, shading="nearest", rasterized=True, zorder=1)
jj, ii = np.nonzero((pG < 0.05) & np.isfinite(rel)); ax.plot(LO[ii], lat[jj], ".", color="k", ms=1.6, zorder=7)
base(ax, ylab=False)
hbar(imt, ax, "Isopycnal tilt trend, σ$_0$ = 27.5 (% per decade)", "both"); S.plabel(ax, "c", "Isopycnal tilt trend, 2004–2025", dy=3)
# d (v3.45): as Figure 3d -- EKE trend between the SAF and PF per 2 deg longitude (cm2 s-2 per decade, calibrated Theil-Sen,
#   dk8) with 95 % band; mean seafloor depth in the band on the right axis
ax = fig.add_subplot(g2[0])
axd = ax.twinx(); axd.fill_between(AL.lon, AL.depth_mean / 1000, 7, color="#E3E3E3", step="mid", lw=0, zorder=0)
axd.set_ylim(7, -3.5); axd.set_yticks([0, 1, 2, 3, 4, 5, 6]); axd.set_ylabel("Depth (km)", color="#777777", fontsize=5.5, labelpad=1.5)
axd.tick_params(axis="y", colors="#777777", labelsize=5); axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
ax.fill_between(AL.lon, AL.lo, AL.hi, color=C["eke"], alpha=0.18, lw=0, zorder=1); ax.plot(AL.lon, AL.slope_dec, color=C["eke"], lw=0.8, zorder=3)
sg = AL.p.round(3) < 0.05
ax.plot(AL.lon[sg], AL.slope_dec[sg], "o", color=C["eke"], ms=1.8, zorder=4); ax.plot(AL.lon[~sg], AL.slope_dec[~sg], "o", mfc="white", color=C["eke"], ms=1.8, mew=0.4, zorder=4)
S.zero(ax); ax.set_xlim(-86, -40); ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels(["80°W", "70°W", "60°W", "50°W", "40°W"])
ax.set_xlabel("Longitude (flow from left to right)"); ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ per decade)"); ax.set_ylim(-60, 110)
for a_, b_ in (PASS,): ax.axvline(a_, color="#BDBDBD", lw=0.5, ls=(0, (2, 2))); ax.axvline(b_, color="#BDBDBD", lw=0.5, ls=(0, (2, 2)))
ax.text(np.mean(PASS), 0.99, "Drake Passage", ha="center", va="top", fontsize=4.8, color="#222222", transform=ax.get_xaxis_transform())
ax.annotate("", xy=(-41, 0.88), xytext=(-47, 0.88), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"))
ax.text(-44, 0.90, "flow", ha="center", va="bottom", fontsize=5, color="#555555", transform=ax.get_xaxis_transform())
S.plabel(ax, "d", "EKE trend between the SAF and PF", dy=4)

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
# e (v3.60): trend of the cross-track velocity along the LM Gould line in passage coordinates (25 km bins, drawn at the mean
#   latitude of the bin): surface geostrophic velocity from altimetry, 2001-2025 (dk5; shading 95 %), and surface baroclinic
#   velocity relative to 1975 dbar from Argo, 2004-2025 (dk45; bars 95 %). Dashed: mean crossings of the SAF and PF levels.
K5 = pd.read_csv(DR / "out" / "dk5" / "dk5_trends.csv"); T39 = pd.read_csv(WO / "dk39" / "dk39_transport_bins.csv"); lat_of = dict(zip(T39.bin, T39.lat))
D45 = pd.read_csv(WO / "dk45" / "dk45_track_trends.csv"); CA = "#8E44AD"; CU = C["tr"]
ax = fig.add_subplot(g2[1]); S.zero(ax); ax.set_xlim(*XL); ax.set_xticks(XT); ax.set_xticklabels(XTL)
ax.axvspan(XL[0], ice_lat, color="#DCEBF7", lw=0, zorder=0)
for fn, fl in FRONTS:
    ax.axvline(fl, color="#1A1A1A", lw=0.7, ls=(0, (4, 2)), zorder=1.5); ax.text(fl, 1.0, fn, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.2, fontweight="bold")
q = K5[(K5.coord == "passage") & (K5.src == "allsat") & (K5["var"] == "u") & (K5.period == "2001-2025")].sort_values("bin").copy(); q["lat"] = q.bin.map(lat_of)
ax.fill_between(q.lat, 10 * q.ci_lo, 10 * q.ci_hi, color=CU, alpha=0.18, lw=0, zorder=1); prof(ax, q.lat, 10 * q.slope, q.p, CU, ms=1.8)
qa = D45[(D45.coord == "passage") & (D45["var"] == "u_surf")].sort_values("bin").copy(); qa["lat"] = qa.bin.map(lat_of)
ax.vlines(qa.lat, qa.ci_lo, qa.ci_hi, color=CA, lw=0.6, zorder=2); prof(ax, qa.lat, qa.slope_dec, qa.p, CA, mk="s", ms=2.2, lw=0.8)
ax.set_ylim(-6, 7); ax.set_ylabel("Cross-track velocity trend\n(cm s$^{-1}$ per decade)")
ax.text(XL[0] + 0.1, 0.04, "seasonally ice-covered", transform=ax.get_xaxis_transform(), fontsize=4.6, color="#4A7BA6", va="bottom")
ax.legend(handles=[Line2D([], [], color=CU, marker="o", ms=1.8, label="surface, altimetry, 2001–2025 (shading: 95 %)"), Line2D([], [], color=CA, marker="s", ms=2.2, lw=0.8, label="surface baroclinic (ref. 1975 dbar), Argo, 2004–2025 (bars: 95 %)"),
                   Line2D([], [], color="#555555", marker="o", mfc="white", ls="", ms=2, label="open: p ≥ 0.05")], loc="upper left", fontsize=4.5)
S.plabel(ax, "e", "Cross-track velocity trend", dy=9)
S.save(fig, FIGD, "Figure2_drake_simple"); plt.close(fig)
print("done")
