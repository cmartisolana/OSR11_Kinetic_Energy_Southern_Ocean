"""Figure 2 (v3.58: maps, along-band EKE trend, velocity trend per streamline pair; option B layout; two periods 2001-2017 / 2018-2025 in d and g, dk37): Drake Passage, rebuilt from the data (replaces the composite of the v2.0 PDF).
Row 1 (maps, 86-40 W): (a) EKE trend 2001-2025 (Product 1, 2x2 pixel average; dk34); (b) occurrence of half-power points
  with the 6 % contours of 2001-2015 and 2016-2025 (dk10). [v2.1, 28 Oct 2026: the former panel c (sigma0 = 27.5 tilt trend)
  and the Argo velocity line were removed; no isopycnal slope in this figure (author request); former d and e are now c and d.]
Row 2: (d) as Fig. 3d: EKE trend between the SAF and PF (cm2 s-2 per decade, dk8) with mean seafloor depth (right);
  (e) annual EKE between the SAF and PF (% of the mean of each bin; dk34).
Row 3, profiles along the LM Gould line (south to the left; upper: mean, lower: trend): (f) cross-track velocity (dk3, dk25),
  (g) surface transport per unit depth accumulated from the north (dk35), (h) EKE (dk35) with the mean depth of sigma0 = 27.3-27.7
  (Argo 2004-2025, dk33) in the upper panel.
v2.1 (29 Oct 2026, author request): (b) 10 % contours drawn as in Fig. 3b (2001-2015 cyan solid, 2016-2025 black dashed)
  with two descriptive metrics of their difference (area overlap, mean latitude shift); (c) mean seafloor depth between the
  fronts at the native 1/8 deg resolution, dark fill; (d) net-transport text removed (not significant); titles aligned.
Usage: python fig2_drake_v361.py <revision dir> <drake dir> <work out dir> <dk10 npz> <figure dir>"""
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
    cax = ax.inset_axes([0.12, -0.17, 0.76, 0.035]); cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend=ext)
    cb.ax.tick_params(labelsize=5, length=2, width=0.4); cb.outline.set_linewidth(0.4); cb.set_label(label, fontsize=5.4)
fig = plt.figure(figsize=(S.W2, 0.80 * S.W2))
gA = fig.add_gridspec(2, 1, height_ratios=[0.9, 0.5], hspace=0.28)
g1 = gA[0].subgridspec(1, 2, wspace=0.10); g2 = gA[1].subgridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.45)
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
LEVB = 10; CEa, CLa = "#00B4E6", "#111111"                      # as Fig. 3b
SM = {k: gaussian_filter(DN[k], 1.0) for k in ("early", "late")}
for key, col, ls, lw in (("early", CEa, "-", 1.3), ("late", CLa, (0, (2.5, 1.2)), 0.9)):
    ax.contour(lc, ac, SM[key], levels=[LEVB], colors=col, linewidths=lw, linestyles=[ls], zorder=4 if key == "early" else 5)
# descriptive difference of the two contours: overlap of the enclosed areas (intersection over union, cos-weighted) and the
#   mean over longitudes of the change in occurrence-weighted latitude inside either area (late minus early; + = northward)
_me, _ml = SM["early"] >= LEVB, SM["late"] >= LEVB; _w = np.cos(np.radians(ac))[:, None] * np.ones_like(_me, float)
IOU = (_w * (_me & _ml)).sum() / (_w * (_me | _ml)).sum(); _U = _me | _ml; _sh = []
for j in range(lc.size):
    m_ = _U[:, j]; de_, dl_ = DN["early"][m_, j], DN["late"][m_, j]
    if m_.any() and de_.sum() > 0 and dl_.sum() > 0: _sh.append((ac[m_] * dl_).sum() / dl_.sum() - (ac[m_] * de_).sum() / de_.sum())
SHIFT = 111.195 * np.mean(_sh); print(f"Fig 2b: 10 % areas IoU = {IOU:.2f}; mean latitude shift = {SHIFT:+.1f} km over {len(_sh)} columns")
base(ax, ylab=False, fronts=False)
ax.legend(handles=[Line2D([], [], color=CEa, lw=1.3, label="2001–2015"), Line2D([], [], color=CLa, ls=(0, (2.5, 1.2)), lw=0.9, label="2016–2025")],
          title=f"{LEVB} % contour", title_fontsize=5.0, fontsize=5.0, loc="upper left", framealpha=0.9, frameon=True, edgecolor="none").set_zorder(10)
ax.text(0.985, 0.03, f"area overlap {100 * IOU:.0f} %\nmean shift {SHIFT:+.1f} km".replace("-", "−"), transform=ax.transAxes, ha="right", va="bottom",
        fontsize=5.0, color="#222222", zorder=12, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
hbar(imd, ax, "Half-power points (% of months per cell)", "max"); S.plabel(ax, "b", "Half-power points, 2001–2025", dy=3)
# d (v3.45): as Figure 3d -- EKE trend between the SAF and PF per 2 deg longitude (cm2 s-2 per decade, calibrated Theil-Sen,
#   dk8) with 95 % band; mean seafloor depth in the band on the right axis
ax = fig.add_subplot(g2[0])
_ab = cl.adt.interp(latitude=bs.latitude, longitude=bs.longitude).values      # mean ADT on the 1/8 deg bathymetry grid
_band = (_ab >= LEVF[0]) & (_ab <= LEVF[1]) & (bs.values < 0)
BDEP = np.nanmean(np.where(_band, -bs.values, np.nan), 0) / 1000; BLON = bs.longitude.values; _k = (BLON >= -86) & (BLON <= -40)
axd = ax.twinx(); axd.fill_between(BLON[_k], BDEP[_k], 7, color="#3A3A3A", lw=0, zorder=0)   # bottom strip, below the EKE data
axd.set_ylim(6.2, -6); axd.set_yticks([3, 4, 5, 6]); axd.set_ylabel("Depth\n(km)", color="#333333", fontsize=5.2)
axd.yaxis.set_label_coords(1.075, 0.13); axd.tick_params(axis="y", colors="#333333", labelsize=5)
axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5); axd.spines["right"].set_bounds(6.2, 2.5)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
ax.fill_between(AL.lon, AL.lo, AL.hi, color=C["eke"], alpha=0.30, lw=0, zorder=1); ax.plot(AL.lon, AL.slope_dec, color=C["eke"], lw=1.0, zorder=3)
sg = AL.p.round(3) < 0.05
ax.plot(AL.lon[sg], AL.slope_dec[sg], "o", color=C["eke"], ms=1.8, zorder=4); ax.plot(AL.lon[~sg], AL.slope_dec[~sg], "o", mfc="white", color=C["eke"], ms=1.8, mew=0.4, zorder=4)
S.zero(ax); ax.set_xlim(-86, -40); ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels(["80°W", "70°W", "60°W", "50°W", "40°W"])
ax.set_xlabel("Longitude (flow from left to right)"); ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ per decade)"); ax.set_ylim(-115, 110); ax.set_yticks([-50, -25, 0, 25, 50, 75, 100])
ax.spines["left"].set_bounds(-50, 110)
for a_, b_ in (PASS,): ax.axvline(a_, color="#BDBDBD", lw=0.5, ls=(0, (2, 2))); ax.axvline(b_, color="#BDBDBD", lw=0.5, ls=(0, (2, 2)))
ax.text(np.mean(PASS), 0.99, "Drake Passage", ha="center", va="top", fontsize=4.8, color="#222222", transform=ax.get_xaxis_transform())
ax.annotate("", xy=(-41, 0.88), xytext=(-47, 0.88), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"))
ax.text(-44, 0.90, "flow", ha="center", va="bottom", fontsize=5, color="#555555", transform=ax.get_xaxis_transform())
S.plabel(ax, "c", "EKE trend between the SAF and PF", dy=9)

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
# e (v3.58): trend of the cross-track surface velocity per pair of SSH streamlines across Drake Passage, 2001-2025 (dk5 streamwise,
#   calibrated Theil-Sen, 95 %); x = ADT at the centre of each pair, south to the left; SAF and PF zones shaded.
K5 = pd.read_csv(DR / "out" / "dk5" / "dk5_trends.csv"); A24 = pd.read_csv(DR / "out" / "dk24" / "dk24_stream_trends.csv")
adt_of = dict(zip(A24.psic.round(3), A24.adt_c)); ED = pd.read_csv(DR / "out" / "dk5" / "dk5_stream_edges.csv"); ed = dict(zip(ED.psi.round(3), ED.adt_edge))
ZONES = [("SAF", ed[-0.6], ed[-0.4]), ("PF", ed[-1.25], ed[-1.0])]
qu = K5[(K5.coord == "stream") & (K5.src == "allsat") & (K5["var"] == "u") & (K5.period == "2001-2025")].copy(); qu["adt"] = qu.bin.round(3).map(adt_of); qu = qu.dropna(subset=["adt"]).sort_values("adt")
ax = fig.add_subplot(g2[1]); CU = C["tr"]
for zn, lo_, hi_ in ZONES:
    ax.axvspan(lo_, hi_, color="#EDEDED", lw=0, zorder=0); ax.text((lo_ + hi_) / 2, 1.0, zn, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.2, fontweight="bold")
S.zero(ax); ax.set_xlim(-1.16, 0.44)
ax.vlines(qu.adt, 10 * qu.ci_lo, 10 * qu.ci_hi, color=CU, lw=0.6, zorder=2); prof(ax, qu.adt, 10 * qu.slope, qu.p, CU, mk="s", ms=2.4, lw=0.9)
ax.set_ylim(-8, 8); ax.set_ylabel("Cross-track velocity trend\n(cm s$^{-1}$ per decade)"); ax.set_xlabel("ADT of the streamline pair (m)   south ←  → north")
# v3.61: Argo surface baroclinic cross-track velocity (relative to 1975 dbar) per streamline pair, 2004-2025 (dk45: dynamic height
#   interpolated onto the track, points assigned to pairs by the mean ADT; only pairs with complete profiles, SAF-PF zone)
ax.legend(handles=[Line2D([], [], color=CU, marker="s", ms=2.4, lw=0.9, label="surface, altimetry, 2001–2025 (bars: 95 %)"),
                   Line2D([], [], color="#555555", marker="s", mfc="white", ls="", ms=2.4, label="open: p ≥ 0.05")], loc="upper left", fontsize=4.5)
S.plabel(ax, "d", "Cross-track velocity trend", dy=9)
S.save(fig, FIGD, "Figure2_drake_simple"); plt.close(fig)
print("done")
