"""Figure 2 (v3.43): Drake Passage, rebuilt from the data (replaces the composite of the v2.0 PDF).
Row 1 (maps, 86-40 W): (a) EKE trend 2001-2025 (Product 1, 2x2 pixel average; dk34); (b) occurrence of half-power points
  with the 6 % contours of 2001-2015 and 2016-2025 (dk10); (c) trend of the isopycnal tilt |grad z| of sigma0 = 27.5
  (Roemmich-Gilson Argo, 1 deg, 2004-2025, % of the mean per decade, cells with mean tilt >= 20 m per 100 km; dots p < 0.05).
Row 2: (d) EKE trend between the SAF and PF per 2 deg longitude (calibrated Theil-Sen, dk8) and seafloor depth;
  (e) annual EKE between the SAF and PF (% of the mean of each bin; dk34).
Row 3, along the LM Gould line (south to the left): (f) mean depth of sigma0 = 27.3, 27.5, 27.7 (upper) and their depth trends
  2004-2025 (lower; least squares, 95 % with n_eff from the lag-1 autocorrelation; dk33); (g) cross-passage velocity:
  altimetry vs in situ ADCP (upper) and velocity trends (lower; dk3, dk25).
Usage: python fig2_drake_v343.py <revision dir> <drake dir> <work out dir> <dk10 npz> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
from scipy.ndimage import gaussian_filter
REV, DR, WO, D10, FIGD = (Path(a) for a in sys.argv[1:6])
sys.path.insert(0, str(REV))
import style_nc as S
from matplotlib.lines import Line2D
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
fig = plt.figure(figsize=(S.W2, 1.12 * S.W2))
gA = fig.add_gridspec(2, 1, height_ratios=[0.62, 1.9], hspace=0.36)
gB = gA[1].subgridspec(2, 1, hspace=0.3)
g1 = gA[0].subgridspec(1, 3, wspace=0.08); g2 = gB[0].subgridspec(1, 2, wspace=0.14); g3 = gB[1].subgridspec(1, 2, wspace=0.14)
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
hbar(imd, ax, "Half-power points (% of months per cell)", "max"); S.plabel(ax, "b", "Occurrence of half-power points", dy=3)
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
hbar(imt, ax, "Isopycnal tilt trend, σ$_0$ = 27.5 (% per decade)", "both"); S.plabel(ax, "c", "Isopycnal tilt trend, 2004–2025 (Argo)", dy=3)
# d along-ACC EKE trend + depth
g2L = g2[0].subgridspec(2, 1, height_ratios=[1, 0.4], hspace=0.08)
ax = fig.add_subplot(g2L[0]); ax.axvspan(*PASS, color="#F1F1F1", lw=0, zorder=0)
ax.fill_between(AL.lon, AL.lo, AL.hi, color=C["eke"], alpha=0.18, lw=0, zorder=1); ax.plot(AL.lon, AL.slope_dec, color=C["eke"], lw=1.0, zorder=3)
sg = AL.p.round(3) < 0.05
ax.plot(AL.lon[sg], AL.slope_dec[sg], "o", color=C["eke"], ms=2.4, zorder=4); ax.plot(AL.lon[~sg], AL.slope_dec[~sg], "o", mfc="white", color=C["eke"], ms=2.4, mew=0.5, zorder=4)
S.zero(ax); ax.set_xlim(-86, -40); ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels([]); ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ dec$^{-1}$)")
ax.text(np.mean(PASS), 1.0, "Drake Passage", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.2, color="#555555")
ax.annotate("", xy=(-41, 0.93), xytext=(-48, 0.93), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"))
ax.text(-44.5, 0.96, "flow", ha="center", va="bottom", fontsize=5, color="#555555", transform=ax.get_xaxis_transform())
S.plabel(ax, "d", "EKE trend between the SAF and PF, and depth", dy=4)
ax = fig.add_subplot(g2L[1]); ax.axvspan(*PASS, color="#F1F1F1", lw=0, zorder=0)
ax.fill_between(AL.lon, -AL.depth_mean / 1000, -6, color="#8A8A8A", step="mid", lw=0, zorder=1)
ax.step(AL.lon, -AL.depth_min / 1000, where="mid", color="#222222", lw=0.8, zorder=2)
ax.text(-85, -2.4, "shallowest depth", fontsize=5, zorder=3); ax.text(-52, -5.0, "mean depth", fontsize=5, color="white", zorder=3)
ax.set_ylim(-5.6, 0); ax.set_yticks([-5, -3, -1]); ax.set_yticklabels(["5", "3", "1"]); ax.set_ylabel("Seafloor\ndepth (km)")
ax.set_xlim(-86, -40); ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels(["80°W", "70°W", "60°W", "50°W", "40°W"])
# e Hovmoller
ax = fig.add_subplot(g2[1])
hl, hy = F34["hov_lon"], F34["hov_years"]
imh = ax.pcolormesh(hl, hy, F34["hov"], cmap="RdBu_r", vmin=-40, vmax=40, shading="nearest", rasterized=True)
for x_ in PASS: ax.axvline(x_, color="#222222", lw=0.6, ls=(0, (1, 1.5)))
ax.text(np.mean(PASS), 1.0, "Drake Passage", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.2, color="#555555")
ax.set_xlim(-86, -40); ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels(["80°W", "70°W", "60°W", "50°W", "40°W"]); ax.set_yticks([2005, 2010, 2015, 2020, 2025])
cax = ax.inset_axes([1.03, 0.1, 0.03, 0.8]); cb = fig.colorbar(imh, cax=cax, extend="both"); cb.ax.tick_params(labelsize=5, length=2, width=0.4)
cb.outline.set_linewidth(0.4); cb.set_label("EKE anomaly (% of mean)", fontsize=5.4)
S.plabel(ax, "e", "Annual EKE between the SAF and PF", dy=4)
# f isopycnals along the LM Gould line
g3L = g3[0].subgridspec(2, 1, height_ratios=[0.8, 1], hspace=0.1); g3R = g3[1].subgridspec(2, 1, height_ratios=[0.8, 1], hspace=0.1)
ice_lat = V.latb[V.months < 300].max() + 0.125; XL = (-63.6, -54.9); XT = [-63, -61, -59, -57, -55]; XTL = ["63°S", "61°S", "59°S", "57°S", "55°S"]
la_ = SEC["tlat"]; yrs3 = SEC["years"]; t3 = yrs3 - yrs3.mean()
top = fig.add_subplot(g3L[0]); bot = fig.add_subplot(g3L[1])
for a_ in (top, bot): a_.axvspan(XL[0], ice_lat, color="#DCEBF7", lw=0, zorder=0); a_.set_xlim(*XL); a_.set_xticks(XT)
for s in (27.3, 27.5, 27.7):
    Zs = SEC[f"p{s}"]; ok = np.isfinite(Zs).all(0); x = la_[ok]; Y = Zs[:, ok]
    top.plot(x, Y.mean(0), color=CS[s], lw=1.1 if s == 27.5 else 0.8, label=f"σ$_0$ = {s}")
    m = Y.mean(0); b3 = (t3[:, None] * (Y - m)).sum(0) / (t3 ** 2).sum(); rr = Y - m - t3[:, None] * b3
    a1 = (rr[1:] * rr[:-1]).sum(0) / (rr * rr).sum(0); n3 = np.clip(len(t3) * (1 - a1) / (1 + a1), 3, len(t3))
    se3 = np.sqrt((rr ** 2).sum(0) / (n3 - 2) / (t3 ** 2).sum()); q = stats.t.ppf(0.975, n3 - 2); p3 = 2 * stats.t.sf(np.abs(b3 / se3), n3 - 2)
    if s == 27.5: bot.fill_between(x, 10 * (b3 - q * se3), 10 * (b3 + q * se3), color=CS[s], alpha=0.18, lw=0, zorder=1)
    bot.plot(x, 10 * b3, color=CS[s], lw=1.1 if s == 27.5 else 0.7, zorder=3)
top.set_ylim(2000, 0); top.set_ylabel("Mean depth\n(dbar)"); top.set_xticklabels([]); top.legend(loc="lower left", fontsize=4.8, ncol=1)
bot.axhline(0, color="#9A9A9A", lw=0.5); bot.set_xticklabels(XTL); bot.set_ylabel("Depth trend\n(dbar per decade)"); bot.set_ylim(-40, 40)
bot.text(XL[0] + 0.1, 0.03, "seasonally ice-covered", transform=bot.get_xaxis_transform(), fontsize=5.0, color="#4A7BA6", va="bottom")
r5 = T33.loc[27.5]
bot.text(0.99, 0.97, f"cross-passage slope (σ$_0$ = 27.5): {r5.slope_trend_pct:+.1f} % per decade ({r5.slope_trend_pct - r5.slope_ci_pct:+.1f} to {r5.slope_trend_pct + r5.slope_ci_pct:+.1f})".replace("-", "−"),
         transform=bot.transAxes, ha="right", va="top", fontsize=4.9, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
S.plabel(top, "f", "Isopycnals along the LM Gould line, 2004–2025 (Argo)", dy=4)
# g velocity
top = fig.add_subplot(g3R[0]); bot = fig.add_subplot(g3R[1])
for a_ in (top, bot):
    a_.axvspan(XL[0], ice_lat, color="#DCEBF7", lw=0, zorder=0); a_.set_xlim(*XL); a_.set_xticks(XT); a_.axhline(0, color="#9A9A9A", lw=0.5, zorder=1)
    a_.spines["right"].set_visible(True); a_.spines["left"].set_visible(False); a_.yaxis.tick_right(); a_.yaxis.set_label_position("right")
top.plot(P3.latb, P3.ualt_mean, color=C["eke"], lw=0.9, marker="o", ms=1.8, label="altimetry (surface)")
top.plot(P3.latb, P3.uin_mean, color=C["obs"], lw=0.9, marker="D", ms=1.8, label="in situ ADCP, 70–118 m")
top.set_ylim(-5, 45); top.set_yticks([0, 20, 40]); top.set_xticklabels([]); top.set_ylabel("Mean cross-track\nvelocity (cm s$^{-1}$)", fontsize=5.8)
top.legend(loc="upper left", fontsize=4.8)
def prof(ax, x, y, p, col, ls="-", ms=2.4, lw=0.9, z=3):
    ax.plot(x, y, color=col, ls=ls, lw=lw, zorder=z); s = np.array([sigp(v) for v in p])
    ax.plot(np.asarray(x)[s], np.asarray(y)[s], "o", color=col, ms=ms, zorder=z + 1); ax.plot(np.asarray(x)[~s], np.asarray(y)[~s], "o", color=col, mfc="white", ms=ms, mew=0.6, zorder=z + 1)
bot.fill_between(V.latb, V.ci_lo_0125, V.ci_hi_0125, color=C["eke"], alpha=0.18, lw=0, zorder=1)
prof(bot, V.latb, V.slope_0125, V.p_0125, C["eke"]); prof(bot, V.latb, V.slope_0519, V.p_0519, C["eke"], ls=(0, (3, 1.5)), ms=1.9, lw=0.7, z=2)
bot.set_ylim(-8, 9); bot.set_xticklabels(XTL); bot.set_ylabel("Velocity trend\n(cm s$^{-1}$ per decade)", fontsize=5.8)
bot.legend(handles=[Line2D([], [], color=C["eke"], marker="o", ms=2.4, label="altimetry, 2001–2025"),
                    Line2D([], [], color=C["eke"], ls=(0, (3, 1.5)), lw=0.7, marker="o", ms=1.9, label="altimetry, Oct 2005–Apr 2019")], loc="upper left", fontsize=4.8)
S.plabel(top, "g", "Cross-passage velocity: agreement and change", dy=4)
S.save(fig, FIGD, "Figure2_drake"); plt.close(fig)
print("done")
