"""Figure 2: Drake Passage (86-40 W), Product 1.
a  EKE trend 2001-2025 per 0.25 deg cell (2 x 2 average of the 1/8 deg grid; least-squares trend of annual means; dk34), mean SAF
   and PF (2001-2025 mean ADT contours, the band of c), 3000 m isobath, LM Gould track (dashed where seasonally ice-covered).
   No field significance is computed for this map.
b  half-power occurrence 2001-2025 (% of months per cell; dk10) with the smoothed 10 % contours of 2001-2015 (cyan) and 2016-2025
   (black dashed). Note: overlap of the two enclosed areas (intersection over union, cos-weighted, computed here) and the
   2016-2025 minus 2001-2015 displacement of the half-power latitude and of the fixed ADT contour through the mean half-power
   point, 70-55 W (x1_drake_period_shift.json, written by 03t_contour_vs_energy_drake.py); the 70-55 W window is marked at the bottom.
c  EKE trend between the mean SAF and PF in 2 deg longitude bins (calibrated Theil-Sen, 95 % band; dk8), mean seafloor depth
   between the fronts at 1/8 deg (dark fill, right axis); dashed lines: the narrow passage, 68-58 W.
d  trend of the surface cross-track velocity per pair of SSH streamlines along the Gould line, 2001-2025 (dk5 streamwise bins,
   calibrated Theil-Sen, 95 % bars), x = mean ADT at the centre of each pair; SAF and PF zones shaded (dk5_stream_edges).
       json: outputs/x1/x1_drake_period_shift.json (03t)

Publication repository: ported from figures_final/scripts/fig2_drake.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/drake/{dk34,dk10,dk8,dk5,dk24}, outputs/d1/d1b_track.csv, outputs/x1/x1_drake_period_shift.json, outputs/base/*
Outputs : rendered/{png,pdf}/Figure2
Run     : python figures/fig2_drake.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from osr11 import style as S, paths as P  # noqa: E402
from matplotlib.colors import TwoSlopeNorm
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else P.FIG
PSJ = P.OUT / "x1" / "x1_drake_period_shift.json"   # 03t_contour_vs_energy_drake.py
plt = S.plt; C = S.C
DR = WO = P.OUT / "drake"
F34 = np.load(WO / "dk34" / "dk34_drake_fields.npz"); DN = np.load(WO / "dk10" / "dk10_hp_density_p2.npz")
AL = pd.read_csv(DR / "dk8" / "dk8_along_acc_fronts_allsat.csv")
trk = pd.read_csv(P.OUT / "d1" / "d1b_track.csv")
box = dict(longitude=slice(-87, -33), latitude=slice(-66, -46))
bs = xr.open_dataset(P.BATHY).elevation.sel(**box)
cl = xr.open_dataset(P.CLIM).sel(**box)
LEV_PF, LEV_SAF = float(F34["lev_pf"]), float(F34["lev_saf"])
EXT = (-86, -40, -65, -50); PASS = (-68, -58); WIN = (-70, -55)
PS = json.load(open(PSJ))["7055"]


def base(ax, ylab=True, fronts=True, track=True):
    S.add_isobath(ax, bs.longitude, bs.latitude, -bs.values)
    ax.contourf(bs.longitude, bs.latitude, bs.values, levels=[0, 1e5], colors=[C["land"]], zorder=6)
    if fronts:
        S.add_fronts(ax, cl.longitude, cl.latitude, cl.adt.values, LEV_PF, LEV_SAF)
        for y_, t_ in ((-55.6, "SAF"), (-60.7, "PF")):
            ax.text(-85.4, y_, t_, fontsize=S.FS_ANN, fontweight="bold", zorder=8)
    if track:
        ax.plot(trk.lon[trk.perm], trk.lat[trk.perm], color="#00A651", lw=1.4, zorder=8, solid_capstyle="butt")
        ax.plot(trk.lon[~trk.perm], trk.lat[~trk.perm], color="#00A651", lw=1.0, ls=(0, (2, 1.2)), zorder=8)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
    ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels([S.lonfmt(v) for v in (-80, -70, -60, -50, -40)])
    ax.set_yticks([-65, -60, -55, -50]); ax.set_yticklabels([S.latfmt(v) for v in (-65, -60, -55, -50)] if ylab else [])
    ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)


H = 146.0                      # 5 Nov 2026: panel e (Qs) removed
fig = S.figure(S.W2, H)
MW, MH = 74.0, 74.0 * 15 * (1 / np.cos(np.radians(57))) / 46          # map aspect of the previous version (1/cos 57 deg)
YT = H - 10.5 - MH
# a ---------------------------------------------------------------- EKE trend map
ax = S.axmm(fig, 11, YT, MW, MH)
tc = plt.get_cmap(S.TREND_CMAP).copy(); tc.set_bad(C["masked"])
im = ax.pcolormesh(F34["slon"], F34["slat"], np.ma.masked_invalid(F34["slope_map"]), cmap=tc,
                   norm=TwoSlopeNorm(0, -S.EKE_TREND_LIM, S.EKE_TREND_LIM), shading="auto", rasterized=True, zorder=1)
base(ax)
for x_, y_, t_ in ((-71.5, -64.2, "Drake Passage"), (-46.5, -61.0, "Scotia Sea")):
    ax.text(x_, y_, t_, fontsize=S.FS_ANN, ha="center", style="italic", color="#333333", zorder=9)
S.letter(fig, 2, H - 5, "a", "EKE trend, 2001–2025")
S.cbar(fig, im, 11 + 0.1 * MW, YT - 11.5, 0.8 * MW, "EKE trend (cm$^2$ s$^{-2}$ per decade)", extend="both", ticks=S.EKE_TREND_TICKS)
# b ---------------------------------------------------------------- composite across the pathway: annual KE and ADT profiles (5 Nov 2026)
# Retained tracks: monthly KE and ADT sampled on each track's longitude within +-3.5 deg of its mean half-power latitude
# (outputs/x2/x2_composite_*.npz, from the monthly fields), averaged over the tracks and over each year 2001-2025; ADT is referenced
# to each track's fixed contour level, so zero marks the fixed ADT contour. Markers: annual means of the monthly half-power latitude
# (triangles, on the KE curves) and of the monthly contour latitude (circles, on the zero line), from the x1 track files.
# 5 Nov 2026: tracks of the longitude-dependent threshold (x6/x7), 86-40 W subset
X7 = P.OUT / "x7"; ZT = np.load(X7 / "x7_tracks.npz"); _m = (ZT["lon"] >= -86) & (ZT["lon"] <= -40)
TD = pd.read_csv(X7 / "x7_trends.csv"); TD = TD[TD.region == "Drake region 86-40W"].set_index("series")
CP = np.load(X7 / "x7_composite_annual_drake.npz")   # annual composites reduced from the monthly ones by scripts/02_circumpolar/02p_pathway_composites.py
HP, LC = ZT["lat_hp"][:, _m], ZT["lat_c"][:, _m]; TR_ = TD; NTR = int(_m.sum())
BX, BW, BY, BH = 103, 63, YT, MH; LETX = 85
KMD = 111.195; yrX = np.repeat(np.arange(2001, 2026), 12); YRS = np.arange(2001, 2026)
xk = CP["offset_deg"] * KMD; ke_y, ad_y = CP["ke_y"], CP["ad_y"]   # annual means over months and tracks: KE (cm2 s-2), ADT minus contour level (cm)
def _ann(a_):
    an = (a_ - np.nanmean(a_, 0)) * KMD; return np.array([np.nanmean(np.nanmean(an[yrX == y], 0)) for y in YRS])
hp_y, ct_y = _ann(HP), _ann(LC)
_e, _l = YRS <= 2017, YRS >= 2018   # 5 Nov 2026: same split as the EKE timing (SI Fig. S3)
rise = np.interp(0, xk, ad_y[_l].mean(0)) - np.interp(0, xk, ad_y[_e].mean(0)); grad = 10 * (np.interp(10, xk, ad_y.mean(0)) - np.interp(-10, xk, ad_y.mean(0))) / 20
print(f"Fig b composite (drake): KE peak {ke_y[_e].max(1).mean():.0f} -> {ke_y[_l].max(1).mean():.0f}; ADT rise 2018-2025 minus 2001-2017 {rise:+.1f} cm; gradient {grad:.2f} cm/10 km; "
      f"half-power period means {hp_y[_e].mean():+.1f}/{hp_y[_l].mean():+.1f}, contour {ct_y[_e].mean():+.1f}/{ct_y[_l].mean():+.1f} km")
from matplotlib.colors import LinearSegmentedColormap as _LSC
YC = _LSC.from_list("years", plt.get_cmap("plasma")(np.linspace(0.0, 0.85, 64))); ycol = lambda y: YC((y - 2001) / 24); XL = 260   # plasma truncated: dark blue (2001) to orange (2025); XL = 260
import matplotlib as _mpl
ax = S.axmm(fig, BX, BY + 0.56 * BH, BW, 0.44 * BH)
for k_, y in enumerate(YRS): ax.plot(xk, ke_y[k_], color=ycol(y), lw=0.7, alpha=0.9, zorder=3)
ax.axvline(0, color="#BDBDBD", lw=0.5, zorder=0)
for k_, y in enumerate(YRS): ax.plot(hp_y[k_], np.interp(hp_y[k_], xk, ke_y[k_]), "v", color=ycol(y), ms=3.0, mec="white", mew=0.3, zorder=5, clip_on=False)
ax.text(0.03, 0.95, "▼ half-power latitude", transform=ax.transAxes, ha="left", va="top", fontsize=S.FS_ANN - 0.5)   # 5 Nov 2026: less text
ax.set_xlim(-XL, XL); S.xticks_minus(ax); ax.tick_params(labelbottom=False); ax.set_ylabel("KE (cm$^2$ s$^{-2}$)", labelpad=2); ax.set_ylim(0, 1.18 * ke_y.max())
cax = fig.add_axes([(BX + BW - 22) / S.W2, (BY + 0.93 * BH) / H, 20 / S.W2, 1.6 / H])
cb = _mpl.colorbar.ColorbarBase(cax, cmap=YC, norm=_mpl.colors.Normalize(2001, 2025), orientation="horizontal", ticks=[2001, 2013, 2025]); cb.ax.tick_params(labelsize=S.FS_ANN - 0.5, length=1.5, pad=1); cb.outline.set_linewidth(0.4)
cax.set_title("annual means", fontsize=S.FS_ANN - 0.5, pad=1.5)
ax2 = S.axmm(fig, BX, BY, BW, 0.50 * BH)
ax2.axhline(0, color="#7F7F7F", lw=0.9, zorder=1); ax2.axvline(0, color="#BDBDBD", lw=0.5, zorder=0)
for k_, y in enumerate(YRS): ax2.plot(xk, ad_y[k_], color=ycol(y), lw=0.7, alpha=0.9, zorder=3)
for k_, y in enumerate(YRS): ax2.plot(ct_y[k_], 0, "o", color=ycol(y), ms=3.0, mec="white", mew=0.3, zorder=5)
ax2.text(0.03, 0.95, "● fixed ADT contour", transform=ax2.transAxes, ha="left", va="top", fontsize=S.FS_ANN - 0.5)
ax2.set_xlim(-XL, XL); ax2.set_ylim(-48, 48); ax2.set_yticks([-40, -20, 0, 20, 40]); S.yticks_minus(ax2); S.xticks_minus(ax2)
ax2.set_ylabel("ADT − contour level (cm)", labelpad=2); ax2.set_xlabel("Distance from the mean half-power point (km, + north)", labelpad=1.5)
rh, rc = TR_.loc["lat_hp"], TR_.loc["lat_contour"]
S.letter(fig, LETX, H - 5, "b", f"KE and ADT across the energetic pathway, 86–40°W ({NTR} tracks)")
# c ---------------------------------------------------------------- EKE trend between the fronts, 2 deg bins
RB, RH = 13.0, 50.0; RT = RB + RH
ax = S.axmm(fig, 16, RB, 60, RH)
_ab = cl.adt.interp(latitude=bs.latitude, longitude=bs.longitude).values
_band = (_ab >= LEV_PF) & (_ab <= LEV_SAF) & (bs.values < 0)
with np.errstate(all="ignore"): BDEP = np.nanmean(np.where(_band, -bs.values, np.nan), 0) / 1000
BLON = bs.longitude.values; _k = (BLON >= -86) & (BLON <= -40)
axd = ax.twinx(); axd.fill_between(BLON[_k], BDEP[_k], 7, color=C["land"], lw=0, zorder=0)
axd.set_ylim(6.2, -6); axd.set_yticks([3, 4, 5, 6]); axd.set_ylabel("Depth (km)", color="#333333")
axd.yaxis.set_label_coords(1.09, 0.16); axd.tick_params(axis="y", colors="#333333")
axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5); axd.spines["right"].set_bounds(6.2, 2.5)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
ax.vlines(AL.lon, AL.lo, AL.hi, color=C["eke"], lw=0.6, alpha=0.8, zorder=2)   # 5 Nov 2026: error bars instead of shading
S.profile(ax, AL.lon, AL.slope_dec, AL.p, C["eke"])
S.zero(ax); ax.set_xlim(-86, -40); ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels([S.lonfmt(v) for v in (-80, -70, -60, -50, -40)])
ax.set_xlabel("Longitude (flow from left to right)"); ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ per decade)")
ax.set_ylim(-115, 112); ax.set_yticks([-50, -25, 0, 25, 50, 75, 100]); ax.spines["left"].set_bounds(-50, 112); S.yticks_minus(ax)
for x_ in PASS: ax.axvline(x_, color="#BDBDBD", lw=0.5, ls=(0, (2, 2)))
ax.text(np.mean(PASS), 0.995, "Drake Passage", ha="center", va="top", fontsize=S.FS_ANN, transform=ax.get_xaxis_transform())
S.flow_arrow(ax, -47.5, -41.5, 0.86)
S.letter(fig, 2, RT + 5, "c", "EKE trend between the SAF and PF (2° bins)")
# d ---------------------------------------------------------------- cross-track velocity trend per streamline pair
K5 = pd.read_csv(WO / "dk5" / "dk5_trends.csv"); A24 = pd.read_csv(WO / "dk24" / "dk24_stream_trends.csv")
adt_of = dict(zip(A24.psic.round(3), A24.adt_c)); ED = pd.read_csv(WO / "dk5" / "dk5_stream_edges.csv"); ed = dict(zip(ED.psi.round(3), ED.adt_edge))
ZONES = [("SAF", ed[-0.6], ed[-0.4]), ("PF", ed[-1.25], ed[-1.0])]
qu = K5[(K5.coord == "stream") & (K5.src == "allsat") & (K5["var"] == "u") & (K5.period == "2001-2025")].copy()
qu["adt"] = qu.bin.round(3).map(adt_of); qu = qu.dropna(subset=["adt"]).sort_values("adt")
ax = S.axmm(fig, 104, RB, 62, RH); CU = C["vel"]
for zn, lo_, hi_ in ZONES:
    ax.axvspan(lo_, hi_, color="#EDEDED", lw=0, zorder=0)
    ax.text((lo_ + hi_) / 2, 1.0, zn, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=S.FS_ANN, fontweight="bold")
S.zero(ax); ax.set_xlim(-1.16, 0.44)
_nar = ED.adt_edge.values[(ED.psi.values <= -1.45)]; ax.axvspan(_nar.min(), _nar.max(), color="#E3E3E3", lw=0, zorder=0)   # 0.05 m bins (4 Nov 2026)
ax.text(_nar.mean(), -7.6, "0.05 m bins", ha="center", va="bottom", fontsize=S.FS_ANN, color=C["grey"])
ax.vlines(qu.adt, 10 * qu.ci_lo, 10 * qu.ci_hi, color=CU, lw=0.6, zorder=2)
S.profile(ax, qu.adt, 10 * qu.slope, qu.p, CU, mk="s", ms=2.4)
ax.set_ylim(-8, 8); ax.set_yticks([-8, -6, -4, -2, 0, 2, 4, 6, 8]); S.yticks_minus(ax); S.xticks_minus(ax)
ax.set_ylabel("Cross-track velocity trend\n(cm s$^{-1}$ per decade)")
ax.set_xlabel("Mean ADT at the bin centre (m)\nsouth ←         → north", linespacing=1.3)
S.letter(fig, 85, RT + 5, "d", "Cross-track velocity trend in ADT bins, 2001–2025")

S.save(fig, FIGD, "Figure2"); plt.close(fig)
