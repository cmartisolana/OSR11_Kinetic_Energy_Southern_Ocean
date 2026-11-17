"""Figure 3: circumpolar EKE and the energetic pathway, Product 1, 2001-2025.
a  EKE trend (cm2 s-2 per decade, 1 deg cells; a1_maps_bc), no significance dots (10 % FDR dots removed 10 Nov 2026); mean SAF and PF, 3000 m isobath, sea ice.
b  (5 Nov 2026) circumpolar-mean latitude anomaly of the half-power point and of the fixed ADT contour through its mean position
   (900 tracks, x1) with their trends; inset: the two trends with 95 % intervals. The former occurrence map is SI Fig. S7b.
c  ACC band (mean ADT -0.6 to +0.2 m): EKE and KE, monthly (thin) and 12-month running means (thick), 2001-2025 trends with 95 %
   wedges (revised calibrated test, retest_summary allsat|ACC_*|2001-2025; solid: p < 0.05), values at the right.
d  EKE trend between the mean SAF and PF in 2 deg longitude bins (dk12, bars 95 %), mean seafloor depth between the fronts
   (1/8 deg, dark fill, right axis); keys of the bathymetric features as on the maps; light shading: the window of Fig. 2c.
Rim keys on the maps (style.KEYS) mark the same longitudes as the keys of d.

Publication repository: ported from figures_final/scripts/fig3_circumpolar.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/a1/a1_maps_bc.nc, outputs/drake/{dk14,dk12}, outputs/a5/a5_series_allsat.csv, outputs/retest/retest_summary.csv,
          outputs/base/*; Natural Earth land (cartopy)
Outputs : rendered/{png,pdf}/Figure3
Run     : python figures/fig3_circumpolar.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from osr11 import style as S, paths as P  # noqa: E402
import cartopy.crs as ccrs, cartopy.feature as cfeature
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else P.FIG
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
OUT = P.OUT
bat = xr.open_dataset(P.BATHY).elevation; bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
clim = xr.open_dataset(P.CLIM)
ice_ds = xr.open_dataset(P.ICE_MASK); ice_g = ice_ds.ice_mask.values > 0.5
lat8, lon8 = ice_ds.latitude.values, ice_ds.longitude.values
RT = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")
M = xr.open_dataset(P.OUT / "a1" / "a1_maps_bc.nc")
DN = np.load(P.OUT / "drake" / "dk14" / "dk14_hp_density_circ.npz")
DK = P.OUT / "drake" / "dk12"; TB2 = pd.read_csv(DK / "dk12_band_trends_allsat.csv").sort_values("lon")
Z = np.load(DK / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z["lev_saf"]), float(Z["lev_pf"])


def base_map(ax, fronts=True):
    S.polar(ax, ccrs, cfeature)
    S.add_isobath(ax, bc.longitude, bc.latitude, -bc.values, transform=PC)
    if fronts: S.add_fronts(ax, clim.longitude, clim.latitude, clim.adt.values, LEV_PF, LEV_SAF, transform=PC)
    S.add_ice(ax, ccrs, lon8, lat8, ice_g)


H = 184.0; D = 68.0; YM = 108.0   # 5 Nov 2026: panel e moved to SI Fig. S8; a and b enlarged
fig = S.figure(S.W2, H)
# a ---------------------------------------------------------------- EKE trend map
ax = S.axmm(fig, 8, YM, D, D, projection=ccrs.SouthPolarStereo()); base_map(ax)
tc = plt.get_cmap(S.TREND_CMAP).copy(); tc.set_bad(C["masked"])
im = ax.pcolormesh(S.edges(M.lon.values), S.edges(M.lat.values), np.ma.masked_invalid(10 * M["EKE_slope"].values), transform=PC, cmap=tc,
                   vmin=-S.EKE_TREND_LIM, vmax=S.EKE_TREND_LIM, shading="flat", rasterized=True, zorder=0)
jj, ii = np.nonzero(M["EKE_fdr"].values == 1); print("Fig 3a FDR cells:", len(ii))
# S.sig_dots(...) removed 10 Nov 2026 at the author's request (no 10 % FDR dots in Fig. 3a)
S.letter(fig, 2, H - 5, "a", "EKE trend, 2001–2025")
S.cbar(fig, im, 8 + 0.1 * D, YM - 10, 0.8 * D, "EKE trend (cm$^2$ s$^{-2}$ per decade)", extend="both", ticks=S.EKE_TREND_TICKS)
# b ---------------------------------------------------------------- composite across the pathway: annual KE and ADT profiles (5 Nov 2026)
# Retained tracks: monthly KE and ADT sampled on each track's longitude within +-3.5 deg of its mean half-power latitude
# (outputs/x2/x2_composite_*.npz, from the monthly fields), averaged over the tracks and over each year 2001-2025; ADT is referenced
# to each track's fixed contour level, so zero marks the fixed ADT contour. Markers: annual means of the monthly half-power latitude
# (triangles, on the KE curves) and of the monthly contour latitude (circles, on the zero line), from the x1 track files.
# 5 Nov 2026: tracks of the longitude-dependent threshold (x6/x7)
X7 = OUT / "x7"; ZT = np.load(X7 / "x7_tracks.npz"); TRX = pd.read_csv(X7 / "x7_trends.csv"); TRX = TRX[TRX.region == "circumpolar"].set_index("series")
CP = np.load(X7 / "x7_composite_annual_circ.npz")   # annual composites reduced from the monthly ones by scripts/02_circumpolar/02p_pathway_composites.py
HP, LC = ZT["lat_hp"], ZT["lat_c"]; TR_ = TRX; NTR = ZT["lat_hp"].shape[1]
BX, BW, BY, BH = 100, 64, YM, D; LETX = 86
KMD = 111.195; yrX = np.repeat(np.arange(2001, 2026), 12); YRS = np.arange(2001, 2026)
xk = CP["offset_deg"] * KMD; ke_y, ad_y = CP["ke_y"], CP["ad_y"]   # annual means over months and tracks: KE (cm2 s-2), ADT minus contour level (cm)
def _ann(a_):
    an = (a_ - np.nanmean(a_, 0)) * KMD; return np.array([np.nanmean(np.nanmean(an[yrX == y], 0)) for y in YRS])
hp_y, ct_y = _ann(HP), _ann(LC)
_z = lambda v: f"{(0.0 if abs(v) < 0.5 else v):+.0f}"   # avoid "-0"
_e, _l = YRS <= 2017, YRS >= 2018   # 5 Nov 2026: same split as the EKE timing (SI Fig. S3)
rise = np.interp(0, xk, ad_y[_l].mean(0)) - np.interp(0, xk, ad_y[_e].mean(0)); grad = 10 * (np.interp(10, xk, ad_y.mean(0)) - np.interp(-10, xk, ad_y.mean(0))) / 20
print(f"Fig b composite (circ): KE peak {ke_y[_e].max(1).mean():.0f} -> {ke_y[_l].max(1).mean():.0f}; ADT rise 2018-2025 minus 2001-2017 {rise:+.1f} cm; gradient {grad:.2f} cm/10 km; "
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
S.letter(fig, LETX, H - 5, "b", f"KE and ADT across the energetic pathway, circumpolar ({NTR} tracks)")
fig.legend(handles=S.map_legend_handles(fdr=False, occ=False), loc="center", bbox_to_anchor=((8 + D / 2) / S.W2, (YM - 21) / H), ncol=3,
           handlelength=1.6, columnspacing=1.1, fontsize=S.FS_ANN)
# c ---------------------------------------------------------------- ACC-band series
A5 = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv"); tA = A5.x.values; LABC = {"KE": 300, "EKE": 252}
ax = S.axmm(fig, 16, 56, 144, 22)
for v, lab, col in (("ACC_EKE", "EKE", C["eke"]), ("ACC_KE", "KE", C["ke"])):
    y = A5[v].values; r = RT.loc[f"allsat|{v}|2001-2025"]
    ax.plot(tA, y, color=col, lw=S.LW["raw"], alpha=0.45, zorder=2)
    ax.plot(tA, S.runfull(y), color=col, lw=S.LW["smooth"] + 0.1, zorder=3)
    S.trend_line(ax, 2001, 2026, tA.mean(), y.mean(), r.slope, r.ci_lo, r.ci_hi, r.p_new, col)
    _m = 10 * r.slope / (r.pct_per_decade / 100)
    ax.text(2026.6, LABC[lab], S.um(f"{lab}\n{r.pct_per_decade:+.1f} % per decade\n({1000 * r.ci_lo / _m:+.1f} to {1000 * r.ci_hi / _m:+.1f}), {S.fp(r.p_new)}"),
            color=col, fontsize=S.FS_ANN, ha="left", va="center", linespacing=1.1)
S.ts_axis(ax, 2000.8, 2026.2, 2032.6); ax.set_ylabel("cm$^2$ s$^{-2}$")
_all = np.r_[A5.ACC_EKE.values, A5.ACC_KE.values]; _lo, _hi = np.nanmin(_all), np.nanmax(_all); _pad = 0.06 * (_hi - _lo)
ax.set_ylim(_lo - _pad, _hi + _pad); ax.set_yticks(np.arange(np.ceil(_lo / 25) * 25, _hi + _pad, 25))
S.hgrid(ax, ax.get_yticks(), [2000.8, 2026.2])
S.letter(fig, 2, 82, "c", "EKE and KE averaged over the ACC band")
# d ---------------------------------------------------------------- EKE trends by longitude between the fronts
ax = S.axmm(fig, 16, 10, 144, 34); axd = ax.twinx()
axd.axvspan(-86, -40, color="#EEEEEE", lw=0, zorder=-1)
_ab = clim.adt.interp(latitude=bat.latitude, longitude=bat.longitude).values
_band = (_ab >= LEV_PF) & (_ab <= LEV_SAF) & (bat.values < 0)
with np.errstate(all="ignore"): BDEP = np.nanmean(np.where(_band, -bat.values, np.nan), 0) / 1000
axd.fill_between(bat.longitude.values, BDEP, 7, color=C["land"], lw=0, zorder=0)
axd.set_ylim(6.2, -6.8); axd.set_yticks([2, 3, 4, 5, 6]); axd.set_ylabel("Depth (km)", color="#333333")
axd.yaxis.set_label_coords(1.045, 0.17); axd.tick_params(axis="y", colors="#333333")
axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5); axd.spines["right"].set_bounds(6.2, 1.2)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
T2 = TB2.dropna(subset=["slope_dec"])
ax.vlines(T2.lon, T2.lo, T2.hi, color=C["eke"], lw=0.6, alpha=0.8, zorder=2)   # 5 Nov 2026: error bars instead of shading
S.profile(ax, T2.lon, T2.slope_dec, T2.p, C["eke"], ms=1.9, lw=S.LW["profile"])
S.zero(ax); S.lon_axis(ax); ax.set_ylim(-380, 222); ax.set_yticks([-100, 0, 100, 200]); ax.spines["left"].set_bounds(-110, 222); S.yticks_minus(ax)
ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ per decade)")
S.keys_on_axis(ax); S.flow_arrow(ax, 166, 179, 0.90)
S.letter(fig, 2, 48, "d", "EKE trend between the SAF and PF (2° bins); mean seafloor depth (bottom)")
S.save(fig, FIGD, "Figure3"); plt.close(fig)
print("values: KE", RT.loc["allsat|ACC_KE|2001-2025", ["pct_per_decade", "p_new"]].tolist(), "EKE", RT.loc["allsat|ACC_EKE|2001-2025", ["pct_per_decade", "p_new"]].tolist(),
      "positive bins", int((T2.slope_dec > 0).sum()), "of", len(T2))
