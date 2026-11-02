"""Figure 3: circumpolar EKE and the energetic pathway, Product 1, 2001-2025.
a  EKE trend (cm2 s-2 per decade, 1 deg cells; a1_maps_bc), white dots: 10 % FDR; mean SAF and PF, 3000 m isobath, sea ice.
b  half-power occurrence 2001-2025 (dk14) with the smoothed 10 % contours of 2001-2015 (cyan) and 2016-2025 (black dashed);
   note: overlap of the two enclosed areas (intersection over union, cos-weighted, computed here).
c  ACC band (mean ADT -0.6 to +0.2 m): EKE and KE, monthly (thin) and 12-month running means (thick), 2001-2025 trends with 95 %
   wedges (revised calibrated test, retest_summary allsat|ACC_*|2001-2025; solid: p < 0.05), values at the right.
d  EKE trend between the mean SAF and PF in 2 deg longitude bins (dk12, shading 95 %), mean seafloor depth between the fronts
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
from scipy.ndimage import gaussian_filter
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


H = 220.0; D = 70.0; YM = 139.0
fig = S.figure(S.W2, H)
# a ---------------------------------------------------------------- EKE trend map
ax = S.axmm(fig, 8, YM, D, D, projection=ccrs.SouthPolarStereo()); base_map(ax)
tc = plt.get_cmap(S.TREND_CMAP).copy(); tc.set_bad(C["masked"])
im = ax.pcolormesh(S.edges(M.lon.values), S.edges(M.lat.values), np.ma.masked_invalid(10 * M["EKE_slope"].values), transform=PC, cmap=tc,
                   vmin=-S.EKE_TREND_LIM, vmax=S.EKE_TREND_LIM, shading="flat", rasterized=True, zorder=0)
jj, ii = np.nonzero(M["EKE_fdr"].values == 1); print("Fig 3a FDR cells:", len(ii))
S.sig_dots(ax, M.lon.values[ii], M.lat.values[jj], transform=PC)
S.letter(fig, 2, H - 5, "a", "EKE trend, 2001–2025")
S.cbar(fig, im, 8 + 0.1 * D, YM - 10, 0.8 * D, "EKE trend (cm$^2$ s$^{-2}$ per decade)", extend="both", ticks=S.EKE_TREND_TICKS)
# b ---------------------------------------------------------------- half-power occurrence
ax = S.axmm(fig, 93, YM, D, D, projection=ccrs.SouthPolarStereo()); base_map(ax, fronts=False)
dd = np.where(DN["all"] > 0, DN["all"], np.nan)
lc = 0.5 * (DN["lon_edges"][1:] + DN["lon_edges"][:-1]); ac = 0.5 * (DN["lat_edges"][1:] + DN["lat_edges"][:-1])
imd = ax.pcolormesh(DN["lon_edges"], DN["lat_edges"], dd, cmap=S.OCC_CMAP, vmin=0, vmax=S.OCC_LIM, shading="flat", transform=PC,
                    rasterized=True, zorder=0)
SM = {k: gaussian_filter(DN[k], 1.0) for k in ("early", "late")}
for key, col, ls, lw in (("early", C["early"], "-", S.LW["early"]), ("late", C["late"], (0, (2.5, 1.2)), S.LW["late"])):
    ax.contour(lc, ac, SM[key], levels=[S.LEVB], colors=col, linewidths=lw, linestyles=[ls], transform=PC, zorder=4 if key == "early" else 5)
_me, _ml = SM["early"] >= S.LEVB, SM["late"] >= S.LEVB; _w = np.cos(np.radians(ac))[:, None] * np.ones_like(_me, float)
IOU = (_w * (_me & _ml)).sum() / (_w * (_me | _ml)).sum(); print(f"Fig 3b: 10 % areas IoU = {IOU:.3f}")
ax.text(1.0, -0.02, f"area overlap\n{100 * IOU:.0f} %", transform=ax.transAxes, ha="right", va="bottom", fontsize=S.FS_ANN, linespacing=1.2)
S.letter(fig, 86, H - 5, "b", "Half-power occurrence, 2001–2025")
S.cbar(fig, imd, 93 + 0.1 * D, YM - 10, 0.8 * D, "Half-power points (% of months per cell)", extend="max", ticks=S.OCC_TICKS)
# shared map legend, one row below the colourbars
fig.legend(handles=S.map_legend_handles(fdr=True, occ=True), loc="center", bbox_to_anchor=(0.5, 117.0 / H), ncol=6,
           handlelength=1.6, columnspacing=1.1, fontsize=S.FS_ANN)
# c ---------------------------------------------------------------- ACC-band series
A5 = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv"); tA = A5.x.values; LABC = {"KE": 300, "EKE": 252}
ax = S.axmm(fig, 16, 70, 144, 36)
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
S.letter(fig, 2, 111, "c", "EKE and KE averaged over the ACC band")
# d ---------------------------------------------------------------- trends by longitude between the fronts
ax = S.axmm(fig, 16, 13, 144, 44); axd = ax.twinx()
axd.axvspan(-86, -40, color="#EEEEEE", lw=0, zorder=-1)                   # window of Fig. 2c, behind the depth fill
_ab = clim.adt.interp(latitude=bat.latitude, longitude=bat.longitude).values
_band = (_ab >= LEV_PF) & (_ab <= LEV_SAF) & (bat.values < 0)
with np.errstate(all="ignore"): BDEP = np.nanmean(np.where(_band, -bat.values, np.nan), 0) / 1000
axd.fill_between(bat.longitude.values, BDEP, 7, color=C["land"], lw=0, zorder=0)
axd.set_ylim(6.2, -6.8); axd.set_yticks([2, 3, 4, 5, 6]); axd.set_ylabel("Depth (km)", color="#333333")
axd.yaxis.set_label_coords(1.045, 0.17); axd.tick_params(axis="y", colors="#333333")
axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5); axd.spines["right"].set_bounds(6.2, 1.2)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
T2 = TB2.dropna(subset=["slope_dec"])
ax.fill_between(T2.lon, T2.lo, T2.hi, color=C["eke"], alpha=0.30, lw=0, zorder=1)
S.profile(ax, T2.lon, T2.slope_dec, T2.p, C["eke"], ms=1.9, lw=S.LW["profile"])
S.zero(ax); S.lon_axis(ax)
ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ per decade)"); ax.set_ylim(-240, 222); ax.set_yticks([-100, -50, 0, 50, 100, 150, 200])
ax.spines["left"].set_bounds(-100, 222); S.yticks_minus(ax)
S.keys_on_axis(ax); S.flow_arrow(ax, 166, 179, 0.80)
S.letter(fig, 2, 62, "d", "EKE trend between the SAF and PF (2° bins); mean seafloor depth (bottom)")
S.save(fig, FIGD, "Figure3"); plt.close(fig)
print("values: KE", RT.loc["allsat|ACC_KE|2001-2025", ["pct_per_decade", "p_new"]].tolist(), "EKE", RT.loc["allsat|ACC_EKE|2001-2025", ["pct_per_decade", "p_new"]].tolist(),
      "positive bins", int((T2.slope_dec > 0).sum()), "of", len(T2))
