"""Figure 4: wind stress, EKE and the isopycnal slope (eddy-saturation test; Thompson and Naveira Garabato, 2014).
a  eastward wind-stress trend 2001-2025 (Product 2, 1 deg; 10-3 N m-2 per decade; f3w_wind_map.nc), white dots: 10 % FDR.
b  trend of the signed cross-stream slope of sigma0 = 27.3 kg m-3 (Roemmich-Gilson Argo, monthly, 2004-2025, % of the local mean
   per decade; calibrated monthly test, f4m_map.nc); domain ACC band + 100 km, ice removed; white dots: 10 % FDR. Colour scale
   saturated at +-8 %.
c  wind stress, EKE and slope over the common cells (ACC band + 100 km), anomalies in % of the 2004-2025 mean, 12-month running
   means, with the 2004-2025 trends and 95 % wedges (f4e D1; solid p < 0.05, dashed otherwise), values at the right.
d  fractional-anomaly regression coefficients of KE, EKE and slope on the previous-year wind-stress anomaly (responses 2005-2025,
   wind 2004-2024; f4s_sens.py), bootstrap 95 % bars; filled p < 0.05, open otherwise; grey reference lines at 1 and 0.7.
e  correlation of detrended annual KE, EKE and slope (2004-2025) with the wind stress L = 0-3 yr earlier (f4e D2); markers filled
   where pointwise p < 0.05; grey: pointwise 95 % surrogate envelope (EKE); legend: family-wise p over the four lags.

Publication repository: ported from figures_final/scripts/fig4_tilt.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/f4/f4_tilt_maps.nc, outputs/f4m/f4m_map.nc, outputs/f3w/f3w_wind_map.nc, outputs/f4e/*, outputs/f4s/*,
          outputs/drake/dk12/dk12_band_series_allsat.npz, outputs/base/*; Natural Earth land
Outputs : rendered/{png,pdf}/Figure4
Run     : python figures/fig4_wind_slope.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from osr11 import style as S, paths as P  # noqa: E402
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else P.FIG
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
OUT = P.OUT
M = xr.open_dataset(OUT / "f4" / "f4_tilt_maps.nc"); MM_ = xr.open_dataset(OUT / "f4m" / "f4m_map.nc")
MW = xr.open_dataset(OUT / "f3w" / "f3w_wind_map.nc")
clim = xr.open_dataset(P.CLIM)
ice_ds = xr.open_dataset(P.ICE_MASK); ice_g = ice_ds.ice_mask.values > 0.5
lat8, lon8 = ice_ds.latitude.values, ice_ds.longitude.values
bat = xr.open_dataset(P.BATHY).elevation; bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
Z = np.load(P.OUT / "drake" / "dk12" / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z["lev_saf"]), float(Z["lev_pf"])
ES = pd.read_csv(OUT / "f4e" / "f4e_series.csv"); ER = pd.read_csv(OUT / "f4e" / "f4e_results.csv")
D1 = ER[ER.test == "D1"].set_index("series"); D2 = ER[ER.test == "D2"]
SR = pd.read_csv(OUT / "f4s" / "f4s_results.csv").set_index("series")
CS, CW, CE = C["slope"], C["wind"], C["eke"]


def base_map(ax):
    S.polar(ax, ccrs, cfeature)
    S.add_isobath(ax, bc.longitude, bc.latitude, -bc.values, transform=PC)
    S.add_fronts(ax, clim.longitude, clim.latitude, clim.adt.values, LEV_PF, LEV_SAF, transform=PC)
    S.add_ice(ax, ccrs, lon8, lat8, ice_g)
    ax.contour(M.lon.values, M.lat.values, M["buffer100"].values.astype(float), levels=[0.5], colors=C["buffer"],
               linewidths=S.LW["buffer"], linestyles=[(0, (2, 1.5))], transform=PC, zorder=5)


def trend_map(ax, lon, lat, field, fdr, vmax):
    base_map(ax)
    tc = plt.get_cmap(S.TREND_CMAP).copy(); tc.set_bad(C["masked"])
    im = ax.pcolormesh(S.edges(lon), S.edges(lat), np.ma.masked_invalid(field), transform=PC, cmap=tc, vmin=-vmax, vmax=vmax,
                       shading="flat", rasterized=True, zorder=0)
    jj, ii = np.nonzero(fdr); S.sig_dots(ax, lon[ii], lat[jj], transform=PC)
    return im, len(ii)


H = 226.0; D = 70.0; YM = 145.0
fig = S.figure(S.W2, H)
# a ---------------------------------------------------------------- wind-stress trend
ax = S.axmm(fig, 8, YM, D, D, projection=ccrs.SouthPolarStereo())
imw, nw = trend_map(ax, MW.lon.values, MW.lat.values, 1000 * MW["wind_trend_dec"].values, MW["wind_fdr"].values > 0, 12)
S.letter(fig, 2, H - 5, "a", "Eastward wind-stress trend, 2001–2025")
S.cbar(fig, imw, 8 + 0.1 * D, YM - 10, 0.8 * D, "Wind-stress trend (10$^{-3}$ N m$^{-2}$ per decade)", extend="both", ticks=[-12, -6, 0, 6, 12])
# b ---------------------------------------------------------------- slope trend
ax = S.axmm(fig, 93, YM, D, D, projection=ccrs.SouthPolarStereo())
ims, ns = trend_map(ax, M.lon.values, M.lat.values, MM_["slope_pct"].values, MM_["slope_fdr"].values > 0, 8)
sp = MM_["slope_pct"].values; fd = MM_["slope_fdr"].values > 0
print(f"Fig 4a FDR cells {nw}; Fig 4b FDR cells {ns} (positive {int((sp[fd] > 0).sum())}, negative {int((sp[fd] < 0).sum())})")
S.letter(fig, 86, H - 5, "b", "Isopycnal-slope trend, σ$_0$ = 27.3 kg m$^{-3}$, 2004–2025")
S.cbar(fig, ims, 93 + 0.1 * D, YM - 10, 0.8 * D, "Slope trend (% of mean per decade)", extend="both", ticks=[-8, -4, 0, 4, 8])
fig.legend(handles=S.map_legend_handles(fdr=True, buffer=True), loc="center", bbox_to_anchor=(0.5, 123.0 / H), ncol=5,
           handlelength=1.8, columnspacing=1.4, fontsize=S.FS_ANN)
# c ---------------------------------------------------------------- common-cell series
ax = S.axmm(fig, 16, 76, 144, 36)
xe = ES.x.values; c4 = ES.year.values >= 2004; xc = xe[c4].mean()
LABY = {"wind": 2.6, "EKE": 12.0, "slope": -6.6}
for k, col, nm, z in (("wind", CW, "wind stress", 2), ("EKE", CE, "EKE", 3), ("slope", CS, "isopycnal slope", 4)):
    v = ES[k].values; mu = np.nanmean(v[c4]); pa = 100 * (v / mu - 1); r = D1.loc[k]
    ax.plot(xe[c4], pd.Series(pa).rolling(12, center=True).mean().values[c4], color=col, lw=S.LW["smooth"] + (0.3 if k == "slope" else 0),
            alpha=0.95, zorder=z)
    S.trend_line(ax, 2004, 2026, xc, 0.0, r.pct / 10, r.pct_lo / 10, r.pct_hi / 10, r.p, col, z=5)
    ax.text(2026.6, LABY[k], S.um(f"{nm}\n{r.pct:+.1f} % per decade\n({r.pct_lo:+.1f} to {r.pct_hi:+.1f}), {S.fp(r.p)}"), color=col,
            fontsize=S.FS_ANN, ha="left", va="center", linespacing=1.1)
S.zero(ax, lim=[2003.8, 2026.2]); S.ts_axis(ax, 2003.8, 2026.2, 2031.4)
ax.set_ylabel("Anomaly\n(% of 2004–2025 mean)"); ax.set_ylim(-15, 19); ax.set_yticks([-10, -5, 0, 5, 10, 15]); S.yticks_minus(ax)
S.hgrid(ax, (-10, -5, 5, 10, 15), [2003.8, 2026.2])
S.letter(fig, 2, 117, "c", "Wind stress, EKE and isopycnal slope over the ACC band + 100 km")
# d ---------------------------------------------------------------- regression coefficients
COLS = {"KE": C["ke"], "EKE": CE, "slope": CS}
ax = S.axmm(fig, 16, 13, 50, 46)
for q, k in enumerate(("KE", "EKE", "slope")):
    r = SR.loc[k]; col = COLS[k]
    ax.bar(q, r.elasticity, 0.6, color=col if S.sigp(r.p) else "white", edgecolor=col, lw=0.9, zorder=3)
    ax.plot([q, q], [r.el_lo, r.el_hi], color="#333333", lw=0.7, zorder=4)
    ax.text(q, max(r.el_hi, r.elasticity) + 0.05, S.um(f"{r.elasticity:+.2f}"), ha="center", va="bottom", fontsize=S.FS_ANN)
for yv, ls, c_ in ((1.0, (0, (3, 1.6)), "#777777"), (0.7, (0, (1, 1.2)), "#9A9A9A")):
    ax.axhline(yv, color=c_, lw=0.6, ls=ls, zorder=0)
    ax.text(2.68, yv, f"{yv:g}", ha="left", va="center", fontsize=S.FS_ANN, color=c_)
S.zero(ax)
ax.set_xticks(range(3)); ax.set_xticklabels(["KE", "EKE", "isopycnal\nslope"]); ax.tick_params(axis="x", length=0); ax.set_xlim(-0.6, 2.6)
ax.set_ylabel("Regression coefficient\n(fractional anomalies)"); ax.set_ylim(-0.25, 1.8); ax.set_yticks([0, 0.5, 1.0, 1.5]); S.yticks_minus(ax)
S.letter(fig, 2, 64, "d", "Sensitivity to the wind (lag 1 yr), 2005–2025")
# e ---------------------------------------------------------------- lagged correlations
ax = S.axmm(fig, 101, 13, 65, 46); hl = []
for k, mk, nm in (("KE", "s", "KE"), ("EKE", "o", "EKE"), ("slope", "D", "isopycnal slope")):
    q = D2[D2.series == k].sort_values("lag"); col = COLS[k]
    if k == "EKE": ax.fill_between(q.lag, -q.env95, q.env95, color="#E6E6E6", lw=0, zorder=0)
    ax.plot(q.lag, q.r, color=col, lw=0.9, zorder=3)
    S.profile(ax, q.lag, q.r, q.p, col, mk=mk, ms=3.0, lw=0)
    hl.append(Line2D([], [], color=col, lw=0.9, marker=mk, ms=3.0, label=f"{nm} (family-wise {S.fp(q.p_fw.iloc[0])})"))
hl.append(Patch(facecolor="#E6E6E6", label="95 % surrogate envelope"))
S.zero(ax); ax.set_xticks([0, 1, 2, 3]); ax.set_xlabel("Lag (yr, wind leading)"); ax.set_xlim(-0.2, 3.2)
ax.set_ylabel("r with wind stress\n(detrended annual)"); ax.set_ylim(-0.6, 1.05); ax.set_yticks([-0.5, 0, 0.5, 1.0]); S.yticks_minus(ax)
ax.legend(handles=hl, loc="upper right", fontsize=S.FS_ANN, handlelength=1.8)
S.letter(fig, 86, 64, "e", "Lagged correlation with the wind, 2004–2025")
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("D1:", D1.loc[["wind", "EKE", "slope"], ["pct", "pct_lo", "pct_hi", "p"]].round(3).to_dict("index"))
print("f4s:", SR[["elasticity", "el_lo", "el_hi", "p"]].round(3).to_dict("index"))
