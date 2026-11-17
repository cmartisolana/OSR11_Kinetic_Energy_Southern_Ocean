"""Figure 4: wind stress, EKE and the isopycnal slope (eddy-saturation test; Thompson and Naveira Garabato, 2014).
a  eastward wind-stress trend 2001-2025 (Product 2, 1 deg; 10-3 N m-2 per decade; f3w_wind_map.nc), no significance dots (removed 4 Nov 2026).
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
# 9 Nov 2026: panel d = power-law exponents (slope of ln X on ln tau of the previous year), response 2004-2025 (x9), with the
# simulation-calibrated 95 % intervals of x10 (trend-retained and detrended); p-values from x9 (phase-randomised wind surrogates)
_x9 = pd.read_csv(OUT / "x9" / "x9_sensitivity_log.csv").set_index("series")
X8 = pd.read_csv(OUT / "x10" / "x10_intervals.csv").set_index("series"); X8["p"] = _x9["p"]
SR = X8.rename(columns={"b": "elasticity", "b_lo": "el_lo", "b_hi": "el_hi"})
CS, CW, CE = C["slope"], C["wind"], C["eke"]


def base_map(ax):
    S.polar(ax, ccrs, cfeature)
    S.add_isobath(ax, bc.longitude, bc.latitude, -bc.values, transform=PC)
    S.add_fronts(ax, clim.longitude, clim.latitude, clim.adt.values, LEV_PF, LEV_SAF, transform=PC)
    S.add_ice(ax, ccrs, lon8, lat8, ice_g)
    ax.contour(M.lon.values, M.lat.values, M["buffer100"].values.astype(float), levels=[0.5], colors=C["buffer"],
               linewidths=S.LW["buffer"], linestyles=[(0, (2, 1.5))], transform=PC, zorder=5)


def trend_map(ax, lon, lat, field, fdr, vmax, bad=C["masked"], cmap=S.TREND_CMAP):
    base_map(ax)
    tc = plt.get_cmap(cmap).copy(); tc.set_bad(bad)
    field = np.where(np.isfinite(field), field, np.nan)                     # masked cells drawn in C["masked"] grey (4 Nov 2026)
    im = ax.pcolormesh(S.edges(lon), S.edges(lat), np.ma.masked_invalid(field), transform=PC, cmap=tc, vmin=-vmax, vmax=vmax,
                       shading="flat", rasterized=True, zorder=0)
    jj, ii = np.nonzero(fdr); S.sig_dots(ax, lon[ii], lat[jj], transform=PC)
    return im, len(ii)


H = 230.0; D = 66.0; YM = 150.0
fig = S.figure(S.W2, H)
# a ---------------------------------------------------------------- wind-stress trend
ax = S.axmm(fig, 8, YM, D, D, projection=ccrs.SouthPolarStereo())
imw, nw = trend_map(ax, MW.lon.values, MW.lat.values, 1000 * MW["wind_trend_dec"].values, np.zeros(MW["wind_fdr"].shape, bool), 12, bad="#DCDCDC", cmap="BrBG_r")   # 5 Nov 2026: own colormap (EKE RdBu_r, wind BrBG_r, slope PuOr)   # 4 Nov 2026: FDR dots removed from (a) at the author's request; (b) keeps them
S.letter(fig, 2, H - 5, "a", "Eastward wind-stress trend, 2001–2025")
S.cbar(fig, imw, 8 + 0.1 * D, YM - 10, 0.8 * D, "Wind-stress trend (10$^{-3}$ N m$^{-2}$ per decade)", extend="both", ticks=[-12, -6, 0, 6, 12])
# b ---------------------------------------------------------------- slope trend
ax = S.axmm(fig, 93, YM, D, D, projection=ccrs.SouthPolarStereo())
ims, ns = trend_map(ax, M.lon.values, M.lat.values, MM_["slope_pct"].values, MM_["slope_fdr"].values > 0, 8, cmap="PuOr")
sp = MM_["slope_pct"].values; fd = MM_["slope_fdr"].values > 0
print(f"Fig 4a FDR cells {nw}; Fig 4b FDR cells {ns} (positive {int((sp[fd] > 0).sum())}, negative {int((sp[fd] < 0).sum())})")
S.letter(fig, 86, H - 5, "b", "Isopycnal-slope trend, σ$_0$ = 27.3 kg m$^{-3}$, 2004–2025")
S.cbar(fig, ims, 93 + 0.1 * D, YM - 10, 0.8 * D, "Slope trend (% of mean per decade)", extend="both", ticks=[-8, -4, 0, 4, 8])
fig.legend(handles=S.map_legend_handles(fdr=False, buffer=True), loc="center", bbox_to_anchor=(0.5, 127.5 / H), ncol=4,
           handlelength=1.8, columnspacing=1.4, fontsize=S.FS_ANN)
# c ---------------------------------------------------------------- common-cell series, three strips (4 Nov 2026)
xe = ES.x.values; c4 = ES.year.values >= 2004; xc = xe[c4].mean()
STRIP = 13.0; Y0 = 73.0; GAP = 2.0
for i, (k, col, nm, yl) in enumerate((("wind", CW, "wind stress", (-14, 14)), ("EKE", CE, "EKE", (-16, 20)), ("slope", CS, "isopycnal slope", (-4, 4)))):
    ax = S.axmm(fig, 16, Y0 + (2 - i) * (STRIP + GAP), 128, STRIP)
    v = ES[k].values; mu = np.nanmean(v[c4]); pa = 100 * (v / mu - 1); r = D1.loc[k]
    ax.plot(xe[c4], pa[c4], color=col, lw=S.LW["raw"], alpha=0.25, zorder=2)
    ax.plot(xe[c4], pd.Series(pa).rolling(12, center=True).mean().values[c4], color=col, lw=S.LW["smooth"] + 0.1, zorder=3)
    S.trend_line(ax, 2004, 2026, xc, 0.0, r.pct / 10, r.pct_lo / 10, r.pct_hi / 10, r.p, col, z=5)
    ax.text(2026.6, 0.0, S.um(f"{nm}\n{r.pct:+.1f} % per decade\n({r.pct_lo:+.1f} to {r.pct_hi:+.1f}), {S.fp(r.p)}" + ("" if S.sigp(r.p) else " (n.s.)")),
            color=col, fontsize=S.FS_ANN, ha="left", va="center", linespacing=1.1)
    S.zero(ax, lim=[2003.8, 2026.2]); S.ts_axis(ax, 2003.8, 2026.2, 2031.4); ax.set_ylim(*yl); S.yticks_minus(ax)
    if i < 2: ax.set_xticklabels([])
    if i == 1: ax.set_ylabel("Anomaly (% of 2004–2025 mean)")
S.letter(fig, 2, Y0 + 3 * STRIP + 2 * GAP + 3, "c", "Wind stress, EKE and isopycnal slope over the ACC band + 100 km (shading: 95 % trend interval)")
# d ---------------------------------------------------------------- regression coefficients
COLS = {"KE": C["ke"], "EKE": CE, "slope": CS}
ax = S.axmm(fig, 16, 13, 50, 46)
# detrended coefficients (4 Nov 2026): same 21 years and normalisation as f4s, both series detrended (OLS), 2000 pairs bootstraps
DET = {k: (X8.loc[k, "b_detr"], X8.loc[k, "bd_lo"], X8.loc[k, "bd_hi"]) for k in ("KE", "EKE", "slope")}
for q, k in enumerate(("KE", "EKE", "slope")):
    r = SR.loc[k]; col = COLS[k]; xb = q - 0.13
    ax.bar(xb, r.elasticity, 0.42, color=col if S.sigp(r.p) else "white", edgecolor=col, lw=0.9, zorder=3)
    ax.plot([xb, xb], [r.el_lo, r.el_hi], color="#333333", lw=0.7, zorder=4)
    ax.text(xb, max(r.el_hi, r.elasticity) + 0.05, S.um(f"{r.elasticity:+.2f}"), ha="center", va="bottom", fontsize=S.FS_ANN, zorder=2.5, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    b, lo, hi = DET[k]; xd = q + 0.25
    ax.plot([xd, xd], [lo, hi], color=col, lw=0.8, zorder=4)
    ax.plot(xd, b, "D", color=col, mfc="white", mec=col, mew=0.8, ms=3.2, zorder=5)
    ax.text(xd, hi + 0.05, S.um(f"{b:+.2f}"), ha="center", va="bottom", fontsize=S.FS_ANN, color=col, zorder=2.5, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
print("detrended:", {k: tuple(round(v, 2) for v in DET[k]) for k in DET})
ax.legend(handles=[Patch(facecolor="#BBBBBB", edgecolor="#555555", label="trends retained"),
                   Line2D([], [], color="#555555", marker="D", mfc="white", ms=3.2, lw=0.8, label="both series detrended")],
          loc="upper right", fontsize=S.FS_ANN, handlelength=1.4, bbox_to_anchor=(1.0, 1.04))
for yv, ls, c_ in ((1.0, (0, (3, 1.6)), "#777777"), (0.7, (0, (1, 1.2)), "#9A9A9A")):
    ax.axhline(yv, color=c_, lw=0.6, ls=ls, zorder=0)
    ax.text(2.72, yv, f"{yv:g}", ha="left", va="center", fontsize=S.FS_ANN, color=c_)
S.zero(ax)
ax.set_xticks(range(3)); ax.set_xticklabels(["KE", "EKE", "isopycnal\nslope"]); ax.tick_params(axis="x", length=0); ax.set_xlim(-0.6, 2.6)
ax.set_ylabel("Exponent $b$ of $X \\propto \\tau^{\\,b}$\n(slope of ln $X$ on ln $\\tau$)"); ax.set_ylim(-0.25, 2.0); ax.set_yticks([0, 0.5, 1.0, 1.5]); S.yticks_minus(ax)
S.letter(fig, 2, 64, "d", "Sensitivity to the wind (lag 1 yr), 2004–2025")
# e ---------------------------------------------------------------- lagged correlations, monthly resolution (x12, 9 Nov 2026)
X12 = pd.read_csv(OUT / "x12" / "x12_lagcorr.csv"); X12s = pd.read_csv(OUT / "x12" / "x12_summary.csv").set_index("series")
ax = S.axmm(fig, 101, 13, 65, 46); hl = []
for k, ls, nm in (("KE", "-", "KE"), ("EKE", "-", "EKE"), ("slope", "-", "isopycnal slope")):
    q = X12[X12.series == k].sort_values("lag_months"); col = COLS[k]; x = q.lag_months.values
    if k == "EKE": ax.fill_between(x, -q.env95, q.env95, color="#E6E6E6", lw=0, zorder=0)
    ax.plot(x, q.r, color=col, lw=0.6, zorder=3)
    sig = q.p_pointwise.values < 0.05; rr = np.where(sig, q.r.values, np.nan); ax.plot(x, rr, color=col, lw=1.8, zorder=4, solid_capstyle="butt")
    sm = X12s.loc[k]; ax.plot(sm.lag_max, sm.r_max, "o", color=col, mfc=col if sm.p_fw < 0.05 else "white", mec=col, mew=0.8, ms=3.4, zorder=5)
    hl.append(Line2D([], [], color=col, lw=1.2, label=f"{nm} (family-wise {S.fp(sm.p_fw)})"))
hl.append(Line2D([], [], color="#777777", lw=1.8, label="thick: pointwise p < 0.05")); hl.append(Patch(facecolor="#E6E6E6", label="95 % surrogate envelope"))
S.zero(ax); ax.set_xticks([0, 6, 12, 18, 24, 30, 36]); ax.set_xlabel("Lag (months, wind leading)"); ax.set_xlim(0, 36)
ax.set_ylabel("r with wind stress\n(detrended, 12-month running means)"); ax.set_ylim(-0.6, 1.3); ax.set_yticks([-0.5, 0, 0.5, 1.0]); S.yticks_minus(ax)
ax.legend(handles=hl, loc="upper right", fontsize=S.FS_ANN - 0.5, handlelength=1.8)
S.letter(fig, 86, 64, "e", "Lagged correlation with the wind, 2004–2025")
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("D1:", D1.loc[["wind", "EKE", "slope"], ["pct", "pct_lo", "pct_hi", "p"]].round(3).to_dict("index"))
print("x8:", SR[["elasticity", "el_lo", "el_hi", "p"]].round(3).to_dict("index"))
