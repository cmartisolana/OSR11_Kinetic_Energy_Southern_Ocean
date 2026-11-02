"""Figure 4 (v3.33): EKE and the wind, in three forms of earlier studies (dk24), laid out as Fig. 3 (maps left).
a  Trend of the eastward wind stress, 2001-2025 (least squares on annual means, % of the cell mean per decade), with the mean
   SAF and PF (cf. the EKE trend map of Fig. 3a; Martinez-Moreno et al. 2021).
b  Composite (Meredith and Hogg 2006): EKE anomaly one year after the 6 years of strongest minus the 6 years of weakest
   ACC-band wind stress (annual, detrended; % of the cell mean); dots: p < 0.05 (Welch t test).
c  EKE trend against wind-stress trend in twelve 30 deg sectors of the ACC band (Zhang et al. 2021), 95 % intervals.
d  Regional lag maps (Hogg et al. 2022): correlation of the annual EKE of four regions (outlined) with the wind stress one
   year earlier on 3 deg cells (both detrended).
e  ACC-band composite difference (strong minus weak years) of EKE at lags 0-2 yr, for the wind stress and the SAM.
Usage: python fig4_v333.py <revision dir> <drake out dir> <dk24 out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, DO, D24, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
Z = np.load(D24 / "dk24_patterns.npz", allow_pickle=True); ST = pd.read_csv(D24 / "dk24_sector_trends.csv")
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
ice_ds = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc")
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); an = A.groupby("year").mean(numeric_only=True); yrs = an.index.values
CTAU, CSAM = "#555555", "#7B3294"


def fronts(ax):
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="#1A1A1A", linewidths=0.5, transform=PC, zorder=6)


fig = plt.figure(figsize=(S.W2, 0.98 * S.W2))
gR = fig.add_gridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.24)
g1 = gR[0].subgridspec(2, 1, hspace=0.16)
gC = gR[1].subgridspec(3, 1, height_ratios=[1.0, 1.25, 0.75], hspace=0.42)
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
L2, A2 = np.meshgrid(Z["lon"], Z["lat"])
# a wind trend map
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
im = ax.pcolormesh(L2, A2, Z["tau_trend_pct"], transform=PC, cmap=cm_, vmin=-20, vmax=20, shading="nearest", rasterized=True, zorder=0)
fronts(ax); add_ice_mask(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_ds.ice_mask.values > 0.5)
S.plabel(ax, "a", "Eastward wind-stress trend, 2001–2025")
S.hcbar(fig, im, ax, "wind-stress trend (% per decade)", shrink=0.7, extend="both")
# b composite
ax = fig.add_subplot(g1[1], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
F22 = np.load(D24.parent / "dk22" / "dk22_fields.npz"); EKa = F22["EK"].reshape(25, 12, 30, 360).mean(1)
EK3 = EKa.reshape(25, 10, 3, 120, 3).mean(axis=(2, 4)); tt = yrs - yrs.mean()
Ep3 = 100 * (EK3 / EK3.mean(0) - 1); Ep3 = Ep3 - Ep3.mean(0) - tt[:, None, None] * ((tt[:, None, None] * (Ep3 - Ep3.mean(0))).sum(0) / (tt ** 2).sum())
hi_ = np.searchsorted(yrs, Z["years_tau_strong"]); lo_ = np.searchsorted(yrs, Z["years_tau_weak"]); h_, l_ = hi_[hi_ + 1 < 25] + 1, lo_[lo_ + 1 < 25] + 1
d3 = Ep3[h_].mean(0) - Ep3[l_].mean(0); p3 = stats.ttest_ind(Ep3[h_], Ep3[l_], axis=0, equal_var=False).pvalue
la3, lo3 = Z["lat"].reshape(-1, 3).mean(1), Z["lon"].reshape(-1, 3).mean(1); L3, A3 = np.meshgrid(lo3, la3)
im = ax.pcolormesh(L3, A3, np.ma.masked_invalid(d3), transform=PC, cmap=cm_, vmin=-20, vmax=20, shading="nearest", rasterized=True, zorder=0)
jj, ii = np.nonzero(p3 < 0.05); ax.scatter(lo3[ii], la3[jj], s=0.8, c="k", lw=0, transform=PC, zorder=7)
fronts(ax)
S.plabel(ax, "b", "EKE one year after strong minus weak wind years, 3° cells")
S.hcbar(fig, im, ax, "EKE anomaly difference (% of mean)", shrink=0.7, extend="both")
# c sector scatter
ax = fig.add_subplot(gC[0]); S.zero(ax); S.zero(ax, "x")
E_ = ST[ST["var"] == "EKE"].set_index("lon"); T_ = ST[ST["var"] == "tau"].set_index("lon")
ax.errorbar(T_.pct, E_.pct, xerr=[T_.pct - T_.lo, T_.hi - T_.pct], yerr=[E_.pct - E_.lo, E_.hi - E_.pct], fmt="none", ecolor="#CFCFCF", lw=0.5, zorder=1)
sg = (E_.p < 0.05).values
ax.scatter(T_.pct[sg], E_.pct[sg], s=14, color=C["eke"], zorder=3, lw=0); ax.scatter(T_.pct[~sg], E_.pct[~sg], s=14, facecolor="white", edgecolor=C["eke"], lw=0.7, zorder=3)
lab = lambda c: f"{abs(c):.0f}°{'W' if c < 0 else 'E'}"
OFF = {-105: (3, 2), -75: (3, -5), -135: (3, -6), 135: (3, 2), 165: (3, -6), -165: (3, 2), 15: (3, -4), 45: (3, -5), 75: (-14, -6), 105: (-6, 4), -15: (-5, 4), -45: (3, 4)}
for c_ in E_.index: ax.annotate(lab(c_), (T_.pct[c_], E_.pct[c_]), xytext=OFF.get(int(c_), (3, 2)), textcoords="offset points", fontsize=4.2, color="#444444")
ax.text(0.98, 0.03, "labels: sector centre longitude", transform=ax.transAxes, ha="right", fontsize=4.3, color="#777777")
r = np.corrcoef(T_.pct, E_.pct)[0, 1]; sp = stats.spearmanr(T_.pct, E_.pct)
ax.text(0.02, 0.97, f"r = {r:.2f} (12 sectors; Spearman {sp.correlation:.2f}, p = {sp.pvalue:.2f})", transform=ax.transAxes, fontsize=4.9, va="top")
ax.set_xlabel("wind-stress trend (% per decade)", labelpad=1); ax.set_ylabel("EKE trend (% per decade)")
S.plabel(ax, "c", "EKE and wind-stress trends in 30° sectors of the ACC band", dy=4)
# d regional lag maps
gd = gC[1].subgridspec(2, 2, wspace=0.0, hspace=0.14)
REG = {"SE Pacific (130-80W)": (-130, -80), "Kerguelen lee (70-100E)": (70, 100), "SE Indian-Macquarie (110-170E)": (110, 170), "Atlantic (50W-20E)": (-50, 20)}
TT = {"SE Pacific (130-80W)": "SE Pacific", "Kerguelen lee (70-100E)": "Kerguelen lee", "SE Indian-Macquarie (110-170E)": "SE Indian–Macquarie", "Atlantic (50W-20E)": "Atlantic"}
L4, A4 = np.meshgrid(Z["lon4"], Z["lat4"])
for k, (nm, (a_, b_)) in enumerate(REG.items()):
    ax = fig.add_subplot(gd[k // 2, k % 2], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature, lons=False)
    imr = ax.pcolormesh(L4, A4, Z[f"reg_{nm}_lag1"], transform=PC, cmap=cm_, vmin=-0.7, vmax=0.7, shading="nearest", rasterized=True, zorder=0)
    fronts(ax)
    xs = np.linspace(a_, b_, 40)
    for la in (-36.5,): ax.plot(xs, np.full_like(xs, la), color="#111111", lw=2.0, transform=PC, zorder=8, solid_capstyle="butt")
    ax.set_title(TT[nm], fontsize=5.2, pad=1.5)
    if k == 0: S.plabel(ax, "d", "Correlation of regional EKE with wind stress one year earlier", dy=12)
bb = gC[1].get_position(fig); cax = fig.add_axes([bb.x1 + 0.005, bb.y0 + 0.02, 0.008, bb.height - 0.05])
cb = fig.colorbar(imr, cax=cax, orientation="vertical", extend="both")
cb.set_label("correlation (annual, detrended)", fontsize=5.2); cb.ax.tick_params(labelsize=5, width=0.4, length=2); cb.outline.set_linewidth(0.4)
# e composite lags, ACC band
ax = fig.add_subplot(gC[2]); S.zero(ax)
ea = 100 * (an.ACC_EKE / an.ACC_EKE.mean() - 1).values; t = yrs - yrs.mean(); ea = ea - np.polyval(np.polyfit(t, ea, 1), t)
w_ = 0.36
for k, (nm, col, lab_) in enumerate((("tau", CTAU, "wind stress"), ("sam", CSAM, "SAM"))):
    hi = np.searchsorted(yrs, Z[f"years_{nm}_strong"]); lo = np.searchsorted(yrs, Z[f"years_{nm}_weak"])
    for L in (0, 1, 2):
        h, l = hi[hi + L < 25] + L, lo[lo + L < 25] + L
        dd = ea[h].mean() - ea[l].mean(); p = stats.ttest_ind(ea[h], ea[l], equal_var=False).pvalue
        se = np.sqrt(ea[h].var(ddof=1) / len(h) + ea[l].var(ddof=1) / len(l)); x = L + (k - 0.5) * w_
        ax.bar(x, dd, w_ * 0.9, color=col if p < 0.05 else "white", edgecolor=col, lw=0.8, zorder=3)
        ax.vlines(x, dd - 1.96 * se, dd + 1.96 * se, color="#333333", lw=0.5, zorder=4)
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["same year", "1 yr later", "2 yr later"], fontsize=5.3); ax.tick_params(axis="x", length=0)
ax.set_ylabel("EKE difference (%)")
ax.legend(handles=[Patch(facecolor=CTAU, label="strong − weak wind years"), Patch(facecolor=CSAM, label="high − low SAM years"),
                   Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="lower left", ncol=3, fontsize=4.7)
S.plabel(ax, "e", "ACC-band EKE after strong and weak years (6 each)", dy=4)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
