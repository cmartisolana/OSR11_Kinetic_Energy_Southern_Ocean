"""Figure (draft, v2): EKE and the history of the westerly forcing, in the layout of Figure 3 (maps left, series right).
a  Correlation of annual EKE (3 deg cells, % of mean, not detrended) with the SAM of the same year, 2001-2025.
b  Same with the SAM history, exponential memory T = 15 yr. Dots: r above the 95th percentile of 2000 phase-randomised SAM
   surrogates passed through the same kernel (dk30). Mean SAF and PF.
c  Annual SAM 1957-2025 and its history for T = 2 and 15 yr and the Hogg et al. (2015) response (dk29).
d  ACC-band EKE anomaly and its fits to the SAM of the same year, the Hogg et al. response and the T = 15 yr history.
e  Correlation of EKE with the SAM history against T: ACC band (black) and sectors (thin), 95 % of surrogates for the
   ACC band (grey), and the fraction of the ACC-band EKE trend reproduced (green); diamonds: Hogg et al. kernel.
Usage: python fig_memory_v2.py <revision dir> <drake out dir> <work out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, DO, WO, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
Z = np.load(WO / "dk29" / "dk29_series.npz", allow_pickle=True); D = pd.read_csv(WO / "dk29" / "dk29_memory.csv")
M = np.load(WO / "dk30" / "dk30_maps.npz"); MS = pd.read_csv(WO / "dk30" / "dk30_summary.csv").set_index("T")
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
ice_ds = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc")
bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation.coarsen(latitude=4, longitude=4, boundary="trim").mean()
CSAM = "#7B3294"; CH = "#D55E00"
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
fig = plt.figure(figsize=(S.W2, 0.98 * S.W2))
gR = fig.add_gridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.24)
g1 = gR[0].subgridspec(2, 1, hspace=0.16)
gC = gR[1].subgridspec(3, 1, height_ratios=[1.0, 1.0, 1.0], hspace=0.48)
L3, A3 = np.meshgrid(M["lon3"], M["lat3"])
for k, (T, lab, ttl) in enumerate((("0", "a", "EKE vs SAM of the same year"), ("15", "b", "EKE vs SAM history (memory T = 15 yr)"))):
    ax = fig.add_subplot(g1[k], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
    im = ax.pcolormesh(L3, A3, M[f"r_{T}"], transform=PC, cmap=cm_, vmin=-0.8, vmax=0.8, shading="nearest", rasterized=True, zorder=0)
    ax.contour(bat.longitude, bat.latitude, -bat, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="#1A1A1A", linewidths=0.6, transform=PC, zorder=6)
    add_ice_mask(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_ds.ice_mask.values > 0.5)
    jj, ii = np.nonzero(M[f"p_{T}"] < 0.05); ax.scatter(M["lon3"][ii], M["lat3"][jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
    q = MS.loc[T]
    ax.text(0.99, -0.03, f"median r = {q.median_r:.2f}\ncells p < 0.05: {100 * q.frac_cells_p05:.0f} %", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=4.8, color="#222222", zorder=12, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    S.plabel(ax, lab, ttl)
    if k == 1: S.hcbar(fig, im, ax, "correlation with EKE, 2001–2025 (dots: p < 0.05 vs SAM surrogates)", shrink=0.8, extend="both")
# c SAM since 1957
ax = fig.add_subplot(gC[0]); S.zero(ax)
ya = Z["yrs_all"] + 0.5; sa = Z["sam_annual"]
ax.axvspan(2001, 2026, color="#EFEFEF", lw=0, zorder=0); ax.text(2013.5, 2.25, "altimetry", ha="center", fontsize=4.8, color="#777777")
ax.bar(ya, sa, 0.8, color=np.where(sa > 0, "#C9A8D6", "#D9D9D9"), lw=0, zorder=1, label="SAM, annual")
ax.plot(ya, Z["sam_full_T2"], color=CSAM, lw=0.8, ls=(0, (3, 1.6)), zorder=3, label="history, T = 2 yr")
ax.plot(ya, Z["sam_full_T15"], color=CSAM, lw=1.3, zorder=3, label="history, T = 15 yr")
ax.plot(ya, Z["sam_full_Hogg"], color=CH, lw=1.0, zorder=3, label="Hogg et al. (2015) response")
ax.set_xlim(1957, 2026); ax.set_ylim(-2.6, 2.6); ax.set_ylabel("SAM index")
ax.legend(loc="lower right", ncol=2, fontsize=4.6)
S.plabel(ax, "c", "Westerly forcing since 1957 (SAM) and its history", dy=4)
# d ACC-band EKE and fits
ax = fig.add_subplot(gC[1]); S.zero(ax)
yy = Z["yrs"] + 0.5
ax.plot(yy, Z["ACC_band|obs"], "-o", color=C["eke"], lw=1.0, ms=2.0, zorder=4, label="ACC-band EKE")
q = D[D.region == "ACC band"].set_index("T")
for key, col, ls, lab in (("0", "#999999", (0, (3, 1.6)), "SAM same year"), ("Hogg", CH, "-", "Hogg et al. response"), ("15", CSAM, "-", "history, T = 15 yr")):
    ax.plot(yy, Z[f"ACC_band|{key}"], color=col, ls=ls, lw=0.9, zorder=3, label=f"fit: {lab} (r = {q.loc[key, 'r']:.2f})")
S.ts_style(ax, (2000.6, 2025.6)); ax.set_ylabel("EKE anomaly (%)")
ax.legend(loc="upper left", ncol=2, fontsize=4.5)
S.plabel(ax, "d", "ACC-band EKE and fits to the forcing history", dy=4)
# e r vs T
ax = fig.add_subplot(gC[2])
def tn(d): d = d[d["T"] != "Hogg"].copy(); d["Tn"] = d["T"].astype(float); d.loc[d.Tn == 0, "Tn"] = 0.25; return d
qa = tn(D[D.region == "ACC band"])
ax.fill_between(qa.Tn, 0, qa.r95, color="#EAEAEA", lw=0, zorder=0)
for reg, ls in (("Pacific", "-"), ("Indian", (0, (3, 1.4))), ("Atlantic", (0, (1, 1)))):
    qs = tn(D[D.region == reg]); ax.plot(qs.Tn, qs.r, color="#777777", lw=0.6, ls=ls, zorder=2)
ax.plot(qa.Tn, qa.r, "-o", color="#111111", ms=2.0, lw=0.9, zorder=3)
ax.plot(qa.Tn, qa.frac, "-s", color=C["eke"], ms=1.8, lw=0.9, zorder=3)
hq = D[(D.region == "ACC band") & (D["T"] == "Hogg")].iloc[0]
ax.plot([45], [hq.r], "D", color=CH, ms=3.0, zorder=4, clip_on=False); ax.plot([45], [hq.frac], "D", color=CH, mfc="white", ms=3.0, zorder=4, clip_on=False)
ax.set_xscale("log"); ax.set_xlim(0.2, 55); ax.set_xticks([0.25, 1, 2, 5, 10, 20, 30, 45]); ax.set_xticklabels(["0", "1", "2", "5", "10", "20", "30", "Hogg"])
ax.axvline(37, color="#BDBDBD", lw=0.5)
ax.set_xlabel("memory time scale T (yr)"); ax.set_ylim(0, 1.32); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_ylabel("correlation / fraction")
ax.axhline(1, color="#BDBDBD", lw=0.5, ls=(0, (2, 2)))
ax.legend(handles=[Line2D([], [], color="#111111", marker="o", ms=2, lw=0.9, label="r, ACC band"),
                   Patch(facecolor="#EAEAEA", label="95 % of SAM surrogates"),
                   Line2D([], [], color="#777777", lw=0.6, label="r, Pacific / Indian / Atlantic (solid / dashed / dotted)"),
                   Line2D([], [], color=C["eke"], marker="s", ms=2, lw=0.9, label="fraction of ACC-band EKE trend reproduced")],
          loc="upper left", fontsize=4.4, ncol=2)
S.plabel(ax, "e", "How long a memory is needed?", dy=4)
S.save(fig, FIGD, "FigureX_forcing_memory"); plt.close(fig)
print("done")
