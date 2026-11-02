"""Figure 4 (v3.32): EKE and its atmospheric forcing, in the forms of earlier studies (Meredith and Hogg 2006: lagged EKE
response to the SAM; Hogg et al. 2015, Zhang et al. 2021: lagged correlation of sector EKE with zonal wind stress) and in
the layout of Fig. 3 (maps left; series, sector lags and the profile along the ACC right). Product 1 unless stated;
wind stress = eastward stress of Product 2; SAM = Marshall index. All correlations: annual means 2001-2025, both series
linearly detrended; positive lag: EKE lags the forcing; p from a t test with effective sample size (dk22, dk23).
a  r(local wind stress, EKE one year later), 3 deg cells (dk22).
b  r(SAM, EKE two years later), 3 deg cells (dk23). Lags in a and b: those of the largest |r| of the ACC-band series in d.
c  Annual anomalies of ACC-band EKE (Products 1 and 4) with the ACC-band wind stress (bars) and the SAM.
d  Lagged correlation of ACC-band and sector EKE with the wind stress over the same band (left) and with the SAM (right);
   grey band: |r| below the p = 0.05 level for n_eff = 24 (lag 1).
e  Correlation along the ACC (6 deg bins, between the mean SAF and PF) of EKE with the local wind stress one year earlier
   and with the SAM two years earlier; grey: |r| below the p = 0.05 level.
Usage: python fig4_v332.py <revision dir> <drake out dir> <dk23 out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, DO, D23, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
CTAU, CSAM = "#555555", "#7B3294"
SEC = [("ACC", "ACC band", "#111111", "-"), ("ACC two-satellite", "ACC band, two-satellite", "#9A9A9A", (0, (3, 1.6))),
       ("Atlantic", "Atlantic", "#E69F00", "-"), ("Indian", "Indian", "#56B4E9", "-"), ("Pacific", "Pacific", "#CC79A7", "-")]
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); an = A.groupby("year").mean(numeric_only=True); yrs = an.index.values + 0.5
SL = pd.read_csv(D23 / "dk23_sector_lagcorr.csv"); PR = pd.read_csv(D23 / "dk23_alongACC.csv")
Zw = np.load(DO / "dk22" / "dk22_corr_3deg.npz"); Zs = np.load(D23 / "dk23_sam_maps_3deg.npz")
LW, LS_ = 1, 2


def sigp(p): return np.isfinite(p) and round(float(p), 3) < 0.05


fig = plt.figure(figsize=(S.W2, 0.98 * S.W2))
gR = fig.add_gridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.24)
g1 = gR[0].subgridspec(2, 1, hspace=0.16)
gC = gR[1].subgridspec(3, 1, height_ratios=[0.9, 1.0, 1.0], hspace=0.48)
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
for k, (Z, rk, pk, lab, ttl) in enumerate(((Zw, f"r_local_lag{LW}", f"p_local_lag{LW}", "a", "EKE vs local wind stress one year earlier"),
                                          (Zs, f"r_sam_lag{LS_}", f"p_sam_lag{LS_}", "b", "EKE vs SAM two years earlier"))):
    ax = fig.add_subplot(g1[k], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
    L2, A2 = np.meshgrid(Z["lon"], Z["lat"])
    im = ax.pcolormesh(L2, A2, np.ma.masked_invalid(Z[rk]), transform=PC, cmap=cm_, vmin=-0.7, vmax=0.7, shading="nearest", rasterized=True, zorder=0)
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="#1A1A1A", linewidths=0.5, transform=PC, zorder=6)
    jj, ii = np.nonzero(Z[pk] < 0.05); ax.scatter(Z["lon"][ii], Z["lat"][jj], s=0.8, c="k", lw=0, transform=PC, zorder=7)
    S.plabel(ax, lab, ttl)
    S.hcbar(fig, im, ax, "correlation (annual means, detrended)", shrink=0.7, extend="both")
# c series
ax = fig.add_subplot(gC[0]); S.zero(ax)
pct = lambda c: 100 * (an[c] / an[c].mean() - 1)
ax.bar(yrs, pct("tau_accmask_ACC"), 0.75, color="#D9D9D9", edgecolor="none", zorder=1, label="wind stress (left)")
ax.plot(yrs, pct("ACC_EKE"), "-o", color=C["eke"], lw=0.8, ms=1.8, zorder=3, label="EKE (left)")
ax.plot(yrs, pct("ACC2s_EKE"), "--", color="#9A9A9A", lw=0.8, zorder=3, label="EKE, two-satellite (left)")
S.ts_style(ax, (2000.6, 2025.6)); ax.set_ylim(-15, 20); ax.set_ylabel("anomaly (%)")
axs = ax.twinx(); axs.spines["right"].set_visible(True); axs.spines["right"].set_linewidth(0.5)
axs.plot(yrs, an.SAM.values, "-", color=CSAM, lw=0.8, label="SAM (right)"); axs.set_ylim(-3.5, 4.0)
axs.set_ylabel("SAM index", color=CSAM, fontsize=5.5); axs.tick_params(axis="y", colors=CSAM, labelsize=5)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axs.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper left", ncol=2, fontsize=4.8)
S.plabel(ax, "c", "ACC-band EKE, wind stress and SAM, annual means")
# d lagged correlation by sector
gd = gC[1].subgridspec(1, 2, wspace=0.08)
for k, fn in enumerate(("wind stress", "SAM")):
    ax = fig.add_subplot(gd[k]); S.zero(ax); ax.axvline(0, color="#BDBDBD", lw=0.5)
    rc = SL[(SL.forcing == fn) & (SL.sector == "ACC") & (SL.lag == 1)].r_crit.iloc[0]
    ax.axhspan(-rc, rc, color="#EFEFEF", lw=0, zorder=0)
    for key, lab, col, ls in SEC:
        q = SL[(SL.forcing == fn) & (SL.sector == key)].sort_values("lag")
        ax.plot(q.lag, q.r, color=col, ls=ls, lw=0.9 if key == "ACC" else 0.7, zorder=3, label=lab)
        sg = q.p < 0.05; ax.plot(q.lag[sg], q.r[sg], "o", color=col, ms=2.4, zorder=4); ax.plot(q.lag[~sg], q.r[~sg], "o", color=col, mfc="white", ms=2.0, mew=0.5, zorder=4)
    ax.set_xlim(-2.3, 5.3); ax.set_xticks(range(-2, 6)); ax.set_ylim(-0.75, 0.8); ax.set_xlabel("lag (years; EKE lags)", labelpad=1)
    ax.set_title(f"with {fn}", fontsize=5.6, pad=2)
    if k == 0:
        ax.set_ylabel("correlation"); S.plabel(ax, "d", "Lagged correlation of EKE, by sector", dy=10)
    else:
        ax.set_yticklabels([])
        ax.legend(handles=[Line2D([], [], color=c_, ls=l_, lw=0.8, label=lb) for _, lb, c_, l_ in SEC], loc="upper left", fontsize=4.4, ncol=1, handlelength=1.6)
# e along the ACC
ax = fig.add_subplot(gC[2]); S.zero(ax)
qt = PR[PR.lag == LW].sort_values("lon"); qs = PR[PR.lag == LS_].sort_values("lon")
rc = float(np.nanmedian(qt.rcrit)); ax.axhspan(-rc, rc, color="#EFEFEF", lw=0, zorder=0)
for q, col, cr, pc, lab in ((qt, CTAU, "r_tau", "p_tau", "local wind stress, EKE 1 yr later"), (qs, CSAM, "r_sam", "p_sam", "SAM, EKE 2 yr later")):
    ax.plot(q.lon, q[cr], color=col, lw=0.8, zorder=3, label=lab)
    sg = q[pc] < 0.05; ax.plot(q.lon[sg], q[cr][sg], "o", color=col, ms=2.2, zorder=4); ax.plot(q.lon[~sg], q[cr][~sg], "o", color=col, mfc="white", ms=1.8, mew=0.45, zorder=4)
TICK = np.arange(-180, 181, 60)
ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels([("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in TICK])
ax.set_ylim(-0.8, 0.9); ax.set_ylabel("correlation"); ax.set_xlabel("Longitude (flow from left to right)")
for x_, s_ in ((-142, "PAR"), (-65, "Drake"), (-40, "Scotia"), (37, "SWIR"), (78, "Kerguelen"), (130, "SEIR"), (160, "MR")):
    ax.text(x_, 0.99, s_, ha="center", va="top", fontsize=4.8, color="#222222", transform=ax.get_xaxis_transform())
ax.legend(loc="lower left", ncol=2, fontsize=4.8)
S.plabel(ax, "e", "EKE correlation along the ACC, between the SAF and PF (6° bins)", dy=4)
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
