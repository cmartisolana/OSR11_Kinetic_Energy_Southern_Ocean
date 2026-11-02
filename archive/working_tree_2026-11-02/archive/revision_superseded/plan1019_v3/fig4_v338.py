"""Figure 4 (v3.38): what can the wind, the SAM and ENSO account for in the 2001-2025 EKE increase? (dk24 wind trend, dk28)
a  Trend of the eastward wind stress (% per decade), mean SAF and PF.
b  ACC band and sectors: observed EKE trend and the trends congruent with the local wind stress, the SAM and Nino 3.4
   (interannual sensitivity, EKE lagging by 1 yr, times the driver trend), with 95 % intervals; diamonds: upper bound for the
   wind if EKE is proportional to the wind stress (equilibrium eddy-saturated models).
c  Wind-congruent EKE trend per 3 deg cell (local wind, lag 1 yr).
d  Residual: observed EKE trend minus the wind-, SAM- and ENSO-congruent trends (3 deg cells).
Usage: python fig4_v338.py <revision dir> <drake out dir> <work out dir> <figure dir>"""
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
Z24 = np.load(WO / "dk24" / "dk24_patterns.npz", allow_pickle=True); M = np.load(WO / "dk28" / "dk28_maps.npz"); D = pd.read_csv(WO / "dk28" / "dk28_attribution.csv")
clim = xr.open_dataset(REV / "outputs" / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
ice_ds = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc")
cm_ = plt.get_cmap("RdBu_r").copy(); cm_.set_bad("#F2F2F2")
CW, CSAM, CONI = "#555555", "#7B3294", "#2C7FB8"


def fronts(ax):
    ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="#1A1A1A", linewidths=0.5, transform=PC, zorder=6)


def polar_map(ax, lon, lat, f, vmin, vmax, lab, ttl, cbl, ice=False):
    S.polar(ax, ccrs, cfeature); L, A_ = np.meshgrid(lon, lat)
    im = ax.pcolormesh(L, A_, np.ma.masked_invalid(f), transform=PC, cmap=cm_, vmin=vmin, vmax=vmax, shading="nearest", rasterized=True, zorder=0)
    fronts(ax)
    if ice: add_ice_mask(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_ds.ice_mask.values > 0.5)
    S.plabel(ax, lab, ttl); S.hcbar(fig, im, ax, cbl, shrink=0.7, extend="both")


fig = plt.figure(figsize=(S.W2, 1.0 * S.W2))
g = fig.add_gridspec(2, 2, width_ratios=[0.85, 1.15], height_ratios=[1, 1], wspace=0.22, hspace=0.18)
polar_map(fig.add_subplot(g[0, 0], projection=ccrs.SouthPolarStereo()), Z24["lon"], Z24["lat"], Z24["tau_trend_pct"], -20, 20, "a",
          "Eastward wind-stress trend, 2001–2025", "wind-stress trend (% per decade)", ice=True)
# b bars
ax = fig.add_subplot(g[0, 1]); S.zero(ax)
REGS = ["ACC band", "Pacific", "Indian", "Atlantic"]; w_ = 0.19
SER = [("observed", None, C["eke"]), ("wind stress", 1, CW), ("SAM", 1, CSAM), ("Nino 3.4", 1, CONI)]
for gi, rg in enumerate(REGS):
    q = D[D.region == rg]
    for k, (nm, L, col) in enumerate(SER):
        x = gi + (k - 1.5) * w_
        if nm == "observed":
            r = q[(q.driver == "wind stress") & (q.lag == 1)].iloc[0]; v, lo, hi, filled = r.obs_trend, r.obs_lo, r.obs_hi, True
        else:
            r = q[(q.driver == nm) & (q.lag == L)].iloc[0]; v, lo, hi = r.congruent, r.cong_lo, r.cong_hi; filled = r.p < 0.05
        ax.bar(x, v, w_ * 0.9, color=col if filled else "white", edgecolor=col, lw=0.8, zorder=3)
        ax.vlines(x, lo, hi, color="#333333", lw=0.5, zorder=4)
    eq = q[q.driver.str.startswith("wind stress, equilibrium")].iloc[0]
    ax.plot(gi - 0.5 * w_, eq.congruent, "D", ms=3.2, color=CW, mfc="white", mew=0.8, zorder=5)
ax.set_xticks(range(len(REGS))); ax.set_xticklabels(REGS, fontsize=5.6); ax.tick_params(axis="x", length=0)
ax.set_ylabel("EKE trend (% per decade)"); ax.set_ylim(-4, 21)
ax.legend(handles=[Patch(facecolor=C["eke"], label="observed EKE trend"), Patch(facecolor=CW, label="wind-congruent"),
                   Patch(facecolor=CSAM, label="SAM-congruent"), Patch(facecolor=CONI, label="Niño 3.4-congruent"),
                   Line2D([], [], ls="", marker="D", color=CW, mfc="white", ms=3, label="wind, if EKE ~ τ (equilibrium)"),
                   Patch(facecolor="white", edgecolor="#555555", label="open: sensitivity p ≥ 0.05")], loc="upper right", ncol=2, fontsize=4.6)
S.plabel(ax, "b", "Observed EKE trend and the part congruent with each driver", dy=4)
polar_map(fig.add_subplot(g[1, 0], projection=ccrs.SouthPolarStereo()), M["lon3"], M["lat3"], np.where(M["accmask3"], M["cong_wind"], np.nan), -20, 20, "c",
          "Wind-congruent EKE trend", "EKE trend (% per decade)")
polar_map(fig.add_subplot(g[1, 1], projection=ccrs.SouthPolarStereo()), M["lon3"], M["lat3"], np.where(M["accmask3"], M["resid"], np.nan), -20, 20, "d",
          "EKE trend not accounted for by wind, SAM and Niño 3.4", "EKE trend (% per decade)")
S.save(fig, FIGD, "Figure4"); plt.close(fig)
print("done")
