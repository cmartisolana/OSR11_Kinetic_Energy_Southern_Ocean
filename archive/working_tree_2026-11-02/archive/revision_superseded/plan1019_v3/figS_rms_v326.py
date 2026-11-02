"""SI figure: RMS of the month-to-month velocity anomalies inside the ACC band (B5e).
a  2001-2025 RMS of the north-south anomaly v' (cm/s) on 1 deg cells; b  trend of the monthly RMS v' (OLS on the
deseasonalised monthly series, cm/s per decade; slope only, as Fig. 3a); c  trends of the band RMS v' and RMS u' in 5 deg
bins (% per decade, calibrated test, 95 % intervals; filled p < 0.05, ringed FDR), with the upstream / ridge /
downstream regions of the five meanders. Equal relative trends of u' and v' = eddies stronger without turning north-south.
Usage: python figS_rms_v326.py <revision dir> <b5 out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, B5, FIGD = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
src = open(Path(__file__).parent / "b5_meander_shape.py").read(); exec(src[src.index("MEANDERS = {"):src.index("files = [")])
SHORT = {"Kerguelen Plateau": "KP", "Southeast Indian Ridge": "SEIR", "Macquarie Ridge": "MR", "Pacific–Antarctic Ridge": "PAR", "Drake Passage–Scotia Sea": "DP"}
M = np.load(B5 / "b5e_map.npz"); D = pd.read_csv(B5 / "b5e_trends.csv"); B = D[D.lon.notna()]
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask
fig = plt.figure(figsize=(S.W2, 0.98 * S.W2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.72], hspace=0.22, wspace=0.12)
for k, (fld, cm, vmin, vmax, ttl, lab) in enumerate(((M["mean"], "viridis", 0, 25, "RMS north–south velocity anomaly, 2001–2025", "RMS v′ (cm s$^{-1}$)"),
                                                     (M["slope"], "RdBu_r", -2, 2, "Trend of RMS v′, 2001–2025", "cm s$^{-1}$ per decade"))):
    ax = fig.add_subplot(gs[0, k], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature, extent_lat=-70)
    im = ax.pcolormesh(M["lon"], M["lat"], fld, transform=PC, cmap=cm, vmin=vmin, vmax=vmax, shading="nearest", rasterized=True, zorder=0)
    ax.contour(mean_adt.longitude, mean_adt.latitude, mean_adt.values, levels=[-0.6, 0.2], colors="#555555", linewidths=0.35, transform=PC, zorder=4)
    add_ice_mask(ax, ccrs, ice.longitude.values, ice.latitude.values, ice.values > 0.5)
    for nm, (u_, r_, d_) in MEANDERS.items():
        ax.text(0.5 * (r_[0] + d_[1]), -39, SHORT[nm], transform=PC, fontsize=5.2, fontweight="bold", ha="center", va="center", zorder=9,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
    S.plabel(ax, "ab"[k], ttl); S.hcbar(fig, im, ax, lab, shrink=0.7, extend="max" if k == 0 else "both")
ax = fig.add_subplot(gs[1, :])
for nm, (u_, r_, d_) in MEANDERS.items():
    for (a, b), col in ((u_, "#56B4E9"), (r_, "#555555"), (d_, "#E69F00")): ax.plot([a, b], [1.0, 1.0], color=col, lw=2.2, solid_capstyle="butt", transform=ax.get_xaxis_transform(), clip_on=False)
    ax.text(np.mean(r_), 1.04, SHORT[nm], transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.0, fontweight="bold")
S.zero(ax)
for v, col, mk, dx in (("rms_v", "#CC79A7", "o", -0.8), ("rms_u", "#7F7F7F", "s", 0.8)):
    d = B[(B["var"] == v) & B.p.notna()].sort_values("lon"); x = d.lon.values + dx; sg = (d.p.round(3) < 0.05).values
    ax.vlines(x, d.pct_lo, d.pct_hi, color=col, lw=0.6); ax.plot(x[sg], d.pct.values[sg], mk, color=col, ms=2.5, zorder=4)
    ax.plot(x[~sg], d.pct.values[~sg], mk, mfc="white", color=col, ms=2.5, mew=0.6, zorder=4)
    fd = d.fdr10.astype(bool).values; ax.plot(x[fd], d.pct.values[fd], "o", mfc="none", mec="#111111", ms=5.5, mew=0.7, zorder=5)
ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 30)); ax.set_ylabel("% per decade")
ax.set_xticklabels([("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in np.arange(-180, 181, 30)]); ax.set_xlabel("Longitude (flow from left to right)")
ax.legend(handles=[Line2D([], [], ls="-", marker="o", color="#CC79A7", ms=2.5, lw=0.6, label="RMS v′ (north–south)"),
                   Line2D([], [], ls="-", marker="s", color="#7F7F7F", ms=2.5, lw=0.6, label="RMS u′ (east–west)")], loc="upper left", fontsize=5.0, ncol=2)
S.plabel(ax, "c", "Trends of the RMS velocity anomalies inside the ACC, by longitude", dy=14)
S.save(fig, FIGD, "FigureS_rms_velocity"); print("done")
