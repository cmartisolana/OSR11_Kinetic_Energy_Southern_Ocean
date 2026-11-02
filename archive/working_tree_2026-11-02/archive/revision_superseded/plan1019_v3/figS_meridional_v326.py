"""SI figure: north-south meandering inside the ACC band, trends from monthly anomalies (B5d; supersedes the B5c annual version).
d  trend of the north-south share of the month-to-month eddy energy (B5d).
a  2001-2025 mean meridional velocity v (cm/s, + north) on 1 deg cells inside the ACC band (mean ADT -0.6..+0.2 m, ice-free;
   cells with standing branches |v| >= 2 cm/s); b  trend of the standing north-south flow s = v' sign(v_mean) from monthly
   anomalies (OLS on the deseasonalised monthly series, cm/s per decade; red: standing branches strengthen); c  trend of the band-mean s in
   5 deg bins (cm/s per decade, calibrated test, 95 % intervals) and the change expected if the flow only accelerated along
   a fixed path (mean s x jet-speed trend of b3, which is measured on the jet-core path, not the band); bars: upstream / ridge / downstream regions of the five meanders.
Usage: python figS_meridional_v326.py <revision dir> <b5 out dir> <figure dir>"""
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
M = np.load(B5 / "b5d_map.npz"); DT = pd.read_csv(B5 / "b5d_trends.csv"); B = DT[DT.lon.notna()]
P = B[B["var"] == "S"].sort_values("lon"); Qp = B[B["var"] == "Q"].sort_values("lon")
SP = pd.read_csv(REV / "outputs" / "b3" / "b3_trends_speed.csv")
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask
fig = plt.figure(figsize=(S.W2, 1.18 * S.W2))
gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 0.62, 0.62], hspace=0.3, wspace=0.12)
for k, (fld, vmax, ttl, lab) in enumerate(((M["vmean"], 20, "Mean meridional velocity inside the ACC, 2001–2025", "v (cm s$^{-1}$, + north)"),
                                           (M["slope"], 3, "Trend of the standing north–south flow (monthly), 2001–2025", "cm s$^{-1}$ per decade (+: branches strengthen)"))):
    ax = fig.add_subplot(gs[0, k], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature, extent_lat=-70)
    im = ax.pcolormesh(M["lon"], M["lat"], fld, transform=PC, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="nearest", rasterized=True, zorder=0)
    ax.contour(mean_adt.longitude, mean_adt.latitude, mean_adt.values, levels=[-0.6, 0.2], colors="#555555", linewidths=0.35, transform=PC, zorder=4)
    add_ice_mask(ax, ccrs, ice.longitude.values, ice.latitude.values, ice.values > 0.5)
    for nm, (u_, r_, d_) in MEANDERS.items():
        ax.text(0.5 * (r_[0] + d_[1]), -39, SHORT[nm], transform=PC, fontsize=5.2, fontweight="bold", ha="center", va="center", zorder=9,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
    S.plabel(ax, "ab"[k], ttl); S.hcbar(fig, im, ax, lab, shrink=0.7, extend="both")
ax = fig.add_subplot(gs[1, :]); axc = ax
for nm, (u_, r_, d_) in MEANDERS.items():
    for (a, b), col in ((u_, "#56B4E9"), (r_, "#555555"), (d_, "#E69F00")): ax.plot([a, b], [1.0, 1.0], color=col, lw=2.2, solid_capstyle="butt", transform=ax.get_xaxis_transform(), clip_on=False)
    ax.text(np.mean(r_), 1.04, SHORT[nm], transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.0, fontweight="bold")
S.zero(ax)
def pts(d, y, lo, hi, col, mk, dx):
    ax = plt.gca(); sg = (d.p.round(3) < 0.05).values; x = d.lon.values + dx
    ax.vlines(x, lo, hi, color=col, lw=0.6); ax.plot(x[sg], y[sg], mk, color=col, ms=2.6, zorder=4); ax.plot(x[~sg], y[~sg], mk, mfc="white", color=col, ms=2.6, mew=0.6, zorder=4)
d = P.dropna(subset=["p"]); pts(d, d.slope_dec.values, d.lo.values, d.hi.values, "#CC79A7", "o", 0.0)
sp = SP.set_index("lon").tr.reindex(d.lon.values).values                       # jet-speed trend, % per decade (b3)
mc = pd.read_csv(B5 / "b5c_profile.csv").set_index("lon")["mean"].reindex(d.lon.values).values   # mean standing strength |v_mean| (cm/s)
exp_ = mc * sp / 100; ok = np.isfinite(exp_) & ~((d.lon.values > -60) & (d.lon.values < -50))
ax.step(d.lon.values[ok], exp_[ok], where="mid", color=C["tr"], lw=0.8, zorder=3)
ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 30)); ax.set_ylabel("cm s$^{-1}$ per decade")
lo_, hi_ = np.nanmin(d.lo), np.nanmax(d.hi); ax.set_ylim(lo_ - 0.1 * (hi_ - lo_), hi_ + 0.1 * (hi_ - lo_))
ax.set_xticklabels([("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in np.arange(-180, 181, 30)])
ax.legend(handles=[Line2D([], [], ls="-", marker="o", color="#CC79A7", ms=2.6, lw=0.6, label="standing north–south flow, band mean (5° bins)"),
                   Line2D([], [], ls="-", color=C["tr"], lw=0.8, label="expected if the flow only sped up (mean × jet-speed trend)")],
          loc="upper left", fontsize=5.0, ncol=2)
S.plabel(ax, "c", "Standing north–south flow inside the ACC, trend by longitude", dy=14)
ax = fig.add_subplot(gs[2, :]); S.zero(ax)
d = Qp.dropna(subset=["p"]); pts(d, d.pct.values, d.pct_lo.values, d.pct_hi.values, "#555555", "D", 0.0)
ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 30)); ax.set_ylabel("% per decade")
ax.set_xticklabels([("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in np.arange(-180, 181, 30)]); ax.set_xlabel("Longitude (flow from left to right)")
S.plabel(ax, "d", "North–south share of the month-to-month eddy energy, v′²/(u′²+v′²), trend by longitude", dy=4)
S.save(fig, FIGD, "FigureS_meridional_flow"); print("done")
