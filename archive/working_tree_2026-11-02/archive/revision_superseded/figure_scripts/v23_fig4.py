"""v2.3 Figure 4: along-stream structure of the trends on the jet-core path (b3_alongstream.py, b3_trends_fast.py).
x: longitude along the ACC (flow from left to right), 5-degree bins. Grey bands: bins where the mean path crosses
water shallower than 3000 m (topographic obstacles). Hatched: path in the Malvinas Current (not ACC).
(a) mean EKE within 150 km of the path (log) and depth under the mean path;
(b) EKE trend, (c) jet speed trend, (d) path sinuosity trend (all % per decade), (e) path latitude trend (km per decade,
+ north); 95 % CI; filled p < 0.05; ringed: significant after FDR (10 %).
Usage: python v23_fig4.py
"""
import numpy as np, pandas as pd
from figstyle import *  # noqa
FIG = OUT / "figures_v23"; FIG.mkdir(exist_ok=True)
D = OUT / "b3"
C = pd.read_csv(D / "b3_climatology.csv")
T = {v: pd.read_csv(D / f"b3_trends_{v}.csv") for v in ("eke", "speed", "sinu", "lat") if (D / f"b3_trends_{v}.csv").exists()}
SHALLOW = C.lon[C.depth_min < 3000].values
MALV = C.lon[(C.lon > -60) & (C.lon < -50)].values
LAB = [(-135, "Pacific–Antarctic\nRidge"), (-62, "Drake\nPassage"), (40, "SW Indian R.\n/ Crozet"), (78, "Kerguelen\nPlateau"),
       (133, "SE Indian\nRidge"), (163, "Macquarie\nRidge")]
fig, axs = plt.subplots(6, 1, figsize=(11.0, 13.0), sharex=True, gridspec_kw=dict(height_ratios=[0.9, 0.45, 1, 1, 1, 1], hspace=0.14))
def bands(ax):
    for l in SHALLOW: ax.axvspan(l - 2.5, l + 2.5, color="#E3E6EA", lw=0, zorder=0)
    for l in MALV: ax.axvspan(l - 2.5, l + 2.5, facecolor="none", edgecolor="#B0B6BD", hatch="////", lw=0, zorder=0)
ax = axs[0]; bands(ax)
ax.semilogy(C.lon, C.eke_mean.where(~C.lon.between(-60, -50)), "-", color=STYLE["eke"], lw=2.0)
ax.set_ylabel("Mean EKE\n(cm$^2$ s$^{-2}$)", fontsize=9); base_timeseries_style(ax); ax.set_ylim(40, 1500)
for x_, s in LAB: ax.text(x_, 1350, s, ha="center", va="top", fontsize=7.5, color="#4A4A4A")
ax.set_title("a) Mean EKE within 150 km of the jet-core path", loc="left", fontsize=10.5, fontweight="bold")
ax = axs[1]; bands(ax)
ax.fill_between(C.lon, 6000, C.depth_min, step="mid", color="#6B7280", lw=0); ax.set_ylim(6000, 0); ax.set_yticks([0, 3000, 6000])
ax.set_ylabel("Depth (m)", fontsize=9); base_timeseries_style(ax)
ax.set_title("b) Shallowest depth under the mean path", loc="left", fontsize=10.5, fontweight="bold")
spec = [("eke", "c) EKE within 150 km of the path", "% dec$^{-1}$", STYLE["eke"]),
        ("speed", "d) Jet speed on the core streamline", "% dec$^{-1}$", "#1D6FA3"),
        ("sinu", "e) Path sinuosity (meander curvature)", "% dec$^{-1}$", "#6A4C93"),
        ("lat", "f) Path latitude (+ north)", "km dec$^{-1}$", STYLE["lat"])]
for ax, (v, ttl, unit, col) in zip(axs[2:], spec):
    bands(ax); ax.axhline(0, color="#495057", lw=0.8)
    if v in T:
        d = T[v].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]; sig = d.p < 0.05
        ax.vlines(d.lon, d.lo, d.hi, color=col, lw=1.1, alpha=0.8)
        ax.plot(d.lon[sig], d.tr[sig], "o", color=col, ms=4.5, zorder=4)
        ax.plot(d.lon[~sig], d.tr[~sig], "o", mfc="white", color=col, ms=4.5, zorder=4)
        f = d.fdr10.astype(bool)
        ax.plot(d.lon[f], d.tr[f], "o", mfc="none", mec="#111111", ms=9, mew=1.2, zorder=5)
        n = len(d); ax.text(0.995, 0.04, f"p < 0.05: {int(sig.sum())} of {n} bins; FDR: {int(f.sum())}", transform=ax.transAxes,
                            ha="right", va="bottom", fontsize=8, color="#4A4A4A")
        lo, hi = np.nanpercentile(np.r_[d.lo, d.hi], [2, 98]); pad = 0.1 * (hi - lo); ax.set_ylim(lo - pad, hi + pad)
    ax.set_ylabel(unit, fontsize=9); base_timeseries_style(ax)
    ax.set_title(ttl, loc="left", fontsize=10.5, fontweight="bold")
axs[-1].set_xlim(-180, 180); axs[-1].set_xticks(np.arange(-180, 181, 30))
axs[-1].set_xticklabels([f"{abs(v)}°{'W' if v < 0 else 'E' if v > 0 else ''}" if abs(v) != 180 else "180°" for v in np.arange(-180, 181, 30)])
axs[-1].set_xlabel("Longitude (flow from left to right)")
fig.savefig(FIG / "Figure4.png"); fig.savefig(FIG / "Figure4.pdf"); print("saved")
