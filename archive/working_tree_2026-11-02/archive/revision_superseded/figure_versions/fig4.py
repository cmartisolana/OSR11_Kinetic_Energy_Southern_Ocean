"""Figure 4 (v2.5): along-stream structure of the trends on the jet-core path (b3_alongstream.py, b3_trends_fast.py).
x: longitude along the ACC (flow from left to right), 5-degree bins. Grey bands: bins where the mean path crosses water
shallower than 3000 m (topographic obstacles). Hatched: path in the Malvinas Current (excluded).
(a) mean EKE within 150 km of the path (log); (b) shallowest depth under the mean path;
(c) EKE trend, (d) jet speed trend, (e) path sinuosity trend (% per decade), (f) path latitude trend (km per decade, + north);
95 % CI; filled p < 0.05; ringed: significant after FDR (10 %).
Usage: python fig4.py
"""
import numpy as np, pandas as pd
from style_nc import *  # noqa

FIGD = OUT / "figures_final"
D = OUT / "b3"
CL = pd.read_csv(D / "b3_climatology.csv")
T = {v: pd.read_csv(D / f"b3_trends_{v}.csv") for v in ("eke", "speed", "sinu", "lat")}
SHALLOW = CL.lon[CL.depth_min < 3000].values
MALV = CL.lon[(CL.lon > -60) & (CL.lon < -50)].values
LAB = [(-135, "Pacific–Antarctic\nRidge"), (-65, "Drake\nPassage"), (40, "SW Indian R.\n/ Crozet"), (78, "Kerguelen\nPlateau"),
       (133, "SE Indian\nRidge"), (163, "Macquarie\nRidge")]
fig, axs = plt.subplots(6, 1, figsize=(W2, 6.0), sharex=True,
                        gridspec_kw=dict(height_ratios=[0.95, 0.42, 1, 1, 1, 1], hspace=0.5))


def bands(ax):
    for l in SHALLOW: ax.axvspan(l - 2.5, l + 2.5, color="#EDEDED", lw=0, zorder=0)
    for l in MALV: ax.axvspan(l - 2.5, l + 2.5, facecolor="none", edgecolor="#BDBDBD", hatch="////", lw=0, zorder=0)


ax = axs[0]; bands(ax)
ax.semilogy(CL.lon, CL.eke_mean.where(~CL.lon.between(-60, -50)), "-", color=C["eke"], lw=1.0)
ax.set_ylabel("cm$^2$ s$^{-2}$"); ax.set_ylim(40, 2500)
for x_, s in LAB: ax.text(x_, 2300, s, ha="center", va="top", fontsize=5.3, color="#555555", linespacing=0.95)
plabel(ax, "a", "Mean EKE within 150 km of the jet-core path", dx=-28)
ax = axs[1]; bands(ax)
ax.fill_between(CL.lon, 6000, CL.depth_min, step="mid", color="#8A8A8A", lw=0); ax.set_ylim(6000, 0); ax.set_yticks([0, 3000, 6000])
ax.set_ylabel("Depth (m)")
plabel(ax, "b", "Shallowest depth under the mean path", dx=-28)
spec = [("eke", "c", "EKE within 150 km of the path", "% per decade", C["eke"]),
        ("speed", "d", "Jet speed along the path", "% per decade", C["tr"]),
        ("sinu", "e", "Path sinuosity (meander curvature)", "% per decade", C["sin"]),
        ("lat", "f", "Path latitude (+ north)", "km per decade", C["pos"])]
for ax, (v, l, ttl, unit, col) in zip(axs[2:], spec):
    bands(ax); zero(ax)
    d = T[v].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]; sig = d.p < 0.05
    ax.vlines(d.lon, d.lo, d.hi, color=col, lw=0.6)
    ax.plot(d.lon[sig], d.tr[sig], "o", color=col, ms=2.6, zorder=4)
    ax.plot(d.lon[~sig], d.tr[~sig], "o", mfc="white", color=col, ms=2.6, mew=0.6, zorder=4)
    f = d.fdr10.astype(bool) if "fdr10" in d else np.zeros(len(d), bool)
    ax.plot(d.lon[f], d.tr[f], "o", mfc="none", mec="#111111", ms=5.5, mew=0.7, zorder=5)
    lo, hi = np.nanmin(d.lo), np.nanmax(d.hi); pad = 0.06 * (hi - lo); ax.set_ylim(lo - pad, hi + pad)
    ax.set_ylabel(unit)
    plabel(ax, l, ttl, dx=-28)
    print(v, "pos", int((d.tr > 0).sum()), "of", len(d), "sig", int(sig.sum()), "fdr", int(f.sum()), "median", round(float(d.tr.median()), 2))
axs[-1].set_xlim(-180, 180); axs[-1].set_xticks(np.arange(-180, 181, 30))
axs[-1].set_xticklabels([f"{abs(v)}°{'W' if v < 0 else 'E' if v > 0 else ''}" if abs(v) != 180 else "180°" for v in np.arange(-180, 181, 30)])
axs[-1].set_xlabel("Longitude (flow from left to right)")
save(fig, FIGD, "Figure4")
