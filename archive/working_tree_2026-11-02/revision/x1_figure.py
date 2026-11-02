"""x1 figure (exploratory, 29 Oct 2026): contour-vs-energy test on the 900 retained tracks.
(a) circumpolar-mean anomaly of the half-power latitude and of the latitude of the fixed ADT contour through the mean
    half-power point (km, 12-month running means; straight lines: 2001-2025 Theil-Sen trends, calibrated test).
(b) ADT at the half-power point and ACC-band mean ADT (cm anomaly): the energy peak drifts up the height field at the
    rate of regional sea-level rise, i.e. it stays in place while the height contours move south underneath it.
Usage: python x1_figure.py <revision dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, FIGD = (Path(a) for a in sys.argv[1:3]); sys.path.insert(0, str(REV))
import style_nc as S
plt = S.plt; C = S.C; KM = 111.195
Z = np.load(REV / "outputs" / "x1" / "x1_tracks.npz"); TR = pd.read_csv(REV / "outputs" / "x1" / "x1_circumpolar_trends.csv").set_index("series")
band = 100 * np.load(REV / "outputs" / "x1" / "x1_band_adt.npy")
t = 2001 + (np.arange(300) + 0.5) / 12
def series(a): return np.nanmean(a - np.nanmean(a, 0), 1)
def fp(p): return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"
def rm(y): return pd.Series(y).rolling(12, center=True).mean().values
fig, axs = plt.subplots(2, 1, figsize=(S.W1 * 1.35, 0.95 * S.W1), sharex=True, gridspec_kw=dict(hspace=0.32))
CH, CC, CA = C["pos"], "#7F7F7F", C["tr"]
# a
ax = axs[0]
for key, col, nm, laby in (("lat_hp", CH, "energy pathway\n(half-power latitude)", 7), ("lat_contour", CC, "fixed ADT contour\nthrough the mean\nhalf-power point", -17)):
    y = series(Z["lat_c" if key == "lat_contour" else key]) * KM; r = TR.loc[key]
    ax.plot(t, y, color=col, lw=0.3, alpha=0.4); ax.plot(t, rm(y), color=col, lw=1.1)
    xx = np.array([2001, 2026]); ax.plot(xx, r.slope / 10 * (xx - t.mean()), color=col, lw=0.9, ls="-" if r.p < 0.05 else (0, (3, 1.6)))
    ax.text(2026.4, laby, f"{nm}\n{r.slope:+.1f} km per decade\n({r.lo:+.1f} to {r.hi:+.1f}), {fp(r.p)}".replace("-", "−"), color=col, fontsize=4.9, ha="left", va="center", linespacing=1.0)
ax.set_ylabel("Latitude anomaly (km)\nnorth ↑"); ax.set_ylim(-36, 28); S.zero(ax)
S.plabel(ax, "a", "Position of the energy pathway and of the height contour, 900 retained longitudes", dy=4)
# b
ax = axs[1]
y = series(Z["adt_hp"]) * 100; r = TR.loc["adt_hp"]
ax.plot(t, y, color=CA, lw=0.3, alpha=0.4); ax.plot(t, rm(y), color=CA, lw=1.1)
xx = np.array([2001, 2026]); ax.plot(xx, r.slope / 10 * (xx - t.mean()), color=CA, lw=0.9)
ax.text(2026.4, 5.5, f"ADT at the half-power point\n{r.slope:+.1f} cm per decade\n({r.lo:+.1f} to {r.hi:+.1f}), {fp(r.p)}".replace("-", "−"), color=CA, fontsize=4.9, ha="left", va="center", linespacing=1.0)
yb = band - band.mean(); ax.plot(t, rm(yb), color="#222222", lw=0.8, ls=(0, (2.5, 1.2)))
ax.text(2026.4, -5.0, "ACC-band mean ADT\n(12-month running mean)", color="#222222", fontsize=4.9, ha="left", va="center", linespacing=1.0)
ax.set_ylabel("ADT anomaly (cm)"); ax.set_ylim(-9, 9); S.zero(ax)
S.plabel(ax, "b", "Height of the energy pathway", dy=4)
for ax in axs:
    ax.set_xlim(2000.8, 2032.5); ax.set_xticks([2005, 2010, 2015, 2020, 2025]); ax.spines["bottom"].set_bounds(2000.8, 2026.2)
S.save(fig, FIGD, "FigureX_contour_vs_energy"); print("done")
