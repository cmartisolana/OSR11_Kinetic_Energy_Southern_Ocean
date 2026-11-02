"""SI Figure S7: energy pathway versus height contour on the 900 retained half-power tracks, Product 1, 2001-2025.
a  circumpolar-mean anomaly of the half-power latitude (purple) and of the latitude of the fixed ADT contour through the mean
   half-power point (grey), km (thin: monthly; thick: 12-month running means; straight lines: Theil-Sen trends, solid p < 0.05,
   dashed otherwise; trends, 95 % intervals and p from the calibrated test, revision/outputs/x1/x1_circumpolar_trends.csv).
b  ADT at the half-power point (blue, same conventions) and the ACC-band mean ADT (black dashed, 12-month running mean), cm anomaly.
Inputs: outputs/x1/x1_tracks.npz and x1_band_adt.npy (02j).
Port of revision/x1_figure.py (FigureX_contour_vs_energy).

Publication repository: ported from figures_final/scripts/figS7_contour_vs_energy.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/x1/x1_tracks.npz, x1_band_adt.npy (02j), x1_circumpolar_trends.csv (02k)
Outputs : rendered/{png,pdf}/FigureS7
Run     : python figures/figS7_contour_vs_energy.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd
from osr11 import style as S, paths as P  # noqa: E402
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else P.FIG; X1 = P.OUT / "x1"   # 02j/02k
plt = S.plt; C = S.C; KM = 111.195
Z = np.load(X1 / "x1_tracks.npz"); TR = pd.read_csv(P.OUT / "x1" / "x1_circumpolar_trends.csv").set_index("series")
band = 100 * np.load(X1 / "x1_band_adt.npy")
t = 2001 + (np.arange(300) + 0.5) / 12; tc = t.mean()
def series(a): return np.nanmean(a - np.nanmean(a, 0), 1)

H = 94.0
fig = S.figure(S.W2, H)
CH, CC, CA = C["pos"], "#7F7F7F", C["vel"]
# a
ax = S.axmm(fig, 20, 51, 112, 36)
for key, arr, col, nm, laby in (("lat_hp", "lat_hp", CH, "half-power latitude", 9), ("lat_contour", "lat_c", CC, "fixed ADT contour through\nthe mean half-power point", -18)):
    y = series(Z[arr]) * KM; r = TR.loc[key]
    ax.plot(t, y, color=col, lw=S.LW["raw"], alpha=0.5); ax.plot(t, S.runfull(y), color=col, lw=S.LW["smooth"] + 0.1)
    S.trend_line(ax, 2001, 2026, tc, 0.0, r.slope / 10, r.lo / 10, r.hi / 10, r.p, col, wedge=False)
    ax.text(2026.8, laby, S.um(f"{nm}\n{r.slope:+.1f} km per decade\n({r.lo:+.1f} to {r.hi:+.1f}), {S.fp(r.p)}"), color=col,
            fontsize=S.FS_ANN, ha="left", va="center", linespacing=1.1, clip_on=False)
ax.set_ylabel("Latitude anomaly (km,\npositive northward)"); ax.set_ylim(-42, 32); S.zero(ax); S.yticks_minus(ax)
S.ts_axis(ax, 2000.8, 2026.2)
S.letter(fig, 2, H - 4, "a", "Position of the energy pathway and of the height contour, 900 retained tracks")
# b
ax = S.axmm(fig, 20, 9, 112, 32)
y = series(Z["adt_hp"]) * 100; r = TR.loc["adt_hp"]
ax.plot(t, y, color=CA, lw=S.LW["raw"], alpha=0.5); ax.plot(t, S.runfull(y), color=CA, lw=S.LW["smooth"] + 0.1)
S.trend_line(ax, 2001, 2026, tc, 0.0, r.slope / 10, r.lo / 10, r.hi / 10, r.p, CA, wedge=False)
ax.text(2026.8, 6.0, S.um(f"ADT at the half-power point\n{r.slope:+.1f} cm per decade\n({r.lo:+.1f} to {r.hi:+.1f}), {S.fp(r.p)}"), color=CA,
        fontsize=S.FS_ANN, ha="left", va="center", linespacing=1.1, clip_on=False)
yb = band - band.mean(); ax.plot(t, S.runfull(yb), color="#222222", lw=0.8, ls=(0, (2.5, 1.2)))
ax.text(2026.8, -1.5, "ACC-band mean ADT\n(12-month running mean)", color="#222222", fontsize=S.FS_ANN, ha="left", va="center", linespacing=1.1, clip_on=False)
ax.set_ylabel("ADT anomaly (cm)"); ax.set_ylim(-9, 9); S.zero(ax); S.yticks_minus(ax)
S.ts_axis(ax, 2000.8, 2026.2)
S.letter(fig, 2, 45, "b", "Height of the energy pathway")
S.save(fig, FIGD, "FigureS7"); plt.close(fig)
from scipy import stats
print("S7 uncalibrated Theil-Sen check (km or cm per decade):", {k: round(10 * stats.theilslopes(series(Z[k]) * (100 if k == "adt_hp" else KM), t)[0], 2) for k in ("lat_hp", "lat_c", "adt_hp")},
      "band ADT", round(1000 * stats.theilslopes(band / 100, t)[0], 2), "cm per decade")
