"""Scale-aware KE trends (dk20): ACC-band KE of the monthly-mean surface geostrophic velocity by length-scale band
(coarse-graining), Products 1 and 4, 2001-2025. a absolute trend per band (cm2 s-2 per decade), b relative trend
(% of the band mean per decade), c mean KE per band. 95 % intervals; open: p >= 0.05.
Usage: python figS_scale_v331.py <revision dir> <drake out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, DO, FIGD = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV))
import style_nc as S
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
plt = S.plt; C = S.C
D = pd.concat([pd.read_csv(DO / "dk20" / f"dk20_trends_{p}.csv") for p in ("allsat", "twosat")])
D = D[(D.field == "KE") & (D.band_km != "all")]
B = ["0-100", "100-200", "200-300", "300-500", "500-1000", ">1000"]
LAB = ["<100", "100–200", "200–300", "300–500", "500–1000", ">1000"]
fig, ax = plt.subplots(1, 3, figsize=(S.W2, 0.36 * S.W2), gridspec_kw=dict(width_ratios=[1, 1, 0.75], wspace=0.38))
for a, col, ylab, lims in ((ax[0], "abs_dec", "KE trend (cm$^2$ s$^{-2}$ per decade)", (-0.5, 8.5)), (ax[1], "pct", "KE trend (% per decade)", (0, 11))):
    for k, (p, c, dx) in enumerate((("allsat", C["ke"], -0.18), ("twosat", C["two"], 0.18))):
        q = D[D["product"] == p].set_index("band_km").loc[B]; x = np.arange(len(B)) + dx
        lo, hi = (q.abs_lo, q.abs_hi) if col == "abs_dec" else (q.pct_lo, q.pct_hi)
        sg = q.p.round(3) < 0.05
        a.bar(x[sg], q[col][sg], 0.34, color=c, zorder=3); a.bar(x[~sg], q[col][~sg], 0.34, color="white", edgecolor=c, lw=0.8, zorder=3)
        a.vlines(x, lo, hi, color="#333333", lw=0.5, zorder=4)
    a.axvspan(-0.5, 0.5, color="#F1F1F1", zorder=0)
    a.set_xticks(range(len(B))); a.set_xticklabels(LAB, fontsize=5.2, rotation=30); a.set_xlabel("length scale (km)")
    a.set_ylabel(ylab); a.set_ylim(*lims); S.zero(a)
ax[0].text(0, 8.3, "below effective\nresolution", ha="center", va="top", fontsize=4.6, color="#777777")
ax[0].legend(handles=[Patch(facecolor=C["ke"], label="all-satellite (Product 1)"), Patch(facecolor=C["two"], label="two-satellite (Product 4)"),
                      Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper right", fontsize=4.9)
S.plabel(ax[0], "a", "Trend by scale band"); S.plabel(ax[1], "b", "Relative trend")
q = D[D["product"] == "allsat"].set_index("band_km").loc[B]
ax[2].bar(range(len(B)), q["mean"], 0.6, color="#BDBDBD"); ax[2].set_xticks(range(len(B))); ax[2].set_xticklabels(LAB, fontsize=5.2, rotation=30)
ax[2].set_ylabel("mean KE (cm$^2$ s$^{-2}$)"); ax[2].set_xlabel("length scale (km)"); S.plabel(ax[2], "c", "Mean, 2001–2025")
S.save(fig, FIGD, "FigureS_scale_KE"); plt.close(fig)
