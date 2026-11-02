"""Scale-band trends at the five standing meanders (dk22), Product 1: (a) eddy part of KE at 100-300 km (eddy dominated),
(b) KE above 1000 km (jets, mean-flow dominated), in the upstream, ridge and downstream sectors and away from the meanders.
% of the band mean per decade, 95 % intervals; filled: p < 0.05; ringed: significant after FDR control (10 %).
Usage: python figS_scale_meanders_v331.py <revision dir> <drake out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, DO, FIGD = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV))
import style_nc as S
from matplotlib.patches import Patch
plt = S.plt; C = S.C
D = pd.read_csv(DO / "dk22" / "dk22_trends.csv")
REG = ["Kerguelen Plateau", "Southeast Indian Ridge", "Macquarie Ridge", "Pacific–Antarctic Ridge", "Drake Passage–Scotia Sea"]
LAB = ["Kerguelen\nPlateau", "Southeast\nIndian Ridge", "Macquarie\nRidge", "Pacific–Antarctic\nRidge", "Drake Passage–\nScotia Sea", "away from\nmeanders"]
SEC = [("upstream", "#BFD7EA"), ("ridge", "#6A9FCB"), ("downstream", "#1F4E79")]
fig, ax = plt.subplots(2, 1, figsize=(S.W2, 0.62 * S.W2), sharex=True, gridspec_kw=dict(hspace=0.22))
for a, (fld, band, title, col_away) in zip(ax, (("EKEm", "100-300", "Eddy part of KE, 100–300 km (eddy dominated)", C["eke"]),
                                               ("KE", ">1000", "KE above 1000 km (jets, mean-flow dominated)", C["ke"]))):
    q = D[(D.field == fld) & (D.band_km == band)]
    for gi, rg in enumerate(REG):
        for si, (sec, c) in enumerate(SEC):
            r = q[(q.region == rg) & (q.sector == sec)].iloc[0]; x = gi + (si - 1) * 0.26; sg = round(r.p, 3) < 0.05
            a.bar(x, r.pct, 0.24, color=c if sg else "white", edgecolor=c, lw=0.9, zorder=3)
            a.vlines(x, r.pct_lo, r.pct_hi, color="#333333", lw=0.5, zorder=4)
            if r.fdr10: a.plot(x, r.pct_hi + 1.2, "v", color="#111111", ms=2.5, zorder=5)
    r = q[q.region.str.startswith("Away")].iloc[0]; sg = round(r.p, 3) < 0.05
    a.bar(5, r.pct, 0.3, color="#9A9A9A" if sg else "white", edgecolor="#9A9A9A", lw=0.9, zorder=3); a.vlines(5, r.pct_lo, r.pct_hi, color="#333333", lw=0.5, zorder=4)
    if r.fdr10: a.plot(5, r.pct_hi + 1.2, "v", color="#111111", ms=2.5, zorder=5)
    a.axvline(4.5, color="#9A9A9A", lw=0.5, ls=(0, (2, 2))); S.zero(a); a.set_ylabel("trend (% per decade)"); a.set_ylim(-20, 32)
    S.plabel(a, "a" if fld == "EKEm" else "b", title)
ax[1].set_xticks(range(6)); ax[1].set_xticklabels(LAB, fontsize=5.5); ax[1].tick_params(axis="x", length=0)
ax[0].legend(handles=[Patch(facecolor=c, edgecolor=c, label=s) for s, c in SEC] + [Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05"),
             plt.Line2D([], [], ls="", marker="v", color="#111111", ms=3, label="significant after FDR")], loc="upper left", ncol=5, fontsize=5)
S.save(fig, FIGD, "FigureS_scale_meanders"); plt.close(fig)
