"""Quick-look figure for E5 (not a manuscript figure). Usage: python e5_quicklook.py <out_dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = Path(sys.argv[1])
q1 = pd.read_csv(out / "e5_q1_longitude_composite.csv"); q1b = pd.read_csv(out / "e5_q1_age_composite.csv")
q2 = pd.read_csv(out / "e5_q2_by_class.csv"); q3 = pd.read_csv(out / "e5_q3_annual.csv", header=[0, 1], index_col=0)
feats = ["KP", "SEIR", "MR", "PAR", "DP"]
fig = plt.figure(figsize=(13, 10))
gs = fig.add_gridspec(3, 5, height_ratios=[1.1, 1, 1])
for k, f in enumerate(feats):
    ax = fig.add_subplot(gs[0, k])
    for per, c, lab in (("early", "0.6", "2001–2011"), ("late", "C3", "2012–2023")):
        d = q1[(q1.feature == f) & (q1.period == per)]
        ax.plot(d.bin + 0.5, d.amp_med * 100, color=c, label=lab)
        if per == "late": ax.fill_between(d.bin + 0.5, d.amp_q25 * 100, d.amp_q75 * 100, color=c, alpha=0.15, lw=0)
    ax.axvspan(-10, 0, color="0.9", zorder=0); ax.set_title(f); ax.set_xlabel("lon east of ridge edge (°)")
    if k == 0: ax.set_ylabel("eddy amplitude, median (cm)"); ax.legend(frameon=False, fontsize=8)
# Q1b age
ax = fig.add_subplot(gs[1, 0:2])
for c, col in zip(("upstream", "ridge", "downstream", "away"), ("C0", "C1", "C3", "0.4")):
    d = q1b[q1b.cls == c]; ax.plot(d.age, d.amp_med * 100, color=col, label=c)
ax.set_xlabel("days since birth"); ax.set_ylabel("amplitude, median (cm)"); ax.set_title("by birth class (survivors only)"); ax.legend(frameon=False, fontsize=8)
# Q2 bars
ax = fig.add_subplot(gs[1, 2:5])
d = q2[(q2.feature == "all") & (q2.pol == "both")].set_index("cls").loc[["upstream", "ridge", "downstream", "away"]]
x = np.arange(4); w = 0.26
ax.bar(x - w, d.amp_max_med * 100, w, label="max amplitude (cm)"); ax.bar(x, d.lifetime_med, w, label="lifetime, median (d)")
ax.bar(x + w, d.dist_med / 10, w, label="distance, median (10 km)")
ax.set_xticks(x, d.index); ax.legend(frameon=False, fontsize=8); ax.set_title("cohort statistics by birth class (n = %s)" % ", ".join(str(int(n)) for n in d.n))
# Q3 series
for k, (v, lab) in enumerate((("lifetime", "mean lifetime (d)"), ("amp_max", "mean max amplitude (m)"), ("dist", "mean distance (km)"))):
    ax = fig.add_subplot(gs[2, k * 5 // 3:(k + 1) * 5 // 3] if k < 2 else gs[2, 3:5])
    for reg, col in (("MEANDER", "C3"), ("AWAY", "0.4")):
        ax.plot(q3.index, q3[(reg, v)], "o-", color=col, ms=3, label=reg)
    ax.set_title(lab); ax.legend(frameon=False, fontsize=8)
fig.suptitle("E5 quick-look: tracked eddies at the five standing meanders (META4.0, 2001–2023)")
fig.tight_layout(); fig.savefig(out / "e5_quicklook.png", dpi=130)
print("saved")
