"""Figure (draft, v1): EKE and the history of the westerly forcing (dk29).
a  Annual SAM index 1957-2025 (bars) and its exponentially weighted history for memory time scales T = 2 and 15 yr, and
   the two-timescale response of Hogg et al. (2015); grey: altimetry period 2001-2025.
b  ACC-band EKE anomaly 2001-2025 (annual, % of mean) and its least-squares fit to the SAM of the same year (T = 0), to the
   Hogg et al. kernel and to the SAM history with T = 15 yr.
c  Against memory time scale T: correlation of ACC-band EKE with the SAM history (black; grey band: 95th percentile of 2000
   phase-randomised SAM surrogates) and fraction of the observed 2001-2025 EKE trend reproduced by the fit (green).
Usage: python fig_memory_v1.py <revision dir> <dk29 out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, D29, FIGD = (Path(a) for a in sys.argv[1:4])
sys.path.insert(0, str(REV))
import style_nc as S
from matplotlib.lines import Line2D
plt = S.plt; C = S.C
Z = np.load(D29 / "dk29_series.npz", allow_pickle=True); D = pd.read_csv(D29 / "dk29_memory.csv")
CSAM = "#7B3294"; CH = "#D55E00"
fig = plt.figure(figsize=(S.W2, 0.62 * S.W2))
g = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1, 1], wspace=0.28, hspace=0.45)
# a
ax = fig.add_subplot(g[0, 0]); S.zero(ax)
ya = Z["yrs_all"] + 0.5; sa = Z["sam_annual"]
ax.axvspan(2001, 2026, color="#EFEFEF", lw=0, zorder=0); ax.text(2013.5, 2.3, "altimetry record", ha="center", fontsize=5, color="#777777")
ax.bar(ya, sa, 0.8, color=np.where(sa > 0, "#C9A8D6", "#D9D9D9"), lw=0, zorder=1, label="SAM, annual")
ax.plot(ya, Z["sam_full_T2"], color=CSAM, lw=0.9, ls=(0, (3, 1.6)), zorder=3, label="SAM history, T = 2 yr")
ax.plot(ya, Z["sam_full_T15"], color=CSAM, lw=1.3, zorder=3, label="SAM history, T = 15 yr")
ax.plot(ya, Z["sam_full_Hogg"], color=CH, lw=1.0, zorder=3, label="Hogg et al. (2015) response")
ax.set_xlim(1957, 2026); ax.set_ylim(-2.6, 2.6); ax.set_ylabel("SAM index")
ax.legend(loc="lower right", ncol=2, fontsize=4.8)
S.plabel(ax, "a", "Westerly forcing since 1957 (SAM) and its history", dy=4)
# b
ax = fig.add_subplot(g[1, 0]); S.zero(ax)
yy = Z["yrs"] + 0.5
ax.plot(yy, Z["ACC_band|obs"], "-o", color=C["eke"], lw=1.0, ms=2.2, zorder=4, label="ACC-band EKE (observed)")
for key, col, ls, lab in (("0", "#999999", (0, (3, 1.6)), "fit: SAM same year"), ("Hogg", CH, "-", "fit: Hogg et al. (2015) response"),
                          ("15", CSAM, "-", "fit: SAM history, T = 15 yr")):
    ax.plot(yy, Z[f"ACC_band|{key}"], color=col, ls=ls, lw=0.9, zorder=3, label=lab)
S.ts_style(ax, (2000.6, 2025.6)); ax.set_ylabel("EKE anomaly (%)")
ax.legend(loc="upper left", ncol=2, fontsize=4.7)
S.plabel(ax, "b", "ACC-band EKE and fits to the forcing history", dy=4)
# c
ax = fig.add_subplot(g[:, 1])
q = D[(D.region == "ACC band") & (D["T"] != "Hogg")].copy(); q["Tn"] = q["T"].astype(float); q.loc[q.Tn == 0, "Tn"] = 0.25
ax.fill_between(q.Tn, 0, q.r95, color="#EAEAEA", lw=0, zorder=0, label="95 % of SAM surrogates")
ax.plot(q.Tn, q.r, "-o", color="#111111", ms=2.2, lw=0.9, zorder=3, label="correlation with EKE")
ax.plot(q.Tn, q.frac, "-s", color=C["eke"], ms=2.0, lw=0.9, zorder=3, label="fraction of EKE trend reproduced")
hq = D[(D.region == "ACC band") & (D["T"] == "Hogg")].iloc[0]
ax.plot([16], [hq.r], "D", color=CH, ms=3.2, zorder=4); ax.plot([16], [hq.frac], "D", color=CH, mfc="white", ms=3.2, zorder=4)
ax.text(17.5, hq.r, "Hogg et al. (2015)", fontsize=4.6, color=CH, va="center")
ax.set_xscale("log"); ax.set_xticks([0.25, 1, 2, 5, 10, 20, 30]); ax.set_xticklabels(["0", "1", "2", "5", "10", "20", "30"])
ax.set_xlabel("memory time scale T (yr)"); ax.set_ylim(0, 1.05); ax.set_ylabel("correlation  /  fraction")
ax.axhline(1, color="#BDBDBD", lw=0.5, ls=(0, (2, 2)))
ax.legend(loc="lower right", fontsize=4.8)
S.plabel(ax, "c", "How long a memory is needed?", dy=4)
S.save(fig, FIGD, "FigureX_forcing_memory"); plt.close(fig)
print("done")
