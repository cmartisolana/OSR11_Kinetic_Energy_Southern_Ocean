"""SI figure: annual EKE between the mean SAF and PF in 2 deg longitude bins, anomaly in % of each bin's 2001-2025 mean
(dk12; moved from the main circumpolar figure), with the topographic obstacles marked.
Usage: python figS_hovmoller_v327.py <revision dir> <dk12 dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, DK, FIGD = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV)); import style_nc as S
plt = S.plt
HV = pd.read_csv(DK / "dk12_band_annual_allsat.csv", index_col=0); lc = HV.columns.astype(float).values; yy = HV.index.values
OBS = [(-142, "PAR"), (-65, "Drake"), (-40, "Scotia"), (37, "SWIR"), (78, "Kerguelen"), (138, "SEIR"), (160, "Macquarie")]
TICK = np.arange(-180, 181, 30)
fig, ax = plt.subplots(figsize=(S.W2, 0.42 * S.W2))
im = ax.pcolormesh(np.append(lc - 1, lc[-1] + 1), np.append(yy, yy[-1] + 1), HV.values, cmap="RdBu_r", vmin=-40, vmax=40, shading="flat", rasterized=True)
for x_, s_ in OBS:
    ax.axvline(x_, color="#555555", lw=0.3, ls=":"); ax.text(x_, 2026.15, s_, ha="center", va="bottom", fontsize=5.0)
ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels([("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in TICK])
ax.set_xlabel("Longitude (flow from left to right)"); ax.set_ylim(2001, 2026)
ax.set_yticks(np.array([2005, 2010, 2015, 2020, 2025]) + 0.5); ax.set_yticklabels(["2005", "2010", "2015", "2020", "2025"])
cb = fig.colorbar(im, ax=ax, extend="both", pad=0.01, aspect=25); cb.set_label("EKE anomaly (% of mean)", fontsize=6); cb.ax.tick_params(labelsize=5.5)
S.save(fig, FIGD, "FigureS_eke_hovmoller"); print("done")
