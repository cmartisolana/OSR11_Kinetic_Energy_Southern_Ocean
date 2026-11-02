"""Figure for the B1 covariability test (supplementary).

Top: lagged correlation between prepared contour length (150 km smoothing) and downstream
EKE, lags -6..+6 months (positive = length leads); grey band = pre-specified window 0-3 months.
Lags outside 0-3 are shown for context only. Bottom: prepared series (standardised) for each meander.
-> outputs/figures/figS_b1_covariability.png
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b1_covar as C

O = HERE / "outputs"
T = pd.read_csv(O / "b1" / "b1_covar_KERG_PAR_SEIR150.csv")
names = {"PAR": "Pacific–Antarctic Ridge", "SEIR150": "SE Indian Ridge (150°E)", "KERG": "Kerguelen (Fawn Trough)"}
fig = plt.figure(figsize=(10, 8))
for j, n in enumerate(names):
    P = pd.read_csv(O / "b1" / f"b1_series_{n}.csv")
    L = C.prep(P, "core_s150_length_km"); D = C.prep(P, "down_eke")
    lags = np.arange(-6, 7)
    r = [C.lagcorr(L, D, k)[0] if k >= 0 else C.lagcorr(D, L, -k)[0] for k in lags]
    ax = fig.add_axes([0.07 + j * 0.32, 0.6, 0.27, 0.33])
    ax.axvspan(-0.3, 3.3, color="0.9")
    ax.axhline(0, color="0.5", lw=0.6)
    ax.plot(lags, r, "-o", ms=3, color="#1f5fa8")
    sub = T[(T.meander == n) & (T.link == "L2_length->down_eke")]
    ax.set_title(f"{names[n]}\nmax-lag surrogate p = {sub.p_surr_maxlag.iloc[0]:.3f}", fontsize=8)
    ax.set_ylim(-0.35, 0.4); ax.set_xlabel("lag (months; + = length leads)", fontsize=8)
    if j == 0:
        ax.set_ylabel("r (contour length, downstream EKE)", fontsize=8)
    ax.tick_params(labelsize=7)
    a2 = fig.add_axes([0.07, 0.38 - j * 0.14, 0.88, 0.11])
    x = P.x.values[1:]
    z = lambda s: (s - np.nanmean(s)) / np.nanstd(s)
    a2.plot(x, z(L), lw=0.9, color="#1f5fa8", label="contour length")
    a2.plot(x, z(D), lw=0.9, color="#d0731c", label="downstream EKE")
    a2.text(0.005, 0.85, names[n], transform=a2.transAxes, fontsize=7)
    a2.set_ylim(-4, 5); a2.tick_params(labelsize=7)
    if j == 0:
        a2.legend(loc="upper right", fontsize=7, frameon=False, ncol=2)
    if j < 2:
        a2.set_xticklabels([])
fig.text(0.5, 0.005, "prepared series: deseasonalised, detrended, 2-month running mean, standardised", ha="center", fontsize=8)
fig.savefig(O / "figures" / "figS_b1_covariability.png", dpi=200)
print("saved")
