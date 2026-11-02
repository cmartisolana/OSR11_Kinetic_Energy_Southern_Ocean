"""Supplementary figure for A4/A5: ACC-band energy terms in the two altimetric products.

Annual means of the deseasonalised monthly series (a5_series_*.csv), as anomalies
from the 2001-2025 mean, for KE, EKE, its two parts (LF = 0.5|U_m - ubar|^2,
SUB = 0.5 var_intra) and the along-stream speed SPD. Calibrated Theil-Sen trends
(a5_trends_*.csv) are given in each panel.
-> outputs/figures/figS_a5_products.png
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
O = HERE / "outputs"
COL = {"allsat": "#1f5fa8", "twosat": "#d0731c"}
LAB = {"allsat": "all-sat (008_047, 1/8°)", "twosat": "two-sat (008_057, 1/4°)"}
VARS = [("KE", "KE (cm² s⁻²)"), ("EKE", "EKE (cm² s⁻²)"), ("LF", "low-frequency part of EKE (cm² s⁻²)"),
        ("SUB", "sub-monthly part of EKE (cm² s⁻²)"), ("SPD", "along-stream speed (cm s⁻¹)")]
fig, axs = plt.subplots(len(VARS), 1, figsize=(7, 11), sharex=True)
for ax, (v, lab) in zip(axs, VARS):
    txt = []
    for p in ("allsat", "twosat"):
        S = pd.read_csv(O / "a5" / f"a5_series_{p}.csv"); T = pd.read_csv(O / "a5" / f"a5_trends_{p}.csv")
        y = S[f"ACC_{v}"] - S.groupby("month")[f"ACC_{v}"].transform("mean")
        a = y.groupby(S.year).mean()
        ax.plot(a.index + 0.5, a.values, "-o", ms=3, lw=1.5, color=COL[p], label=LAB[p])
        r = T[(T.series == f"ACC_{v}") & (T.period == "2001-2025")].iloc[0]
        txt.append(f"{p}: {r.slope:+.3f} yr⁻¹ [{r.ci_cal_lo:+.3f}, {r.ci_cal_hi:+.3f}], {r.pct_per_decade:+.1f} %/dec")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_ylabel(lab, fontsize=8)
    ax.text(0.01, 0.97, "\n".join(txt), transform=ax.transAxes, va="top", fontsize=7)
    ax.grid(alpha=0.25, lw=0.5)
axs[0].legend(loc="lower right", fontsize=7, frameon=False)
axs[-1].set_xlabel("year")
fig.suptitle("ACC band (mean ADT −0.6 to +0.2 m): annual anomalies, both products", fontsize=9)
fig.tight_layout()
(O / "figures").mkdir(exist_ok=True)
fig.savefig(O / "figures" / "figS_a5_products.png", dpi=200)
print("saved")
