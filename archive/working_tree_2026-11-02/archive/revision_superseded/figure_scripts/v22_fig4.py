"""v2.2 Figure 4: path length -> EKE by 30-degree window (b2_tests.py outputs), one panel.
Black: +0.18 m level (SAF side), same window; grey: -0.20 m level (PF side), same window; red: +0.18 m level, EKE 15-35 deg
downstream. Best of lags 0-3 (sign kept); filled: p < 0.05 (phase surrogates, max over lags).
Usage: python v22_fig4.py
"""
import numpy as np, pandas as pd
from figstyle import *  # noqa
FIG22 = OUT / "figures_v22"; FIG22.mkdir(exist_ok=True)
T = pd.read_csv(OUT / "b2" / "b2_window_tests_v2.csv"); P = pd.read_csv(OUT / "b2_pf" / "b2_window_tests_v2.csv")
def best(T, link):
    rr = T[[f"r_{link}_lag{k}" for k in range(4)]].values
    return rr[np.arange(len(T)), np.nan_to_num(np.abs(rr)).argmax(1)], T[f"p_{link}"].values
fig, b = plt.subplots(figsize=(11.0, 4.2))
for Tx, link, lab, col, z in ((T, "E", "EKE within 150 km of the path, same window (+0.18 m level)", "#1F2933", 3),
                              (P, "E", "same, −0.20 m level (Polar Front side)", "#9AA3AB", 2),
                              (T, "D", "EKE within 150 km of the path, 15–35° downstream (+0.18 m)", "#C0392B", 3)):
    r, p = best(Tx, link); ok = p < 0.05
    b.plot(Tx.centre, r, "-", color=col, lw=1.5, label=lab, zorder=z)
    b.plot(Tx.centre[ok], r[ok], "o", color=col, ms=5.5, zorder=z + 1); b.plot(Tx.centre[~ok], r[~ok], "o", mfc="white", color=col, ms=5.5, zorder=z + 1)
b.axhline(0, color="#495057", lw=0.8); b.set_xlim(-180, 180); b.set_ylim(-0.35, 0.7)
b.set_xticks(np.arange(-180, 181, 30)); b.set_xticklabels([f"{abs(v)}°{'W' if v < 0 else 'E' if v > 0 else ''}" if abs(v) != 180 else "180°" for v in np.arange(-180, 181, 30)])
for name, x in (("Drake", -62), ("Kerguelen", 80), ("SEIR", 150), ("PAR", -120)):
    b.text(x, 0.66, name, ha="center", va="top", fontsize=8, color="#4A4A4A")
b.set_xlabel("Window centre"); b.set_ylabel("r (path length → EKE, best of lags 0–3)")
b.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.93)); base_timeseries_style(b)
b.text(0.99, 0.04, "filled: p < 0.05", transform=b.transAxes, ha="right", va="bottom", fontsize=8, color="#4A4A4A")
fig.savefig(FIG22 / "Figure4.png"); fig.savefig(FIG22 / "Figure4.pdf"); print("saved")
