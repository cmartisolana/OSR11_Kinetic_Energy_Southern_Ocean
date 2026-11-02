"""Supplementary figure for A3 (injection-recovery). Static PNG."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
d = pd.read_csv(HERE / "outputs" / "a3" / "a3_realisations.csv")
s = pd.read_csv(HERE / "outputs" / "a2" / "a2_summary.csv")
obs = s[(s.set == "baseline") & (s.region == "ACC")].iloc[0]
r = d[d.real > 0]
C1, C2, C3, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6})
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))

# (a) position: recovered slope vs injected shift
a = ax[0]
inj = {"S0_null": 0, "S1_shift5km": 5, "S1_shift10km": 10, "S1_shift25km": 25, "S1_shift50km": 50}
for col, c, lab, off in [("lat_hp_slope", C1, "matched fixed sample (revised)", -0.6),
                         ("u_median_lat_slope", C2, "all detections, median (v2 metric)", 0.6)]:
    for sc, km in inj.items():
        y = r[r.scenario == sc][col] * 111.2 * 25
        a.errorbar(km + off, y.median(), yerr=[[y.median() - y.quantile(.025)], [y.quantile(.975) - y.median()]],
                   fmt="o", color=c, ms=4, elinewidth=1.2, label=lab if sc == "S0_null" else None)
a.plot([0, 50], [0, 50], color=MUTED, lw=0.8, ls="--")
a.text(42, 46, "1:1", color=MUTED, fontsize=7.5)
a.set_xlabel("Injected northward shift over 25 yr (km)")
a.set_ylabel("Recovered shift over 25 yr (km)")
a.set_title("(a) Frontal-position recovery (median, 95 % of 41 realisations)", loc="left", fontsize=9)
a.legend(frameon=False, fontsize=7.5, loc="upper left")

# (b) width trends produced by amplitude-only growth vs observed
b = ax[1]
amp = {"S0_null": 0.0, "S3_amp0.15": 0.15, "S3_amp0.30": 0.30, "S3_amp0.60": 0.60}
for col, c, lab, off, o in [("seg_km_slope", C2, "above-threshold width (v2 metric)", -0.012, obs.m_seg_km_slope),
                            ("fwhm_km_slope", C1, "FWHM", 0.0, obs.m_fwhm_km_slope),
                            ("w25_km_slope", C3, "width at 25 % of peak", 0.012, obs.m_w25_km_slope)]:
    for sc, a_ in amp.items():
        y = r[r.scenario == sc][col]
        b.errorbar(a_ + off, y.median(), yerr=[[y.median() - y.quantile(.025)], [y.quantile(.975) - y.median()]],
                   fmt="o", color=c, ms=4, elinewidth=1.2, label=lab if sc == "S0_null" else None)
    b.axhline(o, color=c, lw=1.0, ls=":")
b.axhline(0, color=MUTED, lw=0.8)
b.set_xlabel("Imposed KE amplitude growth, no change in shape (% yr⁻¹)")
b.set_ylabel("Width trend (km yr⁻¹)")
b.set_title("(b) Width trends from amplitude growth alone (dotted: observed)", loc="left", fontsize=9)
b.legend(frameon=False, fontsize=7.5, loc="upper left")
fig.tight_layout()
out = HERE / "outputs" / "figures"; out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "figS_a3_detectability.png", dpi=200)
print("saved")
