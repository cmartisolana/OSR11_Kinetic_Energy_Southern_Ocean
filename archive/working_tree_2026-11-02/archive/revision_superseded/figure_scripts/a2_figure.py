"""Supplementary figure for A2 (half-power validation). Static PNG for review."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R

D = HERE / "outputs" / "a2"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"     # validated categorical slots 1-3
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "lines.linewidth": 2})

s = pd.read_csv(D / "a2_summary.csv")
acc = s[(s.region == "ACC") & (s.ntracks >= 100)].reset_index(drop=True)
lon = pd.read_csv(D / "a2f_lon_trends.csv")
fr = pd.read_csv(D / "detection_frequency_lon.csv")

fig, ax = plt.subplots(2, 2, figsize=(10, 8))

# (a) detection frequency
a = ax[0, 0]
a.plot(fr.lon, 100 * fr.freq, color=C1, lw=1.2)
a.set_xlim(-180, 180); a.set_ylim(0, 105)
a.set_xlabel("Longitude (°E)"); a.set_ylabel("Months with ≥1 ACC envelope (%)")
a.set_title("(a) Detection frequency, baseline 200 cm² s⁻² / 100 km / 25 %", loc="left", fontsize=9)

# (b) half-power latitude trend vs parameter set: matched vs original-style
b = ax[0, 1]
y = np.arange(len(acc))
for off, pre, col, lab in [(-0.15, "m_lat_hp", C1, "matched fixed sample"), (0.15, "u_median_lat", C2, "original style (all detections)")]:
    sl = acc[f"{pre}_slope"] * 111.2 * 10
    lo = acc[f"{pre}_lo"] * 111.2 * 10
    hi = acc[f"{pre}_hi"] * 111.2 * 10
    b.errorbar(sl, y + off, xerr=[sl - lo, hi - sl], fmt="o", ms=4, color=col, elinewidth=1.2, capsize=0, label=lab)
b.axvline(0, color=MUTED, lw=0.8)
b.set_yticks(y); b.set_yticklabels(acc.set); b.invert_yaxis()
b.set_xlabel("ACC half-power latitude trend (km per decade, + = north)")
b.set_title("(b) Position trend across parameter sets (95 % CI)", loc="left", fontsize=9)
b.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=7.5)

# (c) width trend by definition
c = ax[1, 0]
for off, pre, col, lab in [(-0.2, "m_seg_km", C2, "above-threshold width (used in v2)"),
                           (0.0, "m_fwhm_km", C1, "FWHM (amplitude-invariant)"),
                           (0.2, "m_w25_km", C3, "width at 25 % of peak (amplitude-invariant)")]:
    sl = acc[f"{pre}_slope"]; lo = acc[f"{pre}_lo"]; hi = acc[f"{pre}_hi"]
    c.errorbar(sl, y + off, xerr=[sl - lo, hi - sl], fmt="o", ms=4, color=col, elinewidth=1.2, label=lab)
c.axvline(0, color=MUTED, lw=0.8)
c.set_yticks(y); c.set_yticklabels(acc.set); c.invert_yaxis()
c.set_xlabel("ACC envelope-width trend (km yr⁻¹), matched fixed sample")
c.set_title("(c) Width trend depends on the width definition", loc="left", fontsize=9)
c.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=7.5)

# (d) longitude-resolved trends
d = ax[1, 1]
lon["fdr"] = R.bh_fdr(lon.p.values, 0.10)
x = lon.lon0 + 2.5
sl = lon.slope_deg_yr * 111.2 * 10; lo = lon.lo * 111.2 * 10; hi = lon.hi * 111.2 * 10
d.errorbar(x, sl, yerr=[sl - lo, hi - sl], fmt="o", ms=4, color=C1, elinewidth=1.0)
d.axhline(0, color=MUTED, lw=0.8)
d.set_xlim(-180, 180)
d.set_xlabel("Longitude (°E), 5° bins with ≥5 matched tracks")
d.set_ylabel("Half-power latitude trend (km per decade)")
d.set_title(f"(d) Longitude-resolved trends: {int((lon.p<0.05).sum())}/{len(lon)} bins p<0.05, {int(lon.fdr.sum())} after FDR", loc="left", fontsize=9)

fig.tight_layout()
out = HERE / "outputs" / "figures"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "figS_a2_halfpower_validation.png", dpi=200)
print("saved")
