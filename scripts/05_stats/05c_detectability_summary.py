"""Summary of the half-power injection-recovery experiments (written for the publication repository; the working copy held
a3_summary.csv without a script). Reproduces the shipped a3_summary.csv exactly (checked 2 Nov 2026).

Purpose : per scenario and metric: noise-free slope (real 0, '_nf'), median slope over the noisy realisations ('_med') and power
          = fraction of noisy realisations outside the central 95 % (2.5-97.5 percentiles) of the S0 null ('_pow') (Fig. S5c,d).
Inputs  : outputs/a3/a3_realisations.csv (02d)
Outputs : outputs/a3/a3_summary.csv
Run     : python scripts/05_stats/05c_detectability_summary.py      (no raw data needed; < 1 s)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd
from osr11 import paths as PA

FIELDS = ["lat_hp_slope", "fwhm_km_slope", "w25_km_slope", "seg_km_slope", "u_median_lat_slope", "u_seg_km_slope", "u_n_slope"]
d = pd.read_csv(PA.OUT / "a3" / "a3_realisations.csv")
n0 = d[(d.scenario == "S0_null") & (d.real > 0)]
rows = []
for sc, g in d.groupby("scenario"):
    r = {"scenario": sc}
    for f in FIELDS:
        lo, hi = np.percentile(n0[f], [2.5, 97.5]); gg = g[g.real > 0][f]
        r[f + "_nf"] = g[g.real == 0][f].iloc[0]; r[f + "_med"] = gg.median(); r[f + "_pow"] = ((gg < lo) | (gg > hi)).mean()
    rows.append(r)
pd.DataFrame(rows).to_csv(PA.out("a3") / "a3_summary.csv", index=False); print("saved", PA.OUT / "a3" / "a3_summary.csv")
