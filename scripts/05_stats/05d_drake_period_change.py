"""Late-versus-early EKE between the SAF and PF in the three Drake Passage domains (Sect. 3.1 sentence; written for the
publication repository, the numbers were previously aggregated by hand).

Purpose : percentage by which the 2016-2025 mean EKE exceeds the 2001-2015 mean, between the mean SAF and PF, in
            upstream   86-68 W  (2-deg bin centres -85 ... -69)
            passage    68-58 W  (-67 ... -59)
            downstream 58-40 W  (-57 ... -41)
          (domains as SI Table S3). Statistic = 100 * (sum of the bin means in 2016-2025 / sum in 2001-2015 - 1), i.e. the
          ratio of the regional (equal-weight bin) period means; the mean of the per-bin ratios is given for reference.
          Cross-check column anomaly_difference_pct: from the annual anomalies of dk12 (in % of each bin's 2001-2025 mean), the
          2016-2025 minus 2001-2015 difference averaged over the bins of the domain (SI table drake_domains_period_difference).
Inputs  : outputs/drake/dk8/dk8_along_acc_fronts_periods_allsat.csv (03j), outputs/drake/dk12/dk12_band_annual_allsat.csv (02i)
Outputs : outputs/drake/dk8/dk8_period_change_regions.csv
Run     : python scripts/05_stats/05d_drake_period_change.py          (no raw data needed; < 1 s)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd
from osr11 import paths as PA

REG = {"upstream": (-86.0, -68.0), "passage": (-68.0, -58.0), "downstream": (-58.0, -40.0)}
D = pd.read_csv(PA.dk(8) / "dk8_along_acc_fronts_periods_allsat.csv").pivot(index="lon", columns="period", values="mean")
A = pd.read_csv(PA.dk(12) / "dk12_band_annual_allsat.csv", index_col=0); A.index = A.index.astype(int); alon = A.columns.astype(float)
rows = []
for name, (a, b) in REG.items():
    s = D[(D.index > a) & (D.index < b)]
    rows.append(dict(region=name, lon_west=a, lon_east=b, n_bins=len(s), mean_2001_2015=s["2001-2015"].mean(),
                     mean_2016_2025=s["2016-2025"].mean(), pct_change=100 * (s["2016-2025"].sum() / s["2001-2015"].sum() - 1),
                     pct_change_mean_of_bins=float(np.mean(100 * (s["2016-2025"] / s["2001-2015"] - 1))),
                     anomaly_difference_pct=float(A.loc[2016:2025, (alon > a) & (alon < b)].mean().mean() - A.loc[2001:2015, (alon > a) & (alon < b)].mean().mean())))
T = pd.DataFrame(rows); T.to_csv(PA.dk(8) / "dk8_period_change_regions.csv", index=False); print(T.round(1).to_string(index=False))
