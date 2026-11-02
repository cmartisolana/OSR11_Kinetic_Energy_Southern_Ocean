"""F4x: cross-bin correlations between the trends of eastward wind stress, EKE and the sigma0 = 27.3 cross-stream slope (28 Oct 2026),
for Figure 4d. Same 72 bins of 5 deg longitude of the Figure 4 domain (ACC band + 100 km, ice removed):
 slope  f4m_lon5_trends.csv (monthly test, 2004-2025), EKE f4d_eke_bins.csv (Product 1, 2004-2025),
 wind   f4b_bins.csv (Product 2, 2001-2025, full wind record). All in % of the bin mean per decade.
Pearson r for each pair; p from all 71 circular shifts of the second profile (keeps the along-path autocorrelation).
Usage: python f4x_corr.py <f4m out dir> <f4d out dir> <f4b out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
FM, FD, FB = (Path(a) for a in sys.argv[1:4])
S = pd.read_csv(FM / "f4m_lon5_trends.csv")[["lon", "pct", "p"]].rename(columns={"pct": "slope", "p": "slope_p"})
E = pd.read_csv(FD / "f4d_eke_bins.csv")[["lon", "eke_pct", "eke_p"]].rename(columns={"eke_pct": "eke"})
W = pd.read_csv(FB / "f4b_bins.csv")[["lon", "tau_pct", "tau_p"]].rename(columns={"tau_pct": "wind", "tau_p": "wind_p"})
B = S.merge(E, on="lon").merge(W, on="lon").sort_values("lon").reset_index(drop=True)
B["sector"] = np.where((B.lon >= -70) & (B.lon < 20), "Atlantic", np.where((B.lon >= 20) & (B.lon < 147), "Indian", "Pacific"))
B.to_csv(FM / "f4x_bins.csv", index=False); rows = []
for a, b in (("wind", "eke"), ("wind", "slope"), ("eke", "slope")):
    x, y = B[a].values, B[b].values; r = np.corrcoef(x, y)[0, 1]
    null = np.array([np.corrcoef(x, np.roll(y, k))[0, 1] for k in range(1, len(y))]); p = (1 + np.sum(np.abs(null) >= abs(r))) / (1 + null.size)
    rows.append(dict(x=a, y=b, r=r, p=p, n=len(B))); print(f"{a} vs {b}: r = {r:+.2f}, p = {p:.3f}")
pd.DataFrame(rows).to_csv(FM / "f4x_results.csv", index=False)
