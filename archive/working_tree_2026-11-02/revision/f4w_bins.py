"""F4w: cross-bin correlation of the monthly-test sigma0 = 27.3 slope trend (f4m_lon5_trends.csv, 2004-2025) with the eastward
wind-stress trend over the same 5 deg bins of the domain (f4b_bins.csv, 2001-2025); p from all circular shifts (28 Oct 2026).
Usage: python f4w_bins.py <f4m out dir> <f4b out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
FM, FB = (Path(a) for a in sys.argv[1:3])
S = pd.read_csv(FM / "f4m_lon5_trends.csv"); W = pd.read_csv(FB / "f4b_bins.csv")
B = S.merge(W[["lon", "tau_pct", "tau_p"]], on="lon").sort_values("lon"); a, b = B.pct.values, B.tau_pct.values
r = np.corrcoef(a, b)[0, 1]; null = np.array([np.corrcoef(a, np.roll(b, k))[0, 1] for k in range(1, len(b))])
p = (1 + np.sum(np.abs(null) >= abs(r))) / (1 + null.size)
pd.DataFrame([dict(series="cross-bin r(monthly slope trend, wind trend 2001-2025)", r=r, p=p, n=len(B))]).to_csv(FM / "f4w_results.csv", index=False)
print(f"r = {r:+.2f}, p = {p:.3f}, n = {len(B)}")
