"""dk41: EKE averaged along the whole LM Gould line (dk24 section fields, Product 1; ice-covered points masked per month,
line mean over available points weighted by ds), monthly 2001-2025; trend: calibrated Theil-Sen (n_outer 150), % per decade.
Usage: python dk41_line_eke.py <revision dir> <dk24 npz> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, SEC, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
Z = np.load(SEC); u, v, dsx = Z["u"], Z["v"], Z["dsx"]; ice = Z["ice"] > 0.10
ek = (0.5 * ((u - np.nanmean(u, 0)) ** 2 + (v - np.nanmean(v, 0)) ** 2) + 0.5 * Z["vi"]) * 1e4; ek[ice] = np.nan
w = np.where(ice, 0, dsx[None, :]); y = np.nansum(ek * w, 1) / w.sum(1)
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
r = R.trend_bcz(T, y, month=MON, n_outer=150); mu = y.mean()
pd.DataFrame([dict(series="line_mean_EKE", mean=mu, slope=r["slope"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"], p=r["p"], pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu)]).to_csv(OUT / "dk41_line_eke.csv", index=False)
print(round(mu, 1), round(1000 * r["slope"] / mu, 2), round(1000 * r["ci_lo"] / mu, 2), round(1000 * r["ci_hi"] / mu, 2), round(r["p"], 3))
