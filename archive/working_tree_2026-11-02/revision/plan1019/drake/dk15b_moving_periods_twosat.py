"""dk15: moving-period trends of the ACC-band KE, EKE and along-stream speed, after the sensitivity test of
Barcelo-Llull et al. (2025, their Fig. 3): trends over start-2025 for start years 2001-2016 (periods of 25 to 10 years),
Product 4 (two-satellite), monthly series of a5_energy.py, trend_bcz (300 null series), % of the period mean per decade and absolute.
Usage: python dk15_moving_periods.py <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from multiprocessing import Pool
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True); sys.path.insert(0, str(REV)); import rstats as R
A = pd.read_csv(REV / "outputs" / "a5" / "a5_series_twosat.csv")
A["year"] = 2001 + A.i // 12; A["month"] = A.i % 12 + 1; A["x"] = A.year + (A.month - 0.5) / 12


def one(args):
    var, y0 = args; d = A[A.year >= y0]
    r = R.trend_bcz(d.x.values, d[var].values.astype(float), month=d.month.values, n_outer=300); mu = float(d[var].mean())
    return dict(var=var, start=y0, end=2025, n_years=2026 - y0, mean=mu, slope_yr=r["slope"], lo_yr=r["ci_lo"], hi_yr=r["ci_hi"], p=r["p"],
                pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu)


if __name__ == "__main__":
    jobs = [(v, y0) for v in ("ACC_EKE", "ACC_KE") for y0 in range(2001, 2017)]
    with Pool(2) as P: rows = P.map(one, jobs)
    D = pd.DataFrame(rows); D.to_csv(OUT / "dk15_moving_periods_twosat.csv", index=False)
    print(D[["var", "start", "pct", "pct_lo", "pct_hi", "p"]].round(3).to_string(index=False)); print("done")
