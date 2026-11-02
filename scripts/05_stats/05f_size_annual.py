"""Purpose : false-positive rate and CI coverage of trend_bcz for short annual AR(1) series (SI S4: up to 9.5 % for 15 values).
Inputs  : none (synthetic)
Outputs : printed (the shipped logs are outputs/stats/size_bc/ann15.log and ann25.log)
Run     : python scripts/05_stats/05f_size_annual.py 15 0.5 200 > outputs/stats/size_bc/ann15.log
          python scripts/05_stats/05f_size_annual.py 25 0.5 200 > outputs/stats/size_bc/ann25.log
Seeds   : default_rng(n*100 + int(phi*10)); per replicate seed = replicate index.
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, numpy as np
from osr11 import rstats as R
n, phi, reps = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(n * 100 + int(phi * 10)); x = np.arange(n, dtype=float) + 2001
rej = []; cov = []
for k in range(reps):
    y = R._simulate_ar(np.array([phi]), 1.0, n, 1, rng, burn=300)[0]
    o = R.trend_bcz(x, y, monthly=False, n_outer=200, seed=k)
    rej.append(o["p"] < 0.05); cov.append(o["ci_lo"] <= 0 <= o["ci_hi"])
    if (k + 1) % 50 == 0:
        print(n, phi, k + 1, "size", np.mean(rej).round(3), "coverage", np.mean(cov).round(3), flush=True)
