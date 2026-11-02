"""Purpose : false-positive rate of the sieve / calibrated tests on trend-free monthly AR(1) series with a seasonal cycle (SI S4).
Inputs  : none (synthetic)
Outputs : outputs/stats/size_phi<phi>.npy; progress printed (redirect to a log)
Run     : python scripts/05_stats/05e_size_monthly.py 0.8 200 0 > outputs/stats/size_phi0.8.log
Seeds   : default_rng(<seed argument>); per replicate seed = replicate index.
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, numpy as np
from osr11 import rstats as R, paths as PA
phi, reps, seed = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(seed); n = 300
x = 2001 + (np.arange(n) + .5) / 12; mo = np.tile(np.arange(1, 13), 25)
out = []
for r in range(reps):
    e = rng.normal(size=n + 200); y = np.zeros(n + 200)
    for i in range(1, n + 200): y[i] = phi * y[i - 1] + e[i]
    y = y[200:] + 2 * np.sin(2 * np.pi * (mo - 1) / 12)     # trend-free, with a seasonal cycle
    o = R.trend_calibrated(x, y, month=mo, n_outer=40, n_boot=150, seed=r)
    out.append((o["p_boot"], o["p_cal"], o["p_hamed_rao"]))
    a = np.array(out)
    print(r + 1, "sieve", np.mean(a[:, 0] < .05).round(3), "cal", np.mean(a[:, 1] < .05).round(3),
          "hamed-rao", np.mean(a[:, 2] < .05).round(3), flush=True)
np.save(PA.out("stats") / f"size_phi{phi}.npy", np.array(out))
