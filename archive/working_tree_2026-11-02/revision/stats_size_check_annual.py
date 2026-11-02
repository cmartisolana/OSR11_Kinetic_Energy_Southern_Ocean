import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[0])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, numpy as np
sys.path.insert(0, _REV)
import rstats as R
n, phi, reps = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(n * 100 + int(phi * 10)); x = np.arange(n, dtype=float) + 2001
rej = []; cov = []
for k in range(reps):
    y = R._simulate_ar(np.array([phi]), 1.0, n, 1, rng, burn=300)[0]
    o = R.trend_bcz(x, y, monthly=False, n_outer=200, seed=k)
    rej.append(o["p"] < 0.05); cov.append(o["ci_lo"] <= 0 <= o["ci_hi"])
    if (k + 1) % 50 == 0:
        print(n, phi, k + 1, "size", np.mean(rej).round(3), "coverage", np.mean(cov).round(3), flush=True)
