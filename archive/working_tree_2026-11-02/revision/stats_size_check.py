import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[0])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, numpy as np
sys.path.insert(0, _REV)
import rstats as R
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
np.save((_REV_OUT + f"/stats/size_phi{phi}.npy"), np.array(out))
