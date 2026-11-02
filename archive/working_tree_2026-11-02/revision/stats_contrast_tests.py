import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[0])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, numpy as np, pandas as pd
sys.path.insert(0, _REV)
import rstats as R
O = _REV_OUT
rows = []
def test(lab, df, a, b):
    ya, yb = df[a].values, df[b].values
    d = ya / ya.mean() - yb / yb.mean()
    x = (df.year + (df.month - .5) / 12).values
    o = R.trend_calibrated(x, d, month=df.month.values, n_outer=100, n_boot=300)
    rows.append(dict(test=lab, diff_pct_dec=1000 * o["slope"], lo=1000 * o["ci_cal_lo"], hi=1000 * o["ci_cal_hi"], p_cal=o["p_cal"]))
    print(rows[-1], flush=True)
for p in ("allsat", "twosat"):
    S = pd.read_csv(f"{O}/a5/a5_series_{p}.csv")
    test(f"ACC EKE-KE {p}", S, "ACC_EKE", "ACC_KE")
BX = pd.read_csv(f"{O}/a6/a6b_drake_box_series.csv")
test("Drake box EKE-KE ice mask", BX, "ice_any_EKE", "ice_any_KE")
test("Drake box EKE-KE no mask", BX, "nomask_EKE", "nomask_KE")
pd.DataFrame(rows).to_csv((_REV_OUT + "/stats/contrast_tests.csv"), index=False)
