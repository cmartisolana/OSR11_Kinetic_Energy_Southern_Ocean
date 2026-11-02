"""DK23b: interannual sensitivity d ln X / d ln tau in the Drake regions, as Fig. 4d (s4_elasticity.py, 'detrended'):
annual means, ln X and ln tau each detrended with Theil-Sen, OLS slope with HAC (maxlags 1) t-interval; lags 0 and 1 yr
(response at year t+L on tau at year t). Usage: python dk23b_interannual.py <revision dir> <dk1 dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
REV, DK1, OUT = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV)); import rstats as R
A = pd.read_csv(DK1 / "dk1_series.csv").groupby("year").mean(numeric_only=True); xa = A.index.values + 0.5
def detr(y):
    b, a = R.theil_sen(xa, y); return y - (a + b * xa)
def ols(x, y):
    f = sm.OLS(y, sm.add_constant(x)).fit(cov_type="HAC", cov_kwds={"maxlags": 1}); b, se = f.params[1], f.bse[1]
    tq = stats.t.ppf(0.975, len(y) - 2); return b, b - tq * se, b + tq * se, 2 * stats.t.sf(abs(b / se), len(y) - 2), f.rsquared
rows = []
for reg in ("up", "box", "down"):
    lt = detr(np.log(A[f"wind_{reg}_TAUX"].values))
    for v in ("SPD", "KE", "EKE"):
        for prod in ("allsat", "twosat"):
            ly = detr(np.log(A[f"{prod}_{reg}_{v}"].values))
            for L in (0, 1):
                b, lo, hi, p, r2 = ols(lt[:len(lt) - L], ly[L:])
                rows.append(dict(region=reg, var=v, product=prod, lag_yr=L, sens=b, lo=lo, hi=hi, p=p, r2=r2, n=len(ly) - L))
T = pd.DataFrame(rows); T.to_csv(OUT / "dk23b_interannual.csv", index=False); print(T[T["product"] == "allsat"].round(3).to_string()); print(T[T["product"] == "twosat"].round(3).to_string())
