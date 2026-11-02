"""DK23: Drake Passage eddy-saturation diagnostics for Fig. 2e-f, from dk1_series.csv (monthly, 2001-2025).
Trend ratio d ln X / d ln tau with X in {along-stream speed SPD, KE, EKE} (Product 1) and tau the eastward wind stress
averaged over the same longitudes (63-54S), per region (upstream 100-80W, passage box 80-52W, downstream 55-35W).
Method as s4_elasticity.py (Fig. 4d): annual means, trend of ln X and ln tau with rstats.trend_bcz (annual, AR(1)),
interval by Monte Carlo from the two (normal) trend uncertainties, 2.5/97.5 percentiles; seed 3.
Usage: python dk23_drake_sensitivity.py <revision dir> <dk1 dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, DK1, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
D = pd.read_csv(DK1 / "dk1_series.csv")
A = D.groupby("year").mean(numeric_only=True); xa = A.index.values + 0.5
rng = np.random.default_rng(3)
def ln_trend(col):
    o = R.trend_bcz(xa, np.log(A[col].values), monthly=False, n_outer=300)
    return o["slope"], (o["ci_hi"] - o["ci_lo"]) / (2 * 1.96), o["p"]
rows = []
for reg in ("up", "box", "down"):
    tc = f"wind_{reg}_TAUX"; assert (A[tc] > 0).all()
    bt, st, pt = ln_trend(tc)
    for v in ("SPD", "KE", "EKE"):
        for prod in ("allsat", "twosat"):
            bx, sx, px = ln_trend(f"{prod}_{reg}_{v}")
            dt = rng.normal(bt, st, 20000); dx = rng.normal(bx, sx, 20000); r = dx / dt
            rows.append(dict(region=reg, var=v, product=prod, lnX_pct_dec=1000 * bx, lnX_p=px, lntau_pct_dec=1000 * bt,
                             lntau_p=pt, ratio=bx / bt, lo=np.percentile(r, 2.5), hi=np.percentile(r, 97.5),
                             P_tau_trend_le0=float(np.mean(dt <= 0))))
T = pd.DataFrame(rows); T.to_csv(OUT / "dk23_trend_ratios.csv", index=False); print(T.round(3).to_string())
