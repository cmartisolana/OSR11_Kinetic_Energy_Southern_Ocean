"""B2 tests from the saved b2 series (b2_pathfollowing.py), with two corrections decided after seeing that the
first pass was contaminated by construction:
  * "downstream" EKE of window c = EKE in the path band of the window centred at c+30 (which does not overlap
    window c), restricted to longitudes (c+15, c+35]; the first pass used the next window (c+10), which
    shares two thirds of its area with c.
  * births are expressed per 1000 km of path, because a longer path has a larger band and collects more
    births for that reason alone.
The "same-window" link (length -> EKE in the band of the same window) is kept for reference only: a long
path samples a different area and an eddy that bends the contour is itself in the band.
Usage: python b2_tests.py <b2 out dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import b1_covar as BC

out = Path(sys.argv[1]); rng = np.random.default_rng(20261014)
S = pd.read_csv(out / "b2_series.csv"); EB = np.load(out / "b2_eke_bins.npy")       # (36, 300, 360)
CENTRES = np.arange(-180, 180, 10); idx = {c: i for i, c in enumerate(CENTRES)}
BINS = np.arange(-180, 180) + 0.5
NT_E = 2023 * 12 + 8 - 2001 * 12; NS = 1000

def wrap(x): return ((x + 180) % 360) - 180

def test(x, y):
    rr = [BC.lagcorr(x, y, L)[0] for L in range(4)]; obs = max(abs(r) for r in rr)
    s = [max(abs(BC.lagcorr(BC.phase_surr(x, rng), y, L)[0]) for L in range(4)) for _ in range(NS)]
    return rr, obs, (1 + np.sum(np.array(s) >= obs)) / (NS + 1)

rows = []
for c in CENTRES:
    L = S[f"len_{c}"].values; cov = float(np.isfinite(L).mean())
    row = dict(centre=c, coverage=cov)
    if cov < 0.8:
        rows.append(row); continue
    cd = wrap(c + 30); rel = wrap(BINS - c)
    sel = (rel > 15) & (rel <= 35)
    down = np.nanmean(EB[idx[cd]][:, sel], axis=1)
    P = S[["year", "month", "x"]].copy(); P["L"] = L; P["E"] = S[f"eke_{c}"].values; P["D"] = down
    xL = BC.prep(P, "L")
    for name in ("E", "D"):
        if np.isfinite(P[name]).mean() < 0.8:
            continue
        rr, obs, p = test(xL, BC.prep(P, name))
        row.update({f"r_{name}_lag{k}": rr[k] for k in range(4)}); row[f"rmax_{name}"] = obs; row[f"p_{name}"] = p
    Pe = P.iloc[:NT_E].copy(); Pe["B"] = 1000 * S[f"births_{c}"].values[:NT_E] / L[:NT_E]
    if np.isfinite(Pe.B).mean() >= 0.8:
        rr, obs, p = test(BC.prep(Pe, "L"), BC.prep(Pe, "B"))
        row.update({f"r_B_lag{k}": rr[k] for k in range(4)}); row["rmax_B"] = obs; row["p_B"] = p
    rows.append(row)
T = pd.DataFrame(rows); T.to_csv(out / "b2_window_tests_v2.csv", index=False)
pd.set_option("display.width", 250)
print(T[["centre", "coverage", "rmax_E", "p_E", "rmax_D", "p_D", "rmax_B", "p_B"]].round(3).to_string())

pool = []
for link in ("E", "D", "B"):
    sub = T[(T.centre % 30 == 10) & T.get(f"p_{link}", pd.Series(np.nan, index=T.index)).notna()]
    zs = [np.sum(np.arctanh(sub[f"r_{link}_lag{k}"].values) * np.sqrt(290 - 3)) / np.sqrt(len(sub)) for k in range(4)]
    best = sub[[f"r_{link}_lag{k}" for k in range(4)]].values
    pool.append(dict(link=link, n_windows=len(sub), **{f"z_lag{k}": zs[k] for k in range(4)},
                     **{f"p_lag{k}": 2 * stats.norm.sf(abs(zs[k])) for k in range(4)},
                     n_sig=int((sub[f"p_{link}"] < 0.05).sum()), frac_positive=float(np.mean(best[np.arange(len(sub)), np.abs(best).argmax(1)] > 0))))
    print(pool[-1])
pd.DataFrame(pool).to_csv(out / "b2_pooled_v2.csv", index=False)
