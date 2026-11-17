"""Sensitivity of the circumpolar half-power trend to the track-retention threshold (SI Sect. S5, Fig. S5a).

The circumpolar set keeps the ACC tracks matched in at least 90 % of months (900 tracks). This script repeats the
trend for retention thresholds of 75 % (1301 tracks) and 50 % (1695 tracks), and for 90 % as a check, with the same
series definition (mean over tracks of each track's anomaly) and the same calibrated test as 02k. The values were
previously computed by hand (x1_results.md, archived) and are quoted in SI S5.

Purpose : SI S5 retention sensitivity ("−0.8 [−4.1, +2.5] and −0.9 [−3.9, +2.0] km per decade").
Inputs  : outputs/a2/tracks_baseline.npz (02c)
Outputs : outputs/a2/a2_retention_trends.csv (min_cov, n_tracks, slope, lo, hi, p; km per decade, positive northward)
Run     : python scripts/05_stats/05j_retention_trends.py      (no raw data needed; ~3 min)
Seeds   : rstats.trend_bcz(seed=0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd
from osr11 import rstats as R, paths as PA, params as P

T = np.load(PA.OUT / "a2" / "tracks_baseline.npz")
acc, cov, lat_hp = T["acc"].astype(bool), T["coverage"], T["lat_hp"].astype(float)
t = P.Y0 + (np.arange(P.NMONTHS) + 0.5) / 12; month = np.tile(np.arange(1, 13), P.Y1 - P.Y0 + 1)
rows = []
for c in (0.9, 0.75, 0.5):
    m = acc & (cov >= c); a = lat_hp[:, m]; y = np.nanmean(a - np.nanmean(a, 0), 1) * P.KM_PER_DEG
    r = R.trend_bcz(t, y, month=month, monthly=True, n_outer=300, n_bias=100, seed=0)
    rows.append(dict(min_cov=c, n_tracks=int(m.sum()), slope=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"]))
    print(f"coverage >= {c:.2f}: {m.sum():4d} tracks  {10*r['slope']:+.2f} km/decade ({10*r['ci_lo']:+.2f} to {10*r['ci_hi']:+.2f}) p={r['p']:.3f}", flush=True)
pd.DataFrame(rows).to_csv(PA.out("a2") / "a2_retention_trends.csv", index=False)
