"""A2e: detection on annual-mean and non-overlapping 3-year-mean KE (2001-2024
for 3-yr blocks: 8 blocks, as in Chambers 2018; 2025 left out), baseline D1
parameters, matched to climatological envelopes, ACC mask."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fronts as F, a2_run as A, rstats as R

ke, lat, lon, acc, t = A.load()
K = np.asarray(ke, np.float64)
p = A.BASE
d = F.meridional_dists(lat)
kemean = np.nanmean(K, 0)
cc, ct = F.detect_field(kemean, lat, d, p["thr"], p["min_w"], p["frac"], A.MAX_ENV)
rows = []
for label, blocks in [("annual", [list(range(12 * i, 12 * i + 12)) for i in range(25)]),
                      ("3yr", [list(range(36 * i, 36 * i + 36)) for i in range(8)])]:
    nb = len(blocks)
    counts = np.zeros((nb, lon.size), np.int16)
    tab = np.full((nb, lon.size, A.MAX_ENV, F.NF), np.nan, np.float32)
    for b, idx in enumerate(blocks):
        c, tb = F.detect_field(np.nanmean(K[idx], 0), lat, d, p["thr"], p["min_w"], p["frac"], A.MAX_ENV)
        counts[b] = c; tab[b] = tb
    arr, meta = A.matched_tracks(counts, tab, cc, ct, lat, lon, acc, min_cov=0.9)
    keep = meta.keep.values
    x = np.array([2001 + (np.mean(idx) + 0.5) / 12 for idx in blocks])
    secs = A.sector(lon[meta.j.values])
    for reg in ["ACC", "Atlantic", "Indian", "Pacific"]:
        sel = keep & ((secs == reg) | (reg == "ACC"))
        for f in ["lat_hp", "seg_km", "fwhm_km", "w25_km", "peak_val"]:
            a = arr[f][:, sel].astype(float)
            y = np.nanmean(a - np.nanmean(a, 0, keepdims=True), 1) * (1e4 if f == "peak_val" else 1)
            if label == "3yr":
                # OLS with t-test on 8 values (as Chambers 2018, independent blocks)
                from scipy import stats
                r = stats.linregress(x, y)
                tq = stats.t.ppf(0.975, len(x) - 2)
                rows.append(dict(avg=label, region=reg, var=f, ntracks=int(sel.sum()), slope=r.slope,
                                 lo=r.slope - tq * r.stderr, hi=r.slope + tq * r.stderr, p=r.pvalue))
            else:
                T = R.trend(x, y, monthly=False, n_boot=1000)
                rows.append(dict(avg=label, region=reg, var=f, ntracks=int(sel.sum()), slope=T.slope,
                                 lo=T.ci_lo, hi=T.ci_hi, p=T.p_boot))
pd.DataFrame(rows).to_csv(HERE / "outputs" / "a2" / "a2e_averaging.csv", index=False)
print(pd.DataFrame(rows).round(5).to_string())
