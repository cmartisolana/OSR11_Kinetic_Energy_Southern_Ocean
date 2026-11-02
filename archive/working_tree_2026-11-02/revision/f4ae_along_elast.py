"""F4ae: elasticity to the wind along the ACC (28 Oct 2026), per 5 deg bin, same definition as the circumpolar f4s_sens.py:
annual anomalies in % of the bin's 2004-2025 mean (trends kept), EKE and slope of year y against the bin's wind stress of year y-1
(2005-2025, 21 points); elasticity = OLS slope (% per %), 95 % from 2000 pairs bootstraps, p from 1000 phase-randomised surrogates
of the bin's annual wind (2004-2025), lagged alike; BH-FDR 10 % across bins. Input: f4v_bin_series.npz.
Usage: python f4ae_along_elast.py <revision dir> <f4v out dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

REV, FV, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R  # noqa: E402
Z = np.load(FV / "f4v_bin_series.npz"); YY = Z["year"]; yrs = np.arange(2004, 2026); lon = Z["lon"]
ann = lambda s: np.array([s[YY == y].mean() for y in yrs])
rng = np.random.default_rng(20261028)
def surr(y):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size)
rows = []
for b, c in enumerate(lon):
    w = ann(Z["wind"][b]); wp = 100 * (w / w.mean() - 1); wx = wp[:-1]; row = dict(lon=c)
    SW = [surr(wp)[:-1] for _ in range(1000)]
    for k in ("eke", "slope"):
        v = ann(Z[k][b]); yp = 100 * (v / v.mean() - 1); yy = yp[1:]
        el = np.polyfit(wx, yy, 1)[0]; r = np.corrcoef(wx, yy)[0, 1]
        bb = [np.polyfit(wx[i], yy[i], 1)[0] for i in (rng.integers(0, 21, 21) for _ in range(2000))]
        null = np.array([np.corrcoef(s_, yy)[0, 1] for s_ in SW]); p = (1 + np.sum(np.abs(null) >= abs(r))) / 1001
        row.update({f"{k}_el": el, f"{k}_lo": np.percentile(bb, 2.5), f"{k}_hi": np.percentile(bb, 97.5), f"{k}_r": r, f"{k}_p": p})
    rows.append(row)
B = pd.DataFrame(rows)
for k in ("eke", "slope"): B[f"{k}_fdr"] = R.bh_fdr(B[f"{k}_p"].values, 0.10)
B.to_csv(OUT / "f4ae_along.csv", index=False)
for k in ("eke", "slope"):
    print(k, "median elasticity", round(B[f"{k}_el"].median(), 2), "positive", int((B[f"{k}_el"] > 0).sum()), "/", len(B),
          "p<0.05", int((B[f"{k}_p"] < 0.05).sum()), "FDR", int(B[f"{k}_fdr"].sum()), flush=True)
