"""F4a: eddy-saturation diagnostics along the ACC (28 Oct 2026), for Figure 4d, e. Inputs: f4v_bin_series.npz (72 bins of 5 deg
longitude of the common cells: ACC band + 100 km, ice removed, slope and EKE defined in all months; monthly 2004-2025 bin means of
eastward wind stress, EKE and the sigma0 = 27.3 slope).
 A1 relative trends per bin, 2004-2025: rstats.trend_bcz on the deseasonalised monthly series (n_outer 150), % of the bin mean per
    decade, for wind and EKE (slope as f4m, recomputed here on the same bins for consistency).
 A2 interannual response per bin: detrended annual means; r(X_y, wind_{y-L}) for L = 0 and 1 yr (lag 1 = the lag of the circumpolar
    maximum, Fig. 4e; stated as such), X = EKE and slope; p from 1000 phase-randomised surrogates of the bin's detrended annual wind;
    BH-FDR 10 % across bins for each X and lag. Elasticity (regression of X'/mean(X) on wind'/mean(wind), lag 1) is also saved.
Usage: python f4a_along.py <revision dir> <f4v out dir> <out dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

REV, FV, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402

Z = np.load(FV / "f4v_bin_series.npz"); lon, x, YY, MM = Z["lon"], Z["x"], Z["year"], Z["month"]
V = {k: Z[k] for k in ("wind", "eke", "slope")}; yrs = np.arange(2004, 2026)
rows = []; rng = np.random.default_rng(20261028)
def detr(y):
    t = np.arange(y.size); b = np.polyfit(t, y, 1); return y - np.polyval(b, t)
def surr(y):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size)
for b, c in enumerate(lon):
    row = dict(lon=c)
    for k, v in V.items():
        s = v[b]; r = R.trend_bcz(x, s, month=MM, n_outer=150); mu = s.mean()
        row.update({f"{k}_pct": 1000 * r["slope"] / mu, f"{k}_lo": 1000 * r["ci_lo"] / mu, f"{k}_hi": 1000 * r["ci_hi"] / mu, f"{k}_p": r["p"]})
    A = {k: np.array([v[b][YY == y].mean() for y in yrs]) for k, v in V.items()}
    wd = detr(A["wind"]); SW = np.array([surr(wd) for _ in range(1000)])
    for k in ("eke", "slope"):
        xd = detr(A[k])
        for L in (0, 1):
            xa, wa = (xd, wd) if L == 0 else (xd[1:], wd[:-1])
            r_ = np.corrcoef(xa, wa)[0, 1]
            null = np.array([np.corrcoef(xa, s_ if L == 0 else s_[:-1])[0, 1] for s_ in SW])
            row[f"r_{k}_lag{L}"] = r_; row[f"p_{k}_lag{L}"] = (1 + np.sum(np.abs(null) >= abs(r_))) / 1001
        row[f"elast_{k}_lag1"] = np.polyfit(wd[:-1] / A["wind"].mean(), xd[1:] / A[k].mean(), 1)[0]
    rows.append(row); print(int(c), {k_: round(v_, 2) for k_, v_ in row.items() if k_.startswith("r_")}, flush=True)
B = pd.DataFrame(rows)
for k in ("eke", "slope"):
    for L in (0, 1): B[f"fdr_{k}_lag{L}"] = R.bh_fdr(B[f"p_{k}_lag{L}"].values, 0.10)
B["sector"] = np.where((B.lon >= -70) & (B.lon < 20), "Atlantic", np.where((B.lon >= 20) & (B.lon < 147), "Indian", "Pacific"))
B.to_csv(OUT / "f4a_along.csv", index=False)
S = []
for k in ("eke", "slope"):
    for L in (0, 1):
        r = B[f"r_{k}_lag{L}"]; S.append(dict(var=k, lag=L, median_r=r.median(), n_pos=int((r > 0).sum()), n=len(B),
                                              n_p05=int((B[f"p_{k}_lag{L}"] < 0.05).sum()), n_fdr=int(B[f"fdr_{k}_lag{L}"].sum())))
        print(S[-1], flush=True)
pd.DataFrame(S).to_csv(OUT / "f4a_summary.csv", index=False)
print("done", flush=True)
