"""False-positive rate of trend_bc (4th arg: single | cal | z = trend_bcz) (bias-corrected sieve) on trend-free series. Usage:
python stats_size_check_bc.py <process> <reps> <out_dir>
process: ar1_0.5 | ar1_0.8 | ar1_0.95 | ar2 | annual_ar1_0.5 | real_ACC_EKE | real_Drake_EKE
'real_*': the AR model (bias-corrected) fitted to the detrended deseasonalised observed series.

Purpose : false-positive rate of the bias-corrected tests (mode single | cal | z = trend_bcz) on trend-free series from AR
          processes, including AR models fitted to the observed ACC and Drake-box EKE (SI S4: 3-5 % for monthly series).
Inputs  : outputs/a5/a5_series_allsat.csv (02a), outputs/a6/a6b_drake_box_series.csv (03u) for the 'real_*' processes
Outputs : outputs/stats/size_bc/size_bc_<mode>_<process>.npy; progress printed (the shipped *.log files)
Run     : python scripts/05_stats/05g_size_bc.py ar1_0.8 300 outputs/stats/size_bc z > outputs/stats/size_bc/z_ar1_0.8.log
          (processes: ar1_0.5 ar1_0.8 ar1_0.95 ar2 annual_ar1_0.5 real_ACC_EKE real_Drake_EKE; env RSTATS_MONTHLY_ORDER=12
          for the ord12_* logs)
Seeds   : CHANGED for the publication repository. The working copy seeded the generator with abs(hash(process)) % 2**32;
          Python string hashing is randomised per interpreter process (PYTHONHASHSEED), so those runs cannot be repeated
          exactly. The seed is now zlib.crc32(process) (stable). A re-run therefore draws different synthetic series from
          the shipped logs and the SI S4 false-positive rates may shift slightly (by Monte Carlo error, about +-1.3
          percentage points at 5 % with 300 replicates).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd
import zlib
from osr11 import rstats as R, paths as PA
proc, reps, outdir = sys.argv[1], int(sys.argv[2]), Path(sys.argv[3]); outdir.mkdir(parents=True, exist_ok=True)
MODE = sys.argv[4] if len(sys.argv) > 4 else "single"   # single | cal | z
CAL = MODE == "cal"
rng = np.random.default_rng(zlib.crc32(proc.encode()))   # stable seed (was abs(hash(proc)) % 2**32, not reproducible)
monthly = not proc.startswith("annual")
n = 300 if monthly else 25
x = 2001 + (np.arange(n) + .5) / 12 if monthly else np.arange(2001, 2001 + n, dtype=float)
mo = np.tile(np.arange(1, 13), 25) if monthly else None
if proc.startswith("real"):
    O = PA.OUT
    y0 = (pd.read_csv(O / "a5" / "a5_series_allsat.csv")["ACC_EKE"] if proc == "real_ACC_EKE"
          else pd.read_csv(O / "a6" / "a6b_drake_box_series.csv")["ice_any_EKE"]).values
    y0 = R.deseasonalise(y0, mo); b, a = R.theil_sen(x, y0); r = y0 - a - b * x
    ph, sg = R._fit_ar(r, 12); phis = R._bias_correct_ar(ph, sg, x, 500, rng)
else:
    phis = {"ar1_0.5": [0.5], "ar1_0.8": [0.8], "ar1_0.95": [0.95], "ar2": [1.2, -0.35], "annual_ar1_0.5": [0.5]}[proc]
    phis, sg = np.array(phis, float), 1.0
print(proc, "true phi sum", round(float(np.sum(phis)), 3), flush=True)
out = []
for k in range(reps):
    y = R._simulate_ar(phis, sg, n, 1, rng, burn=500)[0]
    if monthly:
        y = y + 0.5 * np.std(y) * np.sin(2 * np.pi * (mo - 1) / 12)
    if MODE == "z":
        o = R.trend_bcz(x, y, month=mo, monthly=monthly, n_outer=200, seed=k)
    else:
        o = R.trend_bc(x, y, month=mo, monthly=monthly, n_boot=500 if CAL else 1000, calibrate=CAL, n_outer=40, seed=k)
    if CAL:
        o["p"] = o["p_cal"]
    t = R.trend(x, y, month=mo, monthly=monthly, n_boot=1000, seed=k)
    out.append((o["p"], t.p_boot, o["phi_sum"], o["phi_sum_bc"]))
    a = np.array(out)
    if (k + 1) % 25 == 0:
        print(k + 1, "bc", np.mean(a[:, 0] < .05).round(3), "sieve", np.mean(a[:, 1] < .05).round(3),
              "phi_sum", a[:, 2].mean().round(3), "bc", a[:, 3].mean().round(3), flush=True)
np.save(outdir / f"size_bc_{MODE}_{proc}.npy", np.array(out))
