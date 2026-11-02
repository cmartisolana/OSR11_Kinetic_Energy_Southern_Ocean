"""Section statistics from b1_maps.py path files (annual values, 2001-2025).

Per meander and section (upstream / meander / downstream, as in b1_meander_run.py):
  mean path latitude (annual mean of monthly path latitude, averaged over the section's 0.5° bins)
  within-year SD of path latitude (section mean)
Trends: Theil-Sen on 25 annual values, AR-sieve p (rstats.trend, monthly=False), BH-FDR over all tests.
Epoch difference 2017-25 minus 2001-08.
-> outputs/b1/b1_paths_stats.csv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R
SECS = {"PAR": {"upstream": (-150, -135.13), "meander": (-135.13, -109.38), "downstream": (-109.38, -80)},
        "SEIR150": {"upstream": (125, 138), "meander": (140, 150), "downstream": (150, 158)},
        "KERG": {"upstream": (70, 77), "meander": (76, 84), "downstream": (84, 94)}}
rows = []
yrs = np.arange(2001, 2026)
for n, secs in SECS.items():
    P = xr.open_dataset(HERE / "outputs" / "b1" / f"paths_{n}.nc")
    ann = pd.DataFrame(P.path_lat.values).groupby(np.repeat(yrs, 12)).mean().values
    for s, (a, b) in secs.items():
        m = (P.lon.values >= a) & (P.lon.values <= b)
        for var, Y in (("path_lat_deg", np.nanmean(ann[:, m], 1)), ("within_year_sd_deg", np.nanmean(P.annual_sd.values[:, m], 1))):
            T = R.trend(yrs + 0.5, Y, monthly=False, n_boot=2000)
            e = float(np.mean(Y[yrs >= 2017]) - np.mean(Y[yrs <= 2008]))
            rows.append(dict(meander=n, section=s, var=var, mean=float(np.mean(Y)), interannual_sd=float(np.std(Y, ddof=1)),
                             slope_per_decade=10 * T.slope, ci_lo=10 * T.ci_lo, ci_hi=10 * T.ci_hi, p=T.p_boot,
                             epoch_diff_late_minus_early=e))
D = pd.DataFrame(rows)
D["fdr_sig"] = R.bh_fdr(D.p.values, 0.10)
D.to_csv(HERE / "outputs" / "b1" / "b1_paths_stats.csv", index=False)
pd.set_option("display.width", 200)
print(D.round(3).to_string())
