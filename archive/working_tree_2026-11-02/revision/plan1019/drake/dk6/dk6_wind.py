"""dk6 T5: interannual sensitivity of the passage-mean half-power latitude to eastward wind stress and the shift expected
from the observed wind trend. Annual means, 2001-2025; wind leads by 0-3 years. p from 5000 phase-randomised surrogates of the
half-power series; slope interval from OLS with an AR(1) effective sample size; expected shift = slope x wind trend (Monte Carlo)."""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[3])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import numpy as np, pandas as pd
OUT = "out/dk6"; rng = np.random.default_rng(1)
S1 = pd.read_csv("out/dk1/dk1_series.csv"); TR1 = pd.read_csv("out/dk1/dk1_trends.csv"); TR1 = TR1[TR1.period == "2001-2025"].set_index("name")
REV = _REV
W6 = pd.read_csv(f"{REV}/outputs/a6/a6_wind_stress_series.csv"); RT = pd.read_csv(f"{REV}/outputs/retest/retest_summary.csv").set_index("key")
winds = {"Drake box": (S1.groupby("year").wind_box_TAUX.mean(), TR1.loc["wind_box_TAUX", ["pct", "pct_lo", "pct_hi"]].values),
         "upstream 100-80W": (S1.groupby("year").wind_up_TAUX.mean(), TR1.loc["wind_up_TAUX", ["pct", "pct_lo", "pct_hi"]].values)}
r = RT.loc["tau_Circumpolar_45_60S|2001-2025"]; f = 1000 / (r.slope * 1000 / r.pct_per_decade)
winds["circumpolar 45-60S"] = (W6.groupby("year").Circumpolar_45_60S.mean(), np.array([r.pct_per_decade, r.ci_lo * f, r.ci_hi * f]))
def surr(y, n):
    Y = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, (n, Y.size)); ph[:, 0] = 0
    if y.size % 2 == 0: ph[:, -1] = 0
    return np.fft.irfft(np.abs(Y) * np.exp(1j * ph), n=y.size)
rows = []
for prod in ("allsat", "twosat"):
    S = pd.read_csv(f"{OUT}/dk6_series_{prod}_ice0.csv"); col = "0.5|T1 half-power latitude (km, +N)"
    hp = S.groupby(np.repeat(np.arange(2001, 2026), 12))[col].mean()
    obs = pd.read_csv(f"{OUT}/dk6_trends_{prod}_ice0.csv"); obs = obs[(obs.coverage_min == 0.5) & (obs.series == "T1 half-power latitude (km, +N)")].iloc[0]
    for wn, (w, wt) in winds.items():
        wp = 100 * (w / w.mean() - 1)
        for L in range(4):
            y = hp.values[L:]; x = wp.values[:len(wp) - L]; n = y.size
            b, a = np.polyfit(x, y, 1); rr = np.corrcoef(x, y)[0, 1]
            res = y - (a + b * x); r1 = np.corrcoef(x[1:], x[:-1])[0, 1]; r2 = np.corrcoef(y[1:], y[:-1])[0, 1]
            neff = max(3.0, n * (1 - r1 * r2) / (1 + r1 * r2)); se = np.sqrt(np.sum(res ** 2) / (neff - 2) / np.sum((x - x.mean()) ** 2))
            sr = np.array([np.corrcoef(x, s)[0, 1] for s in surr(y, 5000)]); p = (np.abs(sr) >= abs(rr)).mean()
            bs = rng.normal(b, se, 20000); ws = rng.normal(wt[0], (wt[2] - wt[1]) / 3.92, 20000); e = bs * ws
            rows.append(dict(product=prod, wind=wn, lag_yr=L, n=n, r=rr, p_surr=p, slope_km_per_pct=b, slope_lo=b - 1.96 * se, slope_hi=b + 1.96 * se,
                             wind_trend_pct_dec=wt[0], expected_km_dec=np.median(e), exp_lo=np.percentile(e, 2.5), exp_hi=np.percentile(e, 97.5),
                             observed_km_dec=obs.trend_dec, obs_lo=obs.lo, obs_hi=obs.hi))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/dk6_wind.csv", index=False); pd.set_option("display.width", 250); print(R.round(3).to_string())
