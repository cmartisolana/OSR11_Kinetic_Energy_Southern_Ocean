"""x12 (9 Nov 2026): lagged correlation of KE, EKE and the sigma0 = 27.3 slope with the wind stress at monthly resolution (Fig. 4e),
replacing the annual lags 0-3 yr of 04g (D2).
Series: f4e_series.csv (monthly means over the common cells of the ACC band + 100 km). Each series: anomalies in % of its 2004-2025
mean, mean seasonal cycle removed (2004-2025 climatology for the responses, 2001-2025 for the wind), linear trend removed (responses
over 2004-2025, wind over 2001-2025, as 04g), then a centred 12-month running mean. Response months: the valid running means of
2004-2025 (Jul 2004-Jun 2025, 252 months); wind leading by L = 0-36 months on the same response months (wind valid from Jul 2001).
Significance: 2000 phase-randomised surrogates of the monthly wind anomaly (2001-2025), processed identically; pointwise p and 95 %
envelope at each lag; family-wise p from the maximum |r| over the 37 lags. Lag of the maximum: 95 % range from a moving-block bootstrap
of the response months (24-month blocks, 2000 resamples, wind pairs kept).
Outputs: x12/x12_lagcorr.csv (series, lag_months, r, p_pointwise, env95), x12/x12_summary.csv, x12/x12_log.txt

Origin  : x12_lagcorr_monthly.py (9 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : lagged correlations of Fig. 4e (Sect. 3.3; SI S1.6, S5.2).
Inputs  : outputs/f4e/f4e_series.csv (04g)
Run     : python scripts/04_wind_argo/04k_lag_correlation.py   (no raw data; minutes)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd
from osr11 import paths as PA
OUT = PA.out("x12")
S = pd.read_csv(PA.OUT / "f4e" / "f4e_series.csv").sort_values(["year", "month"]).reset_index(drop=True)
yr, mo = S.year.values, S.month.values; c4 = yr >= 2004
def prep(x, sel):
    x = x.astype(float).copy(); out = np.full(x.size, np.nan)
    a = 100 * (x[sel] / np.nanmean(x[sel]) - 1); m = mo[sel]
    a = a - pd.Series(a).groupby(m).transform("mean").values
    t = np.arange(a.size); a = a - np.polyval(np.polyfit(t, a, 1), t)
    out[sel] = pd.Series(a).rolling(12, center=True).mean().values; return out
W = prep(S.wind.values, np.ones(len(S), bool))
LAGS = np.arange(37)
resp_idx = np.where(c4)[0]
X = {k: prep(S[k].values, c4) for k in ("KE", "EKE", "slope")}
valid = resp_idx[np.isfinite(X["EKE"][resp_idx])]; print("response months:", valid.size, "from", yr[valid[0]], mo[valid[0]], "to", yr[valid[-1]], mo[valid[-1]])
assert np.all(np.isfinite(W[valid - 36]))
def lagr(w, x): return np.array([np.corrcoef(w[valid - L], x[valid])[0, 1] for L in LAGS])
rng = np.random.default_rng(20261109)
raw_w = S.wind.values.astype(float)
def surr(y):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size) + y.mean()
# surrogates of the deseasonalised wind anomaly (seasonal cycle re-added is irrelevant: it is removed again)
wa = 100 * (raw_w / raw_w.mean() - 1); wa = wa - pd.Series(wa).groupby(mo).transform("mean").values
WS = []
for _ in range(2000):
    s_ = surr(wa); t = np.arange(s_.size); s_ = s_ - np.polyval(np.polyfit(t, s_, 1), t); WS.append(pd.Series(s_).rolling(12, center=True).mean().values)
rows, summ, log = [], [], []
P = lambda s: (print(s, flush=True), log.append(s))
for k, x in X.items():
    r = lagr(W, x); nul = np.array([lagr(s_, x) for s_ in WS])
    p_pt = (1 + np.sum(np.abs(nul) >= np.abs(r)[None, :], 0)) / 2001; env = np.percentile(np.abs(nul), 95, axis=0)
    mx = np.abs(nul).max(1); kmax = int(np.argmax(r)); p_fw = (1 + np.sum(mx >= abs(r[kmax]))) / 2001
    # block bootstrap of the lag of the maximum
    nb = int(np.ceil(valid.size / 24)); lm = []
    for _ in range(2000):
        st = rng.integers(0, valid.size - 24 + 1, nb); ii = np.concatenate([np.arange(s0, s0 + 24) for s0 in st])[:valid.size]; v = valid[ii]
        lm.append(LAGS[int(np.argmax([np.corrcoef(W[v - L], x[v])[0, 1] for L in LAGS]))])
    lo, hi = np.percentile(lm, [2.5, 97.5])
    for L in LAGS: rows.append(dict(series=k, lag_months=L, r=r[L], p_pointwise=p_pt[L], env95=env[L]))
    summ.append(dict(series=k, lag_max=LAGS[kmax], r_max=r[kmax], p_pt_max=p_pt[kmax], p_fw=p_fw, lag_max_lo=lo, lag_max_hi=hi, r_lag0=r[0], r_lag12=r[12], r_lag24=r[24], r_lag36=r[36],
                     sig_lags=",".join(str(L) for L in LAGS if p_pt[L] < 0.05)))
    P(f"{k}: max r {r[kmax]:+.2f} at {LAGS[kmax]} months (95 % range of lag {lo:.0f}-{hi:.0f}); pointwise p {p_pt[kmax]:.3f}, family-wise p {p_fw:.3f}; "
      f"r at 0/12/24/36 = {r[0]:+.2f}/{r[12]:+.2f}/{r[24]:+.2f}/{r[36]:+.2f}; pointwise-significant lags: {summ[-1]['sig_lags']}")
pd.DataFrame(rows).to_csv(OUT / "x12_lagcorr.csv", index=False); pd.DataFrame(summ).to_csv(OUT / "x12_summary.csv", index=False)
(OUT / "x12_log.txt").write_text("\n".join(log))
