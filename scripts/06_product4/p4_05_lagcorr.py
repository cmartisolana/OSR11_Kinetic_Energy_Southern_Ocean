"""Step 5: lagged correlation of KE, EKE (and the isopycnal slope) with the wind stress, wind leading 0-36 months (main Fig. 4e),
for Product 1 (reproduction) and Product 4.

Purpose : apply scripts/04_wind_argo/04k_lag_correlation.py (the source of Fig. 4e and of Sect. 3.3: KE r = 0.57 at
          9 months, family-wise p = 0.009; EKE r = 0.47 at 11 months, p = 0.11) unchanged to the Product 4 domain series.
Inputs  : outputs/f4e/f4e_series.csv (Product 1, ACC band + 100 km, the x12 input); outputs/p4/p4_wind_series.csv (Product 4, p4_04);
          outputs/p4/p4_wind_series_p4band.csv (sensitivity domain); outputs/x12/x12_summary.csv (check)
Outputs : outputs/p4/p4_lagcorr.csv          product, series, lag_months, r, p_pointwise, env95
          outputs/p4/p4_lagcorr_summary.csv  product, series, lag_max, r_max, p_pt_max, p_fw, lag_max_lo/hi, r at 0/12/24/36, sig_lags
          outputs/p4/p4_lagcorr_check.csv    Product 1 summary here vs x12_summary.csv
Method  : x12 verbatim (function `run` below = the x12 body): each series in % of its 2004-2025 mean, monthly climatology removed
          (2004-2025 for responses, 2001-2025 for the wind), linear trend removed (responses 2004-2025, wind 2001-2025), centred
          12-month running mean; responses = the valid running means from 2004 (Jul 2004 - Jul 2025, 253 months) correlated with the
          wind L = 0-36 months earlier; 2000 phase-randomised surrogates of the deseasonalised monthly wind (2001-2025) processed
          identically -> pointwise p and 95 % envelope; family-wise p of the maximum r against the surrogate maxima of |r| over the 37
          lags; 95 % range of the lag of the maximum from a moving-block bootstrap (24-month blocks, 2000 resamples).
          RNG: default_rng(20261109) re-created for each product, so the wind surrogates are identical for both products (the wind
          series is identical: same cells).
Run     : python scripts/06_product4/p4_05_lagcorr.py
"""
import p4_common as C
import numpy as np, pandas as pd


def run(S, label):
    S = S.sort_values(["year", "month"]).reset_index(drop=True)
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
    valid = resp_idx[np.isfinite(X["EKE"][resp_idx])]
    print(label, "response months:", valid.size, "from", yr[valid[0]], mo[valid[0]], "to", yr[valid[-1]], mo[valid[-1]])
    assert np.all(np.isfinite(W[valid - 36]))
    def lagr(w, x): return np.array([np.corrcoef(w[valid - L], x[valid])[0, 1] for L in LAGS])
    rng = np.random.default_rng(20261109)
    raw_w = S.wind.values.astype(float)
    def surr(y):
        f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
        if y.size % 2 == 0: ph[-1] = 0
        return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size) + y.mean()
    wa = 100 * (raw_w / raw_w.mean() - 1); wa = wa - pd.Series(wa).groupby(mo).transform("mean").values
    WS = []
    for _ in range(2000):
        s_ = surr(wa); t = np.arange(s_.size); s_ = s_ - np.polyval(np.polyfit(t, s_, 1), t); WS.append(pd.Series(s_).rolling(12, center=True).mean().values)
    rows, summ = [], []
    for k, x in X.items():
        r = lagr(W, x); nul = np.array([lagr(s_, x) for s_ in WS])
        p_pt = (1 + np.sum(np.abs(nul) >= np.abs(r)[None, :], 0)) / 2001; env = np.percentile(np.abs(nul), 95, axis=0)
        mx = np.abs(nul).max(1); kmax = int(np.argmax(r)); p_fw = (1 + np.sum(mx >= abs(r[kmax]))) / 2001
        nb = int(np.ceil(valid.size / 24)); lm = []
        for _ in range(2000):
            st = rng.integers(0, valid.size - 24 + 1, nb); ii = np.concatenate([np.arange(s0, s0 + 24) for s0 in st])[:valid.size]; v = valid[ii]
            lm.append(LAGS[int(np.argmax([np.corrcoef(W[v - L], x[v])[0, 1] for L in LAGS]))])
        lo, hi = np.percentile(lm, [2.5, 97.5])
        for L in LAGS: rows.append(dict(product=label, series=k, lag_months=L, r=r[L], p_pointwise=p_pt[L], env95=env[L]))
        summ.append(dict(product=label, series=k, lag_max=LAGS[kmax], r_max=r[kmax], p_pt_max=p_pt[kmax], p_fw=p_fw, lag_max_lo=lo, lag_max_hi=hi,
                         r_lag0=r[0], r_lag12=r[12], r_lag24=r[24], r_lag36=r[36], sig_lags=",".join(str(L) for L in LAGS if p_pt[L] < 0.05)))
        print(f"{label} {k}: max r {r[kmax]:+.3f} at {LAGS[kmax]} months (lag range {lo:.0f}-{hi:.0f}); pointwise p {p_pt[kmax]:.3f}, "
              f"family-wise p {p_fw:.3f}; significant lags {summ[-1]['sig_lags']}", flush=True)
    return rows, summ


if __name__ == "__main__":
    allr, alls = [], []
    for lab, f in (("P1", C.OUT / "f4e" / "f4e_series.csv"), ("P4", C.HERE / "p4_wind_series.csv"),
                   ("P4 (own-band domain)", C.HERE / "p4_wind_series_p4band.csv")):
        r, s = run(pd.read_csv(f), lab); allr += r; alls += s
    pd.DataFrame(allr).to_csv(C.HERE / "p4_lagcorr.csv", index=False)
    S = pd.DataFrame(alls); S.to_csv(C.HERE / "p4_lagcorr_summary.csv", index=False)
    ref = pd.read_csv(C.OUT / "x12" / "x12_summary.csv").set_index("series")
    chk = []
    for _, q in S[S["product"] == "P1"].iterrows():
        rr = ref.loc[q.series]
        chk.append({"series": q.series, **{f"{c}_here": q[c] for c in ("lag_max", "r_max", "p_fw", "lag_max_lo", "lag_max_hi")},
                    **{f"{c}_x12": rr[c] for c in ("lag_max", "r_max", "p_fw", "lag_max_lo", "lag_max_hi")},
                    "sig_lags_identical": str(q.sig_lags) == str(rr.sig_lags) if isinstance(rr.sig_lags, str) else q.sig_lags == ""})
    pd.DataFrame(chk).to_csv(C.HERE / "p4_lagcorr_check.csv", index=False); print(pd.DataFrame(chk).round(4).to_string(index=False))
