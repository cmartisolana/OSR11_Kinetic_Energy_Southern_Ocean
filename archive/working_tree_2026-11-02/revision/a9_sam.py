"""A9: SAM and the ACC-band KE/EKE/transport-proxy series (Referee A, specific comment).

Specification fixed before results (2026-10-14): monthly series 2001-2025 (all-sat), deseasonalised and Theil-Sen
detrended as in b1_covar.prep (no smoothing, and with the 2-month running mean as a second version); Marshall (2003)
SAM index prepared identically; lag correlation with SAM leading by 0-12 months; significance from 2000 phase-randomised
SAM surrogates, for each lag and for the maximum |r| over lags 0-12 (multiplicity). Also the correlation between the
annual means and of the 12-month running means (the time scale on which the manuscript speaks of 'SAM-like variability').
Usage: python a9_sam.py <out_dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import b1_covar as BC, rstats as R
out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True); rng = np.random.default_rng(20261014)
sam = pd.read_csv(HERE / "data" / "sam_marshall_monthly.txt", sep=r"\s+", comment="#", header=None, na_values="NA")
sam.columns = ["year"] + list(range(1, 13)); sam = sam.melt(id_vars="year", var_name="month", value_name="sam").sort_values(["year", "month"]); sam["month"] = sam.month.astype(int); sam["sam"] = sam.sam.astype(float)
a5 = pd.read_csv(HERE / "outputs" / "a5" / "a5_series_allsat.csv"); a8 = pd.read_csv(HERE / "outputs" / "a8" / "a8_series_allsat.csv")
S = a5.merge(sam, on=["year", "month"]).merge(a8[["year", "month", "ACC_T"]], on=["year", "month"])
S["x"] = S.year + (S.month - 0.5) / 12
assert len(S) == 300
def prep_raw(P, col):
    y = R.deseasonalise(P[col].values.astype(float), P.month.values); ok = np.isfinite(y)
    b, a = R.theil_sen(P.x.values[ok], y[ok]); return y - (a + b * P.x.values)
rows = []
for col in ("ACC_KE", "ACC_EKE", "ACC_SPD", "ACC_T"):
    if col not in S: continue
    for smooth, fx, fy in (("none", prep_raw(S, "sam"), prep_raw(S, col)), ("2-month", BC.prep(S, "sam"), BC.prep(S, col))):
        rr = [BC.lagcorr(fx, fy, L)[0] for L in range(13)]; obs = max(abs(r) for r in rr)
        ps = []; smax = []
        surr = [BC.phase_surr(fx, rng) for _ in range(2000)]
        for L in range(13):
            s = np.array([BC.lagcorr(sx, fy, L)[0] for sx in surr]); ps.append((1 + np.sum(np.abs(s) >= abs(rr[L]))) / 2001)
        smax = np.array([max(abs(BC.lagcorr(sx, fy, L)[0]) for L in range(13)) for sx in surr])
        pmax = (1 + np.sum(smax >= obs)) / 2001
        best = int(np.argmax(np.abs(rr)))
        rows.append(dict(series=col, smoothing=smooth, **{f"r_lag{L}": rr[L] for L in range(13)}, **{f"p_lag{L}": ps[L] for L in range(13)}, best_lag=best, r_best=rr[best], p_max=pmax))
        print(f"{col:8s} {smooth:8s} r(lag0)={rr[0]:+.2f} p={ps[0]:.3f} | best lag {best} r={rr[best]:+.2f} p_max={pmax:.3f}", flush=True)
    # annual means and 12-month running means
    A = S.groupby("year")[[col, "sam"]].mean(); ya = prep_raw(A.assign(x=A.index.values, month=1), col) if False else None
    ra = np.corrcoef(np.diff(A[col].values) if False else A[col].values - np.polyval(np.polyfit(A.index, A[col], 1), A.index), A.sam.values - np.polyval(np.polyfit(A.index, A.sam, 1), A.index))[0, 1]
    xr_ = pd.Series(prep_raw(S, "sam")).rolling(12, center=True).mean().values; yr_ = pd.Series(prep_raw(S, col)).rolling(12, center=True).mean().values
    ok = np.isfinite(xr_) & np.isfinite(yr_); r12 = np.corrcoef(xr_[ok], yr_[ok])[0, 1]
    n_eff = ok.sum() * (1 - BC.r1(xr_[ok]) * BC.r1(yr_[ok])) / (1 + BC.r1(xr_[ok]) * BC.r1(yr_[ok]))
    from scipy import stats
    t = r12 * np.sqrt(max(n_eff - 2, 1) / (1 - r12 ** 2)); p12 = 2 * stats.t.sf(abs(t), max(n_eff - 2, 1))
    rows.append(dict(series=col, smoothing="annual", r_lag0=ra, best_lag=0, r_best=ra)); rows.append(dict(series=col, smoothing="12-month", r_lag0=r12, p_lag0=p12, best_lag=0, r_best=r12, n_eff=n_eff))
    print(f"{col:8s} annual r={ra:+.2f} (n=25) | 12-month running r={r12:+.2f} n_eff={n_eff:.0f} p={p12:.3f}", flush=True)
pd.DataFrame(rows).to_csv(out / "a9_sam_correlations.csv", index=False); print("done")
