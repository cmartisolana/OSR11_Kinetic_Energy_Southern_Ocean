"""x9 (9 Nov 2026): as x8 (power-law exponent b = least-squares slope of ln X on ln tau of the previous year) with the response
period extended to 2004-2025 (22 pairs; wind 2003-2024), the years of Fig. 4c and 4e. Adds a check of the error method:
Monte Carlo coverage of the 95 % intervals (pairs bootstrap, 3-yr block bootstrap, HAC) and size of the surrogate p-value test,
for synthetic responses b*ln tau(y-1) + linear trend + AR(1) noise, with tau from phase-randomised surrogates of the observed ln tau
and trend, noise variance and lag-1 autocorrelation from the observed fit residuals.
Outputs: x9/x9_sensitivity_log.csv, x9/x9_coverage.csv, x9/x9_log.txt

Origin  : x9_wind_sensitivity_log_2004.py (9 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : power-law exponents b of Fig. 4d (Sect. 2.4, 3.3) and their p-values.
Inputs  : outputs/f4e/f4e_series.csv (04g)
Run     : python scripts/04_wind_argo/04i_wind_exponents.py   (no raw data; minutes)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd
from scipy import stats
from osr11 import paths as PA
OUT = PA.out("x9")
S = pd.read_csv(PA.OUT / "f4e" / "f4e_series.csv")
yrs = np.arange(2001, 2026); ann = lambda c: np.array([S[c][S.year == y].mean() for y in yrs])
W = ann("wind"); X = {k: ann(k) for k in ("KE", "EKE", "slope")}
c4 = yrs >= 2004; Y0 = 2004; sel_y = yrs >= Y0; sel_w = (yrs >= Y0 - 1) & (yrs <= 2024)
lw = np.log(W / W.mean()); wx = lw[sel_w]; n = wx.size; t = np.arange(n)
def slope_v(x, y):  # vectorised OLS slope along the last axis
    xc = x - x.mean(-1, keepdims=True); return (xc * (y - y.mean(-1, keepdims=True))).sum(-1) / (xc * xc).sum(-1)
ols = lambda x, y: float(slope_v(x, y))
def detr(y): return y - np.polyval(np.polyfit(np.arange(y.size), y, 1), np.arange(y.size))
def surr(y, rng_):
    f = np.fft.rfft(y - y.mean()); ph = rng_.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size) + y.mean()
def block_idx(n_, L, rng_, m):
    nb = int(np.ceil(n_ / L)); st = rng_.integers(0, n_ - L + 1, (m, nb))
    return (st[:, :, None] + np.arange(L)).reshape(m, -1)[:, :n_]
def hac_ci(x, y, L=2):
    xc = x - x.mean(); b = ols(x, y); e = y - y.mean() - b * xc; u = xc * e; s = np.sum(u * u)
    for l in range(1, L + 1): s += 2 * (1 - l / (L + 1)) * np.sum(u[l:] * u[:-l])
    h = stats.t.ppf(0.975, x.size - 2) * np.sqrt(s) / np.sum(xc * xc); return b - h, b + h
lagw = lambda wf, L: wf[(yrs >= 2004 - L) & (yrs <= 2025 - L)]
NB = 2000
rngA = np.random.default_rng(20261028); SWA = np.array([surr(lw, rngA) for _ in range(NB)])
rngB = np.random.default_rng(20261104); SWB = np.array([surr(lw, rngB) for _ in range(4000)])
rows, log = [], []
P = lambda s: (print(s, flush=True), log.append(s))
for k, v in X.items():
    yy = np.log(v / np.nanmean(v[c4]))[sel_y]
    b = ols(wx, yy); r = np.corrcoef(wx, yy)[0, 1]
    I = rngA.integers(0, n, (NB, n)); bb = slope_v(wx[I], yy[I])
    rs = np.array([np.corrcoef(s_[sel_w], yy)[0, 1] for s_ in SWA]); p = (1 + np.sum(np.abs(rs) >= abs(r))) / (NB + 1)
    yd, wd = detr(yy), detr(wx); bd = ols(wd, yd); rd = np.corrcoef(wd, yd)[0, 1]
    rngD = np.random.default_rng(20261104); I = rngD.integers(0, n, (NB, n)); bbd = slope_v(wd[I], yd[I])
    I3 = block_idx(n, 3, rngB, 4000); bl3 = slope_v(wx[I3], yy[I3]); hl, hh = hac_ci(wx, yy)
    xd_full = detr(np.log(v / np.nanmean(v[c4]))[c4]); bnull = np.empty(4000); rnull = np.empty(4000)
    for i, s_ in enumerate(SWB):
        sd = detr(s_); Ls = (0, 1, 2, 3)[int(np.argmax([abs(np.corrcoef(xd_full, lagw(sd, L))[0, 1]) for L in (0, 1, 2, 3)]))]
        ws = s_[(yrs >= Y0 - Ls) & (yrs <= 2025 - Ls)]; bnull[i] = ols(ws, yy); rnull[i] = np.corrcoef(ws, yy)[0, 1]
    pb_sel = (1 + np.sum(np.abs(bnull) >= abs(b))) / 4001; pr_sel = (1 + np.sum(np.abs(rnull) >= abs(r))) / 4001
    # residual structure of the trend-retained fit with an explicit trend term (for the simulations)
    A = np.c_[np.ones(n), wx, t]; coef, *_ = np.linalg.lstsq(A, yy, rcond=None); e = yy - A @ coef
    phi = np.corrcoef(e[1:], e[:-1])[0, 1]; sig = e.std(ddof=3)
    rows.append(dict(series=k, n=n, b=b, b_lo=np.percentile(bb, 2.5), b_hi=np.percentile(bb, 97.5), r=r, p=p,
                     b_detr=bd, bd_lo=np.percentile(bbd, 2.5), bd_hi=np.percentile(bbd, 97.5), r_detr=rd,
                     blk3_lo=np.percentile(bl3, 2.5), blk3_hi=np.percentile(bl3, 97.5), hac_lo=hl, hac_hi=hh,
                     p_b_lag_reselected=pb_sel, p_r_lag_reselected=pr_sel, shared_trend_frac=(b - bd) / b,
                     resid_ar1=phi, resid_sd=sig, trend_per_yr=coef[2]))
    R_ = rows[-1]
    P(f"{k}: b = {b:.2f} [{R_['b_lo']:.2f}, {R_['b_hi']:.2f}] r {r:.2f} p {p:.3f} | detrended {bd:.2f} [{R_['bd_lo']:.2f}, {R_['bd_hi']:.2f}] r {rd:.2f} | "
      f"block3 [{R_['blk3_lo']:.2f}, {R_['blk3_hi']:.2f}] HAC [{hl:.2f}, {hh:.2f}] | p lag re-selected: b {pb_sel:.3f}, r {pr_sel:.3f} | "
      f"shared-trend {100*(b-bd)/b:.0f} % | residual AR(1) {phi:.2f}")
pd.DataFrame(rows).to_csv(OUT / "x9_sensitivity_log.csv", index=False)
# ---------------- Monte Carlo check of the error method
rng = np.random.default_rng(20261109); NS, NBs = 1000, 1000; cov = []
for R_ in rows:
    for b_true in (R_["b"], 0.0):
        hit = dict(pairs=0, block3=0, hac=0, pairs_detr=0)
        for s in range(NS):
            ws = surr(lw, rng)[sel_w]
            e = np.empty(n); e[0] = rng.normal(0, R_["resid_sd"] / np.sqrt(max(1 - R_["resid_ar1"] ** 2, 1e-3)))
            for i in range(1, n): e[i] = R_["resid_ar1"] * e[i - 1] + rng.normal(0, R_["resid_sd"])
            ys = b_true * ws + R_["trend_per_yr"] * t + e
            I = rng.integers(0, n, (NBs, n)); q = np.percentile(slope_v(ws[I], ys[I]), [2.5, 97.5]); hit["pairs"] += q[0] <= b_true <= q[1]
            I3 = block_idx(n, 3, rng, NBs); q = np.percentile(slope_v(ws[I3], ys[I3]), [2.5, 97.5]); hit["block3"] += q[0] <= b_true <= q[1]
            q = hac_ci(ws, ys); hit["hac"] += q[0] <= b_true <= q[1]
            wsd, ysd = detr(ws), detr(ys); I = rng.integers(0, n, (NBs, n)); q = np.percentile(slope_v(wsd[I], ysd[I]), [2.5, 97.5]); hit["pairs_detr"] += q[0] <= b_true <= q[1]
        row = dict(series=R_["series"], b_true=round(b_true, 3), **{f"coverage_{m}": h / NS for m, h in hit.items()})
        cov.append(row); P(f"  MC {R_['series']} b_true={b_true:.2f}: coverage pairs {row['coverage_pairs']:.3f}, block3 {row['coverage_block3']:.3f}, HAC {row['coverage_hac']:.3f}, pairs (detrended data, true b) {row['coverage_pairs_detr']:.3f}")
pd.DataFrame(cov).to_csv(OUT / "x9_coverage.csv", index=False); (OUT / "x9_log.txt").write_text("\n".join(log))
