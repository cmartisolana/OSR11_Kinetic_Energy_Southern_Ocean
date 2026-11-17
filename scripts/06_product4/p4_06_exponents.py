"""Step 6: power-law exponent b of annual KE and EKE with respect to the previous year's wind stress (main Fig. 4d), trends retained
and both series detrended, for Product 1 (reproduction) and Product 4.

Purpose : apply scripts/04_wind_argo/04i_wind_exponents.py (b, r, p-values) and 04j_exponent_intervals.py
          (95 % intervals), the sources of Fig. 4d and Sect. 3.3 (b = 0.90 EKE, 0.50 KE; detrended 0.55, 0.29), unchanged.
Inputs  : outputs/f4e/f4e_series.csv (P1), outputs/p4/p4_wind_series.csv (P4, p4_04), outputs/p4/p4_wind_series_p4band.csv (sensitivity);
          outputs/x9/x9_sensitivity_log.csv, x10/x10_intervals.csv (check)
Outputs : outputs/p4/p4_exponents.csv        product, series, b, b_lo, b_hi, r, p, b_detr, bd_lo, bd_hi, r_detr, p_b_lag_reselected,
                                     p_r_lag_reselected, shared_trend_frac, p_bd_extra
          outputs/p4/p4_exponents_check.csv  Product 1 values here vs x9/x10
Method  : annual means 2001-2025 of the domain series; response ln(X / mean_2004-2025 X) for 2004-2025, regressor ln(tau / mean_2001-2025
          tau) of the previous year (2003-2024; 22 pairs); b = OLS slope; detrended: both series linearly detrended over the 22 years.
          p (x9): 2000 phase-randomised surrogates of ln tau (default_rng(20261028)), |r| of the surrogate pairs vs the observed |r| at
          the fixed one-year lag; lag-reselected p (x9, SI S5.2): 4000 surrogates (default_rng(20261104)) re-selecting the annual lag
          0-3 yr of maximum detrended |r|. Intervals (x10): 4000 simulations of ln X = a + b ln tau + c t + AR(1) noise (lag-1 coefficient
          bias-corrected phi + (1 + 3 phi)/n) with phase-randomised ln tau, basic interval widened by t(0.975, 19)/z(0.975) = 1.07;
          default_rng(20261109) used sequentially for KE, EKE, slope as in x10 (the x10 nested coverage loop, which follows, is not run).
          p_bd_extra (not in x9): p of the detrended exponent from the same 2000 x9 surrogates, both detrended (as x14 does).
Run     : python scripts/06_product4/p4_06_exponents.py
"""
import p4_common as C
import numpy as np, pandas as pd
from scipy import stats

yrs = np.arange(2001, 2026)


def slope_v(x, y):
    xc = x - x.mean(-1, keepdims=True); return (xc * (y - y.mean(-1, keepdims=True))).sum(-1) / (xc * xc).sum(-1)


def run(S, label):
    ann = lambda c: np.array([S[c][S.year == y].mean() for y in yrs])
    W = ann("wind"); X = {k: ann(k) for k in ("KE", "EKE", "slope")}
    c4 = yrs >= 2004; Y0 = 2004; sel_y = yrs >= Y0; sel_w = (yrs >= Y0 - 1) & (yrs <= 2024)
    lw = np.log(W / W.mean()); wx = lw[sel_w]; n = wx.size; t = np.arange(n)
    ols = lambda x, y: float(slope_v(x, y))
    def detr(y): return y - np.polyval(np.polyfit(np.arange(y.size), y, 1), np.arange(y.size))
    def surr(y, rng_):
        f = np.fft.rfft(y - y.mean()); ph = rng_.uniform(0, 2 * np.pi, f.size); ph[0] = 0
        if y.size % 2 == 0: ph[-1] = 0
        return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size) + y.mean()
    lagw = lambda wf, L: wf[(yrs >= 2004 - L) & (yrs <= 2025 - L)]
    NB = 2000
    rngA = np.random.default_rng(20261028); SWA = np.array([surr(lw, rngA) for _ in range(NB)])
    rngB = np.random.default_rng(20261104); SWB = np.array([surr(lw, rngB) for _ in range(4000)])
    # ---- x9 part
    rows = {}
    for k, v in X.items():
        yy = np.log(v / np.nanmean(v[c4]))[sel_y]
        b = ols(wx, yy); r = np.corrcoef(wx, yy)[0, 1]
        rs = np.array([np.corrcoef(s_[sel_w], yy)[0, 1] for s_ in SWA]); p = (1 + np.sum(np.abs(rs) >= abs(r))) / (NB + 1)
        yd, wd = detr(yy), detr(wx); bd = ols(wd, yd); rd = np.corrcoef(wd, yd)[0, 1]
        bsd = np.array([ols(detr(s_[sel_w]), yd) for s_ in SWA]); p_bd = (1 + np.sum(np.abs(bsd) >= abs(bd))) / (NB + 1)
        xd_full = detr(np.log(v / np.nanmean(v[c4]))[c4]); bnull = np.empty(4000); rnull = np.empty(4000)
        for i, s_ in enumerate(SWB):
            sd = detr(s_); Ls = (0, 1, 2, 3)[int(np.argmax([abs(np.corrcoef(xd_full, lagw(sd, L))[0, 1]) for L in (0, 1, 2, 3)]))]
            ws = s_[(yrs >= Y0 - Ls) & (yrs <= 2025 - Ls)]; bnull[i] = ols(ws, yy); rnull[i] = np.corrcoef(ws, yy)[0, 1]
        rows[k] = dict(product=label, series=k, n=n, b=b, r=r, p=p, b_detr=bd, r_detr=rd,
                       p_b_lag_reselected=(1 + np.sum(np.abs(bnull) >= abs(b))) / 4001,
                       p_r_lag_reselected=(1 + np.sum(np.abs(rnull) >= abs(r))) / 4001, shared_trend_frac=(b - bd) / b, p_bd_extra=p_bd)
    # ---- x10 part (verbatim functions)
    def detr_v(y):
        tc = t - t.mean(); return y - y.mean(-1, keepdims=True) - np.outer(np.atleast_2d(y)[:, :] @ tc / (tc @ tc), tc).reshape(y.shape)
    F0 = np.fft.rfft(lw - lw.mean()); AMP = np.abs(F0)
    def surr_v(m, rng):
        ph = rng.uniform(0, 2 * np.pi, (m, F0.size)); ph[:, 0] = 0
        if lw.size % 2 == 0: ph[:, -1] = 0
        return np.fft.irfft(AMP * np.exp(1j * ph), n=lw.size, axis=-1)[:, sel_w] + lw.mean()
    def ar1_v(m, phi, sig, rng):
        e = np.empty((m, n)); e[:, 0] = rng.normal(0, sig / np.sqrt(1 - phi ** 2), m)
        for i in range(1, n): e[:, i] = phi * e[:, i - 1] + rng.normal(0, sig, m)
        return e
    def fit(ws, ys):
        A = np.c_[np.ones(n), ws, t]; coef, *_ = np.linalg.lstsq(A, ys, rcond=None); e = ys - A @ coef
        phi = np.corrcoef(e[1:], e[:-1])[0, 1]; sig = e.std(ddof=3) * np.sqrt(1 - min(phi, 0.95) ** 2)
        return coef, phi, sig
    TF = stats.t.ppf(0.975, n - 3) / stats.norm.ppf(0.975)
    def interval(ws, ys, rng, m=4000):
        coef, phi, sig = fit(ws, ys); phic = float(np.clip(phi + (1 + 3 * phi) / n, -0.9, 0.95))
        b = float(slope_v(ws, ys)); bd = float(slope_v(detr_v(ws), detr_v(ys)))
        sw = surr_v(m, rng); e = ar1_v(m, phic, sig, rng)
        ysim = b * sw + coef[2] * t + e
        bs = slope_v(sw, ysim); lo, hi = b - TF * np.percentile(bs - b, [97.5, 2.5])
        ysd = bd * sw + coef[2] * t + e
        bsd = slope_v(detr_v(sw), detr_v(ysd)); lod, hid = bd - TF * np.percentile(bsd - bd, [97.5, 2.5])
        return b, lo, hi, bd, lod, hid, phi, phic
    rng = np.random.default_rng(20261109)
    for k, v in X.items():
        yy = np.log(v / np.nanmean(v[c4]))[sel_y]
        b, lo, hi, bd, lod, hid, phi, phic = interval(wx, yy, rng)
        rows[k].update(b_lo=lo, b_hi=hi, bd_lo=lod, bd_hi=hid, resid_ar1=phi, resid_ar1_corrected=phic)
        q = rows[k]
        print(f"{label} {k}: b {q['b']:.3f} [{lo:.2f}, {hi:.2f}] r {q['r']:.2f} p {q['p']:.3f} | detrended {bd:.3f} [{lod:.2f}, {hid:.2f}] "
              f"r {q['r_detr']:.2f} p_extra {q['p_bd_extra']:.3f} | lag-reselected p {q['p_b_lag_reselected']:.3f}", flush=True)
    return list(rows.values())


if __name__ == "__main__":
    rows = []
    for lab, f in (("P1", C.OUT / "f4e" / "f4e_series.csv"), ("P4", C.HERE / "p4_wind_series.csv"),
                   ("P4 (own-band domain)", C.HERE / "p4_wind_series_p4band.csv")):
        rows += run(pd.read_csv(f), lab)
    cols = ["product", "series", "n", "b", "b_lo", "b_hi", "r", "p", "b_detr", "bd_lo", "bd_hi", "r_detr", "p_bd_extra",
            "p_b_lag_reselected", "p_r_lag_reselected", "shared_trend_frac", "resid_ar1", "resid_ar1_corrected"]
    D = pd.DataFrame(rows)[cols]; D.to_csv(C.HERE / "p4_exponents.csv", index=False)
    x9 = pd.read_csv(C.OUT / "x9" / "x9_sensitivity_log.csv").set_index("series")
    x10 = pd.read_csv(C.OUT / "x10" / "x10_intervals.csv").set_index("series")
    chk = []
    for _, q in D[D["product"] == "P1"].iterrows():
        a, c = x9.loc[q.series], x10.loc[q.series]
        chk.append(dict(series=q.series, b=q.b, b_x9=a.b, p=q.p, p_x9=a.p, b_detr=q.b_detr, b_detr_x9=a.b_detr,
                        p_b_lag_reselected=q.p_b_lag_reselected, p_b_lag_reselected_x9=a.p_b_lag_reselected,
                        b_lo=q.b_lo, b_lo_x10=c.b_lo, b_hi=q.b_hi, b_hi_x10=c.b_hi, bd_lo=q.bd_lo, bd_lo_x10=c.bd_lo, bd_hi=q.bd_hi, bd_hi_x10=c.bd_hi))
    chk = pd.DataFrame(chk); chk.to_csv(C.HERE / "p4_exponents_check.csv", index=False); print(chk.round(4).T.to_string())
