"""3. 'Eddy-saturation diagram': elasticities d ln X / d ln tau from annual means.
X in {EKE, KE, SPD (along-stream speed), T (surface transport proxy)}, circumpolar and sectors; tau = sector-matched
45-60S band mean ('z4560') or ACC-mask mean ('accmask'). Lags 0, 1, 2 yr fixed a priori (s3 peaks for EKE fall at
~0.5-1 yr, T at 0); response at year t+L regressed on tau at year t.
(a) 'raw'      : ln X on ln tau (trend + interannual), OLS slope, HAC (Newey-West, 1 lag) SE, t(n-2) 95% CI
(b) 'detrended': both Theil-Sen-detrended first (interannual sensitivity only)
(c) 'trend ratio': (trend of ln X) / (trend of ln tau), annual trend_bcz on each; CI by Monte Carlo from the two
    (assumed normal) trend uncertainties, reported as 2.5/97.5 percentiles plus P(tau trend <= 0).
"""
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
from common import OUT, R

D = pd.read_csv(OUT / "s3_all_series.csv")
A = D.groupby("year").mean(numeric_only=True)
xa = A.index.values + 0.5
rng = np.random.default_rng(3)


def detr(y):
    b, a = R.theil_sen(xa, y); return y - (a + b * xa)


def ols(x, y):
    X = sm.add_constant(x); f = sm.OLS(y, X).fit(cov_type="HAC", cov_kwds={"maxlags": 1})
    b, se = f.params[1], f.bse[1]; tq = stats.t.ppf(0.975, len(y) - 2)
    return b, b - tq * se, b + tq * se, 2 * stats.t.sf(abs(b / se), len(y) - 2), f.rsquared


trend_cache = {}


def ln_trend(col):
    if col not in trend_cache:
        o = R.trend_bcz(xa, np.log(A[col].values), monthly=False, n_outer=300)
        trend_cache[col] = (o["slope"], (o["ci_hi"] - o["ci_lo"]) / (2 * 1.96), o["p"])
    return trend_cache[col]


rows = []
for sec in ("ACC", "Atlantic", "Indian", "Pacific"):
    for tdef in ("z4560", "accmask"):
        tc = f"tau_{tdef}_{sec}"; lt = np.log(A[tc].values)
        bt, st, pt = ln_trend(tc)
        for v in ("EKE", "KE", "SPD", "T"):
            col = f"{sec}_{v}"; ly = np.log(A[col].values)
            bx, sx, px = ln_trend(col)
            dt = rng.normal(bt, st, 20000); dx = rng.normal(bx, sx, 20000); ratio = dx / dt
            base = dict(sector=sec, tau=tdef, var=v, trend_lnX_pct_dec=1000 * bx, trend_lnX_p=px,
                        trend_lntau_pct_dec=1000 * bt, trend_lntau_p=pt,
                        trend_ratio=bx / bt, trend_ratio_lo=np.percentile(ratio, 2.5), trend_ratio_hi=np.percentile(ratio, 97.5),
                        P_tau_trend_le0=float(np.mean(dt <= 0)))
            for L in (0, 1, 2):
                xr_, yr_ = lt[:len(lt) - L], ly[L:]
                for ver in ("raw", "detrended"):
                    if ver == "detrended":
                        x_, y_ = detr(lt)[:len(lt) - L], detr(ly)[L:]
                    else:
                        x_, y_ = xr_, yr_
                    b, lo, hi, p, r2 = ols(x_, y_)
                    rows.append(dict(**base, lag_yr=L, version=ver, elasticity=b, ci_lo=lo, ci_hi=hi, p=p, r2=r2, n=len(y_)))
E = pd.DataFrame(rows); E.to_csv(OUT / "s4_elasticities.csv", index=False)
pd.set_option("display.width", 250)
print(E[E.tau == "accmask"][["sector", "var", "lag_yr", "version", "elasticity", "ci_lo", "ci_hi", "p", "r2"]].round(2).to_string())
TR = E.drop_duplicates(["sector", "tau", "var"])[["sector", "tau", "var", "trend_lnX_pct_dec", "trend_lnX_p", "trend_lntau_pct_dec",
                                                   "trend_lntau_p", "trend_ratio", "trend_ratio_lo", "trend_ratio_hi", "P_tau_trend_le0"]]
TR.to_csv(OUT / "s4_trend_ratio.csv", index=False); print(TR.round(3).to_string())

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for k, sec in enumerate(("ACC", "Atlantic", "Indian", "Pacific")):
    for j, v in enumerate(("EKE", "KE", "SPD", "T")):
        for m, (L, ver) in enumerate(((0, "raw"), (1, "raw"), (0, "detrended"), (1, "detrended"))):
            q = E[(E.sector == sec) & (E["var"] == v) & (E.tau == "accmask") & (E.lag_yr == L) & (E.version == ver)].iloc[0]
            xp = j * 5 + m + k * 0.18
            ax[0].errorbar(xp, q.elasticity, yerr=[[q.elasticity - q.ci_lo], [q.ci_hi - q.elasticity]], fmt="o^sD"[k], color=f"C{m}",
                           ms=4, lw=.8, label=f"{sec}" if (j == 0 and m == 0) else None)
ax[0].axhline(1, color="k", ls="--", lw=.7); ax[0].axhline(0, color="k", lw=.7)
ax[0].set_xticks([1.5, 6.5, 11.5, 16.5]); ax[0].set_xticklabels(["EKE", "KE", "SPD", "T"])
ax[0].set_ylabel("d ln X / d ln τ  (95% CI)"); ax[0].legend(fontsize=7, title="marker = sector", title_fontsize=7)
ax[0].set_title("Annual elasticities vs ACC-mask τ\ncolours: lag0 raw, lag1 raw, lag0 detr., lag1 detr.", fontsize=9)
ax[0].set_ylim(-4, 6)
# scatter diagram circumpolar
lt = np.log(A["tau_accmask_ACC"].values / A["tau_accmask_ACC"].mean())
for v, c in (("EKE", "C3"), ("T", "C0"), ("SPD", "C2")):
    ly = np.log(A[f"ACC_{v}"].values / A[f"ACC_{v}"].mean())
    ax[1].scatter(lt[:-1], ly[1:], color=c, s=14, label=f"{v} (year t+1)")
ax[1].plot([-.1, .1], [-.1, .1], "k--", lw=.7, label="slope 1")
ax[1].axhline(0, color="k", lw=.4); ax[1].axvline(0, color="k", lw=.4)
ax[1].set_xlabel("ln(τ/τ̄), ACC mask, year t"); ax[1].set_ylabel("ln(X/X̄)")
ax[1].set_title("Circumpolar eddy-saturation diagram (annual, raw)", fontsize=9); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(OUT / "s4_elasticity.png", dpi=130)
