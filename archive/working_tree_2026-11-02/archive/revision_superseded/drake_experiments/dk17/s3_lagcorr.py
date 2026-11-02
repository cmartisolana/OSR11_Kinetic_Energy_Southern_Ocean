"""2. Lagged response (Meredith & Hogg 2006; Hogg et al. 2015 style).
Responses: EKE, KE, along-stream speed (a5, all-sat; circumpolar also two-sat) and surface transport proxy T (a8),
circumpolar and per sector. Drivers: tau_x 45-60S (sector-matched), tau_x over the ACC-band mask (sector-matched),
circumpolar tau_x 45-60S, Marshall SAM.
Preparation: deseasonalise (rstats.deseasonalise); detrended (Theil-Sen) or not; then either 12-month centred running
mean of the monthly series (lags -24..+48 months) or annual means (lags -2..+4 yr). Positive lag = ocean lags driver.
Significance: 1000 phase-randomised surrogates of the prepared driver (spectrum preserved); per-lag p and a global p
for the maximum |r| over all lags (accounts for searching over lags).
"""
import numpy as np, pandas as pd
from common import REV, OUT, R, load_sam, time_frame

rng = np.random.default_rng(20261020); NSUR = 1000
a5 = pd.read_csv(REV / "outputs" / "a5" / "a5_series_allsat.csv")
a52 = pd.read_csv(REV / "outputs" / "a5" / "a5_series_twosat.csv")
a8 = pd.read_csv(REV / "outputs" / "a8" / "a8_series_allsat.csv")
W = pd.read_csv(OUT / "s1_wind_series.csv")
D = time_frame().merge(load_sam(), on=["year", "month", "x"])
for s in ("ACC", "Atlantic", "Indian", "Pacific"):
    for v in ("EKE", "KE", "SPD"):
        D[f"{s}_{v}"] = a5[f"{s}_{v}"].values
    D[f"{s}_T"] = a8[f"{s}_T"].values
for v in ("EKE", "KE", "SPD"):
    D[f"ACC2s_{v}"] = a52[f"ACC_{v}"].values
for c in W.columns:
    if c.startswith("tau") or c.startswith("lat_"):
        D[c] = W[c].values
assert len(D) == 300 and D.SAM.notna().all()
D.to_csv(OUT / "s3_all_series.csv", index=False)


def prep(col, detrend, scale):
    y = R.deseasonalise(D[col].values.astype(float), D.month.values); x = D.x.values
    if scale == "annual":
        y = pd.Series(y).groupby(D.year.values).mean().values; x = np.arange(2001, 2026) + 0.5
    if detrend:
        b, a = R.theil_sen(x, y); y = y - (a + b * x)
    if scale == "12mo":
        y = pd.Series(y).rolling(12, center=True).mean().values
    return y


def phase_surr(x, n):
    X = np.fft.rfft(x - x.mean()); m = len(x)
    ph = rng.uniform(0, 2 * np.pi, (n, X.size)); ph[:, 0] = 0
    if m % 2 == 0: ph[:, -1] = 0
    return np.fft.irfft(np.abs(X)[None, :] * np.exp(1j * ph), n=m, axis=1) + x.mean()


def lagcorr_many(Dm, y, lags):
    """Dm (k, n) drivers, y (n): r[k, L] = corr(D[t], y[t+L])."""
    out = np.full((Dm.shape[0], len(lags)), np.nan)
    for j, L in enumerate(lags):
        if L >= 0:
            a, b = Dm[:, :Dm.shape[1] - L], y[L:]
        else:
            a, b = Dm[:, -L:], y[:L]
        a = a - a.mean(1, keepdims=True); b = b - b.mean()
        out[:, j] = (a @ b) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum())
    return out


responses = [f"{s}_{v}" for s in ("ACC", "Atlantic", "Indian", "Pacific") for v in ("EKE", "KE", "SPD", "T")] + \
            [f"ACC2s_{v}" for v in ("EKE", "KE", "SPD")]
rows, curves = [], []
for scale, lags in (("12mo", np.arange(-24, 49)), ("annual", np.arange(-2, 5))):
    for detrend in (True, False):
        cache = {}
        for resp in responses:
            sec = resp.split("_")[0].replace("ACC2s", "ACC")
            drivers = {"tau_sector": f"tau_z4560_{sec}", "tau_accmask": f"tau_accmask_{sec}",
                       "tau_circ": "tau_z4560_ACC", "SAM": "SAM"}
            if sec == "ACC":
                drivers.pop("tau_circ")
            y = prep(resp, detrend, scale)
            for dname, dcol in drivers.items():
                key = dcol
                if key not in cache:
                    d = prep(dcol, detrend, scale); ok = np.isfinite(d)
                    cache[key] = (d, ok, phase_surr(d[ok], NSUR))
                d, ok, Sg = cache[key]
                yy = y[ok]; dd = d[ok]
                r = lagcorr_many(dd[None, :], yy, lags)[0]
                rs = lagcorr_many(Sg, yy, lags)
                p_lag = (1 + (np.abs(rs) >= np.abs(r)[None, :]).sum(0)) / (1 + NSUR)
                gmax = np.abs(rs).max(1)
                j = int(np.nanargmax(np.abs(r)))
                jp = int(np.nanargmax(np.where(lags >= 0, np.abs(r), -1)))   # best non-negative (ocean lags) lag
                p_glob = (1 + (gmax >= np.abs(r[j])).sum()) / (1 + NSUR)
                p_glob_pos = (1 + (np.abs(rs[:, lags >= 0]).max(1) >= np.abs(r[jp])).sum()) / (1 + NSUR)
                unit = 1 if scale == "12mo" else 12
                rows.append(dict(scale=scale, detrended=detrend, response=resp, driver=dname,
                                 r0=r[lags == 0][0], p0=p_lag[lags == 0][0],
                                 peak_r=r[j], peak_lag_months=lags[j] * unit, p_peak_local=p_lag[j], p_peak_global=p_glob,
                                 peakpos_r=r[jp], peakpos_lag_months=lags[jp] * unit, p_peakpos_local=p_lag[jp],
                                 p_peakpos_global=p_glob_pos, r_crit95_global=np.quantile(gmax, 0.95)))
                for L, rv, pv, q in zip(lags, r, p_lag, np.quantile(np.abs(rs), 0.95, axis=0)):
                    curves.append(dict(scale=scale, detrended=detrend, response=resp, driver=dname,
                                       lag_months=L * unit, r=rv, p=pv, r95_pointwise=q))
T = pd.DataFrame(rows); T.to_csv(OUT / "s3_lagcorr_summary.csv", index=False)
C = pd.DataFrame(curves); C.to_csv(OUT / "s3_lagcorr_curves.csv", index=False)
pd.set_option("display.width", 250)
print(T[(T.scale == "12mo") & T.detrended].round(3).to_string())
print(T[(T.scale == "annual") & T.detrended].round(3).to_string())

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
for scale in ("12mo", "annual"):
    fig, ax = plt.subplots(4, 4, figsize=(15, 11), sharex=True, sharey=True)
    for i, s in enumerate(("ACC", "Atlantic", "Indian", "Pacific")):
        for j, v in enumerate(("EKE", "KE", "SPD", "T")):
            a = ax[i, j]
            for dname, c in (("tau_sector", "C0"), ("tau_accmask", "C2"), ("SAM", "C3")):
                for det, ls in ((True, "-"), (False, ":")):
                    q = C[(C.scale == scale) & (C.detrended == det) & (C.response == f"{s}_{v}") & (C.driver == dname)]
                    a.plot(q.lag_months, q.r, color=c, ls=ls, lw=1.2 if det else .9,
                           label=f"{dname}{'' if det else ' (not detr.)'}")
                    if det and dname == "tau_sector":
                        a.fill_between(q.lag_months, -q.r95_pointwise, q.r95_pointwise, color="0.85", zorder=0)
            g = T[(T.scale == scale) & T.detrended & (T.response == f"{s}_{v}") & (T.driver == "tau_sector")].r_crit95_global.iloc[0]
            a.axhline(g, color="k", lw=.5, ls="--"); a.axhline(-g, color="k", lw=.5, ls="--")
            a.axhline(0, color="k", lw=.4); a.axvline(0, color="k", lw=.4)
            a.set_title(f"{s} {v}", fontsize=9)
    ax[0, 0].legend(fontsize=6)
    for a in ax[-1]: a.set_xlabel("lag (months; + = ocean lags driver)")
    for a in ax[:, 0]: a.set_ylabel("r")
    fig.suptitle(f"Lag correlation, {'12-month running mean' if scale=='12mo' else 'annual means'}; grey = pointwise 95% "
                 f"(phase-randomised τ, detrended); dashed = 95% for max|r| over lags")
    fig.tight_layout(); fig.savefig(OUT / f"s3_lagcorr_{scale}.png", dpi=110)
