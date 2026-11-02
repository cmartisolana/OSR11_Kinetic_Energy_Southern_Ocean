"""dk19: relation between intensification and displacement (Fig. 4e, f), Product 1 and Product 2.
(1) Monthly series: circumpolar half-power latitude anomaly over the 900 tracked longitudes (a2 baseline, km, + north) and
    latitude of the zonal-mean westerly maximum (dk17/s1 'lat_tau_zmax', km); trends (trend_bcz, 300 null series); lagged
    correlations (deseasonalised, detrended 12-month running means, lags -24..48 months, 1000 phase-randomised surrogates
    of the driver, global p over non-negative lags) of the half-power latitude with the ACC-band wind stress, the
    westerly-maximum latitude and the surface transport proxy.
(2) Across 5 deg longitude bins: trend of the jet speed along the path (b3) against the trend of the half-power latitude
    (retest 'lonbin|', 33 bins) and of the jet-path latitude (b3, 70 bins, Malvinas excluded); Pearson r, p from 2000
    circular shifts of one profile along longitude.
Usage: python dk19_speed_displacement.py <revision dir> <drake out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, DO = Path(sys.argv[1]), Path(sys.argv[2]); OUT = DO / "dk19"; OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
KM = 111.195; rng = np.random.default_rng(20261021)
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25); YR = np.repeat(np.arange(2001, 2026), 12)
TB = np.load(REV / "outputs" / "a2" / "tracks_baseline.npz"); keep = TB["keep"].astype(bool); a = TB["lat_hp"][:, keep]
hp = KM * np.nanmean(a - np.nanmean(a, 0, keepdims=True), 1)
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv")
W = pd.read_csv(DO / "dk17" / "s1_wind_series.csv")
wl = KM * (W.lat_tau_zmax.values - W.lat_tau_zmax.mean())
S = pd.DataFrame(dict(x=T, year=YR, month=MON, hp_km=hp, westerly_lat_km=wl, tau=A.tau_accmask_ACC.values, T=A.ACC_T.values))
S.to_csv(OUT / "dk19_series.csv", index=False)
rows = []
for c in ("hp_km", "westerly_lat_km"):
    r = R.trend_bcz(T, S[c].values, month=MON, n_outer=300)
    rows.append(dict(series=c, km_dec=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"]))
TR = pd.DataFrame(rows); TR.to_csv(OUT / "dk19_trends.csv", index=False); print(TR.round(3).to_string(index=False))


def prep(y):
    y = R.deseasonalise(np.asarray(y, float), MON); b, a_ = R.theil_sen(T, y); y = y - (a_ + b * T)
    return pd.Series(y).rolling(12, center=True).mean().values


def surr(x, n):
    X = np.fft.rfft(x - x.mean()); ph = rng.uniform(0, 2 * np.pi, (n, X.size)); ph[:, 0] = 0
    if x.size % 2 == 0: ph[:, -1] = 0
    return np.fft.irfft(np.abs(X)[None, :] * np.exp(1j * ph), n=x.size, axis=1) + x.mean()


def lagc(D, y, lags):
    out = np.full((D.shape[0], len(lags)), np.nan)
    for j, L in enumerate(lags):
        a_, b_ = (D[:, :D.shape[1] - L], y[L:]) if L >= 0 else (D[:, -L:], y[:L])
        a_ = a_ - a_.mean(1, keepdims=True); b_ = b_ - b_.mean(); out[:, j] = a_ @ b_ / np.sqrt((a_ ** 2).sum(1) * (b_ ** 2).sum())
    return out


lags = np.arange(-24, 49); y = prep(S.hp_km); curves = []; summ = []
for dn in ("tau", "westerly_lat_km", "T"):
    d = prep(S[dn]); ok = np.isfinite(d) & np.isfinite(y); dd, yy = d[ok], y[ok]
    r = lagc(dd[None, :], yy, lags)[0]; rs = lagc(surr(dd, 1000), yy, lags)
    pos = lags >= 0; j = int(np.nanargmax(np.where(pos, np.abs(r), -1)))
    pg = (1 + (np.abs(rs[:, pos]).max(1) >= abs(r[j])).sum()) / 1001
    summ.append(dict(driver=dn, r0=r[lags == 0][0], peak_r=r[j], peak_lag=lags[j], p_global=pg, r_crit=np.quantile(np.abs(rs[:, pos]).max(1), 0.95)))
    for L, rv, q in zip(lags, r, np.quantile(np.abs(rs), 0.95, 0)): curves.append(dict(driver=dn, lag=L, r=rv, r95=q))
pd.DataFrame(curves).to_csv(OUT / "dk19_lagcorr.csv", index=False)
SM = pd.DataFrame(summ); SM.to_csv(OUT / "dk19_lagcorr_summary.csv", index=False); print(SM.round(3).to_string(index=False))
# (2) across bins
sp = pd.read_csv(REV / "outputs" / "b3" / "b3_trends_speed.csv").set_index("lon")
la = pd.read_csv(REV / "outputs" / "b3" / "b3_trends_lat.csv").set_index("lon")
RT = pd.read_csv(REV / "outputs" / "retest" / "retest_summary.csv").set_index("key")
LB = RT[RT.index.str.startswith("lonbin|")].copy(); LB["lon"] = LB.index.str.split("|").str[1].astype(float) + 2.5
LB["hp_km_dec"] = LB.slope * KM * 10; LB["hp_lo"] = LB.ci_lo * KM * 10; LB["hp_hi"] = LB.ci_hi * KM * 10; LB = LB.set_index("lon")
cen = sp.index.values; excl = (cen > -60) & (cen < -50)


def circ(x, y, n=2000):
    ok = np.isfinite(x) & np.isfinite(y); r = np.corrcoef(x[ok], y[ok])[0, 1]; rr = []
    for _ in range(n):
        k = rng.integers(1, x.size); ys = np.roll(y, k); o2 = np.isfinite(x) & np.isfinite(ys)
        if o2.sum() > 5: rr.append(np.corrcoef(x[o2], ys[o2])[0, 1])
    rr = np.abs(np.array(rr)); return r, (1 + (rr >= abs(r)).sum()) / (1 + rr.size), int(ok.sum())


spd = sp.tr.reindex(cen).values.copy(); spd[excl] = np.nan
hpb = LB.hp_km_dec.reindex(cen).values
jl = la.tr.reindex(cen).values.copy(); jl[excl] = np.nan
B = pd.DataFrame(dict(lon=cen, speed_pct_dec=spd, speed_lo=sp.lo.reindex(cen).values, speed_hi=sp.hi.reindex(cen).values, speed_p=sp.p.reindex(cen).values,
                      hp_km_dec=hpb, hp_lo=LB.hp_lo.reindex(cen).values, hp_hi=LB.hp_hi.reindex(cen).values, hp_p=LB.p_new.reindex(cen).values,
                      path_km_dec=jl, path_p=la.p.reindex(cen).values))
B.to_csv(OUT / "dk19_bins.csv", index=False)
c1 = circ(spd, hpb); c2 = circ(spd, jl)
pd.DataFrame([dict(pair="speed vs half-power latitude", r=c1[0], p=c1[1], n=c1[2]), dict(pair="speed vs path latitude", r=c2[0], p=c2[1], n=c2[2])]).to_csv(OUT / "dk19_bin_corr.csv", index=False)
print("speed vs hp", np.round(c1, 3), "speed vs path", np.round(c2, 3))
