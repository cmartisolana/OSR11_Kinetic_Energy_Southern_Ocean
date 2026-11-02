"""dk29: does the 2001-2025 EKE record follow the integrated history of the westerly forcing? (memory test; eddy saturation)
The wind-stress record used elsewhere starts in 2001, so the forcing history is represented by the SAM (Marshall 2003
station index, monthly from 1957: newsam.1957.2007 for 1957-2000, the current BAS table for 2001-2025), which over
2001-2025 correlates with the ACC-band eastward wind stress (dk17/s3; printed below).
Forcing memory: X_T(y) = sum_{k>=0} h_T(k) SAM(y-k), annual means, with an exponential kernel h_T(k) ~ exp(-k/T) normalised to
unit sum (T = 0: same year), for T = 0..30 yr; and the two-timescale step response of Hogg et al. (2015), fast ~2 yr and slow
~15 yr with the fast amplitude twice the slow (impulse response h = (2/3) e^{-k/2}/2 + (1/3) e^{-k/15}/15, normalised).
EKE: ACC-band and sector EKE (% of the 2001-2025 mean), annual, NOT detrended (the trend is what is to be explained).
Fit EKE = a + b X over 2001-2025 (least squares): r, fraction of the observed EKE trend reproduced (trend of the fit /
observed trend), and the implied response. Null: 2000 phase-randomised surrogates of the 1957-2025 monthly SAM (same
spectrum), passed through the same kernels; 95th percentile of r for each T and of max_T r (multiplicity over T).
Usage: python dk29_forcing_memory.py <drake out dir> <sam 1957-2000 txt> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
DO, SAMOLD, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(29)
old = pd.read_csv(SAMOLD, sep=r"\s+", comment="#", header=None); old.columns = ["year"] + list(range(1, 13))
old = old.melt(id_vars="year", var_name="month", value_name="SAM")
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv")
new = A[["year", "month", "SAM"]]
S = pd.concat([old, new]).sort_values(["year", "month"]).reset_index(drop=True)
assert len(S) == 69 * 12 and S.SAM.notna().all()
an_ = A.groupby("year").mean(numeric_only=True)
print("ACC-band wind stress vs SAM, annual 2001-2025: r =", np.corrcoef(an_.tau_accmask_ACC, an_.SAM)[0, 1].round(2),
      "| detrended:", np.corrcoef(*(v - np.polyval(np.polyfit(np.arange(25), v, 1), np.arange(25)) for v in (an_.tau_accmask_ACC.values, an_.SAM.values)))[0, 1].round(2))
yrs_all = np.arange(1957, 2026); i0 = np.where(yrs_all == 2001)[0][0]
def annual(m): return m.reshape(-1, 12).mean(1)
def kern(T, n=69):
    k = np.arange(n)
    if T == "Hogg": h = (2 / 3) * np.exp(-k / 2.0) / 2.0 + (1 / 3) * np.exp(-k / 15.0) / 15.0
    elif T == 0: h = (k == 0).astype(float)
    else: h = np.exp(-k / T)
    return h / h.sum()
def filt(sa, h):   # causal, annual series 1957-2025 -> 2001-2025
    return np.array([np.dot(h[:i + 1], sa[i::-1]) for i in range(i0, len(sa))])
TS = [0, 0.5, 1, 2, 3, 4, 5, 7, 10, 12, 15, 20, 25, 30]
yrs = np.arange(2001, 2026); t = yrs - yrs.mean()
SER = {"ACC band": "ACC_EKE", "Pacific": "Pacific_EKE", "Indian": "Indian_EKE", "Atlantic": "Atlantic_EKE"}
sam_a = annual(S.SAM.values)
# surrogates of the monthly SAM
def phase_rand(x):
    X = np.fft.rfft(x - x.mean()); ph = np.exp(1j * rng.uniform(0, 2 * np.pi, X.size)); ph[0] = 1
    if x.size % 2 == 0: ph[-1] = np.sign(np.cos(np.angle(ph[-1])))
    return np.fft.irfft(X * ph, n=x.size) + x.mean()
NS = 2000; SUR = np.array([annual(phase_rand(S.SAM.values)) for _ in range(NS)])
rows, fits = [], {}
for reg, col in SER.items():
    e = A.groupby("year")[col].mean().values; ep = 100 * (e / e.mean() - 1)
    obs_tr = 10 * np.polyfit(t, ep, 1)[0]
    rmax_null = np.full(NS, -np.inf)
    for T in TS + ["Hogg"]:
        h = kern(T); x = filt(sam_a, h)
        b, a = np.polyfit(x, ep, 1); fit = a + b * x; r = np.corrcoef(x, ep)[0, 1]
        fit_tr = 10 * np.polyfit(t, fit, 1)[0]
        xs = np.array([filt(s, h) for s in SUR]); xs_c = xs - xs.mean(1, keepdims=True); ec = ep - ep.mean()
        rs = (xs_c @ ec) / np.sqrt((xs_c ** 2).sum(1) * (ec ** 2).sum())
        if T != "Hogg": rmax_null = np.maximum(rmax_null, rs)
        rows.append(dict(region=reg, T=str(T), r=r, r95=np.quantile(rs, 0.95), p=np.mean(rs >= r), b=b, obs_trend=obs_tr, fit_trend=fit_tr,
                         frac=fit_tr / obs_tr, sam_filtered_trend=10 * np.polyfit(t, x, 1)[0]))
        fits[f"{reg}|{T}"] = fit
    D_ = pd.DataFrame([r_ for r_ in rows if r_["region"] == reg and r_["T"] not in ("Hogg",)])
    best = D_.loc[D_.r.idxmax()]
    for r_ in rows:
        if r_["region"] == reg: r_["r95_maxT"] = np.quantile(rmax_null, 0.95); r_["p_maxT_best"] = np.mean(rmax_null >= best.r)
    fits[f"{reg}|obs"] = ep
D = pd.DataFrame(rows); D.to_csv(OUT / "dk29_memory.csv", index=False)
pd.set_option("display.width", 200); print(D.round(3).to_string(index=False))
np.savez(OUT / "dk29_series.npz", yrs=yrs, yrs_all=yrs_all, sam_annual=sam_a,
         sam_T={str(T): filt(sam_a, kern(T)) for T in (0, 2, 5, 15, "Hogg")}, **{k.replace(" ", "_"): v for k, v in fits.items()},
         sam_T2=filt(sam_a, kern(2)), sam_T5=filt(sam_a, kern(5)), sam_T15=filt(sam_a, kern(15)), sam_Hogg=filt(sam_a, kern("Hogg")),
         sam_ma_full={str(T): np.array([np.dot(kern(T)[:i + 1], sam_a[i::-1]) for i in range(len(sam_a))]) for T in (2, 15)},
         sam_full_T2=np.array([np.dot(kern(2)[:i + 1], sam_a[i::-1]) for i in range(len(sam_a))]),
         sam_full_T15=np.array([np.dot(kern(15)[:i + 1], sam_a[i::-1]) for i in range(len(sam_a))]),
         sam_full_Hogg=np.array([np.dot(kern("Hogg")[:i + 1], sam_a[i::-1]) for i in range(len(sam_a))]))
