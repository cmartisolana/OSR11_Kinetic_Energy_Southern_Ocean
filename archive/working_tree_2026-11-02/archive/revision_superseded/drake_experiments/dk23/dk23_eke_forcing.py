"""dk23: EKE against the wind stress and the SAM, in the forms used by earlier studies (Meredith and Hogg 2006: lagged
response of circumpolar EKE to the SAM; Hogg et al. 2015 and Zhang et al. 2021: lagged correlation of sector EKE with
zonal wind stress) and in the forms of Fig. 3 (maps; profile along the ACC).
All series: annual means 2001-2025, linearly detrended. Correlations r(lag) = corr(forcing[t], EKE[t + lag]), lag in years
(positive: EKE lags). Significance: two-sided t test with effective sample size n_eff = n (1 - a1 b1)/(1 + a1 b1)
(lag-1 autocorrelations; Bretherton et al. 1999).
(1) Sector series (dk17/s3): EKE of ACC, Atlantic, Indian, Pacific band vs the wind stress over the same sector band
    ('tau_accmask_*') and vs the SAM (Marshall index), lags -2..+5 yr.
(2) Maps: 3 deg cells (dk22 fields): EKE vs the SAM at lags 0..3 yr (EKE vs local wind: dk22).
(3) Along the ACC: EKE between the mean SAF and PF in 6 deg longitude bins (from the 2 deg bins of dk12) vs the eastward
    wind stress averaged over the same bins between the fronts, and vs the SAM, lags 0..3 yr.
Usage: python dk23_eke_forcing.py <revision dir> <drake out dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
REV, DO, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
LAGS = range(-2, 6)


def detr(a):
    a = np.asarray(a, float); t = np.arange(a.shape[0]) - (a.shape[0] - 1) / 2
    tt = t.reshape((-1,) + (1,) * (a.ndim - 1)); m = a.mean(0)
    b = (tt * (a - m)).sum(0) / (t ** 2).sum()
    return a - m - tt * b


def lagcorr(x, y, lag):
    """x forcing, y response (time first); positive lag: y lags x."""
    if lag >= 0: xa, ya = x[:len(x) - lag], y[lag:]
    else: xa, ya = x[-lag:], y[:len(y) + lag]
    xa = xa - xa.mean(0); ya = ya - ya.mean(0); n = xa.shape[0]
    r = (xa * ya).sum(0) / np.sqrt((xa * xa).sum(0) * (ya * ya).sum(0))
    a1 = (xa[1:] * xa[:-1]).sum(0) / (xa * xa).sum(0); b1 = (ya[1:] * ya[:-1]).sum(0) / (ya * ya).sum(0)
    ne = np.clip(n * (1 - a1 * b1) / (1 + a1 * b1), 3, n)
    p = 2 * stats.t.sf(np.abs(r * np.sqrt((ne - 2) / np.maximum(1 - r ** 2, 1e-9))), ne - 2)
    rc = stats.t.ppf(0.975, ne - 2) / np.sqrt(ne - 2 + stats.t.ppf(0.975, ne - 2) ** 2)   # |r| for p = 0.05
    return r, p, ne, rc


A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); an = A.groupby("year").mean(numeric_only=True)
sam = detr(an.SAM.values)
rows = []
for sec in ("ACC", "Atlantic", "Indian", "Pacific"):
    e = detr(an[f"{sec}_EKE"].values); tau = detr(an[f"tau_accmask_{sec}"].values)
    for fn, x in (("wind stress", tau), ("SAM", sam)):
        for L in LAGS:
            r, p, ne, rc = lagcorr(x, e, L); rows.append(dict(sector=sec, forcing=fn, lag=L, r=r, p=p, n_eff=ne, r_crit=rc))
    e2 = detr(an.ACC2s_EKE.values) if sec == "ACC" else None
    if e2 is not None:
        for fn, x in (("wind stress", tau), ("SAM", sam)):
            for L in LAGS:
                r, p, ne, rc = lagcorr(x, e2, L); rows.append(dict(sector="ACC two-satellite", forcing=fn, lag=L, r=r, p=p, n_eff=ne, r_crit=rc))
S = pd.DataFrame(rows); S.to_csv(OUT / "dk23_sector_lagcorr.csv", index=False)
# trends of the forcing (context)
print(S.pivot_table(index=["sector", "forcing"], columns="lag", values="r").round(2).to_string())
# (2) SAM maps
Z = np.load(DO.parent.parent.parent / "plan1019" / "drake" / "out" / "dk22" / "dk22_fields.npz") if False else np.load(OUT.parent / "dk22" / "dk22_fields.npz")
EK = Z["EK"]; CG = 3
cg = lambda a: a.reshape(a.shape[0], a.shape[1] // CG, CG, a.shape[2] // CG, CG).mean(axis=(2, 4))
E3 = detr(cg(EK).reshape(25, 12, *cg(EK).shape[1:]).mean(1))
maps = {}
for L in range(0, 4):
    r, p, ne, rc = lagcorr(sam[:, None, None] * np.ones(E3.shape), E3, L); maps[f"r_sam_lag{L}"] = r; maps[f"p_sam_lag{L}"] = p
    print("SAM map lag", L, "median r", np.nanmedian(r).round(2), "frac p<0.05", np.round(np.nanmean(p[np.isfinite(p)] < 0.05), 3))
np.savez(OUT / "dk23_sam_maps_3deg.npz", lat=Z["lat"].reshape(-1, CG).mean(1), lon=Z["lon"].reshape(-1, CG).mean(1), **maps)
# (3) along the ACC
B = np.load(DO / "dk12" / "dk12_band_series_allsat.npz"); lon2 = B["lon"]; eke2 = B["eke"]
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt
W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress.sel(latitude=slice(-66, -34))
mw = madt.interp(latitude=W.latitude, longitude=W.longitude).values
band = (mw >= float(B["lev_pf"])) & (mw <= float(B["lev_saf"]))
wl = W.longitude.values; wa = np.cos(np.radians(W.latitude.values))[:, None] * band
TW = W.values  # (300, lat, lon)
ed = np.arange(-180, 181, 6); nb = ed.size - 1
tauB = np.full((300, nb), np.nan); ekeB = np.full((300, nb), np.nan)
for b in range(nb):
    m = (wl >= ed[b]) & (wl < ed[b + 1]); w_ = wa[:, m]
    if w_.sum() > 0: tauB[:, b] = np.nansum(TW[:, :, m] * w_, axis=(1, 2)) / w_.sum()
    m2 = (lon2 >= ed[b]) & (lon2 < ed[b + 1])
    if m2.any(): ekeB[:, b] = np.nanmean(eke2[:, m2], 1)
ann = lambda a: a.reshape(25, 12, -1).mean(1)
tA, eA = detr(ann(tauB)), detr(ann(ekeB))
prow = []
for L in range(0, 4):
    r, p, ne, rc = lagcorr(tA, eA, L); rs, ps, _, rcs = lagcorr(sam[:, None] * np.ones(eA.shape), eA, L)
    for b in range(nb):
        prow.append(dict(lon=0.5 * (ed[b] + ed[b + 1]), lag=L, r_tau=r[b], p_tau=p[b], r_sam=rs[b], p_sam=ps[b], rcrit=rc[b]))
P = pd.DataFrame(prow); P.to_csv(OUT / "dk23_alongACC.csv", index=False)
for L in range(4):
    q = P[P.lag == L]; print("along-ACC lag", L, "tau: median r", q.r_tau.median().round(2), "n p<0.05", (q.p_tau < 0.05).sum(), "/", q.r_tau.notna().sum(),
                              "| SAM: median r", q.r_sam.median().round(2), "n p<0.05", (q.p_sam < 0.05).sum())
