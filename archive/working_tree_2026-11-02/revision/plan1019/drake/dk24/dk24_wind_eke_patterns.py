"""dk24: EKE and wind stress in three forms of earlier studies (Product 1 EKE; Product 2 eastward wind stress; 1 deg fields of
dk22; 2001-2025).
(1) Trend patterns (Martinez-Moreno et al. 2021; Zhang et al. 2021): least-squares trend of annual means on 1 deg cells for
    the wind stress (as the EKE map of Fig. 3a); and, in twelve 30 deg sectors (Zhang et al. 2021), the trends of EKE and
    wind stress averaged over the ACC band (mean ADT -0.6..+0.2 m) with rstats.trend_bcz (% of the mean per decade).
(2) Composites (Meredith and Hogg 2006): detrended annual EKE anomalies (% of the cell mean) averaged over the 6 years of
    strongest minus the 6 years of weakest ACC-band wind stress (detrended), at lags 0, 1, 2 yr (EKE lagging); also for the
    SAM. Significance per cell: two-sample t test (Welch).
(3) Regional lag maps (Hogg et al. 2022): annual EKE averaged over the ACC band in four regions, correlated with the wind
    stress on 3 deg cells at lags 0 and 1 yr (both detrended), to see whether the EKE of a region follows local or remote wind.
Usage: python dk24_wind_eke_patterns.py <revision dir> <drake out dir> <dk22 out dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import stats
from multiprocessing import Pool
REV, DO, D22, OUT = (Path(a) for a in sys.argv[1:5]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
Z = np.load(D22 / "dk22_fields.npz"); EK, TAU, lat, lon = Z["EK"], Z["TAU"], Z["lat"], Z["lon"]
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt
m1 = madt.coarsen(latitude=8, longitude=8, boundary="trim").mean().values
acc = (m1 >= -0.6) & (m1 <= 0.2) & np.isfinite(EK).all(0)
w = np.cos(np.radians(lat))[:, None] * np.ones(lon.size)[None, :]
ann = lambda a: a.reshape(25, 12, *a.shape[1:]).mean(1)
yrs = np.arange(2001, 2026); t = yrs - yrs.mean()


def detr(a):
    a = np.asarray(a, float); tt = t.reshape((-1,) + (1,) * (a.ndim - 1)); m = a.mean(0)
    return a - m - tt * ((tt * (a - m)).sum(0) / (t ** 2).sum())


EA, TA = ann(EK), ann(TAU)
# (1a) wind trend map (per decade, N m-2 and %)
sl = (t[:, None, None] * (TA - TA.mean(0))).sum(0) / (t ** 2).sum() * 10
out = dict(lat=lat, lon=lon, tau_trend=sl, tau_trend_pct=100 * sl / TA.mean(0), tau_mean=TA.mean(0))
# (1b) sector trends
X = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
ed = np.arange(-180, 181, 30); jobs = []
for k in range(12):
    m = acc & ((lon >= ed[k]) & (lon < ed[k + 1]))[None, :]
    if m.sum() < 5: continue
    for v, F in (("EKE", EK), ("tau", TAU)):
        y = np.nansum(F[:, m] * w[m], 1) / w[m].sum(); jobs.append((ed[k] + 15, v, y))


def tr(j):
    c, v, y = j; r = R.trend_bcz(X, y, month=MON, n_outer=300); mu = y.mean()
    return dict(lon=c, var=v, pct=1000 * r["slope"] / mu, lo=1000 * r["ci_lo"] / mu, hi=1000 * r["ci_hi"] / mu, p=r["p"])


# (2) composites
A = pd.read_csv(DO / "dk17" / "s3_all_series.csv"); an = A.groupby("year").mean(numeric_only=True)
Ep = detr(100 * (EA / EA.mean(0) - 1))
for nm, idx in (("tau", detr(an.tau_accmask_ACC.values)), ("sam", detr(an.SAM.values))):
    o = np.argsort(idx); lo_, hi_ = o[:6], o[-6:]
    out[f"years_{nm}_strong"] = yrs[hi_]; out[f"years_{nm}_weak"] = yrs[lo_]
    for L in (0, 1, 2):
        hs, ls_ = hi_[hi_ + L < 25] + L, lo_[lo_ + L < 25] + L
        d = Ep[hs].mean(0) - Ep[ls_].mean(0); tt_, p = stats.ttest_ind(Ep[hs], Ep[ls_], axis=0, equal_var=False)
        out[f"comp_{nm}_lag{L}"] = d; out[f"comp_{nm}_lag{L}_p"] = p
        ok = acc
        print(f"composite {nm} lag {L}: ACC-band mean diff {np.nansum(d[ok] * w[ok]) / w[ok].sum():+.1f} %, frac p<0.05 {np.nanmean(p[ok] < 0.05):.3f}",
              "strong", yrs[hi_], "weak", yrs[lo_])
# (3) regional lag maps on 3 deg wind cells (Hogg et al. 2022 used 4 deg)
REG = {"SE Pacific (130-80W)": (-130, -80), "Kerguelen lee (70-100E)": (70, 100), "SE Indian-Macquarie (110-170E)": (110, 170), "Atlantic (50W-20E)": (-50, 20)}
T4 = TA.reshape(25, lat.size // 3, 3, lon.size // 3, 3).mean(axis=(2, 4)); T4d = detr(T4)
for nm, (a_, b_) in REG.items():
    m = acc & ((lon >= a_) & (lon < b_))[None, :]
    e = detr(np.nansum(EA[:, m] * w[m], 1) / w[m].sum())
    for L in (0, 1):
        x = T4d[:25 - L] if L else T4d; y = e[L:]
        xa = x - x.mean(0); ya = (y - y.mean())[:, None, None]
        r = (xa * ya).sum(0) / np.sqrt((xa ** 2).sum(0) * (ya ** 2).sum())
        out[f"reg_{nm}_lag{L}"] = r
        print(f"{nm} lag {L}: max r {np.nanmax(r):.2f} at", np.unravel_index(np.nanargmax(r), r.shape))
out["lat4"] = lat.reshape(-1, 3).mean(1); out["lon4"] = lon.reshape(-1, 3).mean(1); out["regnames"] = np.array(list(REG))
out["acc"] = acc
np.savez(OUT / "dk24_patterns.npz", **out)
if __name__ == "__main__":
    with Pool(2) as P: rows = P.map(tr, jobs)
    D = pd.DataFrame(rows); D.to_csv(OUT / "dk24_sector_trends.csv", index=False)
    W_ = D.pivot(index="lon", columns="var", values="pct"); print(W_.round(1).to_string())
    print("r(sector EKE trend, wind trend) =", np.corrcoef(W_.EKE, W_.tau)[0, 1].round(2))
