"""B5c: trends of the standing north-south flow inside the ACC band (mean ADT -0.6..+0.2 m, ice-free; Sect. 2.2).
Standing meridional flow s = v * sign(v_mean): the annual-mean meridional velocity projected on the sign of the 2001-2025
mean; s > 0 strengthens the standing northward/southward branches, s < 0 weakens them. Signed, so eddies average out
instead of being rectified as in |v|. Cells with |v_mean| < 2 cm/s (no standing branch) are excluded.
Map: OLS slope of the 25 annual values on 1 deg cells (cm/s per decade; as the Fig. 2a map). Profile: cos-lat weighted
band mean of s in 5 deg longitude bins, trend_bcz on annual values (AR(1), calibrated, 300 null series), FDR 10 %.
Usage: python b5c_meridional_trends.py <revision dir> <b5 out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, B5 = Path(sys.argv[1]), Path(sys.argv[2]); sys.path.insert(0, str(REV)); import rstats as R
with xr.open_dataset(REV / "data" / "monthly_allsat" / "2001" / "01.nc") as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
V = np.load(B5 / "b5_annual_fields.npz")["V"].astype(float) * 100                       # cm/s
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
acc = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice & np.all(np.isfinite(V), 0)
vm = V.mean(0); keep = acc & (np.abs(vm) >= 2.0)
Ssd = V * np.sign(vm)[None]
X = np.arange(2001, 2026) + 0.5; xa = X - X.mean()
# 1 deg map (8 x 8 cells of 1/8 deg)
ny, nx = (lat.size // 8) * 8, lon.size
blk = lambda a: a[..., :ny, :].reshape(a.shape[:-2] + (ny // 8, 8, nx // 8, 8))
w = np.cos(np.radians(lat))[:, None] * keep
num = np.nansum(blk(Ssd * w[None]), axis=(-3, -1)); den = blk(w).sum(axis=(-3, -1))
S1 = np.where(den[None] > 0.25 * 64 * np.cos(np.radians(55)), num / np.where(den > 0, den, np.nan)[None], np.nan)
slope = 10 * np.tensordot(xa, S1 - np.nanmean(S1, 0), axes=(0, 0)) / (xa ** 2).sum()
vm1 = np.nansum(blk(vm * w), axis=(-3, -1)) / np.where(den > 0, den, np.nan)
np.savez(B5 / "b5c_map.npz", lat=lat[:ny].reshape(-1, 8).mean(1), lon=lon.reshape(-1, 8).mean(1), slope=slope, s_mean=np.nanmean(S1, 0), vmean=vm1)
# 5 deg profile
EDGES = np.arange(-180, 181, 5.0); rows = []
for a in EDGES[:-1]:
    m = (lon >= a) & (lon < a + 5); ww = w[:, m]
    if ww.sum() == 0: rows.append(dict(lon=a + 2.5)); continue
    y = np.nansum(Ssd[:, :, m] * ww[None], axis=(1, 2)) / ww.sum()
    r = R.trend_bcz(X, y, monthly=False, n_outer=300); mu = y.mean()
    rows.append(dict(lon=a + 2.5, mean=mu, slope_dec=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"],
                     pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, n_pix=int((ww > 0).sum())))
D = pd.DataFrame(rows); ok = D.p.notna(); pv = D.p[ok].values; o = np.argsort(pv); q = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
q = np.minimum.accumulate(q[::-1])[::-1]; f = np.zeros(len(pv), bool); f[o] = q <= 0.10; D["fdr10"] = False; D.loc[ok, "fdr10"] = f
D.to_csv(B5 / "b5c_profile.csv", index=False)
print("bins", int(ok.sum()), "positive", int((D.slope_dec > 0).sum()), "p<0.05", int((D.p.round(3) < 0.05).sum()), "FDR", int(D.fdr10.sum()))
print(D[D.p.round(3) < 0.05][["lon", "mean", "pct", "p"]].round(3).to_string(index=False)); print("done")
