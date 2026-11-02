"""1. Eastward wind-stress series (CMEMS WIND L4 MY monthly, ERA5 bias-corrected with scatterometers).
- zonal-band mean 45-60S (cos-lat weights, ocean points valid in all 300 months), circumpolar and sectors
- mean over the ACC-band mask (1/8 deg mask block-averaged 2x2 to the 1/4 deg wind grid, fraction >= 0.5)
- zonal-mean stress maximum and its latitude (westerly-jet proxy), zonal-mean curl profile
  (circumpolar zonal-mean curl = -(1/(a cos phi)) d(taux_bar cos phi)/d phi, exact without tau_y up to land gaps)
Trends: rstats.trend_bcz on monthly series (deseasonalised), n_outer=300.
"""
import numpy as np, pandas as pd, xarray as xr
from common import REV, OUT, R, acc_mask_eighth, sector_of, time_frame, pct_dec

W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress
lat = W.latitude.values.astype(float); lon = W.longitude.values.astype(float)
T = W.values.astype(np.float64)                         # (300, 160, 1440)
valid = np.isfinite(T).all(axis=0)
coslat = np.cos(np.deg2rad(lat))[:, None] * np.ones(lon.size)[None, :]
sec = np.broadcast_to(sector_of(lon)[None, :], valid.shape)

# ACC mask -> 1/4 deg
m8, la8, lo8, _ = acc_mask_eighth()
m4 = m8.reshape(120, 2, 1440, 2).mean(axis=(1, 3)) >= 0.5
la4 = la8.reshape(120, 2).mean(1); lo4 = lo8.reshape(1440, 2).mean(1)
assert np.allclose(lo4, lon)
iy0 = int(np.argmin(np.abs(lat - la4[0])))
assert np.allclose(lat[iy0:iy0 + 120], la4)
accw = np.zeros(valid.shape, bool); accw[iy0:iy0 + 120] = m4
accw &= valid

band = (lat >= -60) & (lat <= -45)
S = time_frame()
for s in ("ACC", "Atlantic", "Indian", "Pacific"):
    selsec = np.ones_like(valid) if s == "ACC" else (sec == s)
    for tag, msk in (("z4560", valid & band[:, None] & selsec), ("accmask", accw & selsec)):
        w = coslat * msk
        S[f"tau_{tag}_{s}"] = np.nansum(np.where(msk, T, 0) * w, axis=(1, 2)) / w.sum()

# zonal-mean profile (ocean points valid all months), max and latitude of max, curl
zw = coslat * valid
tz = np.nansum(np.where(valid, T, 0) * zw, axis=2) / zw.sum(1)        # (300, nlat)
sl = (lat >= -65) & (lat <= -35) & (zw.sum(1) > 0) & np.isfinite(tz).all(0)
sm = (lat >= -62) & (lat <= -38)
ii = np.nanargmax(np.where(sm[sl][None, :], tz[:, sl], -np.inf), axis=1); latsl = lat[sl]; tzs = tz[:, sl]
# parabolic refinement
k = np.clip(ii, 1, sl.sum() - 2); r = np.arange(300)
y0, y1, y2 = tzs[r, k - 1], tzs[r, k], tzs[r, k + 1]
d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
S["tau_zmax"] = y1 - 0.25 * (y0 - y2) * d
S["lat_tau_zmax"] = latsl[k] + d * 0.25
a = 6.371e6; phi = np.deg2rad(lat)
curl = -np.gradient(tz * np.cos(phi)[None, :], phi, axis=1) / (a * np.cos(phi)[None, :])   # N m-3
S.to_csv(OUT / "s1_wind_series.csv", index=False)
pd.DataFrame(tz, columns=[f"{v:.3f}" for v in lat]).assign(x=S.x).to_csv(OUT / "s1_wind_zonalmean_profile.csv", index=False)

# trends
rows = []
cols = [c for c in S.columns if c.startswith("tau_") or c.startswith("lat_")]
for c in cols:
    o = R.trend_bcz(S.x.values, S[c].values, month=S.month.values, n_outer=300)
    mean = S[c].mean()
    ann = S.groupby("year")[c].mean()
    b_ann, _ = R.theil_sen(ann.index.values + 0.5, ann.values)
    rows.append(dict(series=c, mean=mean, slope_per_yr=o["slope"], ci_lo=o["ci_lo"], ci_hi=o["ci_hi"], p=o["p"],
                     pct_dec=np.nan if c.startswith("lat_") else pct_dec(o["slope"], mean),
                     pct_dec_lo=np.nan if c.startswith("lat_") else pct_dec(o["ci_lo"], mean),
                     pct_dec_hi=np.nan if c.startswith("lat_") else pct_dec(o["ci_hi"], mean),
                     annual_TS_slope=b_ann))
    print(rows[-1])
TR = pd.DataFrame(rows); TR.to_csv(OUT / "s1_wind_trends.csv", index=False)
S.groupby("year")[cols].mean().to_csv(OUT / "s1_wind_annual.csv")

# trend profiles of zonal-mean tau and curl (annual means, Theil-Sen + bcz annual)
ann_tz = pd.DataFrame(tz).groupby(S.year.values).mean().values
ann_curl = pd.DataFrame(curl).groupby(S.year.values).mean().values
xa = np.arange(2001, 2026) + 0.5
prof = []
for j in np.where(sl)[0]:
    bt, _ = R.theil_sen(xa, ann_tz[:, j]); bc, _ = R.theil_sen(xa, ann_curl[:, j])
    prof.append(dict(lat=lat[j], tau_mean=ann_tz[:, j].mean(), tau_trend_dec=10 * bt,
                     curl_mean=ann_curl[:, j].mean(), curl_trend_dec=10 * bc))
P = pd.DataFrame(prof); P.to_csv(OUT / "s1_wind_lat_profile.csv", index=False)

# figure
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
A = S.groupby("year").mean()
for s, c in zip(("ACC", "Atlantic", "Indian", "Pacific"), ("k", "C0", "C1", "C2")):
    ax[0, 0].plot(A.index, A[f"tau_z4560_{s}"], color=c, label=f"{s} 45-60S")
    ax[0, 0].plot(A.index, A[f"tau_accmask_{s}"], color=c, ls="--")
ax[0, 0].set_ylabel(r"$\tau_x$ (N m$^{-2}$), annual"); ax[0, 0].legend(fontsize=7, title="solid 45-60S, dashed ACC mask", title_fontsize=7)
ax[0, 0].set_title("Eastward wind stress")
ax[0, 1].plot(A.index, A.lat_tau_zmax, "k"); ax[0, 1].set_ylabel("lat of zonal-mean max (deg)")
ax2 = ax[0, 1].twinx(); ax2.plot(A.index, A.tau_zmax, "C3"); ax2.set_ylabel("max zonal-mean τx", color="C3")
ax[0, 1].set_title("Westerly maximum: position (black), strength (red)")
ax[1, 0].plot(P.lat, P.tau_mean, "k", label="mean"); ax[1, 0].set_xlabel("latitude"); ax[1, 0].set_ylabel("τx (N m⁻²)")
ax3 = ax[1, 0].twinx(); ax3.plot(P.lat, P.tau_trend_dec, "C3"); ax3.axhline(0, color="C3", lw=.5); ax3.set_ylabel("trend (N m⁻² dec⁻¹)", color="C3")
ax[1, 0].set_title("Zonal-mean τx: mean and Theil-Sen trend")
ax[1, 1].plot(P.lat, P.curl_mean * 1e7, "k"); ax[1, 1].set_ylabel("curl (1e-7 N m⁻³)"); ax[1, 1].set_xlabel("latitude")
ax4 = ax[1, 1].twinx(); ax4.plot(P.lat, P.curl_trend_dec * 1e7, "C3"); ax4.axhline(0, color="C3", lw=.5); ax4.set_ylabel("trend (1e-7 N m⁻³ dec⁻¹)", color="C3")
ax[1, 1].set_title("Zonal-mean curl: mean and trend")
fig.tight_layout(); fig.savefig(OUT / "s1_wind.png", dpi=130)
