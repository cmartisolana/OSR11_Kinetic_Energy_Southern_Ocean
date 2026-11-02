"""A8: surface geostrophic transport proxy across the ACC band, by longitude, sector and circumpolar (Referee A).

Specification fixed on 2026-10-14 before any trend was computed.
  T(lon, t) = integral over the ACC band of the monthly-mean eastward geostrophic velocity u_mean dy  [m2 s-1]
  Band: fixed set of grid points per longitude, mean ADT in [-0.6, 0.2] m (2001-2025 climatology), ice-free in every
  month (same mask as the KE/EKE band). A fixed set is used for the same reason as in the Drake proxy: a band whose
  extent varies produces a spurious trend. No depth integration: a surface transport proxy per unit depth.
  Also T_abs = integral of |u_mean| (sign-independent) as a check on westward recirculations inside the band.
  Series: circumpolar mean of T over longitude; sector means (Atlantic 290-20E, Indian 20-147E, Pacific 147-290E,
  as a5_energy.py); trends 2001-2025 with the revised test (rstats.trend_bcz), in % of the mean per decade.
  Longitude-resolved trends in 5-degree bins for the Supplement.
Usage: python a8_surface_transport.py <revision/data dir> <repo outputs dir> <out_dir> [twosat]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import rstats as R

data, repo, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
prod = sys.argv[4] if len(sys.argv) > 4 else "allsat"
files = [data / f"monthly_{prod}" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)}); t["x"] = t.year + (t.month - 0.5) / 12
with xr.open_dataset(files[0]) as d0:
    lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
mean_adt = xr.open_dataset(repo / "mean_adt.nc").mean_adt
ice = xr.open_dataset(repo / "sea_ice_mask.nc").ice_mask
if prod != "allsat":                                     # two-sat is on a 1/4 deg grid: interpolate the masks
    mean_adt = mean_adt.interp(latitude=lat, longitude=lon); ice = ice.interp(latitude=lat, longitude=lon)
band = np.isfinite(mean_adt.values) & (mean_adt.values >= -0.6) & (mean_adt.values <= 0.2) & ~(ice.values > 0.5)
dy = np.deg2rad(np.abs(np.diff(lat)).mean()) * 6371e3
T = np.full((len(files), lon.size), np.nan); TA = np.full_like(T, np.nan)
for k, f in enumerate(files):
    with xr.open_dataset(f) as ds:
        u = ds.u_mean.values.astype(np.float64)
    u = np.where(band, u, np.nan)
    T[k] = np.nansum(u, 0) * dy; TA[k] = np.nansum(np.abs(u), 0) * dy
    T[k, band.sum(0) == 0] = np.nan; TA[k, band.sum(0) == 0] = np.nan
lon360 = lon % 360
sec = np.where((lon360 >= 20) & (lon360 < 147), "Indian", np.where((lon360 >= 147) & (lon360 < 290), "Pacific", "Atlantic"))
S = t.copy()
S["ACC_T"] = np.nanmean(T, 1); S["ACC_Tabs"] = np.nanmean(TA, 1)
for s in ("Atlantic", "Indian", "Pacific"):
    S[f"{s}_T"] = np.nanmean(T[:, sec == s], 1); S[f"{s}_Tabs"] = np.nanmean(TA[:, sec == s], 1)
S.to_csv(out / f"a8_series_{prod}.csv", index=False)
rows = []
for c in [c for c in S.columns if c.endswith("_T") or c.endswith("_Tabs")]:
    y = S[c].values; o = R.trend_bcz(S.x.values, y, month=S.month.values, n_outer=300); m = float(np.nanmean(y))
    rows.append(dict(series=c, product=prod, mean_m2s=m, slope=o["slope"], ci_lo=o["ci_lo"], ci_hi=o["ci_hi"],
                     pct_per_decade=1000 * o["slope"] / m, pct_lo=1000 * o["ci_lo"] / m, pct_hi=1000 * o["ci_hi"] / m, p=o["p"]))
    print(f"{prod} {c:14s} mean {m:9.0f} m2/s  trend {rows[-1]['pct_per_decade']:+6.2f} %/dec [{rows[-1]['pct_lo']:+.2f}, {rows[-1]['pct_hi']:+.2f}]  p={o['p']:.3f}", flush=True)
pd.DataFrame(rows).to_csv(out / f"a8_trends_{prod}.csv", index=False)
# 5-degree bins
bins = np.arange(-180, 181, 5); lb = []
for b0, b1 in zip(bins[:-1], bins[1:]):
    sel = (lon >= b0) & (lon < b1)
    y = np.nanmean(T[:, sel], 1)
    if np.isfinite(y).mean() < 0.9:
        continue
    o = R.trend_bcz(S.x.values, y, month=S.month.values, n_outer=300); m = float(np.nanmean(y))
    lb.append(dict(lon=(b0 + b1) / 2, mean_m2s=m, pct_per_decade=1000 * o["slope"] / m, pct_lo=1000 * o["ci_lo"] / m, pct_hi=1000 * o["ci_hi"] / m, p=o["p"]))
pd.DataFrame(lb).to_csv(out / f"a8_lon_trends_{prod}.csv", index=False)
print("done")
