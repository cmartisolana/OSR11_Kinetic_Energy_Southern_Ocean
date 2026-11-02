"""E5: nonlinearity ratio U/c of tracked eddies (AVISO META4.0 DT), ACC band, 2001-01 .. 2023-08.
U = rotational speed (speed_average: mean speed along the contour of maximum speed), c = translation speed from the
displacement between observations of the same track LAG days apart (great-circle distance / LAG days; LAG = 7 by
default, which filters the day-to-day jitter of the detected centre; LAG = 1 gives c about twice as large); the ratio
is assigned to the first of the two observations. Monthly ACC-band series: median U/c, fraction of observations with
U/c > 1 (nonlinear, Chelton et al. 2011), mean c and mean U; trends with the revised test (rstats.trend_bcz).
Usage: python eddies_nonlinearity.py <subset_dir> <repo outputs dir> <out_dir> [LAG days]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, netCDF4 as nc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

sub, repo, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
T0, T1 = pd.Timestamp("2001-01-01"), pd.Timestamp("2023-09-01")
RE = 6371.0e3
LAG = int(sys.argv[4]) if len(sys.argv) > 4 else 7


def load(pol):
    ds = nc.Dataset(sub / f"eddies_so_{pol}_all.nc")
    d = {v: np.asarray(ds.variables[v][:]) for v in ("track", "time", "latitude", "longitude", "speed_average")}
    df = pd.DataFrame(d); df = df[df.track > 0]
    df["time"] = pd.Timestamp("1950-01-01") + pd.to_timedelta(df.time, unit="D")
    df = df[(df.time >= T0) & (df.time < T1)].copy(); df["pol"] = pol
    return df


E = pd.concat([load("anticyclonic"), load("cyclonic")], ignore_index=True)
E = E.sort_values(["pol", "track", "time"]).reset_index(drop=True)
print("tracked observations", len(E), flush=True)
# displacement to the next observation of the same track, 1 day apart
same = (E.pol.values[:-LAG] == E.pol.values[LAG:]) & (E.track.values[:-LAG] == E.track.values[LAG:])
dt = (E.time.values[LAG:] - E.time.values[:-LAG]) / np.timedelta64(1, "D")
ok = same & (dt == LAG)
la0, la1 = np.radians(E.latitude.values[:-LAG]), np.radians(E.latitude.values[LAG:])
dlo = np.radians(E.longitude.values[LAG:] - E.longitude.values[:-LAG])
a = np.sin((la1 - la0) / 2) ** 2 + np.cos(la0) * np.cos(la1) * np.sin(dlo / 2) ** 2
dist = 2 * RE * np.arcsin(np.sqrt(a))                      # m per day
c = np.full(len(E), np.nan); c[:-LAG][ok] = dist[ok] / (86400.0 * LAG)
E["c"] = c; E["U"] = E.speed_average; E["ratio"] = E.U / E.c
E = E[np.isfinite(E.ratio) & (E.c > 0)].copy()
print("pairs", len(E), flush=True)

adt = xr.open_dataset(repo / "mean_adt.nc").mean_adt
ice = xr.open_dataset(repo / "sea_ice_mask.nc").ice_mask.values > 0.5
acc = np.isfinite(adt.values) & (adt.values >= -0.6) & (adt.values <= 0.2) & ~ice
lat, lon = adt.latitude.values, adt.longitude.values
lo = ((E.longitude.values + 180) % 360) - 180
i = np.clip(np.searchsorted(lat, E.latitude.values), 0, lat.size - 1); j = np.clip(np.searchsorted(lon, lo), 0, lon.size - 1)
A = E[acc[i, j]].copy(); A["ym"] = A.time.dt.year * 12 + A.time.dt.month - 1
print("ACC-band pairs", len(A), "overall median U/c", round(float(A.ratio.median()), 2),
      "fraction U/c>1", round(float((A.ratio > 1).mean()), 3), flush=True)

ym = np.arange(2001 * 12, 2023 * 12 + 8)
months = pd.DataFrame({"year": ym // 12, "month": ym % 12 + 1}); months["x"] = months.year + (months.month - 0.5) / 12
g = A.groupby("ym")
S = pd.concat([months, pd.DataFrame({"ratio_median": g.ratio.median(), "frac_nonlinear": g.ratio.apply(lambda r: np.mean(r > 1)),
                                     "c_mean": g.c.mean(), "U_mean": g.U.mean(), "n": g.size()}).reindex(ym).reset_index(drop=True)], axis=1)
S.to_csv(out / f"e5_nonlinearity_series_lag{LAG}.csv", index=False)
rows = []
for col in ("ratio_median", "frac_nonlinear", "c_mean", "U_mean"):
    y = S[col].values; o = R.trend_bcz(S.x.values, y, month=S.month.values, n_outer=300); m = float(np.nanmean(y))
    rows.append(dict(series=col, mean=m, slope=o["slope"], ci_lo=o["ci_lo"], ci_hi=o["ci_hi"], pct_per_decade=1000 * o["slope"] / m,
                     pct_lo=1000 * o["ci_lo"] / m, pct_hi=1000 * o["ci_hi"] / m, p=o["p"]))
    print(col, "mean", round(m, 4), "trend", round(rows[-1]["pct_per_decade"], 2), "%/dec [", round(rows[-1]["pct_lo"], 2), round(rows[-1]["pct_hi"], 2), "] p", round(o["p"], 3), flush=True)
pd.DataFrame(rows).to_csv(out / f"e5_nonlinearity_trends_lag{LAG}.csv", index=False)
print("done")
