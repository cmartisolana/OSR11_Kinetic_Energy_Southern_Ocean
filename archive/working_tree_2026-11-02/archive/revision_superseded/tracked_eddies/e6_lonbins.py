"""E6: tracked-eddy distributions and trends by longitude in the ACC band (for a candidate Figure 4, 26 Oct 2026).

Atlas META4.0 DT all-sat subset, tracked observations (track > 0), ACC band (mean ADT -0.6..+0.2 m, ice-free),
2001-01 .. 2023-08 (last full month). Per-observation variables: amplitude (m), speed radius (km), rotational speed
(m/s); per-track: lifetime (d) and propagation distance (km) for tracks born in the band (observation_number = 0,
births <= 2022-08-31).
  (1) Distributions: histograms of the four observation/track variables in 2001-2011 and 2012-2023 (density),
      with medians.
  (2) Longitude bins of 10 deg (36): monthly series of eddies present per day, mean amplitude, mean speed radius,
      mean rotational speed and sum(amplitude^2) per day; trend_bcz on the deseasonalised monthly series
      (both noise models, n_outer = NOUT), % of the bin mean per decade; BH-FDR 10 % across bins per variable.
      ACC-band totals from E2 are repeated here for reference (same test).
Usage: python e6_lonbins.py <subset_dir> <repo outputs dir> <revision dir> <out_dir> [NOUT]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, netCDF4 as nc

sub, repo, rev, out = (Path(a) for a in sys.argv[1:5]); out.mkdir(parents=True, exist_ok=True)
NOUT = int(sys.argv[5]) if len(sys.argv) > 5 else 150
sys.path.insert(0, str(rev))
import rstats as R  # noqa: E402

T0, T1 = pd.Timestamp("2001-01-01"), pd.Timestamp("2023-09-01")
BIRTH_MAX = pd.Timestamp("2022-08-31")


def load(pol):
    ds = nc.Dataset(sub / f"eddies_so_{pol}_all.nc")
    cols = ("track", "observation_number", "time", "latitude", "longitude", "amplitude", "speed_radius", "speed_average")
    df = pd.DataFrame({v: np.asarray(ds.variables[v][:]) for v in cols})
    df["time"] = pd.Timestamp("1950-01-01") + pd.to_timedelta(df.time.astype(int), unit="D")
    df["lon"] = ((df.longitude + 180) % 360) - 180
    df = df[(df.track > 0) & (df.time >= T0) & (df.time < T1)].copy(); df["pol"] = pol
    return df


E = pd.concat([load("anticyclonic"), load("cyclonic")], ignore_index=True)
E["tid"] = E.pol.str[0] + E.track.astype(str)
adt = xr.open_dataset(repo / "mean_adt.nc").mean_adt
ice = xr.open_dataset(repo / "sea_ice_mask.nc").ice_mask.values > 0.5
acc = np.isfinite(adt.values) & (adt.values >= -0.6) & (adt.values <= 0.2) & ~ice
lat, lon = adt.latitude.values, adt.longitude.values
i = np.clip(np.searchsorted(lat, E.latitude.values), 0, lat.size - 1); j = np.clip(np.searchsorted(lon, E.lon.values), 0, lon.size - 1)
E["acc"] = acc[i, j]
A = E[E.acc].copy()
A["ym"] = A.time.dt.year * 12 + A.time.dt.month - 1
A["bin"] = np.floor((A.lon + 180) / 10).astype(int).clip(0, 35)
A["period"] = np.where(A.time.dt.year <= 2011, "early", "late")
print("band observations", len(A), flush=True)

# per-track (born in band)
E = E.sort_values(["tid", "time"])
g = E.groupby("tid", sort=False)
dlon = g.lon.diff().values; dlat = g.latitude.diff().values
dlon = np.where(dlon > 180, dlon - 360, np.where(dlon < -180, dlon + 360, dlon))
E["step_km"] = np.nan_to_num(np.hypot(dlat * 111.2, dlon * 111.2 * np.cos(np.radians(E.latitude.values))))
first = g.head(1).set_index("tid"); last = g.tail(1).set_index("tid")
T = pd.DataFrame({"born": first.observation_number.values == 0, "acc0": first.acc.values, "t0": first.time,
                  "lifetime_d": (last.time - first.time).dt.days + 1, "dist_km": E.groupby("tid", sort=False).step_km.sum(),
                  "lon0": first.lon})
T = T[T.born & T.acc0 & (T.t0 <= BIRTH_MAX)].copy()
T["period"] = np.where(T.t0.dt.year <= 2011, "early", "late")
T["bin"] = np.floor((T.lon0 + 180) / 10).astype(int).clip(0, 35)

# (1) distributions
edges = {"amplitude": np.arange(0, 0.405, 0.01), "speed_radius": np.arange(20e3, 150e3 + 1, 2.5e3),
         "speed_average": np.arange(0, 0.81, 0.02), "lifetime_d": np.arange(0, 365, 7), "dist_km": np.arange(0, 2000, 40)}
rows = []
for v, e in edges.items():
    src = A if v in A else T
    for per in ("early", "late"):
        x = src.loc[src.period == per, v].values.astype(float)
        h, _ = np.histogram(x, bins=e, density=True)
        for k in range(len(h)):
            rows.append(dict(var=v, period=per, left=e[k], right=e[k + 1], density=h[k]))
        rows.append(dict(var=v, period=per, left=np.nan, right=np.nan, density=np.nan, n=len(x), median=np.median(x),
                         mean=x.mean(), q90=np.quantile(x, .9)))
pd.DataFrame(rows).to_csv(out / "e6_distributions.csv", index=False)
print(pd.DataFrame(rows).dropna(subset=["n"])[["var", "period", "n", "median", "mean", "q90"]].to_string(), flush=True)

# (2) longitude-bin monthly series and trends
ym = np.arange(2001 * 12, 2023 * 12 + 8)
months = pd.DataFrame({"year": ym // 12, "month": ym % 12 + 1}); months["x"] = months.year + (months.month - 0.5) / 12
months["ndays"] = pd.to_datetime(months[["year", "month"]].assign(day=1)).dt.days_in_month.values
def series(d):
    gg = d.groupby("ym")
    s = pd.DataFrame({"n_obs": gg.size(), "amp": gg.amplitude.mean(), "rad": gg.speed_radius.mean(), "spd": gg.speed_average.mean(),
                      "amp2": gg.amplitude.apply(lambda a: np.sum(a ** 2))}).reindex(ym)
    s["present_per_day"] = s.n_obs.values / months.ndays.values; s["amp2_per_day"] = s.amp2.values / months.ndays.values
    return s[["present_per_day", "amp", "rad", "spd", "amp2_per_day"]].reset_index(drop=True)
VARS = ["present_per_day", "amp", "rad", "spd", "amp2_per_day"]
allser = []; rows = []
def tr(label, b, s):
    for v in VARS:
        y = s[v].values.astype(float); m = float(np.nanmean(y))
        o = R.trend_bcz(months.x.values, y, month=months.month.values, n_outer=NOUT)
        rows.append(dict(region=label, bin=b, lon_c=(-175 + 10 * b) if b >= 0 else np.nan, var=v, mean=m,
                         pct=1000 * o["slope"] / m, pct_lo=1000 * o["ci_lo"] / m, pct_hi=1000 * o["ci_hi"] / m, p=o["p"], n_valid=int(np.isfinite(y).sum())))
    print(label, b, " ".join(f"{r['var']}:{r['pct']:+.1f}(p{r['p']:.2f})" for r in rows[-5:]), flush=True)
s = series(A); s.insert(0, "bin", -1); allser.append(s); tr("ACC", -1, s)
for b in range(36):
    s = series(A[A.bin == b]); s.insert(0, "bin", b); allser.append(s); tr("bin", b, s)
S = pd.concat(allser, ignore_index=True); S["x"] = np.tile(months.x.values, 37); S.to_csv(out / "e6_lonbin_series.csv", index=False)
Rz = pd.DataFrame(rows)
Rz["fdr10"] = False
for v in VARS:
    m = (Rz.region == "bin") & (Rz["var"] == v)
    Rz.loc[m, "fdr10"] = R.bh_fdr(Rz.loc[m, "p"].values, 0.10)
Rz.to_csv(out / "e6_lonbin_trends.csv", index=False)
print("done", flush=True)
