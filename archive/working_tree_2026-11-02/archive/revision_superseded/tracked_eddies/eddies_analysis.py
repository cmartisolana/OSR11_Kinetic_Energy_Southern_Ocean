"""E1-E4: AVISO META4.0 DT eddy atlas (2001-01 .. 2023-08) against the KE/EKE and meander results.

Pre-specified (2026-10-13), before any result was looked at:
  E1  Generation sites: density of track births (first observation of each track) and of all observations
      on 1° cells, south of 35°S; compared with the meander sections and downstream boxes of Fig. 4.
  E2  What the EKE trend is made of (ACC band, 2001-01..2023-08, monthly series, revised trend test):
      eddies present per day, tracks born per month, mean amplitude, mean speed radius, mean rotational
      speed, and the "eddy kinetic energy proxy" sum(amplitude^2) per day; alongside the gridded EKE and KE
      (all-sat, a5) over the same months.
  E3  Meander wakes vs rest of the band: relative trends of count and amplitude in the union of the four
      downstream boxes (a5_energy BOXES) against the rest of the ACC band (difference of normalised series).
  E4  Shedding vs path distortion: monthly births in the meander section + downstream box of each meander
      (b1 definitions) against the path length (b1 series, same preparation as b1_covar), lags -3..6,
      surrogate p for the maximum |r| over lags 0-3 (path length leading).
Only tracked observations (track > 0) are used for track statistics; all observations for presence/density.
Usage: python eddies_analysis.py <subset_dir> <repo outputs dir> <out_dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, netCDF4 as nc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402
import b1_covar as BC  # noqa: E402

sub, repo, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
T0, T1 = pd.Timestamp("2001-01-01"), pd.Timestamp("2023-09-01")     # last full month: 2023-08


def load(pol):
    ds = nc.Dataset(sub / f"eddies_so_{pol}_all.nc")
    d = {v: np.asarray(ds.variables[v][:]) for v in ("track", "time", "latitude", "longitude", "amplitude",
                                                     "speed_radius", "speed_average")}
    df = pd.DataFrame(d); df["time"] = pd.Timestamp("1950-01-01") + pd.to_timedelta(df.time, unit="D")
    df["lon"] = ((df.longitude + 180) % 360) - 180
    df = df[(df.time >= T0) & (df.time < T1)].copy(); df["pol"] = pol
    return df


E = pd.concat([load("anticyclonic"), load("cyclonic")], ignore_index=True)
E["ym"] = E.time.dt.year * 12 + E.time.dt.month - 1
print("observations", len(E), flush=True)
# births
tracked = E[E.track > 0]
first = tracked.sort_values("time").groupby(["pol", "track"], sort=False).head(1).copy()
print("tracks", len(first), flush=True)

# ACC band mask at 1/8 deg (mean ADT -0.6..0.2, ice-free in all months)
adt = xr.open_dataset(repo / "mean_adt.nc").mean_adt
ice = xr.open_dataset(repo / "sea_ice_mask.nc").ice_mask.values > 0.5
acc = np.isfinite(adt.values) & (adt.values >= -0.6) & (adt.values <= 0.2) & ~ice
lat, lon = adt.latitude.values, adt.longitude.values
def in_acc(df):
    i = np.clip(np.searchsorted(lat, df.latitude.values), 0, lat.size - 1); j = np.clip(np.searchsorted(lon, df.lon.values), 0, lon.size - 1)
    return acc[i, j]
E["acc"] = in_acc(E); first["acc"] = in_acc(first)

# ---- E1 density maps (1 deg)
def density(df):
    H, _, _ = np.histogram2d(df.latitude, df.lon, bins=[np.arange(-80, -34.9, 1), np.arange(-180, 180.1, 1)])
    return H
nyears = (T1 - T0).days / 365.25
xr.Dataset({"births_per_year": (("lat", "lon"), density(first) / nyears),
            "obs_per_day": (("lat", "lon"), density(E) / (T1 - T0).days)},
           coords={"lat": np.arange(-79.5, -34.9, 1), "lon": np.arange(-179.5, 180, 1)}).to_netcdf(out / "e1_density.nc")

# ---- E2 ACC-band monthly series
ym = np.arange(2001 * 12, 2023 * 12 + 8)
months = pd.DataFrame({"year": ym // 12, "month": ym % 12 + 1}); months["x"] = months.year + (months.month - 0.5) / 12
months["ndays"] = pd.to_datetime(months[["year", "month"]].assign(day=1)).dt.days_in_month
A = E[E.acc]; Fb = first[first.acc]
def series(dfobs, dffirst, label):
    g = dfobs.groupby("ym")
    s = pd.DataFrame({"n_obs": g.size(), "amp": g.amplitude.mean(), "rad": g.speed_radius.mean(), "spd": g.speed_average.mean(),
                      "amp2": g.amplitude.apply(lambda a: np.sum(a ** 2))}).reindex(ym)
    s["births"] = dffirst.groupby("ym").size().reindex(ym).fillna(0).values
    s["present_per_day"] = s.n_obs.values / months.ndays.values; s["amp2_per_day"] = s.amp2.values / months.ndays.values
    s.columns = [f"{label}_{c}" for c in s.columns]
    return s.reset_index(drop=True)
S = pd.concat([months, series(A, Fb, "ACC")], axis=1)
# a5 EKE/KE over the same months
a5 = pd.read_csv(HERE / "outputs" / "a5" / "a5_series_allsat.csv")
a5 = a5[(a5.year * 12 + a5.month - 1 >= ym[0]) & (a5.year * 12 + a5.month - 1 <= ym[-1])].reset_index(drop=True)
S["ACC_EKE_allsat"] = a5.ACC_EKE.values; S["ACC_KE_allsat"] = a5.ACC_KE.values; S["ACC_SUB_allsat"] = a5.ACC_SUB.values

# ---- E3 wakes vs rest
BOXES = {"Kerguelen": (75, 95), "SEIR_Macquarie": (145, 175), "PAR": (-150, -125), "Drake": (-55, -35)}
wake = np.zeros(len(E), bool); wakef = np.zeros(len(first), bool)
for lo0, lo1 in BOXES.values():
    wake |= (E.lon.values >= lo0) & (E.lon.values < lo1); wakef |= (first.lon.values >= lo0) & (first.lon.values < lo1)
S = pd.concat([S, series(E[E.acc & wake], first[first.acc & wakef], "WAKE"), series(E[E.acc & ~wake], first[first.acc & ~wakef], "REST")], axis=1)
S.to_csv(out / "e2_series.csv", index=False)

rows = []
def tr(name, y, scale=1.0):
    o = R.trend_bcz(S.x.values, y, month=S.month.values, n_outer=300)
    m = float(np.nanmean(y))
    rows.append(dict(series=name, mean=m, slope=o["slope"] * scale, ci_lo=o["ci_lo"] * scale, ci_hi=o["ci_hi"] * scale,
                     pct_per_decade=1000 * o["slope"] / m if m else np.nan, p=o["p"]))
    print(name, round(rows[-1]["pct_per_decade"], 2), "%/dec p", round(o["p"], 3), flush=True)
for c in ["ACC_present_per_day", "ACC_births", "ACC_amp", "ACC_rad", "ACC_spd", "ACC_amp2_per_day",
          "ACC_EKE_allsat", "ACC_KE_allsat", "ACC_SUB_allsat"]:
    tr(c, S[c].values)
for v in ("present_per_day", "births", "amp", "amp2_per_day"):
    a, b = S[f"WAKE_{v}"].values, S[f"REST_{v}"].values
    d = a / np.nanmean(a) - b / np.nanmean(b)
    o = R.trend_bcz(S.x.values, d, month=S.month.values, n_outer=300)
    rows.append(dict(series=f"WAKE-REST_{v}", mean=np.nan, slope=1000 * o["slope"], ci_lo=1000 * o["ci_lo"], ci_hi=1000 * o["ci_hi"],
                     pct_per_decade=1000 * o["slope"], p=o["p"]))
    print(f"WAKE-REST_{v}", round(1000 * o["slope"], 2), "points/dec p", round(o["p"], 3), flush=True)
pd.DataFrame(rows).to_csv(out / "e2_e3_trends.csv", index=False)

# ---- E4 shedding vs path length
MEA = {"PAR": dict(mea=(-135.13, -109.38), mlat=(-62, -45), down=(-109.38, -80), band=(-60, -48)),
       "SEIR150": dict(mea=(140, 150), mlat=(-60, -46), down=(150, 158), band=(-60, -47)),
       "KERG": dict(mea=(76, 84), mlat=(-58, -48), down=(84, 94), band=(-57, -48))}
rng = np.random.default_rng(20261013)
e4 = []
def lc(x, y, L): return BC.lagcorr(x, y, L)[0] if L >= 0 else BC.lagcorr(y, x, -L)[0]
for name, c in MEA.items():
    inm = (first.lon >= c["mea"][0]) & (first.lon < c["mea"][1]) & (first.latitude >= c["mlat"][0]) & (first.latitude <= c["mlat"][1])
    ind = (first.lon >= c["down"][0]) & (first.lon < c["down"][1]) & (first.latitude >= c["band"][0]) & (first.latitude <= c["band"][1])
    b = pd.DataFrame({"year": months.year, "month": months.month, "x": months.x})
    b["births_mea"] = first[inm].groupby("ym").size().reindex(ym).fillna(0).values
    b["births_down"] = first[ind].groupby("ym").size().reindex(ym).fillna(0).values
    b["births_both"] = b.births_mea + b.births_down
    P = pd.read_csv(HERE / "outputs" / "b1" / f"b1_series_{name}.csv")
    P = P[(P.year * 12 + P.month - 1 >= ym[0]) & (P.year * 12 + P.month - 1 <= ym[-1])].reset_index(drop=True)
    L = BC.prep(P, "core_s150_length_km")
    for col in ("births_mea", "births_down", "births_both"):
        y = BC.prep(b, col)
        rr = {k: lc(L, y, k) for k in range(-3, 7)}
        obs = max(abs(rr[k]) for k in range(4))
        s = [max(abs(lc(BC.phase_surr(L, rng), y, k)) for k in range(4)) for _ in range(1000)]
        p = (1 + np.sum(np.array(s) >= obs)) / 1001
        e4.append(dict(meander=name, births=col, mean_per_month=float(b[col].mean()), **{f"r_lag{k}": rr[k] for k in rr}, max_r_0_3=obs, p_surr=p))
        print(name, col, "mean/month", round(b[col].mean(), 1), {k: round(v, 2) for k, v in rr.items()}, "p", round(p, 3), flush=True)
    b.to_csv(out / f"e4_births_{name}.csv", index=False)
pd.DataFrame(e4).to_csv(out / "e4_shedding_vs_length.csv", index=False)
print("done")
