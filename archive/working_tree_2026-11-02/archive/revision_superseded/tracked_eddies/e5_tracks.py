"""E5: track-following eddies at the five standing meanders (specification in e5_spec.md, fixed before any result).

Atlas subset META4.0 DT all-sat (south of 35 S, 2001-01 .. 2023-09-08), tracked observations only (track > 0).
Regions = dk13 longitude ranges (upstream / ridge / downstream) of the five meanders; a track is assigned by the
position of its first observation (observation_number == 0, a true birth inside the subset) inside the ACC band
(mean ADT -0.6..+0.2 m, ice-free in every month). Births in the band outside all ranges = "away".
  Q1 composites of amplitude / speed radius / rotational speed against longitude relative to the eastern edge of
     each ridge (-10..+25 deg, 1 deg bins; all band observations) and against days since birth (0-120 d) by birth class.
  Q2 lifetime, along-track propagation distance and net zonal displacement by birth class (births <= 2022-08-31).
  Q3 annual birth-year series 2001-2022 of mean lifetime, mean maximum amplitude and mean propagation distance for
     MEANDER (ridge + downstream, five features pooled) and AWAY; the test is the trend of MEANDER/mean - AWAY/mean
     (rstats.trend_bcz, annual, AR(1)); BH-FDR 10 % over the three contrasts. Per feature: descriptive.
Usage: python e5_tracks.py <subset_dir> <repo outputs dir> <revision dir> <out_dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, netCDF4 as nc

sub, repo, rev, out = (Path(a) for a in sys.argv[1:5]); out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(rev))
import rstats as R  # noqa: E402

T0, TEND = pd.Timestamp("2001-01-01"), pd.Timestamp("2023-09-08")
BIRTH_MAX = pd.Timestamp("2022-08-31")
MEANDERS = {"KP": ((55, 65), (65, 80), (80, 95)), "SEIR": ((120, 130), (130, 140), (140, 150)),
            "MR": ((150, 155), (155, 160), (160, 175)), "PAR": ((-155, -145), (-145, -115), (-115, -100)),
            "DP": ((-80, -70), (-70, -60), (-50, -35))}
CLASSES = ("upstream", "ridge", "downstream")


def load(pol):
    ds = nc.Dataset(sub / f"eddies_so_{pol}_all.nc")
    cols = ("track", "observation_number", "time", "latitude", "longitude", "amplitude", "speed_radius", "speed_average")
    df = pd.DataFrame({v: np.asarray(ds.variables[v][:]) for v in cols})
    df["time"] = pd.Timestamp("1950-01-01") + pd.to_timedelta(df.time.astype(int), unit="D")
    df["lon"] = ((df.longitude + 180) % 360) - 180
    df = df[df.track > 0].copy(); df["pol"] = pol
    return df


E = pd.concat([load("anticyclonic"), load("cyclonic")], ignore_index=True)
E["tid"] = E.pol.str[0] + E.track.astype(str)
E = E.sort_values(["tid", "time"]).reset_index(drop=True)
print("tracked observations", len(E), "tracks", E.tid.nunique(), flush=True)

# ACC band mask (as E1-E4)
adt = xr.open_dataset(repo / "mean_adt.nc").mean_adt
ice = xr.open_dataset(repo / "sea_ice_mask.nc").ice_mask.values > 0.5
acc = np.isfinite(adt.values) & (adt.values >= -0.6) & (adt.values <= 0.2) & ~ice
lat, lon = adt.latitude.values, adt.longitude.values
def in_acc(la, lo):
    i = np.clip(np.searchsorted(lat, la), 0, lat.size - 1); j = np.clip(np.searchsorted(lon, lo), 0, lon.size - 1)
    return acc[i, j]
E["acc"] = in_acc(E.latitude.values, E.lon.values)

# per-track summary
g = E.groupby("tid", sort=False)
dlon = g.lon.diff().values; dlat = g.latitude.diff().values
dlon = np.where(dlon > 180, dlon - 360, np.where(dlon < -180, dlon + 360, dlon))
step = np.hypot(dlat * 111.2, dlon * 111.2 * np.cos(np.radians(E.latitude.values)))
E["step_km"] = np.nan_to_num(step)
first = g.head(1).set_index("tid"); last = g.tail(1).set_index("tid")
T = pd.DataFrame({"pol": first.pol, "born_in_subset": first.observation_number.values == 0,
                  "t0": first.time, "t1": last.time, "lon0": first.lon, "lat0": first.latitude, "lon1": last.lon,
                  "acc0": first.acc.values, "n_obs": g.size(),
                  "lifetime_d": (last.time - first.time).dt.days + 1, "dist_km": g.step_km.sum(),
                  "amp_max": g.amplitude.max(), "amp_mean": g.amplitude.mean(), "rad_mean": g.speed_radius.mean(),
                  "spd_mean": g.speed_average.mean()})
dz = (T.lon1 - T.lon0).values; dz = np.where(dz > 180, dz - 360, np.where(dz < -180, dz + 360, dz))
T["dlon_deg"] = dz
T["censored"] = T.t1 >= TEND - pd.Timedelta(days=1)
T["feature"] = ""; T["cls"] = ""
for f, rngs in MEANDERS.items():
    for c, (a, b) in zip(CLASSES, rngs):
        m = (T.lon0 >= a) & (T.lon0 < b)
        T.loc[m, "feature"] = f; T.loc[m, "cls"] = c
T.loc[T.cls == "", "cls"] = "away"
T["cohort_ok"] = T.born_in_subset & T.acc0 & (T.t0 <= BIRTH_MAX)
T.to_csv(out / "e5_tracks.csv")
print("births in band", int((T.born_in_subset & T.acc0).sum()), "cohort (<= 2022-08)", int(T.cohort_ok.sum()),
      "censored among cohort", int((T.cohort_ok & T.censored).sum()), flush=True)

# ---- Q2 lifetime / distance / displacement by class (and by feature)
rows = []
def q2(sel, label, feature="all"):
    d = T[sel]
    for pol in ("both", "anticyclonic", "cyclonic"):
        dd = d if pol == "both" else d[d.pol == pol]
        if len(dd) < 20: continue
        rows.append(dict(feature=feature, cls=label, pol=pol, n=len(dd), n_censored=int(dd.censored.sum()),
                         lifetime_med=dd.lifetime_d.median(), lifetime_mean=dd.lifetime_d.mean(), lifetime_q75=dd.lifetime_d.quantile(.75),
                         frac_ge90d=(dd.lifetime_d >= 90).mean(), dist_med=dd.dist_km.median(), dist_mean=dd.dist_km.mean(),
                         dlon_med=dd.dlon_deg.median(), frac_eastward=(dd.dlon_deg > 0).mean(),
                         amp_max_med=dd.amp_max.median(), amp_mean_med=dd.amp_mean.median(), rad_med=dd.rad_mean.median(),
                         spd_med=dd.spd_mean.median(), births_per_year=len(dd) / ((BIRTH_MAX - T0).days / 365.25)))
ok = T.cohort_ok
for c in CLASSES + ("away",):
    q2(ok & (T.cls == c), c)
q2(ok & T.cls.isin(["ridge", "downstream"]), "meander")
for f in MEANDERS:
    for c in CLASSES:
        q2(ok & (T.feature == f) & (T.cls == c), c, f)
Q2 = pd.DataFrame(rows); Q2.to_csv(out / "e5_q2_by_class.csv", index=False)
print(Q2[Q2.pol == "both"][["feature", "cls", "n", "lifetime_med", "frac_ge90d", "dist_med", "dlon_med", "amp_max_med", "rad_med"]].to_string(), flush=True)

# ---- Q1a composites against longitude relative to the eastern edge of each ridge (all band observations)
A = E[E.acc]
rows = []
for f, (up, ridge, down) in MEANDERS.items():
    rel = A.lon.values - ridge[1]; rel = np.where(rel < -180, rel + 360, np.where(rel > 180, rel - 360, rel))
    for per, (y0, y1) in {"all": (2001, 2023), "early": (2001, 2011), "late": (2012, 2023)}.items():
        m = (rel >= -10) & (rel < 25) & (A.time.dt.year.values >= y0) & (A.time.dt.year.values <= y1)
        d = A[m].assign(rel=rel[m], bin=np.floor(rel[m]).astype(int))
        gg = d.groupby("bin")
        q = pd.DataFrame({"n": gg.size(), "amp_med": gg.amplitude.median(), "amp_q25": gg.amplitude.quantile(.25),
                          "amp_q75": gg.amplitude.quantile(.75), "rad_med": gg.speed_radius.median(), "spd_med": gg.speed_average.median(),
                          "amp2_sum_per_day": gg.amplitude.apply(lambda a: np.sum(a ** 2)) / ((pd.Timestamp(f"{min(y1,2023)}-09-01") if y1 >= 2023 else pd.Timestamp(f"{y1+1}-01-01")) - pd.Timestamp(f"{y0}-01-01")).days})
        q = q.reset_index().assign(feature=f, period=per)
        rows.append(q)
Q1a = pd.concat(rows, ignore_index=True); Q1a.to_csv(out / "e5_q1_longitude_composite.csv", index=False)

# ---- Q1b amplitude against days since birth, by birth class (survivor-weighted)
E2 = E.merge(T[["cls", "feature", "cohort_ok", "t0"]], left_on="tid", right_index=True)
E2 = E2[E2.cohort_ok]
E2["age"] = (E2.time - E2.t0).dt.days
rows = []
for c in CLASSES + ("away",):
    d = E2[(E2.cls == c) & (E2.age <= 120)]
    gg = d.groupby("age")
    rows.append(pd.DataFrame({"n": gg.size(), "amp_med": gg.amplitude.median(), "amp_mean": gg.amplitude.mean(),
                              "rad_med": gg.speed_radius.median(), "spd_med": gg.speed_average.median()}).reset_index().assign(cls=c))
Q1b = pd.concat(rows, ignore_index=True); Q1b.to_csv(out / "e5_q1_age_composite.csv", index=False)

# ---- Q3 annual birth-year series and the contrast test
T["year"] = T.t0.dt.year
C = T[ok]
def annual(sel):
    d = C[sel].groupby("year").agg(n=("lifetime_d", "size"), lifetime=("lifetime_d", "mean"), amp_max=("amp_max", "mean"),
                                   dist=("dist_km", "mean"), rad=("rad_mean", "mean"))
    return d.reindex(range(2001, 2023))
S = {"MEANDER": annual(C.cls.isin(["ridge", "downstream"])), "AWAY": annual(C.cls == "away")}
for f in MEANDERS:
    S[f] = annual((C.feature == f) & C.cls.isin(["ridge", "downstream"]))
ser = pd.concat({k: v for k, v in S.items()}, axis=1); ser.to_csv(out / "e5_q3_annual.csv")
x = np.arange(2001, 2023) + 0.5
rows = []
def tr(name, y, kind):
    o = R.trend_bcz(x, y, monthly=False, n_outer=300)
    m = float(np.nanmean(y))
    rows.append(dict(series=name, kind=kind, mean=m, slope=o["slope"], ci_lo=o["ci_lo"], ci_hi=o["ci_hi"], p=o["p"],
                     pct_per_decade=1000 * o["slope"] / m if kind == "level" else 1000 * o["slope"],
                     pct_lo=1000 * o["ci_lo"] / m if kind == "level" else 1000 * o["ci_lo"],
                     pct_hi=1000 * o["ci_hi"] / m if kind == "level" else 1000 * o["ci_hi"]))
    print(name, kind, round(rows[-1]["pct_per_decade"], 2), "p", round(o["p"], 3), flush=True)
for v in ("lifetime", "amp_max", "dist"):
    a, b = S["MEANDER"][v].values, S["AWAY"][v].values
    tr(f"MEANDER-AWAY_{v}", a / np.nanmean(a) - b / np.nanmean(b), "contrast")
primary = pd.DataFrame(rows)
primary["fdr10"] = R.bh_fdr(primary.p.values, 0.10)
rows = rows  # continue appending secondary
for k in ("MEANDER", "AWAY"):
    for v in ("lifetime", "amp_max", "dist", "rad", "n"):
        tr(f"{k}_{v}", S[k][v].values.astype(float), "level")
for f in MEANDERS:
    for v in ("lifetime", "amp_max", "dist"):
        y = S[f][v].values.astype(float)
        if np.isfinite(y).sum() >= 18:
            tr(f"{f}-AWAY_{v}", y / np.nanmean(y) - S["AWAY"][v].values / np.nanmean(S["AWAY"][v].values), "contrast")
Q3 = pd.DataFrame(rows); Q3["fdr10_primary"] = False
Q3.loc[:2, "fdr10_primary"] = primary.fdr10.values
Q3.to_csv(out / "e5_q3_trends.csv", index=False)
print("done", flush=True)
