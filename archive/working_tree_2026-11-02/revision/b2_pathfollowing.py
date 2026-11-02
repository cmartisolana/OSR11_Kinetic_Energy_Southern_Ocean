"""B2: box-free versions of the meander -> EKE test (three alternatives), all-sat 1/8 deg, 2001-2025.

Specification fixed on 2026-10-14 before any result was looked at.

Common machinery (replaces the hand-drawn lat-lon boxes of b1):
  * 36 overlapping longitude windows, 30 deg wide, centred every 10 deg (centre c: lon in [c-15, c+15]).
    Latitude range of a window = extent of the mean ACC band (mean ADT in [-0.6, 0.2] m, Product 1
    climatology) within the window, padded by 2 deg and clipped to [-64.9, -36].
  * Jet-core path: ONE circumpolar relative-ADT level (ADT minus its monthly mean over the ACC band, which
    removes sea-level rise), chosen once from the 2001-2025 climatology as the level whose contour carries the
    highest mean geostrophic speed around the ACC (candidates -0.5..0.3 m, step 0.02); per window and month the
    longest contour piece spanning >= 90 % of the window width, smoothed along the arc at 150 km
    (b1_meander.path_metrics). [Revised 2026-10-14 after a first pass with a per-window core level picked the
    Brazil-Malvinas confluence and other non-ACC jets in several windows.]
    Path length (km) is the measure of distortion, as in b1.
  * "Path band" = grid points within 150 km of the monthly path (nearest-neighbour distance in local km
    coordinates), excluding points with ice fraction > 0.10 in any month. EKE from the a5 cache.
  * Series preparation and test as in b1_covar.py: monthly climatology removed, Theil-Sen trend removed,
    2-month running mean; lags 0-3 months (x leads y); surrogate p (phase-randomised x, 1000) for the
    maximum |r| over the four lags; windows with path coverage < 80 % of months are reported but not tested.

Alternative 2 (circumpolar sliding window, primary):
  for every window w: r(length_w -> EKE in the band of w) and r(length_w -> EKE in the band of w+1, the
  next window downstream, centre c+10). Curves of r (best lag) against longitude; pooled pre-specified
  test = Stouffer z over the 12 non-overlapping windows (centres -170, -140, ..., 160) for each link,
  with surrogate p for the max over lags.

Alternative 1 (path-following, three meander windows: c = 80 (Kerguelen), 150 (SEIR), -120 (PAR)):
  EKE in the path band binned by 1 deg of longitude (along-stream coordinate), from the window's own path
  and, beyond its eastern edge, from the paths of the next windows; r(length_w -> EKE(lon bin)) at lags 0-3,
  giving the EKE response as a function of distance downstream without a downstream box.

Alternative 3 (tracked eddies, no EKE): births (first observation of META4.0 tracks) within 150 km of the
  monthly path of window w, counted per month; r(length_w -> births_w) at lags 0-3, 2001-01..2023-08.

Usage: OSR11_CACHE=<cache> [B2_LEVEL=<m>] python b2_pathfollowing.py <revision/data dir> <repo outputs dir> <out_dir>
"""
from __future__ import annotations
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[0])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import os, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, netCDF4 as nc
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b1_meander as B  # noqa: E402
import b1_covar as BC  # noqa: E402

data, repo, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
CACHE = Path(os.environ.get("OSR11_CACHE", _CACHE))
WIDTH, STEP, BAND_KM, SMOOTH = 30.0, 10.0, 150.0, 150.0
CENTRES = np.arange(-180, 180, STEP)                     # 36 windows
MEANDERS = {"KERG": 80, "SEIR150": 150, "PAR": -120}
NS = 1000
rng = np.random.default_rng(20261014)
R_KM = 6371.0

lat = np.load(CACHE / "lat.npy"); lon = np.load(CACHE / "lon.npy")
eke = np.load(CACHE / "eke.npy", mmap_mode="r")           # (300, 240, 2880) m2/s2
mean_adt = xr.open_dataset(repo / "mean_adt.nc").mean_adt.values
ice = xr.open_dataset(repo / "sea_ice_mask.nc").ice_mask.values > 0.5
files = [data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
NT = len(files)
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
t["x"] = t.year + (t.month - 0.5) / 12
print("loading ADT ...", flush=True)
ADT = np.empty((NT, lat.size, lon.size), np.float32)
UM = np.zeros((lat.size, lon.size)); VM = np.zeros((lat.size, lon.size))
for k, f in enumerate(files):
    with xr.open_dataset(f) as ds:
        ADT[k] = ds.adt_mean.values; UM += np.nan_to_num(ds.u_mean.values) / NT; VM += np.nan_to_num(ds.v_mean.values) / NT
print("ADT loaded", flush=True)
ACCM = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice
for k in range(NT):
    ADT[k] -= np.nanmean(ADT[k][ACCM])
RELC = np.nanmean(ADT, 0); SPDC = np.hypot(UM, VM)
def level_speed(lev):
    gen = contourpy.contour_generator(lon, lat, np.ma.masked_invalid(RELC), name="serial")
    v, n = 0.0, 0
    for Lp in gen.lines(lev):
        if np.ptp(Lp[:, 0]) < 5:          # skip closed/short pieces (eddies)
            continue
        ii = np.clip(np.searchsorted(lat, Lp[:, 1]), 0, lat.size - 1); jj = np.clip(np.searchsorted(lon, Lp[:, 0]), 0, lon.size - 1)
        sp = SPDC[ii, jj]; v += np.nansum(sp); n += np.isfinite(sp).sum()
    return v / max(n, 1)
import contourpy
LEVELS = np.round(np.arange(-0.5, 0.31, 0.02), 2)
LS = {lev: level_speed(lev) for lev in LEVELS}
CORE = float(os.environ["B2_LEVEL"]) if "B2_LEVEL" in os.environ else max(LS, key=LS.get)   # B2_LEVEL: sensitivity run at a fixed level
print("circumpolar core level (relative ADT, m):", CORE, "mean speed", round(LS[CORE], 3), flush=True)
pd.DataFrame({"level": list(LS), "mean_speed": list(LS.values())}).to_csv(out / "b2_level_speed.csv", index=False)

# eddy births (alternative 3)
def births():
    fr = []
    for pol in ("anticyclonic", "cyclonic"):
        ds = nc.Dataset(data / "aviso_eddies" / "subset" / f"eddies_so_{pol}_all.nc")
        d = pd.DataFrame({v: np.asarray(ds.variables[v][:]) for v in ("track", "time", "latitude", "longitude")})
        d = d[d.track > 0]; d["time"] = pd.Timestamp("1950-01-01") + pd.to_timedelta(d.time, unit="D")
        fr.append(d.sort_values("time").groupby("track", sort=False).head(1))
    f = pd.concat(fr); f = f[(f.time >= "2001-01-01") & (f.time < "2023-09-01")].copy()
    f["lon"] = ((f.longitude + 180) % 360) - 180; f["k"] = (f.time.dt.year - 2001) * 12 + f.time.dt.month - 1
    return f
FIRST = births(); print("births", len(FIRST), flush=True)
NT_E = 2023 * 12 + 8 - 2001 * 12                           # months with eddy data

def window_idx(c):
    """column indices of the window (handles the dateline by wrapping), and the lon values relative to c."""
    rel = ((lon - c + 180) % 360) - 180
    j = np.nonzero(np.abs(rel) <= WIDTH / 2)[0]
    order = np.argsort(rel[j]); j = j[order]
    return j, rel[j] + c                                    # unwrapped longitudes, monotonic

def to_km(lo, la, c, la0):
    return R_KM * np.cos(np.radians(la)) * np.radians(lo - c), R_KM * np.radians(la - la0)

W = {}
for c in CENTRES:
    j, lonw = window_idx(c)
    band = np.isfinite(mean_adt[:, j]) & (mean_adt[:, j] >= -0.6) & (mean_adt[:, j] <= 0.2)
    rows = np.nonzero(band.any(1))[0]
    la0, la1 = max(lat[rows.min()] - 2, -64.9), min(lat[rows.max()] + 2, -36.0)
    i = np.nonzero((lat >= la0) & (lat <= la1))[0]
    latw = lat[i]
    core = CORE
    LON, LAT = np.meshgrid(lonw, latw)
    gx, gy = to_km(LON, LAT, c, la0)
    icew = ice[np.ix_(i, j)]
    W[c] = dict(i=i, j=j, lonw=lonw, latw=latw, core=core, gx=gx, gy=gy, ice=icew, la0=la0,
                w=np.cos(np.radians(LAT)))
    print(f"window {c:5.0f}: lat {la0:.1f}..{la1:.1f} core level {core}", flush=True)

# monthly paths, band EKE, lon-binned band EKE, births
BINS = np.arange(-180, 181, 1.0)
res = {c: dict(length=np.full(NT, np.nan), eke_band=np.full(NT, np.nan), births=np.full(NT_E, np.nan),
               eke_bins=np.full((NT, 360), np.nan)) for c in CENTRES}
for k in range(NT):
    E = np.asarray(eke[k], np.float32) * 1e4
    for c in CENTRES:
        w = W[c]
        if w["core"] is None:
            continue
        i, j = w["i"], w["j"]
        rel = ADT[k][np.ix_(i, j)]
        line = B.spanning_contour(rel, w["latw"], w["lonw"], w["core"])
        if line is None:
            continue
        m = B.path_metrics(line, smooth=SMOOTH)
        if m is None:
            continue
        res[c]["length"][k] = m["length_km"]
        px, py = to_km(line[:, 0], line[:, 1], c, w["la0"])
        # resample path at 10 km for the tree
        s = np.r_[0, np.cumsum(np.hypot(np.diff(px), np.diff(py)))]
        si = np.arange(0, s[-1], 10.0); px, py = np.interp(si, s, px), np.interp(si, s, py)
        tree = cKDTree(np.c_[px, py])
        d, _ = tree.query(np.c_[w["gx"].ravel(), w["gy"].ravel()])
        inb = (d.reshape(w["gx"].shape) <= BAND_KM) & ~w["ice"]
        Ew = E[np.ix_(i, j)]
        ok = inb & np.isfinite(Ew)
        res[c]["eke_band"][k] = np.sum(Ew[ok] * w["w"][ok]) / np.sum(w["w"][ok])
        lw = np.broadcast_to(w["lonw"], Ew.shape)
        bi = np.digitize(((lw[ok] + 180) % 360) - 180, BINS) - 1
        acc = np.bincount(bi, weights=Ew[ok] * w["w"][ok], minlength=360); cnt = np.bincount(bi, weights=w["w"][ok], minlength=360)
        with np.errstate(invalid="ignore", divide="ignore"):
            res[c]["eke_bins"][k] = np.where(cnt > 0, acc / cnt, np.nan)
        if k < NT_E:
            f = FIRST[FIRST.k == k]
            rl = ((f.lon.values - c + 180) % 360) - 180
            sel = np.abs(rl) <= WIDTH / 2
            if sel.any():
                bx, by = to_km(rl[sel] + c, f.latitude.values[sel], c, w["la0"])
                db, _ = tree.query(np.c_[bx, by]); res[c]["births"][k] = int((db <= BAND_KM).sum())
            else:
                res[c]["births"][k] = 0
    if k % 50 == 0:
        print("month", k, flush=True)

# save series
S = t.copy()
for c in CENTRES:
    S[f"len_{c:.0f}"] = res[c]["length"]; S[f"eke_{c:.0f}"] = res[c]["eke_band"]
    S[f"births_{c:.0f}"] = np.r_[res[c]["births"], np.full(NT - NT_E, np.nan)]
S.to_csv(out / "b2_series.csv", index=False)
np.save(out / "b2_eke_bins.npy", np.array([res[c]["eke_bins"] for c in CENTRES]))

# tests
def test(x, y, nsurr=NS):
    rr = [BC.lagcorr(x, y, L)[0] for L in range(4)]
    obs = max(abs(r) for r in rr)
    s = [max(abs(BC.lagcorr(BC.phase_surr(x, rng), y, L)[0]) for L in range(4)) for _ in range(nsurr)]
    return rr, obs, (1 + np.sum(np.array(s) >= obs)) / (nsurr + 1)

rows = []
nxt = {c: ((c + STEP + 180) % 360) - 180 for c in CENTRES}
for c in CENTRES:
    cov = float(np.isfinite(res[c]["length"]).mean())
    row = dict(centre=c, coverage=cov, core=W[c]["core"])
    if cov >= 0.8:
        P = S[["year", "month", "x"]].copy(); P["L"] = res[c]["length"]; P["E"] = res[c]["eke_band"]
        P["En"] = res[nxt[c]]["eke_band"]
        xL = BC.prep(P, "L"); yE = BC.prep(P, "E"); yN = BC.prep(P, "En")
        for name, y in (("same", yE), ("next", yN)):
            if np.isfinite(y).mean() < 0.8:
                continue
            rr, obs, p = test(xL, y)
            row.update({f"r_{name}_lag{L}": rr[L] for L in range(4)}); row[f"rmax_{name}"] = obs; row[f"p_{name}"] = p
        Pe = P.iloc[:NT_E].copy(); Pe["Bi"] = res[c]["births"]
        if np.isfinite(Pe.Bi).mean() >= 0.8:
            xLe = BC.prep(Pe, "L"); yB = BC.prep(Pe, "Bi")
            rr, obs, p = test(xLe, yB)
            row.update({f"r_births_lag{L}": rr[L] for L in range(4)}); row["rmax_births"] = obs; row["p_births"] = p
    rows.append(row); print(row, flush=True)
T = pd.DataFrame(rows); T.to_csv(out / "b2_window_tests.csv", index=False)

# pooled pre-specified test over the 12 non-overlapping windows
from scipy import stats
pool = []
for link in ("same", "next", "births"):
    sub = T[(T.centre % 30 == 10) & T[f"p_{link}"].notna()] if f"p_{link}" in T else T.iloc[0:0]
    # centres -170, -140, ..., 160 are the ones with centre % 30 == 10 (in [-180,180): -170%30=10)
    if len(sub) == 0:
        continue
    zs = []
    for L in range(4):
        r = sub[f"r_{link}_lag{L}"].values; n = 290
        z = np.arctanh(r) * np.sqrt(n - 3); zs.append(z.sum() / np.sqrt(len(z)))
    zmax = max(abs(z) for z in zs); p_lag = [2 * stats.norm.sf(abs(z)) for z in zs]
    pool.append(dict(link=link, n_windows=len(sub), **{f"z_lag{L}": zs[L] for L in range(4)},
                     **{f"p_lag{L}": p_lag[L] for L in range(4)}, frac_positive_best=float(np.mean(sub[[f"r_{link}_lag{L}" for L in range(4)]].values.max(1) > 0))))
    print(pool[-1], flush=True)
pd.DataFrame(pool).to_csv(out / "b2_pooled.csv", index=False)

# alternative 1: lon-binned response for the three meander windows
alt1 = []
for name, c in MEANDERS.items():
    P = S[["year", "month", "x"]].copy(); P["L"] = res[c]["length"]; xL = BC.prep(P, "L")
    order = [c, nxt[c], nxt[nxt[c]]]
    for b in range(360):
        lonb = BINS[b] + 0.5
        rel = ((lonb - c + 180) % 360) - 180
        if rel < -WIDTH / 2 or rel > WIDTH / 2 + 2 * STEP + 5:
            continue
        y = None
        for cc in order:                                   # own path first, then the next windows' paths
            yy = res[cc]["eke_bins"][:, b]
            if np.isfinite(yy).mean() >= 0.8:
                y = yy; src = cc; break
        if y is None:
            continue
        P["Y"] = y; yb = BC.prep(P, "Y")
        rr = [BC.lagcorr(xL, yb, L)[0] for L in range(4)]
        alt1.append(dict(meander=name, centre=c, lon=lonb, rel_lon=rel, path_src=src, **{f"r_lag{L}": rr[L] for L in range(4)}))
pd.DataFrame(alt1).to_csv(out / "b2_alt1_lonbins.csv", index=False)
print("done")
