"""B3: along-stream structure of the trends (replaces the B2 window correlation as the main-text meander analysis).

Specification fixed on 2026-10-15 before any result was computed. Answers Referee C ("estimate the curvature of the
fronts ... and look at how it changes in time at specific meanders; check whether the jets upstream of a standing
meander are accelerating") and Referee A ("relationship between the jets' position and EKE").

Jet path (as B2): monthly contour of ADT relative to its monthly mean over the ACC band (removes sea-level rise),
level +0.18 m (highest circumpolar mean speed; Subantarctic Front side). Subtropical branches excluded by masking
mean ADT > 0.45 m and latitudes north of 40S. All contour pieces spanning more than 5 deg of longitude are kept; each is
resampled every 10 km and smoothed along the arc with a 150 km Gaussian (the scale of b1/b2).
Bins: 72 bins of 5 deg longitude. For each month and bin, from the longest piece covering the bin:
  lat      mean latitude of the smoothed path in the bin                                   -> trend in km/decade
  sinu     sinuosity = path length in the bin / zonal length of the bin at the path latitude -> trend in %/decade
           (curvature measure: a meander that flexes or lengthens raises it)
  speed    mean geostrophic speed along the path in the bin (jet speed on the core streamline) -> %/decade
  eke      EKE (a5 cache) averaged over grid points within 150 km of the monthly path in the bin, ice-free  -> %/decade
Climatology per bin: mean depth under the 2001-2025 mean path (GEBCO), mean EKE.
Trends: rstats.trend_bcz on monthly series (deseasonalised Theil-Sen, bias-corrected AR noise, max of AIC/AR12 p);
FDR (Benjamini-Hochberg, 10 %) across the 72 bins of each variable. Bins with < 80 % monthly coverage are not tested.
Relations (pre-specified, descriptive): spatial correlation across bins between the EKE trend and (i) the speed trend,
(ii) the sinuosity trend, (iii) the latitude trend; significance from 2000 circular shifts of one profile along
longitude (keeps the spatial autocorrelation).
Usage: OSR11_CACHE=<cache> python b3_alongstream.py <revision/data dir> <repo outputs dir> <out_dir> [vars...]
       (vars: subset of lat sinu speed eke for the trend step; series are always computed/cached first)
"""
from __future__ import annotations
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[0])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import os, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, contourpy
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

data, repo, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
VARS = sys.argv[4:] or ["lat", "sinu", "speed", "eke"]
CACHE = Path(os.environ.get("OSR11_CACHE", _CACHE))
LEVEL, SMOOTH, DS, BAND_KM, RE = 0.18, 150.0, 10.0, 150.0, 6371.0
EDGES = np.arange(-180, 181, 5.0); NB = EDGES.size - 1; CEN = EDGES[:-1] + 2.5
lat = np.load(CACHE / "lat.npy"); lon = np.load(CACHE / "lon.npy")
mean_adt = xr.open_dataset(repo / "mean_adt.nc").mean_adt.values
ice = xr.open_dataset(repo / "sea_ice_mask.nc").ice_mask.values > 0.5
ACCM = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice
EXCL = (mean_adt > 0.45) | (lat[:, None] > -40) | ~np.isfinite(mean_adt)
files = [data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
NT = len(files)
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)}); t["x"] = t.year + (t.month - 0.5) / 12

# wrapped grid (+-30 deg) so that contours cross the dateline
PAD = int(30 / np.diff(lon).mean())
lonw = np.r_[lon[-PAD:] - 360, lon, lon[:PAD] + 360]
def wrap(a): return np.c_[a[:, -PAD:], a, a[:, :PAD]]


def pieces(rel):
    gen = contourpy.contour_generator(lonw, lat, np.ma.masked_invalid(wrap(rel)), name="serial")
    out_ = []
    for L in gen.lines(LEVEL):
        if np.ptp(L[:, 0]) < 5: continue
        if L[0, 0] > L[-1, 0]: L = L[::-1]
        # resample + smooth in local km
        la0 = np.radians(L[:, 1]); x = RE * np.cos(la0) * np.radians(L[:, 0]); y = RE * np.radians(L[:, 1])
        s = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
        if s[-1] < 5 * DS: continue
        si = np.arange(0, s[-1], DS)
        lo_i = np.interp(si, s, L[:, 0]); la_i = np.interp(si, s, L[:, 1])
        sg = SMOOTH / DS
        lo_s = gaussian_filter1d(lo_i, sg, mode="nearest"); la_s = gaussian_filter1d(la_i, sg, mode="nearest")
        out_.append((lo_s, la_s))
    return out_


def bin_metrics(pcs, spd, E, tree_pts):
    res = {v: np.full(NB, np.nan) for v in ("lat", "sinu", "speed", "eke")}
    best = np.zeros(NB)
    for lo_s, la_s in pcs:
        lo_m = ((lo_s + 180) % 360) - 180
        # only the copy inside [-180,180) of the wrapped piece is used: keep points with original lon in range
        inside = (lo_s >= -180) & (lo_s < 180)
        span = np.ptp(lo_s)
        b = np.clip(np.digitize(lo_m, EDGES) - 1, 0, NB - 1)
        for k in np.unique(b[inside]):
            sel = inside & (b == k)
            if sel.sum() < 3 or span <= best[k]: continue
            best[k] = span
            la_b = la_s[sel]; lo_b = lo_s[sel]
            seg = np.hypot(RE * np.cos(np.radians(la_b[1:])) * np.radians(np.diff(lo_b)), RE * np.radians(np.diff(la_b)))
            zonal = RE * np.cos(np.radians(la_b.mean())) * np.radians(np.ptp(lo_b))
            res["lat"][k] = la_b.mean()
            res["sinu"][k] = seg.sum() / zonal if np.ptp(lo_b) > 3.5 else np.nan
            ii = np.clip(np.searchsorted(lat, la_b), 0, lat.size - 1); jj = np.clip(np.searchsorted(lon, lo_m[sel]), 0, lon.size - 1)
            res["speed"][k] = np.nanmean(spd[ii, jj])
            # EKE within BAND_KM of the path piece in this bin
            px = RE * np.cos(np.radians(la_b)) * np.radians(lo_b); py = RE * np.radians(la_b)
            res["eke"][k] = band_eke(px, py, lo_b, E, tree_pts, k)
    return res


# grid points for the EKE band: rows -65..-38, all longitudes, ice-free; projected per bin on the fly
GI, GJ = np.nonzero((lat[:, None] <= -38) & np.ones((1, lon.size), bool) & ~ice)
GLAT, GLON = lat[GI], lon[GJ]; GW = np.cos(np.radians(GLAT))
GBIN = np.clip(np.digitize(GLON, EDGES) - 1, 0, NB - 1)
BIN_IDX = [np.nonzero(GBIN == k)[0] for k in range(NB)]


def band_eke(px, py, lo_b, E, _tree, k):
    idx = BIN_IDX[k]
    lo_c = np.mean(lo_b)
    glo = GLON[idx] + 360 * np.round((lo_c - GLON[idx]) / 360)     # same branch as the path copy
    gx = RE * np.cos(np.radians(GLAT[idx])) * np.radians(glo); gy = RE * np.radians(GLAT[idx])
    d, _ = cKDTree(np.c_[px, py]).query(np.c_[gx, gy], distance_upper_bound=BAND_KM + 1)
    ok = np.isfinite(d) & (d <= BAND_KM)
    e = E[GI[idx][ok], GJ[idx][ok]]; w = GW[idx][ok]; f = np.isfinite(e)
    return float(np.sum(e[f] * w[f]) / np.sum(w[f])) if f.sum() > 10 else np.nan


ser_file = out / "b3_bin_series.npz"
if not ser_file.exists():
    eke = np.load(CACHE / "eke.npy", mmap_mode="r")
    S = {v: np.full((NT, NB), np.nan) for v in ("lat", "sinu", "speed", "eke")}
    ADT_sum = np.zeros_like(mean_adt); cnt = 0
    for k, f in enumerate(files):
        with xr.open_dataset(f) as ds:
            A = ds.adt_mean.values.astype(np.float64); U = ds.u_mean.values; V = ds.v_mean.values
        rel = A - np.nanmean(A[ACCM]); rel[EXCL] = np.nan
        ADT_sum += np.nan_to_num(rel); cnt += 1
        spd = np.hypot(U, V)
        E = np.asarray(eke[k], np.float32) * 1e4
        r = bin_metrics(pieces(rel), spd, E, None)
        for v in S: S[v][k] = r[v]
        if k % 25 == 0: print("month", k, {v: int(np.isfinite(S[v][k]).sum()) for v in S}, flush=True)
    np.savez(ser_file, **S)
    # climatology along the mean path
    rel_c = mean_adt - np.nanmean(mean_adt[ACCM]); rel_c[EXCL] = np.nan
    pc = pieces(rel_c)
    bat = xr.open_dataset(data / "bathymetry_so.nc").elevation
    clim = []
    for k in range(NB):
        pts = [(lo, la) for lo_s, la_s in pc for lo, la in zip(lo_s, la_s) if -180 <= lo < 180 and EDGES[k] <= lo < EDGES[k + 1]]
        if not pts: clim.append(dict(lon=CEN[k], depth=np.nan, lat_mean=np.nan)); continue
        lo_, la_ = np.array(pts).T
        dep = -bat.interp(longitude=xr.DataArray(lo_), latitude=xr.DataArray(la_)).values
        clim.append(dict(lon=CEN[k], depth=float(np.nanmean(dep)), depth_min=float(np.nanmin(dep)), lat_mean=float(la_.mean())))
    C = pd.DataFrame(clim)
    C["eke_mean"] = np.nanmean(S["eke"], 0); C["speed_mean"] = np.nanmean(S["speed"], 0); C["sinu_mean"] = np.nanmean(S["sinu"], 0)
    for v in S: C[f"cov_{v}"] = np.isfinite(S[v]).mean(0)
    C.to_csv(out / "b3_climatology.csv", index=False)
    print("series saved", flush=True)
S = dict(np.load(ser_file))

for v in [v for v in VARS if v in S]:
    rows = []
    for k in range(NB):
        y = S[v][:, k]
        if np.isfinite(y).mean() < 0.8:
            rows.append(dict(lon=CEN[k], var=v, cov=float(np.isfinite(y).mean()))); continue
        o = R.trend_bcz(t.x.values, y, month=t.month.values, n_outer=300); m = float(np.nanmean(y))
        f = 111.195 * 10 if v == "lat" else 1000 / m        # km/decade for lat, %/decade otherwise
        rows.append(dict(lon=CEN[k], var=v, cov=float(np.isfinite(y).mean()), mean=m, slope=o["slope"],
                         tr=o["slope"] * f, lo=o["ci_lo"] * f, hi=o["ci_hi"] * f, p=o["p"]))
        print(v, CEN[k], round(rows[-1]["tr"], 2), round(o["p"], 3), flush=True)
    D = pd.DataFrame(rows)
    pv = D.p.values; ok = np.isfinite(pv); q = np.full(len(D), False)
    if ok.any():
        ps = np.sort(pv[ok]); m = ok.sum(); thr = ps[ps <= 0.10 * np.arange(1, m + 1) / m]
        q[ok] = pv[ok] <= (thr.max() if thr.size else -1)
    D["fdr10"] = q
    D.to_csv(out / f"b3_trends_{v}.csv", index=False)
print("done")
