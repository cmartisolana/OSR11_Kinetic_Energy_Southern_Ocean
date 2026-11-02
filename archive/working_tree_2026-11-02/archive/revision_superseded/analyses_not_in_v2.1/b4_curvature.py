"""B4: curvature of the monthly jet-core path (Referee Gutierrez-Villanueva, major comment 3: "estimate the curvature of
the fronts ... and look at how it changes in time at specific meanders"). Same path as B3 (b3_alongstream.py): monthly
+0.18 m contour of ADT relative to its ACC-band mean, pieces spanning > 5 deg, resampled every 10 km and smoothed along
the arc with a 150 km Gaussian; in each 5 deg bin the longest piece covering the bin is used.
Curvature kappa = d(theta)/ds, theta = direction of the path in local km (dx = R cos(lat) dlon, dy = R dlat), ds = 10 km.
Per month and bin: rms_kappa (rad per 1000 km) = RMS of kappa over the bin (meander tightness); turning = integral of
|kappa| ds over the bin (rad). Half-power curvature was tested separately (SI S4, gates G1/G2: usable only at SEIR).
Trends per bin: trend_bcz (300 null series), %/decade, FDR 10 % over bins (bins with >= 80 % monthly coverage).
Usage: python b4_curvature.py <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, contourpy
from scipy.ndimage import gaussian_filter1d
from multiprocessing import Pool
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
LEVEL, SMOOTH, DS, RE = 0.18, 150.0, 10.0, 6371.0
EDGES = np.arange(-180, 181, 5.0); NB = EDGES.size - 1; CEN = EDGES[:-1] + 2.5
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
ACCM = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice
EXCL = (mean_adt > 0.45) | (lat[:, None] > -40) | ~np.isfinite(mean_adt)
PAD = int(30 / np.diff(lon).mean()); lonw = np.r_[lon[-PAD:] - 360, lon, lon[:PAD] + 360]
wrap = lambda a: np.c_[a[:, -PAD:], a, a[:, :PAD]]


def pieces(rel):
    gen = contourpy.contour_generator(lonw, lat, np.ma.masked_invalid(wrap(rel)), name="serial"); out_ = []
    for L in gen.lines(LEVEL):
        if np.ptp(L[:, 0]) < 5: continue
        if L[0, 0] > L[-1, 0]: L = L[::-1]
        la0 = np.radians(L[:, 1]); x = RE * np.cos(la0) * np.radians(L[:, 0]); y = RE * np.radians(L[:, 1])
        s = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
        if s[-1] < 5 * DS: continue
        si = np.arange(0, s[-1], DS); lo_i = np.interp(si, s, L[:, 0]); la_i = np.interp(si, s, L[:, 1])
        sg = SMOOTH / DS
        out_.append((gaussian_filter1d(lo_i, sg, mode="nearest"), gaussian_filter1d(la_i, sg, mode="nearest")))
    return out_


def curv(lo_s, la_s):
    dx = RE * np.cos(np.radians(0.5 * (la_s[1:] + la_s[:-1]))) * np.radians(np.diff(lo_s)); dy = RE * np.radians(np.diff(la_s))
    th = np.unwrap(np.arctan2(dy, dx)); ds = np.hypot(dx, dy)
    k = np.diff(th) / (0.5 * (ds[1:] + ds[:-1]))                                   # rad/km at interior points
    return np.r_[np.nan, k, np.nan], np.r_[ds, np.nan]


def month(f):
    with xr.open_dataset(f) as d: A = d.adt_mean.values.astype(float)
    rel = A - np.nanmean(A[ACCM]); rel[EXCL] = np.nan
    rk = np.full(NB, np.nan); tu = np.full(NB, np.nan); best = np.zeros(NB)
    for lo_s, la_s in pieces(rel):
        kap, ds = curv(lo_s, la_s); lo_m = ((lo_s + 180) % 360) - 180
        inside = (lo_s >= -180) & (lo_s < 180); span = np.ptp(lo_s)
        b = np.clip(np.digitize(lo_m, EDGES) - 1, 0, NB - 1)
        for k in np.unique(b[inside]):
            sel = inside & (b == k) & np.isfinite(kap)
            if sel.sum() < 3 or span <= best[k]: continue
            best[k] = span
            rk[k] = 1000 * np.sqrt(np.mean(kap[sel] ** 2)); tu[k] = np.nansum(np.abs(kap[sel]) * ds[sel])
    return rk, tu


T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)


def trend(args):
    v, k, y = args
    if np.isfinite(y).mean() < 0.8: return dict(lon=CEN[k], var=v, cov=float(np.isfinite(y).mean()))
    o = R.trend_bcz(T, y, month=MON, n_outer=300); m = float(np.nanmean(y)); f = 1000 / m
    return dict(lon=CEN[k], var=v, cov=float(np.isfinite(y).mean()), mean=m, slope=o["slope"], tr=o["slope"] * f,
                lo=o["ci_lo"] * f, hi=o["ci_hi"] * f, p=o["p"])


if __name__ == "__main__":
    sf = OUT / "b4_bin_series.npz"
    if not sf.exists():
        with Pool(2) as P: res = P.map(month, files)
        np.savez(sf, rms_kappa=np.array([r[0] for r in res]), turning=np.array([r[1] for r in res]), lon=CEN)
        print("series saved", flush=True)
    Z = np.load(sf)
    for v in ("rms_kappa",):
        with Pool(2) as P: rows = P.map(trend, [(v, k, Z[v][:, k]) for k in range(NB)])
        D = pd.DataFrame(rows); pv = D.p.values if "p" in D else np.full(len(D), np.nan); ok = np.isfinite(pv); q = np.zeros(len(D), bool)
        if ok.any():
            ps = np.sort(pv[ok]); m = ok.sum(); thr = ps[ps <= 0.10 * np.arange(1, m + 1) / m]
            q[ok] = pv[ok] <= (thr.max() if thr.size else -1)
        D["fdr10"] = q; D.to_csv(OUT / f"b4_trends_{v}.csv", index=False)
        print(v, "bins positive", int((D.tr > 0).sum()), "of", int(D.tr.notna().sum()), "p<0.05", int((D.p.round(3) < 0.05).sum()), "FDR", int(q.sum()), flush=True)
    print("done")
