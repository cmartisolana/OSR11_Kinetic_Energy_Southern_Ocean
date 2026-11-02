#!/usr/bin/env python3
"""Reproduce Liu et al. (2026, JGR Oceans, doi:10.1029/2025JC023527) PAR meander metrics.

Input: monthly detection counts from extract_liu.py (revision/data/liu_par/<grid>/YYYY_MM.nc).
Steps (Liu et al. Sect. 2.2; window from their 04-...m script: 120 days centred on the 15th):
  F_m      4-month occurrence = sum of daily detections over months m-1..m+2 (calendar-month
           approximation of their 120-day window; centring offset ~2 weeks)
  axis     per longitude (150W-80W): latitude (60S-48S) of max F_m; ties -> larger summed
           gradient -> closest to previous month's axis
  width    distance between the nearest zero-occurrence latitudes north and south of the axis (km)
  speed    4-month mean geostrophic speed at the axis cell
  position axis latitude, section mean
  peaks/troughs: southernmost / northernmost axis latitude within Liu Table 1 windows
           (their 3-longitude and 70 %-persistence screening is not reproduced)
  amplitude half the latitude difference Tr2 - Pk3
Trends, three ways:
  (a) Liu-style: Sen slope on the raw monthly series, plain Mann-Kendall p, and their
      "modified MK" (ported from 06-Modified-Mann-Kendall-Trend-Test.m)
  (b) Hamed-Rao (1998) as implemented in pymannkendall
  (c) this study: deseasonalised Theil-Sen, AR-sieve double-bootstrap calibrated p (rstats)
Periods: 1993-01..2023-01 (Liu: 361 months), 1993-2025, 2001-2025.
Usage: python liu_reproduce.py <liu_par dir> <g025|g0125> <out_dir> [--quick]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

SECT = {"Whole": (-150.0, -80.0), "Upper-Ridge": (-150.0, -135.13), "Lower-Ridge": (-135.13, -109.38),
        "Downstream": (-109.38, -80.0), "Flat": (-107.63, -96.38)}
PEAKS = {"Pk1": (-139.63, -138.13), "Pk2": (-129.38, -127.88), "Pk3": (-125.88, -124.88)}
TROUGHS = {"Tr1": (-135.63, -134.38), "Tr2": (-122.63, -121.63), "Tr3": (-116.13, -115.13)}
LIU_TABLE2 = {  # position (deg/dec), width (km/dec), speed (m/s/dec)
    "Whole": (0.03, 1.44, 0.01), "Upper-Ridge": (0.04, 0.02, 0.01), "Lower-Ridge": (0.09, 3.52, 0.02),
    "Downstream": (-0.03, 0.39, 0.01), "Flat": (-0.12, 0.34, 0.01)}
KM_PER_DEG = 111.195


def load(root, grid):
    files = sorted((Path(root) / grid).glob("*_*.nc"))
    C, G, S, N, keys = [], [], [], [], []
    for f in files:
        with xr.open_dataset(f) as d:
            C.append(d["count"].values.astype(float)); G.append(d.gradsum.values); S.append(d.speedsum.values)
            N.append(d.ndays.values.astype(float)); keys.append((int(d.attrs["year"]), int(d.attrs["month"])))
            lat, lon = d.latitude.values, d.longitude.values
    return np.array(C), np.array(G), np.array(S), np.array(N), keys, lat, lon


def metrics(C, G, S, N, keys, lat, lon):
    nt = len(keys)
    core_lat = (lat >= -60) & (lat <= -48)
    jl = np.nonzero((lon >= -150) & (lon <= -80))[0]
    axis = np.full((nt, jl.size), np.nan); width = axis.copy(); speed = axis.copy(); whalf = axis.copy()
    prev = None
    for k in range(1, nt - 2):                       # needs months k-1..k+2
        F = C[k - 1:k + 3].sum(0); Gs = G[k - 1:k + 3].sum(0)
        Sp = S[k - 1:k + 3].sum(0) / np.maximum(N[k - 1:k + 3].sum(0), 1)
        cur = np.full(jl.size, np.nan)
        for q, j in enumerate(jl):
            col = np.where(core_lat, F[:, j], -1)
            mx = col.max()
            if mx <= 0:
                continue
            cand = np.nonzero(col == mx)[0]
            if cand.size > 1:
                g = Gs[cand, j]; cand = cand[g == g.max()]
            if cand.size > 1 and prev is not None and np.isfinite(prev[q]):
                cand = cand[[np.argmin(np.abs(lat[cand] - prev[q]))]]
            i0 = cand[0]
            cur[q] = lat[i0]; speed[k, q] = Sp[i0, j]
            # zero-to-zero width
            north = np.nonzero(F[i0:, j] == 0)[0]; south = np.nonzero(F[:i0 + 1, j][::-1] == 0)[0]
            if north.size and south.size:
                width[k, q] = (lat[i0 + north[0]] - lat[i0 - south[0]]) * KM_PER_DEG
            # alternative (ours, not Liu's): contiguous band with F >= 0.5 F(axis), edges interpolated
            col_f = F[:, j].astype(float); half = 0.5 * col_f[i0]
            n_ = i0
            while n_ + 1 < lat.size and col_f[n_ + 1] >= half:
                n_ += 1
            s_ = i0
            while s_ - 1 >= 0 and col_f[s_ - 1] >= half:
                s_ -= 1
            def edge(a, b):  # interpolate crossing between index a (>=half) and b (<half)
                if b < 0 or b >= lat.size:
                    return lat[a]
                fa, fb = col_f[a], col_f[b]
                return lat[a] + (lat[b] - lat[a]) * (fa - half) / (fa - fb) if fa != fb else lat[a]
            whalf[k, q] = (edge(n_, n_ + 1) - edge(s_, s_ - 1)) * KM_PER_DEG
        axis[k] = cur; prev = cur
    t = pd.DataFrame(keys, columns=["year", "month"])
    out = t.copy()
    L = lon[jl]
    for s, (a, b) in SECT.items():
        m = (L >= a) & (L <= b)
        out[f"{s}_position"] = np.nanmean(axis[:, m], 1)
        out[f"{s}_width"] = np.nanmean(width[:, m], 1)
        out[f"{s}_speed"] = np.nanmean(speed[:, m], 1)
        out[f"{s}_halfwidth"] = np.nanmean(whalf[:, m], 1)
    for nme, (a, b) in PEAKS.items():
        m = (L >= a) & (L <= b); out[nme] = np.nanmin(axis[:, m], 1)
    for nme, (a, b) in TROUGHS.items():
        m = (L >= a) & (L <= b); out[nme] = np.nanmax(axis[:, m], 1)
    out["amplitude_Tr2_Pk3"] = 0.5 * (out.Tr2 - out.Pk3)
    return out, axis, width, speed, L, whalf


# ---------- Liu et al. statistics (ported from their MATLAB) ----------
def mk_liu(x):
    x = x[np.isfinite(x)]; n = len(x)
    S = sum(np.sign(x[k + 1:] - x[k]).sum() for k in range(n - 1))
    var = n * (n - 1) * (2 * n + 5) / 18
    Z = (S - 1) / np.sqrt(var) if S > 0 else ((S + 1) / np.sqrt(var) if S < 0 else 0)
    return float(2 * (1 - norm.cdf(abs(Z)))), S, n


def mmk_liu(x):
    """06-Modified-Mann-Kendall-Trend-Test.m: S from all n, variance from n_eff = n/tau."""
    x = x[np.isfinite(x)]
    p0, S, n = mk_liu(x)
    L = min(24, n // 4); xc = x - x.mean(); c0 = np.dot(xc, xc) / n
    r = np.array([np.dot(xc[:n - k], xc[k:]) / n / c0 for k in range(1, L + 1)])
    tau = 1 + 2 * np.sum((1 - np.arange(1, L + 1) / n) * r)
    ne = n / tau; var = ne * (ne - 1) * (2 * ne + 5) / 18
    Z = (S - 1) / np.sqrt(var) if S > 0 else ((S + 1) / np.sqrt(var) if S < 0 else 0)
    return float(2 * (1 - norm.cdf(abs(Z))))


def sen_liu(x):
    x = x[np.isfinite(x)]; i, j = np.triu_indices(len(x), 1)
    return float(np.median((x[j] - x[i]) / (j - i))) * 120      # per decade (monthly index)


def trends(M, quick=False, only=None):
    import pymannkendall as pmk
    rows = []
    cols = [c for c in M.columns if c not in ("year", "month")]
    x = (M.year + (M.month - 0.5) / 12).values
    periods = {"1993-2023(Liu)": (M.year >= 1993) & ((M.year < 2023) | ((M.year == 2023) & (M.month == 1))),
               "1993-2025": M.year >= 1993, "2001-2025": M.year >= 2001}
    if only:
        periods = {k: v for k, v in periods.items() if k in only}
    for c in cols:
        for pname, sel in periods.items():
            y = M[c].values[sel.values]; xx = x[sel.values]; mo = M.month.values[sel.values]
            ok = np.isfinite(y)
            if ok.sum() < 60:
                continue
            row = dict(series=c, period=pname, n=int(ok.sum()), sen_liu_per_dec=sen_liu(y),
                       p_mk_liu=mk_liu(y)[0], p_mmk_liu=mmk_liu(y),
                       p_hamed_rao=float(pmk.hamed_rao_modification_test(y[ok]).p))
            headline = any(c.endswith(k) for k in ("_position", "_width", "_speed", "_halfwidth")) or c.startswith("amplitude")
            if quick or not headline:
                o = R.trend(xx[ok], y[ok], month=mo[ok], n_boot=200).as_dict(); o["p_cal"] = np.nan
            else:
                o = R.trend_calibrated(xx[ok], y[ok], month=mo[ok], n_outer=100, n_boot=300)
            row.update(slope_deseas_per_dec=10 * o["slope"], ci_lo=10 * o.get("ci_cal_lo", o.get("ci_lo")),
                       ci_hi=10 * o.get("ci_cal_hi", o.get("ci_hi")), p_cal=o["p_cal"], p_sieve=o.get("p_boot", np.nan))
            rows.append(row); print(c, pname, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k not in ("series", "period")}, flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    root, grid, outdir = sys.argv[1], sys.argv[2], Path(sys.argv[3]); quick = "--quick" in sys.argv
    outdir.mkdir(parents=True, exist_ok=True)
    C, G, S, N, keys, lat, lon = load(root, grid)
    M, axis, width, speed, L, whalf = metrics(C, G, S, N, keys, lat, lon)
    M.to_csv(outdir / f"liu_metrics_{grid}.csv", index=False)
    xr.Dataset({"axis_lat": (("time", "lon"), axis), "width_km": (("time", "lon"), width), "speed": (("time", "lon"), speed),
                "halfwidth_km": (("time", "lon"), whalf)},
               coords=dict(time=np.arange(len(keys)), lon=L, year=("time", M.year.values), month=("time", M.month.values))
               ).to_netcdf(outdir / f"liu_axis_{grid}.nc")
    print("means:", M[[c for c in M.columns if c.startswith("Whole") or c.startswith("Lower")]].mean().round(3).to_dict())
    T = trends(M, quick, only=["1993-2023(Liu)", "1993-2025"])
    T.to_csv(outdir / f"liu_trends_{grid}.csv", index=False)
    print("Liu Table 2 (position deg/dec, width km/dec, speed m/s/dec):", LIU_TABLE2)
