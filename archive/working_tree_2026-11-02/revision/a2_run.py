"""A0 + A2: baseline re-detection and half-power validation on the stored
monthly KE fields (all-sat, 2001-2025).

Inputs : outputs/monthly_half-power_points/YYYY/MM.nc  (field `ke`, m2 s-2)
         outputs/mean_adt.nc                            (ACC mask, -0.6..0.2 m)
Outputs: revision/outputs/a2/*.csv, *.npz

Baseline (D1, fixed before any trend was computed): 200 cm2 s-2, 100 km, 25 %.
"""

from __future__ import annotations

import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fronts as F          # noqa: E402
import rstats as R          # noqa: E402

import os
REPO_OUT = Path(os.environ.get("OSR11_OUTPUTS", HERE.parent / "outputs"))
CACHE = Path(os.environ.get("OSR11_CACHE", HERE / "cache"))
OUT = HERE / "outputs" / "a2"
OUT.mkdir(parents=True, exist_ok=True)

BASE = dict(thr=0.02, min_w=100.0, frac=0.25)
SETS = {
    "baseline": BASE,
    "thr150": dict(BASE, thr=0.015), "thr250": dict(BASE, thr=0.025),
    "thr300": dict(BASE, thr=0.03), "thr400": dict(BASE, thr=0.04),
    "frac20": dict(BASE, frac=0.20), "frac33": dict(BASE, frac=0.33), "frac50": dict(BASE, frac=0.50),
    "minw50": dict(BASE, min_w=50.0), "minw200": dict(BASE, min_w=200.0), "minw300": dict(BASE, min_w=300.0),
    "relthr": dict(BASE, thr=0.20, relative=True),   # 20 % of each longitude's climatological max KE
}
MAX_ENV = 8
FI = {k: i for i, k in enumerate(F.FIELDS)}


def sector(lon):
    l = np.mod(lon, 360.0)
    s = np.full(l.shape, "Atlantic", dtype=object)
    s[(l >= 20) & (l < 147)] = "Indian"
    s[(l >= 147) & (l < 290)] = "Pacific"
    return s


def build_cache():
    """Stack the stored monthly KE/EKE fields (300 x 240 x 2880, float32)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    base = REPO_OUT / "monthly_half-power_points"
    ke = np.empty((300, 240, 2880), np.float32); eke = np.empty_like(ke)
    k = 0
    for y in range(2001, 2026):
        for m in range(1, 13):
            with xr.open_dataset(base / f"{y}" / f"{m:02d}.nc") as ds:
                ke[k] = ds.ke.values; eke[k] = ds.eke.values
                if k == 0:
                    np.save(CACHE / "lat.npy", ds.latitude.values); np.save(CACHE / "lon.npy", ds.longitude.values)
            k += 1
    np.save(CACHE / "ke.npy", ke); np.save(CACHE / "eke.npy", eke)


def load():
    if not (CACHE / "ke.npy").exists():
        build_cache()
    ke = np.load(CACHE / "ke.npy", mmap_mode="r")
    lat = np.load(CACHE / "lat.npy").astype(float)
    lon = np.load(CACHE / "lon.npy").astype(float)
    adt = xr.open_dataset(REPO_OUT / "mean_adt.nc")["mean_adt"].values
    acc = np.isfinite(adt) & (adt >= -0.6) & (adt <= 0.2)
    t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12),
                      "month": np.tile(np.arange(1, 13), 25)})
    t["x"] = t.year + (t.month - 0.5) / 12
    return ke, lat, lon, acc, t


def run_set(ke, lat, lon, p, clim_max=None):
    d = F.meridional_dists(lat)
    nt = ke.shape[0]
    counts = np.zeros((nt, lon.size), np.int16)
    tab = np.full((nt, lon.size, MAX_ENV, F.NF), np.nan, np.float32)
    for k in range(nt):
        fld = np.asarray(ke[k], float)
        if p.get("relative"):
            fld = fld / clim_max[None, :]
        c, tb = F.detect_field(fld, lat, d, p["thr"], p["min_w"], p["frac"], MAX_ENV)
        counts[k] = c
        tab[k] = tb
    return counts, tab


def clim_envelopes(kemean, lat, lon, p, clim_max):
    d = F.meridional_dists(lat)
    fld = kemean / clim_max[None, :] if p.get("relative") else kemean
    return F.detect_field(fld, lat, d, p["thr"], p["min_w"], p["frac"], MAX_ENV)


def in_acc(acc, lat, lon, lat_pts, j):
    i = np.clip(np.searchsorted(lat, lat_pts), 0, lat.size - 1)
    return acc[i, j]


def unmatched_series(counts, tab, lat, lon, acc, t):
    """Original-style monthly statistics over all ACC-masked detections."""
    rows = []
    secs = sector(lon)
    for k in range(counts.shape[0]):
        jj, kk = np.nonzero(np.arange(MAX_ENV)[None, :] < counts[k][:, None])
        if jj.size == 0:
            continue
        r = tab[k, jj, kk]
        keep = in_acc(acc, lat, lon, r[:, FI["lat_hp"]], jj)
        r, jj = r[keep], jj[keep]
        for sname in ["ACC", "Atlantic", "Indian", "Pacific"]:
            m = np.ones(jj.size, bool) if sname == "ACC" else (secs[jj] == sname)
            rr = r[m]
            rows.append(dict(i=k, region=sname, n=int(m.sum()),
                             median_lat=np.median(rr[:, FI["lat_hp"]]) if m.any() else np.nan,
                             seg_km=np.mean(rr[:, FI["seg_km"]]) if m.any() else np.nan,
                             fwhm_km=np.nanmean(rr[:, FI["fwhm_km"]]) if m.any() else np.nan,
                             peak=np.mean(rr[:, FI["peak_val"]]) if m.any() else np.nan))
    df = pd.DataFrame(rows)
    return df.merge(t.reset_index().rename(columns={"index": "i"}), on="i")


def matched_tracks(counts, tab, ccounts, ctab, lat, lon, acc, min_cov=0.9):
    """Assign monthly envelopes to climatological envelopes (same longitude)
    whose segment contains the monthly half-power latitude; keep the strongest
    if several. Returns dict of arrays (nt, ntracks) and track metadata."""
    nt = counts.shape[0]
    meta = []
    for j in range(lon.size):
        for c in range(ccounts[j]):
            s = int(ctab[j, c, FI["i_start"]]); e = int(ctab[j, c, FI["i_end"]])
            meta.append((j, c, lat[s], lat[e], ctab[j, c, FI["lat_hp"]]))
    meta = pd.DataFrame(meta, columns=["j", "c", "lat_lo", "lat_hi", "lat_hp_clim"])
    meta["acc"] = in_acc(acc, lat, lon, meta.lat_hp_clim.values, meta.j.values)
    fields = ["lat_hp", "peak_lat", "peak_val", "seg_km", "fwhm_km", "w25_km", "sig_km", "edge"]
    arr = {f: np.full((nt, len(meta)), np.nan, np.float32) for f in fields}
    jidx = meta.j.values
    lo = meta.lat_lo.values[:, None]
    hi = meta.lat_hi.values[:, None]
    fidx = [FI[f] for f in fields]
    rows = np.arange(len(meta))
    for k in range(nt):
        Rk = tab[k, jidx]                                  # (ntracks, MAX_ENV, NF)
        lh = Rk[:, :, FI["lat_hp"]]
        inside = np.isfinite(lh) & (lh >= lo) & (lh <= hi)
        pk = np.where(inside, Rk[:, :, FI["peak_val"]], -np.inf)
        best = np.argmax(pk, axis=1)
        ok = inside.any(axis=1)
        sel = Rk[rows, best][:, fidx]                      # (ntracks, nfields)
        for q, f in enumerate(fields):
            arr[f][k, ok] = sel[ok, q]
    cov = np.isfinite(arr["lat_hp"]).mean(axis=0)
    meta["coverage"] = cov
    meta["keep"] = meta.acc & (cov >= min_cov)
    return arr, meta


def track_series(arr, meta, keep, groups, t):
    """Average per-track anomalies (relative to each track's own mean) by group."""
    out = []
    for gname, gmask in groups.items():
        sel = keep & gmask
        if sel.sum() == 0:
            continue
        rec = {"region": gname, "ntracks": int(sel.sum())}
        for f in ["lat_hp", "peak_lat", "peak_val", "seg_km", "fwhm_km", "w25_km", "sig_km"]:
            a = arr[f][:, sel].astype(float)
            anom = a - np.nanmean(a, axis=0, keepdims=True)
            rec[f] = np.nanmean(anom, axis=1)
            rec[f + "_mean"] = float(np.nanmean(a))
        rec["edge_frac"] = float(np.nanmean(arr["edge"][:, sel]))
        rec["coverage"] = np.isfinite(arr["lat_hp"][:, sel]).mean(axis=1)
        out.append(rec)
    return out


def tr(t, y, nb=600, seed=0):
    T = R.trend(t.x.values, np.asarray(y, float), month=t.month.values, n_boot=nb, seed=seed)
    return T


def main():
    ke, lat, lon, acc, t = load()
    kemean = np.nanmean(np.asarray(ke, np.float64), axis=0)
    clim_max = np.nanmax(kemean, axis=0)
    secs = sector(lon)
    summary, lonrows = [], []
    for name, p in SETS.items():
        counts, tab = run_set(ke, lat, lon, p, clim_max)
        ccounts, ctab = clim_envelopes(kemean, lat, lon, p, clim_max)
        um = unmatched_series(counts, tab, lat, lon, acc, t)
        arr, meta = matched_tracks(counts, tab, ccounts, ctab, lat, lon, acc)
        keep = meta.keep.values
        jsec = secs[meta.j.values]
        groups = {"ACC": np.ones(len(meta), bool), "Atlantic": jsec == "Atlantic",
                  "Indian": jsec == "Indian", "Pacific": jsec == "Pacific"}
        ts = track_series(arr, meta, keep, groups, t)
        np.savez_compressed(OUT / f"tracks_{name}.npz", **{f: arr[f] for f in arr},
                            j=meta.j.values, c=meta.c.values, keep=keep, coverage=meta.coverage.values,
                            lat_hp_clim=meta.lat_hp_clim.values, acc=meta.acc.values)
        um.to_csv(OUT / f"unmatched_{name}.csv", index=False)
        for rec in ts:
            reg = rec["region"]
            u = um[um.region == reg].sort_values("i")
            row = dict(set=name, region=reg, ntracks=rec["ntracks"],
                       n_tracks_all_acc=int(meta.acc.sum()),
                       edge_frac=rec["edge_frac"],
                       mean_cov=float(np.mean(rec["coverage"])))
            for f, unit_scale in [("lat_hp", 1.0), ("peak_lat", 1.0), ("seg_km", 1.0), ("fwhm_km", 1.0),
                                  ("w25_km", 1.0), ("sig_km", 1.0), ("peak_val", 1e4)]:
                T = tr(t, rec[f] * unit_scale, nb=600 if name == "baseline" else 300)
                row.update({f"m_{f}_slope": T.slope, f"m_{f}_lo": T.ci_lo, f"m_{f}_hi": T.ci_hi,
                            f"m_{f}_p": T.p_boot, f"m_{f}_mean": rec[f + "_mean"] * unit_scale})
            for f, sc in [("median_lat", 1.0), ("seg_km", 1.0), ("fwhm_km", 1.0), ("n", 1.0), ("peak", 1e4)]:
                T = tr(t, u[f].values * sc, nb=300)
                row.update({f"u_{f}_slope": T.slope, f"u_{f}_lo": T.ci_lo, f"u_{f}_hi": T.ci_hi,
                            f"u_{f}_p": T.p_boot})
            summary.append(row)
        if name == "baseline":
            # A2f: longitude-resolved (5 deg bins) matched half-power latitude trends
            bins = np.arange(-180, 181, 5)
            lonq = lon[meta.j.values]
            for b0 in bins[:-1]:
                sel = keep & (lonq >= b0) & (lonq < b0 + 5)
                if sel.sum() < 5:
                    continue
                a = arr["lat_hp"][:, sel].astype(float)
                y = np.nanmean(a - np.nanmean(a, 0, keepdims=True), axis=1)
                T = tr(t, y, nb=400)
                pk = arr["peak_lat"][:, sel].astype(float)
                yp = np.nanmean(pk - np.nanmean(pk, 0, keepdims=True), axis=1)
                Tp = tr(t, yp, nb=200)
                lonrows.append(dict(lon0=b0, ntracks=int(sel.sum()), lat_clim=float(np.nanmean(a)),
                                    slope_deg_yr=T.slope, lo=T.ci_lo, hi=T.ci_hi, p=T.p_boot,
                                    p_hamed_rao=T.p_hamed_rao, peak_slope=Tp.slope, peak_p=Tp.p_boot,
                                    std_hp=float(np.nanstd(y)), std_peak=float(np.nanstd(yp))))
            # detection frequency (fraction of months with >=1 ACC detection) per longitude
            freq = np.zeros(lon.size)
            for k in range(counts.shape[0]):
                for j in np.nonzero(counts[k])[0]:
                    r = tab[k, j, :counts[k, j]]
                    if in_acc(acc, lat, lon, r[:, FI["lat_hp"]], np.full(r.shape[0], j)).any():
                        freq[j] += 1
            pd.DataFrame({"lon": lon, "freq": freq / counts.shape[0]}).to_csv(OUT / "detection_frequency_lon.csv", index=False)
            # 2-D density of half-power points (0.5 deg bins)
            H = np.zeros((60, 720))
            for k in range(counts.shape[0]):
                jj, kk = np.nonzero(np.arange(MAX_ENV)[None, :] < counts[k][:, None])
                la = tab[k, jj, kk, FI["lat_hp"]]
                keepa = in_acc(acc, lat, lon, la, jj)
                h, _, _ = np.histogram2d(la[keepa], lon[jj][keepa], bins=[np.arange(-65, -34.9, 0.5), np.arange(-180, 180.01, 0.5)])
                H += h
            np.savez_compressed(OUT / "hp_density.npz", H=H / counts.shape[0])
        print("done", name, flush=True)
        pd.DataFrame(summary).to_csv(OUT / "a2_summary.csv", index=False)
    pd.DataFrame(lonrows).to_csv(OUT / "a2f_lon_trends.csv", index=False)


if __name__ == "__main__":
    main()
