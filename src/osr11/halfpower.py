"""Half-power envelope tracking on the stored monthly KE fields (library part of the former a2_run.py).

Used by 02b_halfpower_sets.py, 02c_halfpower_tracks.py, 02d_halfpower_detectability.py, 02e_halfpower_twosat.py,
02f_eke_maps.py, 02g_eke_maps_bc.py and 05a_retest.py.

Inputs : $OSR11_DATA/monthly_half-power_points/YYYY/MM.nc  (field `ke`, m2 s-2; 01b_halfpower_fields.py)
         outputs/base/mean_adt.nc                          (ACC mask, -0.6..0.2 m)
Cache  : $OSR11_CACHE/{ke,eke,lat,lon}.npy (1.6 GB, built on first use by build_cache())
Outputs of the callers go to outputs/a2/ (OUT below).

Baseline (D1, fixed before any trend was computed): 200 cm2 s-2, 100 km, 25 %.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from osr11 import fronts as F          # noqa: E402
from osr11 import rstats as R          # noqa: E402
from osr11 import paths as P           # noqa: E402

CACHE = P.CACHE
OUT = P.OUT / "a2"                     # created by the scripts that write into it

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
    base = P.HP_POINTS
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
    adt = xr.open_dataset(P.MEAN_ADT)["mean_adt"].values
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




def coarsen(cube):
    """Block-average (nt, ny, nx) 1/8-deg fields to 1 deg (8 x 8 pixels); cells with < 50 % valid pixels -> NaN
    (moved here from the former a1_maps.py so that 02g and 02h share it)."""
    nt, ny, nx = cube.shape
    c = cube.reshape(nt, ny // 8, 8, nx // 8, 8)
    valid = np.isfinite(c).sum(axis=(2, 4))
    m = np.nanmean(c, axis=(2, 4))
    m[valid < 32] = np.nan
    return m
