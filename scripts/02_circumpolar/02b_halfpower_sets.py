"""Half-power envelope detection for the baseline and all sensitivity parameter sets (former a2_run.py, main part).

Purpose : detect half-power envelopes in every monthly KE field, match them to the climatological envelopes
          (900 ACC longitudes for the baseline), and compute the trend of the matched half-power latitude per
          region, per 5-degree longitude bin, plus detection frequency and the half-power density.
Inputs  : $OSR11_DATA/monthly_half-power_points/YYYY/MM.nc (01b) via the cache $OSR11_CACHE; outputs/base/mean_adt.nc
Outputs : outputs/a2/a2_summary.csv, a2f_lon_trends.csv, tracks_<set>.npz, unmatched_<set>.csv,
          detection_frequency_lon.csv, hp_density.npz
Run     : python scripts/02_circumpolar/02b_halfpower_sets.py        (needs raw data; several hours)
Seeds   : trend bootstraps use seed=0 (rstats.trend).


"""
from __future__ import annotations
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np
import pandas as pd

from osr11.halfpower import (OUT, SETS, sector, load, run_set, clim_envelopes, unmatched_series,  # noqa: F401
                             matched_tracks, track_series, tr, in_acc, MAX_ENV, FI)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
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
