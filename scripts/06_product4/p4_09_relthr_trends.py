"""Step 9: trends of the mean half-power latitude and of the fixed-contour latitude, circumpolar, by sector and in 86-40 W, for
Product 1 (relative threshold, x7 tracks) and Product 4 (p4_08 tracks), with the statistic of x7_relthr_products.py; plus the two
validations of the Product 4 pipeline.

Purpose : position trends of panel f of Fig. S3 (km per decade, positive northward).
Inputs  : outputs/x7/x7_tracks.npz (Product 1, 1826 retained relative-threshold tracks: lat_hp, lat_c, lon);
          outputs/x7/x7_trends.csv (Product 1 reference values); outputs/p4/p4_08_tracks_<thr>_<ice>.npz (Product 4, p4_08);
          outputs/a2/a2h_twosat_vs_allsat.csv (Product 4 fixed-threshold reference, 02e)
Outputs : outputs/p4/p4_09_relthr_trends.csv  product, set, region, series, slope, lo, hi, p, n_tracks (km per decade, north +)
          outputs/p4/p4_09_validation.csv     check, here, reference, abs diff (slope, lo, hi, p, n)
Method  : x7 'series': for the tracks of a region, each track's monthly latitude minus its own 2001-2025 mean, averaged over the
          tracks (nanmean) month by month, x 111.195 km per degree; trend = osr11.rstats.trend_bcz(t, y, month=month, monthly=True)
          with its defaults (n_outer = 300, n_bias = 100, seed 0; p = max of AIC-order and AR(12) p, interval = union), x 10 per
          decade - the x7 call; 02e calls the same function with n_outer = 300 explicitly, so both references use one statistic.
          Regions (02e convention, lon in -180..180): Atlantic -70 <= lon < 20, Indian 20 <= lon < 147, Pacific otherwise
          (147 E-70 W); 86-40 W: -86 <= lon <= -40 (x7). Contour series use every retained track (x7 does not apply the
          >= 200-valid-month filter of x6b; the two give the same circumpolar trend, -16.24 km per decade).
          Validation 1 (statistics): the Product 1 circumpolar and 86-40 W rows recomputed here from x7_tracks.npz must equal
          x7_trends.csv. Validation 2 (tracking): p4_08 run with the fixed threshold and the 01c permanent ice mask (the 02e
          inputs) must reproduce the two-sat rows of a2h (track count and position trends).
Run     : python scripts/06_product4/p4_09_relthr_trends.py
Seeds   : rstats.trend_bcz default seed (0).
"""
from multiprocessing import Pool

import numpy as np
import pandas as pd

import p4_common as PC
from osr11 import rstats as R

KM = 111.195
t = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25)


def regions(lon):
    return {"circumpolar": np.ones(lon.size, bool),
            "Pacific": (lon >= 147) | (lon < -70),
            "Atlantic": (lon >= -70) & (lon < 20),
            "Indian": (lon >= 20) & (lon < 147),
            "86-40W": (lon >= -86) & (lon <= -40)}


def series(a):
    return np.nanmean(a - np.nanmean(a, 0), 1)


def job(args):
    prod, setname, reg, key, y, n = args
    r = R.trend_bcz(t, y, month=month, monthly=True)
    return dict(product=prod, set=setname, region=reg, series=key, slope=10 * r["slope"], lo=10 * r["ci_lo"],
                hi=10 * r["ci_hi"], p=r["p"], n_tracks=n)


def jobs_for(prod, setname, lat_hp, lat_c, lon, keys=("lat_hp", "lat_contour")):
    out = []
    for reg, m in regions(lon).items():
        for key, arr in (("lat_hp", lat_hp), ("lat_contour", lat_c)):
            if key in keys:
                out.append((prod, setname, reg, key, series(arr[:, m].astype(float)) * KM, int(m.sum())))
    return out


if __name__ == "__main__":
    X = np.load(PC.OUT / "x7" / "x7_tracks.npz")
    J = jobs_for("P1", "relthr", X["lat_hp"], X["lat_c"], X["lon"])
    for tag in ("relthr_none", "relthr_perm", "baseline_perm"):
        Z = np.load(PC.HERE / f"p4_08_tracks_{tag}.npz"); k = Z["keep"]
        J += jobs_for("P4", tag, Z["lat_hp"][:, k], Z["lat_c"][:, k], Z["lon"][k],
                      keys=("lat_hp",) if tag == "baseline_perm" else ("lat_hp", "lat_contour"))
    with Pool(2) as pool:
        rows = pool.map(job, J)
    df = pd.DataFrame(rows); df.to_csv(PC.HERE / "p4_09_relthr_trends.csv", index=False)
    print(df.round(3).to_string())
    # validations
    V = []
    ref = pd.read_csv(PC.OUT / "x7" / "x7_trends.csv")
    for reg_ref, reg in (("circumpolar", "circumpolar"), ("Drake region 86-40W", "86-40W")):
        for key in ("lat_hp", "lat_contour"):
            a = df[(df["product"] == "P1") & (df.region == reg) & (df.series == key)].iloc[0]
            b = ref[(ref.region == reg_ref) & (ref.series == key)].iloc[0]
            for c in ("slope", "lo", "hi", "p", "n_tracks"):
                V.append(dict(check=f"P1 relthr {reg} {key} {c} vs x7_trends.csv", here=a[c], reference=b[c], abs_diff=abs(a[c] - b[c])))
    h = pd.read_csv(PC.OUT / "a2" / "a2h_twosat_vs_allsat.csv"); h = h[h["product"] == "two-sat"]
    for reg_ref, reg in (("ACC", "circumpolar"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian")):
        a = df[(df["set"] == "baseline_perm") & (df.region == reg) & (df.series == "lat_hp")].iloc[0]
        b = h[h.region == reg_ref].iloc[0]
        for c, cr in (("slope", "pos_km_dec"), ("lo", "pos_lo"), ("hi", "pos_hi"), ("p", "pos_p"), ("n_tracks", "n_lon")):
            V.append(dict(check=f"P4 fixed-threshold {reg} lat_hp {c} vs a2h two-sat", here=a[c], reference=b[cr], abs_diff=abs(a[c] - b[cr])))
    V = pd.DataFrame(V); V.to_csv(PC.HERE / "p4_09_validation.csv", index=False)
    print(V.to_string()); print("max abs diff", V.abs_diff.max())
