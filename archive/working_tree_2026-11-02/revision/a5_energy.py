"""A5a/A5b/A5c and A4: what the KE/EKE trends can and cannot show, in both products.

Reads revision/data/monthly_<product>/YYYY/MM.nc (from extract_monthly.py).
For each month m and pixel, with ubar = 2001-2025 mean velocity of that product:
  KE_m      = 0.5 |U_m|^2                          (KE of the monthly-mean flow, as in v2)
  EKE_m     = 0.5 |U_m - ubar|^2 + 0.5 var_intra   (identical to v2's EKE, exact)
  LF_m      = 0.5 |U_m - ubar|^2                   (low-frequency anomaly energy)
  SUB_m     = 0.5 var_intra                        (sub-monthly variance energy)
  PROJ_m    = U_m . ubar                           (projection on the mean-flow direction)
  SPD_m     = U_m . ubar / |ubar|                  (along-stream speed)
  EKE3_m    = 0.5 |U_m - U_3yr(m)|^2 + 0.5 var_intra   (EKE w.r.t. centred 36-month running mean)
Identity: EKE_m - KE_m = 0.5|ubar|^2 - PROJ_m + SUB_m.
Area means (cos-lat weights) over the ACC band (mean ADT -0.6..+0.2 m of the
same product), the three sectors inside the band, and the 35-65S box.
Also: KE of annual- and 3-year-block-mean velocities.
Usage: python a5_energy.py <data_root> <product> <out_dir>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

root, product, outdir = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
outdir.mkdir(parents=True, exist_ok=True)
files = [root / f"monthly_{product}" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
missing = [f for f in files if not f.exists()]
if missing:
    sys.exit(f"{len(missing)} monthly files missing, e.g. {missing[0]}")

d0 = xr.open_dataset(files[0])
lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float)
nt = len(files)


def read(k):
    with xr.open_dataset(files[k]) as ds:
        return (ds.u_mean.values.astype(np.float64), ds.v_mean.values.astype(np.float64),
                ds.var_intra.values.astype(np.float64), ds.adt_mean.values.astype(np.float64),
                ds.ice_frac.values.astype(np.float64))


# pass 1: means (streaming; memory-light so it also runs on small machines)
shape = (lat.size, lon.size)
su = np.zeros(shape); sv = np.zeros(shape); sa = np.zeros(shape); imax = np.zeros(shape)
nvalid = np.zeros(shape)
for k in range(nt):
    u, v, s2, a, ice = read(k)
    ok = np.isfinite(u) & np.isfinite(v)
    su += np.where(ok, u, 0); sv += np.where(ok, v, 0); nvalid += ok
    sa += np.nan_to_num(a); imax = np.fmax(imax, np.nan_to_num(ice))
ice_any = imax > 0.10                      # documented choice: >10 % of days in any month
valid = (nvalid == nt) & ~ice_any
ub = su / np.maximum(nvalid, 1); vb = sv / np.maximum(nvalid, 1)
adt_mean = sa / nt
acc = valid & (adt_mean >= -0.6) & (adt_mean <= 0.2)
lon360 = np.mod(lon, 360)
secname = np.where((lon360 >= 20) & (lon360 < 147), "Indian", np.where((lon360 >= 147) & (lon360 < 290), "Pacific", "Atlantic"))
S = np.broadcast_to(secname[None, :], acc.shape)
masks = {"ACC": acc, "Atlantic": acc & (S == "Atlantic"), "Indian": acc & (S == "Indian"),
         "Pacific": acc & (S == "Pacific"), "box35_65": valid}
# A5d boxes inside the ACC band, fixed from bathymetry (GEBCO 2026, 1/8 deg block means) and the
# CLIMATOLOGICAL EKE before any box trend was computed (2026-10-08). Longitudes in deg E (-180..180).
#   feature                          upstream (flatter, deep)   downstream (energetic wake)
#   Kerguelen Plateau                50..60                     75..95
#   SE Indian R.-Macquarie-Campbell  110..130                   145..175
#   Pacific-Antarctic Ridge          -175..-160                 -150..-125
#   Drake Passage / Shackleton       -100..-80                  -55..-35
BOXES = {"Kerguelen": ((50, 60), (75, 95)), "SEIR_Macquarie": ((110, 130), (145, 175)),
         "PAR": ((-175, -160), (-150, -125)), "Drake": ((-100, -80), (-55, -35))}
L2 = np.broadcast_to(lon[None, :], acc.shape)
for feat, ((u0, u1), (d0_, d1_)) in BOXES.items():
    masks[f"up_{feat}"] = acc & (L2 >= u0) & (L2 < u1)
    masks[f"down_{feat}"] = acc & (L2 >= d0_) & (L2 < d1_)
w = np.broadcast_to(np.cos(np.radians(lat))[:, None], acc.shape)
spd_b = np.hypot(ub, vb)
HALF = 18
rows = []
win_u = np.zeros(shape); win_v = np.zeros(shape); a_, b_ = 0, 0
blk_sum = {12: [np.zeros(shape), np.zeros(shape)], 36: [np.zeros(shape), np.zeros(shape)]}
blk = []
for k in range(nt):
    a_new, b_new = max(0, k - HALF), min(nt, k + HALF)
    while b_ < b_new:
        uu, vv = read(b_)[:2]; win_u += np.nan_to_num(uu); win_v += np.nan_to_num(vv); b_ += 1
    while a_ < a_new:
        uu, vv = read(a_)[:2]; win_u -= np.nan_to_num(uu); win_v -= np.nan_to_num(vv); a_ += 1
    U3k = win_u / (b_ - a_); V3k = win_v / (b_ - a_)
    Uk, Vk, S2k = read(k)[:3]
    ke = 0.5 * (Uk ** 2 + Vk ** 2)
    lf = 0.5 * ((Uk - ub) ** 2 + (Vk - vb) ** 2)
    sub = 0.5 * S2k
    proj = Uk * ub + Vk * vb
    spd = proj / np.where(spd_b > 0, spd_b, np.nan)
    eke3 = 0.5 * ((Uk - U3k) ** 2 + (Vk - V3k) ** 2) + sub
    rec = {"i": k}
    for name, m in masks.items():
        ww = w[m]
        for vn, A_ in [("KE", ke), ("EKE", lf + sub), ("LF", lf), ("SUB", sub), ("PROJ", proj),
                       ("SPD", spd), ("EKE3", eke3)]:
            aa = A_[m]; ok = np.isfinite(aa)
            rec[f"{name}_{vn}"] = float(np.sum(aa[ok] * ww[ok]) / np.sum(ww[ok]))
    rows.append(rec)
    # block-mean-flow KE (annual, 3-yr)
    for L in (12, 36):
        blk_sum[L][0] += np.nan_to_num(Uk); blk_sum[L][1] += np.nan_to_num(Vk)
        if (k + 1) % L == 0:
            Ub, Vb = blk_sum[L][0] / L, blk_sum[L][1] / L
            keb = 0.5 * (Ub ** 2 + Vb ** 2) * 1e4
            r = {"block": "annual" if L == 12 else "3yr", "x": 2001 + (k + 1 - L / 2) / 12}
            for name, m in masks.items():
                aa = keb[m]; ww = w[m]
                r[f"{name}_KEblock"] = float(np.sum(aa * ww) / np.sum(ww))
            blk.append(r)
            blk_sum[L] = [np.zeros(shape), np.zeros(shape)]
ts = pd.DataFrame(rows)
ts["year"] = np.repeat(np.arange(2001, 2026), 12); ts["month"] = np.tile(np.arange(1, 13), 25)
ts["x"] = ts.year + (ts.month - 0.5) / 12
# convert energies to cm2 s-2, speeds to cm s-1
for c in ts.columns:
    if any(c.endswith(s) for s in ("_KE", "_EKE", "_LF", "_SUB", "_PROJ", "_EKE3")):
        ts[c] *= 1e4
    if c.endswith("_SPD"):
        ts[c] *= 1e2
ts.to_csv(outdir / f"a5_series_{product}.csv", index=False)

pd.DataFrame(blk).to_csv(outdir / f"a5_blockKE_{product}.csv", index=False)

# trends
out = []
for c in [c for c in ts.columns if "_" in c and c.split("_")[-1] in ("KE", "EKE", "LF", "SUB", "PROJ", "SPD", "EKE3")]:
    for per, sel in [("2001-2025", np.ones(nt, bool)), ("2001-2020", ts.year.values <= 2020)]:
        if per == "2001-2020" and not c.startswith("box35_65"):
            continue
        y = ts[c].values[sel]
        reg, var = c.rsplit("_", 1)
        if reg in ("ACC", "Atlantic", "Indian", "Pacific", "box35_65") and var in ("KE", "EKE", "SPD", "EKE3", "LF", "SUB"):
            o = R.trend_calibrated(ts.x.values[sel], y, month=ts.month.values[sel], n_outer=100, n_boot=300)
        else:   # boxes and remaining terms: sieve only (no double bootstrap) to keep run time reasonable
            o = R.trend(ts.x.values[sel], y, month=ts.month.values[sel], n_boot=1000).as_dict()
            o["p_cal"] = np.nan
        o.update(series=c, period=per, mean=float(np.mean(y)), pct_per_decade=1000 * o["slope"] / float(np.mean(y)))
        out.append(o)
        print(product, c, per, round(o["slope"], 4), round(o["p_cal"], 3), flush=True)
B = pd.DataFrame(blk)
for lab in ["annual", "3yr"]:
    bb = B[B.block == lab]
    for c in [c for c in bb.columns if c.endswith("KEblock")]:
        if lab == "3yr":
            # 8 blocks: an AR sieve cannot be fitted; OLS t-test on independent blocks (as a2e, Chambers 2018)
            from scipy import stats
            r = stats.linregress(bb.x.values, bb[c].values); tq = stats.t.ppf(0.975, len(bb) - 2)
            slope, lo, hi, p = r.slope, r.slope - tq * r.stderr, r.slope + tq * r.stderr, r.pvalue
        else:
            T = R.trend(bb.x.values, bb[c].values, monthly=False, n_boot=1000)
            slope, lo, hi, p = T.slope, T.ci_lo, T.ci_hi, T.p_boot
        out.append(dict(series=f"{c}_{lab}", period="2001-2025", slope=slope, ci_lo=lo, ci_hi=hi,
                        p_boot=p, mean=float(bb[c].mean()), pct_per_decade=1000 * slope / float(bb[c].mean())))
pd.DataFrame(out).to_csv(outdir / f"a5_trends_{product}.csv", index=False)
print("done", product)
