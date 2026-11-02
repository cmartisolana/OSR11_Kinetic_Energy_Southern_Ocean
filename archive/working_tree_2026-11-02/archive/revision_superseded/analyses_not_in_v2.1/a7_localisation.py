"""A7: is the ACC energy increase concentrated downstream of the standing meanders?

Question (eddy-saturation localisation, Thompson & Naveira Garabato 2014): if the meanders are where
the mean flow's extra energy goes into eddies, the EKE increase should be proportionally larger
downstream of the major topographic features than over the rest of the ACC band.

Specification (written before computing any of the statistics below; 2026-10-08):
  Regions: the A5d boxes of a5_energy.py (fixed from bathymetry and climatological EKE before any
    box trend was computed), inside the ACC band of each product (mean ADT -0.6..+0.2 m, pixels
    with ice fraction > 0.10 in any month excluded).
      WAKE = union of the four downstream boxes (Kerguelen 75-95E, SEIR-Macquarie 145-175E,
             PAR 150-125W, Drake/Scotia 55-35W)
      REST = ACC band minus WAKE
    Caveat: the A5d box trends were already seen (descriptively) before this test was designed.
  Statistic: relative-trend difference. For region R, n_R(t) = X_R(t) / mean(X_R). The deseasonalised
    Theil-Sen trend of d(t) = n_WAKE(t) - n_REST(t), x 1000, is the difference in % per decade.
    Null: the same proportional change everywhere (d has no trend). Calibrated AR-sieve test
    (rstats.trend_calibrated), as for all other series.
  Primary: X = EKE, WAKE vs REST, 2001-2025, each product separately (one test per product).
  Secondary (BH-FDR, alpha 0.10, within product): X = KE, along-stream speed SPD, sub-monthly SUB,
    low-frequency LF for WAKE vs REST; EKE and KE for each feature, downstream vs upstream box and
    downstream box vs REST.
  Sensitivity: 2001-2017 (before the main divergence of the two products).
  Descriptive: share of the ACC-band area-integrated EKE trend located in WAKE, compared with WAKE's
    share of area and of mean EKE.
Usage: python a7_localisation.py <data_root> <product> <a5 dir> <out_dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

root, product, a5dir, outdir = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
outdir.mkdir(parents=True, exist_ok=True)
FEAT = ["Kerguelen", "SEIR_Macquarie", "PAR", "Drake"]
BOXES = {"Kerguelen": ((50, 60), (75, 95)), "SEIR_Macquarie": ((110, 130), (145, 175)),
         "PAR": ((-175, -160), (-150, -125)), "Drake": ((-100, -80), (-55, -35))}   # = a5_energy.py
VARS = ["EKE", "KE", "SPD", "SUB", "LF"]

# --- mask areas (same masks as a5_energy.py) ---
files = [root / f"monthly_{product}" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
d0 = xr.open_dataset(files[0]); lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float)
shape = (lat.size, lon.size)
nvalid = np.zeros(shape); sa = np.zeros(shape); imax = np.zeros(shape)
for f in files:
    with xr.open_dataset(f) as ds:
        u, v = ds.u_mean.values, ds.v_mean.values
        nvalid += np.isfinite(u) & np.isfinite(v); sa += np.nan_to_num(ds.adt_mean.values)
        imax = np.fmax(imax, np.nan_to_num(ds.ice_frac.values))
valid = (nvalid == len(files)) & ~(imax > 0.10)
adt = sa / len(files)
acc = valid & (adt >= -0.6) & (adt <= 0.2)
w = np.broadcast_to(np.cos(np.radians(lat))[:, None], shape)
L2 = np.broadcast_to(lon[None, :], shape)
A = {"ACC": float(w[acc].sum())}
for f_, ((u0, u1), (d0_, d1_)) in BOXES.items():
    A[f"up_{f_}"] = float(w[acc & (L2 >= u0) & (L2 < u1)].sum())
    A[f"down_{f_}"] = float(w[acc & (L2 >= d0_) & (L2 < d1_)].sum())
A["WAKE"] = sum(A[f"down_{f_}"] for f_ in FEAT); A["REST"] = A["ACC"] - A["WAKE"]
pd.Series(A).to_csv(outdir / f"a7_areas_{product}.csv", header=["area_coslat_units"])
print({k: round(v / A["ACC"], 3) for k, v in A.items()}, flush=True)

# --- region series from the a5 area means (area-weighted algebra) ---
S = pd.read_csv(a5dir / f"a5_series_{product}.csv")
for v in VARS:
    S[f"WAKE_{v}"] = sum(A[f"down_{f_}"] * S[f"down_{f_}_{v}"] for f_ in FEAT) / A["WAKE"]
    S[f"REST_{v}"] = (A["ACC"] * S[f"ACC_{v}"] - A["WAKE"] * S[f"WAKE_{v}"]) / A["REST"]
S.to_csv(outdir / f"a7_series_{product}.csv", index=False)

tests = [("primary", "EKE", "WAKE", "REST")]
tests += [("secondary", v, "WAKE", "REST") for v in VARS[1:]]
for f_ in FEAT:
    for v in ("EKE", "KE"):
        tests += [("secondary", v, f"down_{f_}", f"up_{f_}"), ("secondary", v, f"down_{f_}", "REST")]
rows = []
for per, sel in (("2001-2025", np.ones(len(S), bool)), ("2001-2017", S.year.values <= 2017)):
    for kind, v, a, b in tests:
        if per == "2001-2017" and not (kind == "primary" or (a == "WAKE")):
            continue
        ya, yb = S[f"{a}_{v}"].values[sel], S[f"{b}_{v}"].values[sel]
        d = ya / ya.mean() - yb / yb.mean()
        o = R.trend_calibrated(S.x.values[sel], d, month=S.month.values[sel], n_outer=100, n_boot=300)
        ta = R.trend(S.x.values[sel], ya, month=S.month.values[sel], n_boot=200).as_dict()
        tb = R.trend(S.x.values[sel], yb, month=S.month.values[sel], n_boot=200).as_dict()
        rows.append(dict(product=product, period=per, kind=kind, var=v, region_a=a, region_b=b,
                         pct_dec_a=1000 * ta["slope"] / ya.mean(), pct_dec_b=1000 * tb["slope"] / yb.mean(),
                         diff_pct_dec=1000 * o["slope"], diff_lo=1000 * o["ci_cal_lo"], diff_hi=1000 * o["ci_cal_hi"],
                         p_cal=o["p_cal"]))
        print(per, v, a, "vs", b, {k: round(rows[-1][k], 3) for k in ("pct_dec_a", "pct_dec_b", "diff_pct_dec", "p_cal")}, flush=True)
T = pd.DataFrame(rows)
T["fdr_sig"] = False
for per in T.period.unique():
    m = (T.kind == "secondary") & (T.period == per)
    T.loc[m, "fdr_sig"] = R.bh_fdr(T.loc[m, "p_cal"].values, 0.10)
T.to_csv(outdir / f"a7_tests_{product}.csv", index=False)

# descriptive: share of the ACC-integrated EKE trend in WAKE
x, mo = S.x.values, S.month.values
bw = R.trend(x, S["WAKE_EKE"].values, month=mo, n_boot=200).slope * A["WAKE"]
bt = R.trend(x, S["ACC_EKE"].values, month=mo, n_boot=200).slope * A["ACC"]
share = dict(product=product, area_share=A["WAKE"] / A["ACC"],
             mean_eke_share=S["WAKE_EKE"].mean() * A["WAKE"] / (S["ACC_EKE"].mean() * A["ACC"]),
             trend_share=bw / bt)
pd.DataFrame([share]).to_csv(outdir / f"a7_share_{product}.csv", index=False)
print(share)
