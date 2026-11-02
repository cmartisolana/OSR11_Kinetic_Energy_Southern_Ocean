"""A6b: sensitivity of the Drake-box KE/EKE trends to the pixel mask.

Box 80-52W, 63-54S, cos-lat weights, fixed 2001-2025 reference velocity.
Masks:
  nomask   : every finite pixel in each month (as a6_drake.py)
  ice_any  : drop pixels with monthly ice_frac > 0.10 in any month, and pixels
             not valid in all 300 months (same rule as a5_energy.py)
The stored v2 cube (a1_recompute.py, Drake_box) uses v2's own mask, which removes
more pixels; its numbers are in outputs/a1_ke_eke_trends.csv.

Purpose : Drake-box KE/EKE monthly series with and without the ice mask (inputs of 05a CONTRAST and of 05g real_Drake_EKE).
Inputs  : $OSR11_DATA/monthly_allsat (00a)
Outputs : outputs/a6/a6b_drake_box_series.csv, a6b_drake_box_trends.csv (distributed)
Run     : python scripts/03_drake/03u_drake_box_masks.py
Seeds   : rstats defaults (seed=0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

from osr11 import rstats as R, paths as PA

data, outdir = PA.DATA, PA.out("a6")
outdir.mkdir(parents=True, exist_ok=True)
files = [data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
t["x"] = t.year + (t.month - 0.5) / 12
U, V, S2, ICE = [], [], [], []
for f in files:
    with xr.open_dataset(f) as ds:
        s = ds.sel(longitude=slice(-80, -52), latitude=slice(-63, -54))
        U.append(s.u_mean.values); V.append(s.v_mean.values); S2.append(s.var_intra.values); ICE.append(s.ice_frac.values)
lat = s.latitude.values
U, V, S2, ICE = (np.array(a, dtype=np.float64) for a in (U, V, S2, ICE))
ub, vb = np.nanmean(U, 0), np.nanmean(V, 0)
KE = 0.5 * (U ** 2 + V ** 2); EKE = 0.5 * ((U - ub) ** 2 + (V - vb) ** 2) + 0.5 * S2
w = np.cos(np.radians(lat))[:, None] * np.ones(U.shape[2])
fin_all = np.all(np.isfinite(U) & np.isfinite(V), axis=0)
ice_any = np.nanmax(np.nan_to_num(ICE), axis=0) > 0.10
masks = {"nomask": np.ones_like(fin_all), "ice_any": fin_all & ~ice_any}
rows, out = {}, []
for mn, m in masks.items():
    for vn, F in (("KE", KE), ("EKE", EKE)):
        y = np.array([np.sum((F[k] * w)[m & np.isfinite(F[k])]) / np.sum(w[m & np.isfinite(F[k])]) for k in range(len(files))]) * 1e4
        rows[f"{mn}_{vn}"] = y
        for per, sel in (("2001-2025", t.year <= 2025), ("2001-2024", t.year <= 2024)):
            o = R.trend_calibrated(t.x.values[sel], y[sel], month=t.month.values[sel], n_outer=100, n_boot=300)
            o.update(series=f"Drake_{vn}_{mn}", period=per, mean=float(y[sel].mean()),
                     n_pix=int(m.sum()), pct_per_decade=1000 * o["slope"] / float(y[sel].mean()))
            out.append(o); print(o["series"], per, round(o["slope"], 3), round(o["p_cal"], 3), flush=True)
pd.DataFrame(rows).assign(year=t.year, month=t.month).to_csv(outdir / "a6b_drake_box_series.csv", index=False)
pd.DataFrame(out).to_csv(outdir / "a6b_drake_box_trends.csv", index=False)
print("pixels: total", fin_all.size, "finite-all", int(fin_all.sum()), "ice_any", int(ice_any.sum()), "kept", int(masks["ice_any"].sum()))
