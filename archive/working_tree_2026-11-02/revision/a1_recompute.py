"""A1: recompute the KE/EKE trend numbers quoted in v2 with the single
statistical framework (deseasonalised Theil-Sen; AR-sieve CI/p; double-
bootstrap calibrated p; Hamed-Rao for comparison).

Series are rebuilt from the stored monthly fields (outputs/monthly_half-power_points):
  ACC band  : mean ADT between -0.6 and +0.2 m (as Fig. 1e/f in v2), cos-lat weights
  sectors   : (i) inside the ACC band, (ii) whole 35-65S box (what Fig. 1f sectors use)
  Drake box : 80-52W, 63-54S (as notebook 07)
Periods: 2001-2025 (text) and 2001-2024 (transport window / map period).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import a2_run as A, rstats as R

ke, lat, lon, acc, t = A.load()
eke = np.load(A.CACHE / "eke.npy", mmap_mode="r")
w = np.cos(np.radians(lat))[:, None] * np.ones((1, lon.size))
sec = A.sector(lon)[None, :] * np.ones((lat.size, 1), dtype=object)
L2, O2 = np.meshgrid(lat, lon, indexing="ij")
masks = {
    "ACC_band": acc,
    "Atlantic_band": acc & (sec == "Atlantic"), "Indian_band": acc & (sec == "Indian"), "Pacific_band": acc & (sec == "Pacific"),
    "Atlantic_box": sec == "Atlantic", "Indian_box": sec == "Indian", "Pacific_box": sec == "Pacific",
    "Drake_box": (O2 >= -80) & (O2 <= -52) & (L2 >= -63) & (L2 <= -54),
}
series = {}
for name, m in masks.items():
    for vname, V in [("KE", ke), ("EKE", eke)]:
        y = np.empty(V.shape[0])
        for k in range(V.shape[0]):
            f = np.asarray(V[k], float)
            ok = m & np.isfinite(f)
            y[k] = np.sum(f[ok] * w[ok]) / np.sum(w[ok]) * 1e4      # cm2 s-2
        series[(name, vname)] = y
rows = []
for (name, vname), y in series.items():
    for per, sel in [("2001-2025", np.ones(300, bool)), ("2001-2024", t.year.values <= 2024)]:
        o = R.trend_calibrated(t.x.values[sel], y[sel], month=t.month.values[sel], n_outer=100, n_boot=400)
        o.update(region=name, var=vname, period=per, mean=float(np.mean(y[sel])),
                 pct_per_decade=100 * 10 * o["slope"] / float(np.mean(y[sel])))
        rows.append(o)
        print(name, vname, per, {k: round(v, 4) if isinstance(v, float) else v for k, v in o.items() if k in ("slope", "ci_lo", "ci_hi", "p_boot", "p_cal", "p_hamed_rao")}, flush=True)
out = pd.DataFrame(rows)
out.to_csv(HERE / "outputs" / "a1_ke_eke_trends.csv", index=False)
pd.DataFrame({f"{k[0]}_{k[1]}": v for k, v in series.items()}).assign(x=t.x.values).to_csv(HERE / "outputs" / "a1_series.csv", index=False)
