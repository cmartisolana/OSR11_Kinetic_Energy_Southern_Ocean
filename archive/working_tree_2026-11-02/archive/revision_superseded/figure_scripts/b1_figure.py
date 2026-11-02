"""B1 pilot figure (Pacific-Antarctic Ridge): pathway, curvature, upstream speed, downstream EKE.

(a) climatological core-level pathway (relative ADT contour) on GEBCO bathymetry with the
    Liu et al. (2026) sections; (b-d) annual means of the deseasonalised monthly series
    with calibrated trends and the minimum detectable 25-yr change (MDC).
-> outputs/figures/figS_b1_par_pilot.png
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b1_meander as B

data = Path(sys.argv[1])
O = HERE / "outputs"
P = pd.read_csv(O / "b1" / "b1_pilot_series.csv"); T = pd.read_csv(O / "b1" / "b1_pilot_trends.csv").set_index("series")
REG = dict(lon0=-152.0, lon1=-78.0, lat0=-62.0, lat1=-45.0)
A = []
for y in range(2001, 2026):
    for m in range(1, 13):
        with xr.open_dataset(data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc") as ds:
            s = ds.adt_mean.sel(longitude=slice(REG["lon0"], REG["lon1"]), latitude=slice(REG["lat0"], REG["lat1"]))
            A.append(s.values)
lat, lon = s.latitude.values, s.longitude.values
Ab = np.nanmean(np.array(A), 0)
mbox = dict(lon0=-135.13, lon1=-109.38, lat0=-62.0, lat1=-45.0)
rel, la, lo = B.box_relative(Ab, lat, lon, mbox)
core = -0.18
bat = xr.open_dataset(data / "bathymetry_so.nc").elevation.sel(longitude=slice(-152, -78), latitude=slice(-62, -45))

fig = plt.figure(figsize=(7.5, 9))
ax = fig.add_axes([0.1, 0.68, 0.85, 0.28])
cf = ax.pcolormesh(bat.longitude, bat.latitude, -bat, cmap="Blues", vmin=0, vmax=5500, shading="auto")
plt.colorbar(cf, ax=ax, label="depth (m)", pad=0.01)
for lev, ls in ((core - 0.10, ":"), (core, "-"), (core + 0.10, "--")):
    line = B.spanning_contour(rel, la, lo, lev)
    if line is not None:
        ax.plot(line[:, 0], line[:, 1], "k", ls=ls, lw=1.2)
for x in (-150, -135.13, -109.38, -80):
    ax.axvline(x, color="#b2182b", lw=0.8)
for xm, t in ((-142.5, "upstream"), (-122.3, "meander"), (-94.7, "downstream")):
    ax.text(xm, -46, t, ha="center", fontsize=8, color="#b2182b")
ax.set_ylim(-62, -45); ax.set_title("(a) mean pathway (core level, ±0.10 m) and Liu et al. (2026) sections", fontsize=8)

def panel(pos, col, lab, title):
    a = fig.add_axes(pos)
    y = P[col] - P.groupby("month")[col].transform("mean")
    an = y.groupby(P.year).mean()
    a.plot(an.index + 0.5, an.values, "-o", ms=3, color="#333333")
    r = T.loc[col]
    a.text(0.01, 0.95, f"trend {r.slope:+.3g} yr⁻¹ [{r.ci_cal_lo:+.3g}, {r.ci_cal_hi:+.3g}], p_cal={r.p_cal:.2f}; "
           f"interannual SD {r.interannual_sd:.3g}; MDC(25 yr) {r.mdc_25yr:.3g}", transform=a.transAxes, va="top", fontsize=7)
    a.set_ylabel(lab, fontsize=8); a.set_title(title, fontsize=8); a.grid(alpha=0.25, lw=0.5)
    return a
panel([0.1, 0.47, 0.85, 0.15], "core_s100_int_abs_kappa", "∫|κ| ds (rad)", "(b) pathway curvature, meander section (100 km smoothing)")
panel([0.1, 0.27, 0.85, 0.15], "up_core_speed", "speed (m s⁻¹)", "(c) upstream jet-core speed (max over 60–48°S)")
a = panel([0.1, 0.07, 0.85, 0.15], "down_eke", "EKE (cm² s⁻²)", "(d) downstream EKE")
a.set_xlabel("year")
(O / "figures").mkdir(exist_ok=True)
fig.savefig(O / "figures" / "figS_b1_par_pilot.png", dpi=200)
print("saved")
