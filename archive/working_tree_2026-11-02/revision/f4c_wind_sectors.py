"""F4c: eastward wind-stress trends over the Figure 4 domain by sector, 2001-2025 (28 Oct 2026), for Figure 4e.
Same wind processing as f4b_tilt_wind.py (Product 2 monthly eastward stress -> 1 deg RG cells -> annual means) and the
same cells as the sigma0 = 27.3 slope (ACC band + 100 km, ice removed, slope defined in all years). Trends: rstats.trend_bcz
on annual series (AR(1), 300 null series), % of the mean per decade. Sectors as f4_tilt.py.
Usage: python f4c_wind_sectors.py <revision dir> <f4 out dir> <out dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

REV, F4, OUT = (Path(a) for a in sys.argv[1:4]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402

M = xr.open_dataset(F4 / "f4_tilt_maps.nc"); lat, lon = M.lat.values, M.lon.values
dom = (M["buffer100"].values > 0) & (M["ice"].values == 0) & np.isfinite(M["signed_27p3_pct"].values)
LON, LAT = np.meshgrid(lon, lat); w = np.cos(np.radians(LAT))
W = xr.open_dataset(REV / "data" / "wind_taux.nc").eastward_stress.sel(latitude=slice(lat[0] - 0.5, lat[-1] + 0.5))
Wy = W.groupby("time.year").mean("time"); yrs = Wy.year.values
la_e = np.r_[lat - 0.5, lat[-1] + 0.5]; lo_e = np.r_[lon - 0.5, lon[-1] + 0.5]
il = np.digitize(Wy.latitude.values, la_e) - 1; jl = np.digitize(Wy.longitude.values, lo_e) - 1
A = Wy.values; T = np.full((yrs.size, lat.size, lon.size), np.nan)
for i in range(lat.size):
    sub = A[:, il == i, :]
    if sub.shape[1] == 0: continue
    for j in range(lon.size):
        T[:, i, j] = np.nanmean(sub[:, :, jl == j], axis=(1, 2))
SECT = {"ACC": np.ones_like(dom), "Atlantic": (LON >= -70) & (LON < 20), "Indian": (LON >= 20) & (LON < 147),
        "Pacific": (LON >= 147) | (LON < -70)}
rows, ser = [], {}
for name, ms in SECT.items():
    m = dom & ms
    y = np.array([np.nansum(T[t][m] * w[m]) / np.sum(w[m] * np.isfinite(T[t][m])) for t in range(yrs.size)])
    ser[name] = y
    for y0 in (2001, 2004):
        s = yrs >= y0; r = R.trend_bcz(yrs[s] + 0.5, y[s], monthly=False, n_outer=300); mu = y[s].mean()
        rows.append(dict(region=name, period=f"{y0}-2025", mean=mu, pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu,
                         pct_hi=1000 * r["ci_hi"] / mu, p=r["p"], ncell=int(m.sum())))
        print(name, y0, f"{rows[-1]['pct']:+.2f} [{rows[-1]['pct_lo']:+.2f}, {rows[-1]['pct_hi']:+.2f}] p={r['p']:.3f}", flush=True)
pd.DataFrame(rows).to_csv(OUT / "f4c_wind_trends.csv", index=False)
np.savez(OUT / "f4c_wind_series.npz", years=yrs, **ser)
