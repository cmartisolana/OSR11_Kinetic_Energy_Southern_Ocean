"""Number of common cells of the wind / KE / EKE / slope comparison (Sect. 3.3, Fig. 4c-e; written for the publication
repository so that the cell count quoted in the text has a stored source).

Purpose : reproduce the domain mask of 04g_common_cells.py without recomputing the series: the sigma0 = 27.3 slope domain
          (ACC band + 100 km, ice removed, slope defined in all 264 months; outputs/f4m/f4m_map.nc 'domain', 2835 cells) restricted
          to 1-deg cells with valid KE and EKE in every month (at least one finite, ice-free 1/8-deg pixel per cell and month,
          with the 2001-2025 mean velocity defined) and, when $OSR11_DATA/wind_taux.nc is present, valid wind in every month.
Inputs  : outputs/f4/f4_tilt_maps.nc, outputs/f4m/f4m_map.nc, outputs/base/sea_ice_mask.nc, $OSR11_DATA/monthly_allsat (00a),
          optional $OSR11_DATA/wind_taux.nc (00b)
Outputs : outputs/f4e/f4e_common_cells.csv (n_slope_domain, n_with_ke_eke, n_with_wind or empty if the wind file is absent)
Run     : python scripts/05_stats/05i_common_cells_count.py      (~1 min; needs the monthly files)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd, xarray as xr
from osr11 import paths as PA

M = xr.open_dataset(PA.OUT / "f4" / "f4_tilt_maps.nc"); lat, lon = M.lat.values, M.lon.values
dom = xr.open_dataset(PA.OUT / "f4m" / "f4m_map.nc")["domain"].values > 0
la_e = np.r_[lat - 0.5, lat[-1] + 0.5]; lo_e = np.r_[lon - 0.5, lon[-1] + 0.5]
files = PA.monthly_files("allsat"); assert all(f.exists() for f in files)
with xr.open_dataset(files[0]) as d0: la8, lo8 = d0.latitude.values, d0.longitude.values
ice8 = xr.open_dataset(PA.ICE_MASK).ice_mask.values > 0.5
i8 = np.digitize(la8, la_e) - 1; j8 = np.digitize(lo8, lo_e) - 1
allfin = np.ones((la8.size, lo8.size), bool); fin = []
for f in files:
    with xr.open_dataset(f) as d:
        ok = np.isfinite(d.u_mean.values) & np.isfinite(d.v_mean.values); okv = ok & np.isfinite(d.var_intra.values)
    allfin &= ok; fin.append(okv)


def tocell(b8):
    out = np.zeros((lat.size, lon.size), bool)
    for i in range(lat.size):
        ri = (i8 == i)
        if not ri.any(): continue
        sub = b8[ri, :].any(0)
        for j in range(lon.size):
            out[i, j] = sub[j8 == j].any()
    return out


valid = np.ones((lat.size, lon.size), bool)
for okv in fin:
    valid &= tocell(okv & allfin & ~ice8)
n1 = int((dom & valid).sum())
nw = None
if PA.WIND.exists():
    W = xr.open_dataset(PA.WIND).eastward_stress.sel(latitude=slice(lat[0] - 0.5, lat[-1] + 0.5))
    il = np.digitize(W.latitude.values, la_e) - 1; jl = np.digitize(W.longitude.values, lo_e) - 1; A = np.isfinite(W.values)
    wv = np.zeros((lat.size, lon.size), bool)
    for i in range(lat.size):
        for j in range(lon.size):
            sub = A[:, il == i, :][:, :, jl == j]; wv[i, j] = sub.size > 0 and sub.any(axis=(1, 2)).all()
    nw = int((dom & valid & wv).sum())
T = pd.DataFrame([dict(n_slope_domain=int(dom.sum()), n_with_ke_eke=n1, n_with_wind=nw)])
T.to_csv(PA.out("f4e") / "f4e_common_cells.csv", index=False); print(T.to_string(index=False))
