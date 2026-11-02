"""Monthly KE and EKE fields used by the half-power analysis (replaces the first-submission half_power_EKE.py output).

Purpose : write, for every month of 2001-2025, the two fields that the circumpolar half-power analysis stacks into its cache
          (osr11.halfpower.build_cache): KE of the monthly-mean surface geostrophic velocity and EKE relative to the
          2001-2025 mean of the daily velocities,
              ke  = 0.5 |U_m|^2                                   (m2 s-2)
              eke = 0.5 |U_m - U_lt|^2 + 0.5 var_intra             (m2 s-2)
          on the 1/8 deg grid, NaN on the permanent sea-ice mask.
          The first-submission files (outputs/monthly_half-power_points/YYYY/MM.nc of the working repository, on the author's Mac) were written by
          half_power_EKE.py from the daily fields; that script is archived (it also stored the first-submission front
          variables, which no retained step reads). This script rebuilds the same two fields from the monthly files;
          checked on 2001-01 against the first-submission file: identical NaN pattern, max |difference| 7.0e-8 (ke) and
          6.7e-8 m2 s-2 (eke), i.e. float32 rounding of the monthly files.
Inputs  : $OSR11_DATA/monthly_allsat/YYYY/MM.nc (00a), outputs/base/sea_ice_mask.nc, $OSR11_DATA/preprocessing.nc (u_ltm, v_ltm; 01a)
Outputs : $OSR11_DATA/monthly_half-power_points/YYYY/MM.nc (ke, eke; float64), 300 files, ~3.5 GB
Run     : python scripts/01_base/01b_halfpower_fields.py [--check FILE]   (FILE: a first-submission month to compare with)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import argparse
from pathlib import Path

import numpy as np
import xarray as xr

from osr11 import paths as PA

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("--check", type=Path, default=None, help="compare 2001-01 with this first-submission file and stop")
a = ap.parse_args()

ice = xr.open_dataset(PA.ICE_MASK).ice_mask.values.astype(bool)
pre = xr.open_dataset(PA.PREPROC)
ult, vlt = pre.u_ltm.values.astype(np.float64), pre.v_ltm.values.astype(np.float64)


def fields(f):
    with xr.open_dataset(f) as m:
        u, v = m.u_mean.values.astype(np.float64), m.v_mean.values.astype(np.float64)
        vi = m.var_intra.values.astype(np.float64)
        lat, lon = m.latitude.values, m.longitude.values
    ke = 0.5 * (u ** 2 + v ** 2)
    eke = 0.5 * ((u - ult) ** 2 + (v - vlt) ** 2) + 0.5 * vi
    ke[ice] = np.nan; eke[ice] = np.nan
    return lat, lon, ke, eke


if a.check:
    lat, lon, ke, eke = fields(PA.monthly_files()[0])
    with xr.open_dataset(a.check) as h:
        for k, x in (("ke", ke), ("eke", eke)):
            r = h[k].values
            print(k, "NaN pattern identical:", bool((np.isnan(r) == np.isnan(x)).all()), "max |diff|:", float(np.nanmax(np.abs(r - x))))
    raise SystemExit

for f in PA.monthly_files():
    y, mo = f.parent.name, f.stem
    lat, lon, ke, eke = fields(f)
    out = PA.HP_POINTS / y / f"{mo}.nc"; out.parent.mkdir(parents=True, exist_ok=True)
    xr.Dataset({"ke": (("latitude", "longitude"), ke), "eke": (("latitude", "longitude"), eke)},
               coords={"latitude": lat, "longitude": lon},
               attrs={"description": "KE of the monthly-mean velocity and EKE w.r.t. the 2001-2025 mean (m2 s-2)",
                      "eke_reference": "long-term mean over full study period (preprocessing.nc)",
                      "ice_masking": "permanent (sea_ice_mask.nc)", "year": y, "month": mo,
                      "script": "scripts/01_base/01b_halfpower_fields.py"}).to_netcdf(out)
print("done", PA.HP_POINTS)
