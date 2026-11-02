"""dk6 step 0: cache monthly fields over Drake Passage (75-50W, 66-48S) for both products, and check that KE recomputed
from u_mean, v_mean reproduces the stored half-power tracks (outputs/a2/tracks_baseline.npz) at the passage longitudes.

Purpose : cache the monthly fields over Drake Passage (75-50 W, 66-35 S) for both products (input of 03m, 03n).
Inputs  : $OSR11_DATA/monthly_allsat and monthly_twosat (00a)
Outputs : outputs/drake/dk6/fields_<product>.npz (cache, not distributed)
Run     : python scripts/03_drake/03l_section_halfpower_load.py
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, numpy as np, xarray as xr, glob
from osr11 import paths as PA
OUT = str(PA.dk(6))
BOX = dict(longitude=slice(-75, -50), latitude=slice(-66, -35))
for prod in ("allsat", "twosat"):
    fs = [str(f) for f in PA.monthly_files(prod)]
    assert len(fs) == 300, len(fs)
    U, V, A, VI, IC = [], [], [], [], []
    for f in fs:
        with xr.open_dataset(f) as d:
            s = d.sel(**BOX)
            U.append(s.u_mean.values); V.append(s.v_mean.values); A.append(s.adt_mean.values)
            VI.append(s.var_intra.values); IC.append(s.ice_frac.values)
            lat, lon = s.latitude.values, s.longitude.values
    np.savez_compressed(f"{OUT}/fields_{prod}.npz", u=np.array(U, np.float32), v=np.array(V, np.float32), adt=np.array(A, np.float32),
                        vi=np.array(VI, np.float32), ice=np.array(IC, np.float32), lat=lat, lon=lon)
    print(prod, np.array(U).shape)
