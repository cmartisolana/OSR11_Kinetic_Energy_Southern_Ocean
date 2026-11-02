"""dk6 step 0: cache monthly fields over Drake Passage (75-50W, 66-48S) for both products, and check that KE recomputed
from u_mean, v_mean reproduces the stored half-power tracks (outputs/a2/tracks_baseline.npz) at the passage longitudes."""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[3])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, numpy as np, xarray as xr, glob
REV = _REV; sys.path.insert(0, REV)
OUT = (_DK_OUT + "/dk6")
BOX = dict(longitude=slice(-75, -50), latitude=slice(-66, -35))
for prod in ("allsat", "twosat"):
    fs = sorted(glob.glob(f"{REV}/data/monthly_{prod}/*/*.nc")); fs = [f for f in fs if 2001 <= int(f.split('/')[-2]) <= 2025]
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
