"""dk10: occurrence density of monthly half-power points (all envelopes detected each month, baseline 200 cm2 s-2, 25 %,
100 km, Product 1) in the Drake Passage domain (86-34W), 2001-2025. Points kept where the 2001-2025 mean ADT at the point
is in the ACC band (-0.6..0.2 m). Density: percentage of months with a point in each 0.5 deg lon x 0.25 deg lat cell, per
longitude column (so a stationary front gives a narrow, high band, a meandering front a broad, low one). Also saved for
2001-2012 and 2014-2025 to compare halves."""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[3])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, glob, numpy as np, xarray as xr
sys.path.insert(0, _DK); from dk6_common import F, FI, P, MAXE, detect_all
REV = _REV; OUT = (_DK_OUT + "/dk10")
BOX = dict(longitude=slice(-87, -33), latitude=slice(-66, -35))
fs = sorted(glob.glob(f"{REV}/data/monthly_allsat/*/*.nc")); fs = [f for f in fs if 2001 <= int(f.split('/')[-2]) <= 2025]
KE = []
for f in fs:
    with xr.open_dataset(f) as d:
        s = d.sel(**BOX); KE.append((0.5 * (s.u_mean.values ** 2 + s.v_mean.values ** 2)).astype(np.float32)); lat, lon = s.latitude.values.astype(float), s.longitude.values.astype(float)
C, TB = detect_all(np.array(KE), lat)
adt = xr.open_dataset((_ROOT_OUT + "/mean_adt.nc")).mean_adt.sel(**BOX).interp(latitude=lat, longitude=lon).values
le = np.arange(-86, -33.99, 0.5); ae = np.arange(-64, -45.99, 0.25)
yr = np.repeat(np.arange(2001, 2026), 12)
def dens(mask):
    H = np.zeros((ae.size - 1, le.size - 1))
    for k in np.where(mask)[0]:
        for j in range(lon.size):
            for c in range(C[k, j]):
                la = TB[k, j, c, FI["lat_hp"]]
                if not np.isfinite(la): continue
                i = min(np.searchsorted(lat, la), lat.size - 1)
                if not (-0.6 <= adt[i, j] <= 0.2): continue
                if -62 <= lon[j] < -56 and la > -52: continue      # Malvinas Current (shelf-break branch)
                if not (ae[0] <= la < ae[-1] and le[0] <= lon[j] < le[-1]): continue
                H[np.searchsorted(ae, la) - 1, np.searchsorted(le, lon[j]) - 1] += 1
    ncol = np.array([((lon >= a) & (lon < a + 0.5)).sum() for a in le[:-1]])   # grid longitudes per column
    return 100 * H / (mask.sum() * np.maximum(ncol, 1))[None, :]
np.savez(f"{OUT}/dk10_hp_density.npz", lon_edges=le, lat_edges=ae, all=dens(np.ones(300, bool)),
         early=dens(yr <= 2012), late=dens(yr >= 2014))
print("ok")
