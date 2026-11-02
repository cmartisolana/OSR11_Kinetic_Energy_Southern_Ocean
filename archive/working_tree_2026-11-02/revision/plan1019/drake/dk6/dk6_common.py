"""Shared helpers for the dk6 Drake Passage half-power tests."""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[3])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, numpy as np, xarray as xr
REV = _REV; sys.path.insert(0, REV)
import fronts as F, rstats
OUT = (_DK_OUT + "/dk6")
FI = {k: i for i, k in enumerate(F.FIELDS)}
P = dict(thr=0.02, min_w=100.0, frac=0.25); MAXE = 8
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25); YR = np.repeat(np.arange(2001, 2026), 12)

def load(prod, icemask):
    z = np.load(f"{OUT}/fields_{prod}.npz")
    u, v = z["u"].astype(float), z["v"].astype(float)
    ke = 0.5 * (u ** 2 + v ** 2)
    ub, vb = np.nanmean(u, 0), np.nanmean(v, 0)
    eke = 0.5 * ((u - ub) ** 2 + (v - vb) ** 2) + 0.5 * z["vi"]
    adt = z["adt"].astype(float)
    if icemask:
        m = z["ice"] > 0.10
        ke[m] = np.nan; eke[m] = np.nan
    return dict(ke=ke, eke=eke, adt=adt, ice=z["ice"], lat=z["lat"].astype(float), lon=z["lon"].astype(float))

def detect_all(ke, lat):
    d = F.meridional_dists(lat); nt, _, nx = ke.shape
    C = np.zeros((nt, nx), int); TB = np.full((nt, nx, MAXE, F.NF), np.nan)
    for k in range(nt):
        C[k], TB[k] = F.detect_field(np.ascontiguousarray(ke[k]), lat, d, P["thr"], P["min_w"], P["frac"], MAXE)
    return C, TB

def match(C, TB, cc, ctb, lat):
    """as a2_run.matched_tracks: monthly envelope whose half-power lies inside the climatological segment; strongest if several"""
    meta = []
    for j in range(cc.size):
        for c in range(cc[j]):
            s, e = int(ctb[j, c, FI["i_start"]]), int(ctb[j, c, FI["i_end"]])
            meta.append((j, c, lat[s], lat[e], ctb[j, c, FI["lat_hp"]], ctb[j, c, FI["lat_s"]], ctb[j, c, FI["lat_n"]]))
    meta = np.array(meta)
    nt = C.shape[0]; out = np.full((nt, len(meta), F.NF), np.nan)
    for k in range(nt):
        for q, (j, c, lo, hi, *_r) in enumerate(meta):
            R = TB[k, int(j)]; lh = R[:, FI["lat_hp"]]
            ins = np.isfinite(lh) & (lh >= lo) & (lh <= hi)
            if ins.any():
                b = np.argmax(np.where(ins, R[:, FI["peak_val"]], -np.inf)); out[k, q] = R[b]
    return meta, out

def tr(y, monthly=True, x=None):
    x = T if x is None else x
    r = rstats.trend_bcz(x, y, month=MON if monthly else None, monthly=monthly)
    return r
