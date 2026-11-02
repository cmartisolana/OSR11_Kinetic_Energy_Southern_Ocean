"""B5b: change of the mean jet-core path at the five meanders, 2013-2025 vs 2001-2012 (replaces the max-min amplitude of
B5, which switches between near-equal crests and mixes in the large-scale tilt of the path).
Mean path: +0.18 m relative-ADT contour of the period-mean ADT (as B3-B5, 150 km smoothing). Latitude interpolated on a
0.25 deg longitude grid from the upstream start to the downstream end of each meander (where the path is single-valued in
longitude; not possible at Drake Passage-Scotia Sea, where the path turns north with the Malvinas Current and back).
Statistic: RMS latitude difference (km) between the two period paths; null: 500 random splits of the 25 years into 12 and
13 years (no temporal order), p = fraction of splits with RMS >= observed. Also the mean signed displacement (km, + north), and a
curve distance valid for every meander: symmetric mean nearest-point distance between the two paths in the box (Drake:
south of 52 S only, excluding the Malvinas Current loop), with the same null.
Usage: python b5b_path_change.py <revision dir> <b5 out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, contourpy
from scipy.ndimage import gaussian_filter1d
REV, B5 = Path(sys.argv[1]), Path(sys.argv[2])
src = open(Path(__file__).parent / "b5_meander_shape.py").read()
exec(src[src.index("MEANDERS = {"):src.index("files = [")])
with xr.open_dataset(REV / "data" / "monthly_allsat" / "2001" / "01.nc") as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
A = np.load(B5 / "b5_annual_fields.npz")["A"]
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
ACCM = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice
EXCL = (mean_adt > 0.45) | (lat[:, None] > -40) | ~np.isfinite(mean_adt)
exec(src[src.index("PAD = int(30"):src.index("def amp(")])


def latline(P, a, b):
    g = np.arange(a, b + 1e-6, 0.25)
    for lo_s, la_s in P:
        s = (lo_s >= a - 1) & (lo_s <= b + 1)
        if s.sum() < 3 or lo_s[s].min() > a or lo_s[s].max() < b: continue
        x = lo_s[s]
        if np.any(np.diff(x) <= 0): return None                 # folds back in longitude
        return np.interp(g, x, la_s[s])
    return None


from scipy.spatial import cKDTree
LATMAX = {"Drake Passage–Scotia Sea": -52.0}                  # exclude the Malvinas Current loop north of 52 S
def pts(P, a, b, lmax=-40.0):                              # path points inside the meander box, local km
    X = [np.c_[lo_s, la_s][(lo_s >= a) & (lo_s <= b) & (la_s <= lmax)] for lo_s, la_s in P]; X = np.vstack([x for x in X if len(x)]) if X else np.empty((0, 2))
    return np.c_[111.195 * np.cos(np.radians(52)) * X[:, 0], 111.195 * X[:, 1]]
def cdist(P1, P2, a, b, lmax=-40.0):                       # symmetric mean nearest-point distance between the two paths (km)
    x1, x2 = pts(P1, a, b, lmax), pts(P2, a, b, lmax)
    if len(x1) < 3 or len(x2) < 3: return np.nan
    return 0.5 * (cKDTree(x2).query(x1)[0].mean() + cKDTree(x1).query(x2)[0].mean())
rng = np.random.default_rng(20261020); NULL = 500
splits = [rng.permutation(25) for _ in range(NULL)]
cache = {}
def lines_for(idx):
    P = path(np.nanmean(A[np.sort(idx)], 0)); return P, {nm: latline(P, u_[0], d_[1]) for nm, (u_, r_, d_) in MEANDERS.items()}
(PE, obs_e), (PL, obs_l) = lines_for(np.arange(12)), lines_for(np.arange(12, 25))
null = {nm: [] for nm in MEANDERS}; nulld = {nm: [] for nm in MEANDERS}
for k, sp in enumerate(splits):
    (P1, e), (P2, l) = lines_for(sp[:12]), lines_for(sp[12:])
    for nm, (u_, r_, d_) in MEANDERS.items():
        if e[nm] is not None and l[nm] is not None: null[nm].append(111.195 * np.sqrt(np.mean((l[nm] - e[nm]) ** 2)))
        nulld[nm].append(cdist(P1, P2, u_[0], d_[1], LATMAX.get(nm, -40.0)))
    if k % 100 == 0: print("split", k, flush=True)
rows = []
for nm, (u_, r_, d_) in MEANDERS.items():
    od = cdist(PE, PL, u_[0], d_[1], LATMAX.get(nm, -40.0)); nd = np.array(nulld[nm]); nd = nd[np.isfinite(nd)]
    dist = dict(curve_dist_km=od, curve_null_median_km=float(np.median(nd)), curve_null_95_km=float(np.percentile(nd, 95)), curve_p=float((nd >= od).mean()))
    if obs_e[nm] is None or obs_l[nm] is None:
        rows.append(dict(meander=nm, rms_km=np.nan, note="path not single-valued in longitude", **dist)); continue
    dlat = 111.195 * (obs_l[nm] - obs_e[nm]); rms = float(np.sqrt(np.mean(dlat ** 2))); nl = np.array(null[nm])
    rows.append(dict(meander=nm, rms_km=rms, mean_shift_km=float(dlat.mean()), max_abs_km=float(np.abs(dlat).max()),
                     null_median_km=float(np.median(nl)), null_95_km=float(np.percentile(nl, 95)), p=float((nl >= rms).mean()), n_null=nl.size, **dist))
D = pd.DataFrame(rows); D.to_csv(B5 / "b5b_path_change.csv", index=False); print(D.round(3).to_string(index=False))
np.savez(B5 / "b5b_lines.npz", **{f"{nm}|early": v for nm, v in obs_e.items() if v is not None}, **{f"{nm}|late": v for nm, v in obs_l.items() if v is not None})
