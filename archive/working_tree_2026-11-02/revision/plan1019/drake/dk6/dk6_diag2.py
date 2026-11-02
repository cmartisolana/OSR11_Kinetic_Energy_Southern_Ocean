import numpy as np, xarray as xr
from dk6_common import *
TBs = np.load(f"{REV}/outputs/a2/tracks_baseline.npz"); glon = xr.open_dataset(f"{REV}/outputs/clim_allsat.nc").longitude.values
z = np.load(f"{OUT}/fields_allsat.npz"); lat = z["lat"].astype(float); lon = z["lon"].astype(float)
ke0 = 0.5 * (z["u"].astype(float) ** 2 + z["v"].astype(float) ** 2); ice = z["ice"]
kec_ref = np.nanmean(ke0, 0)
for name, mask in (("none", None), ("ice>0", ice > 0), ("ice>0.1", ice > 0.1), ("ice>0.5", ice > 0.5)):
    ke = ke0.copy()
    if mask is not None: ke[mask] = np.nan
    C, TB = detect_all(ke, lat)
    cc, ctb = F.detect_field(np.ascontiguousarray(np.nanmean(ke, 0)), lat, F.meridional_dists(lat), P["thr"], P["min_w"], P["frac"], MAXE)
    meta, arr = match(C, TB, cc, ctb, lat)
    tot = same = 0
    for q, (j, c, lo, hi, hpc, ls, ln) in enumerate(meta):
        if not (-70 <= lon[int(j)] < -55 and hpc < -52): continue
        js = np.argmin(abs(glon - lon[int(j)])); cand = np.where((TBs["j"] == js) & (abs(TBs["lat_hp_clim"] - hpc) < 0.05))[0]
        if not len(cand): continue
        a = TBs["lat_hp"][:, cand[0]]; b = arr[:, q, FI["lat_hp"]]
        ok = np.isfinite(a) | np.isfinite(b); tot += ok.sum(); same += (np.abs(np.nan_to_num(a, nan=99) - np.nan_to_num(b, nan=99)) < 1e-3)[ok].sum()
    print(name, "identical months fraction", same / tot)
