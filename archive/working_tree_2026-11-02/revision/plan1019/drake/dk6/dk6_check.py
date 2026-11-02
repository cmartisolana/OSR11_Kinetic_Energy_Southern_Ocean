import numpy as np, xarray as xr
from dk6_common import *
TBs = np.load(f"{REV}/outputs/a2/tracks_baseline.npz"); glon = xr.open_dataset(f"{REV}/outputs/clim_allsat.nc").longitude.values
for icem in (False, True):
    D = load("allsat", icem); C, TB = detect_all(D["ke"], D["lat"])
    cc, ctb = F.detect_field(np.nanmean(D["ke"], 0), D["lat"], F.meridional_dists(D["lat"]), P["thr"], P["min_w"], P["frac"], MAXE)
    meta, arr = match(C, TB, cc, ctb, D["lat"])
    # compare with stored tracks at the same longitude and climatological half-power
    diffs = []; n = 0
    for q, (j, c, lo, hi, hpc, ls, ln) in enumerate(meta):
        lo_ = D["lon"][int(j)]; js = np.argmin(abs(glon - lo_))
        cand = np.where((TBs["j"] == js) & (abs(TBs["lat_hp_clim"] - hpc) < 0.05))[0]
        if len(cand) and hpc < -52:
            a = TBs["lat_hp"][:, cand[0]]; b = arr[:, q, FI["lat_hp"]]; ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() > 10: diffs.append(np.abs(a[ok] - b[ok]).max()); n += 1
    print("icemask", icem, "matched tracks compared", n, "max |diff| deg: median", np.median(diffs) if diffs else None, "max", max(diffs) if diffs else None)
    md=[]; fe=[]
    for q, (j, c, lo, hi, hpc, ls, ln) in enumerate(meta):
        js = np.argmin(abs(glon - D["lon"][int(j)]))
        cand = np.where((TBs["j"] == js) & (abs(TBs["lat_hp_clim"] - hpc) < 0.05))[0]
        if len(cand) and hpc < -52:
            a = TBs["lat_hp"][:, cand[0]]; b = arr[:, q, FI["lat_hp"]]; ok = np.isfinite(a) & np.isfinite(b)
            md.append(np.median(np.abs(a[ok]-b[ok]))); fe.append((np.abs(a[ok]-b[ok])<1e-3).mean())
    print("  median |diff| per track (median over tracks)", np.median(md), "fraction of months identical", np.median(fe))
