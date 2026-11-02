"""H0: Antarctic-slope sea-level index after Hogg et al. (2015) / Hughes et al. (2014) -- TOOLING, needs new data.

Hogg et al. (2015, JGR Oceans 120, 257-267) averaged AVISO gridded sea level over the Antarctic continental shelf and
slope "as defined in Hughes et al. (2014)", removed a bias and the annual + semiannual cycles, extracted the
long-period part h_f by fitting sines/cosines with periods of 10, 20 and 40 yr, and converted with the frequency-
dependent factors of Hughes et al. (2014):  T = -3.69 (h - h_f) - 1.11 h_f   [Sv, h in cm]
(signs as read from the extracted text of Hogg et al. -- CHECK AGAINST THE PDF; the factors are model-derived
(NEMO ORCA025/ORCA12) and called "preliminary" by Hogg et al.; no numerical transport trend is reported there).
Hughes et al. (2014) define the slope region from model correlations (< -0.5), average bottom pressure in 50 depth bins
of 100 m (0-5000 m) and then average the bins; correlations are strongest at 1-3 km depth.

Altimetry-only approximation used here (all choices are arguments, fixed before looking at any index):
  region: ocean points south of LAT_MAX (default 60 S) with depth DEP0..DEP1 (default 1000-3000 m) in the bathymetry,
          RESTRICTED to the continental slope contiguous with Antarctica: points shallower than DEP1 (incl. land and
          shelf) are labelled into connected regions (8-connectivity, periodic in longitude) and only the region(s)
          touching the southernmost grid row (the continent) are kept, so that isolated rises, plateaus and ridge
          crests (Maud Rise, southern Kerguelen Plateau, Pacific-Antarctic Ridge) are excluded (review v3, M6a).
          Set env H0_CONTIGUOUS=0 to reproduce the unrestricted selection as a sensitivity test;
          the kept/removed point counts are printed so the mask can be checked on a map before use.
  per month: ADT anomaly (minus each point's 2001-2025 mean); average within 100-m depth bins (cos-lat weights), then
  equal-weight mean of the bins (as Hughes et al.); ALSO the same in 12 longitude bands of 30 deg;
  coverage: fraction of region points with valid ADT and ice_frac <= 0.10 in that month (sea ice!); months with
  coverage < MINCOV flagged and excluded from the index.
Outputs: index h (cm), coverage, h_f, and the converted series T_hogg (Sv) LABELLED as model-calibrated.
Data needed (not on disk yet; the existing monthly fields stop at 65 S):
  python revision/extract_monthly.py --product allsat --lat-min -80 --lat-max -58 --out-dir revision/data/monthly_allsat_south
  python revision/extract_monthly.py --product twosat --start-year 1993 --lat-min -80 --lat-max -58 --out-dir revision/data/monthly_twosat_south
  bathymetry to 80 S: python h0_prepare_bathymetry_80S.py <global ETOPO/GEBCO file>  -> revision/data/bathymetry_80s.nc
Usage: python h0_antarctic_slope_index.py <monthly dir> <bathymetry nc> <out csv> [y0 y1 LAT_MAX DEP0 DEP1 MINCOV]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

mdir, bat_f, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
y0, y1 = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (2001, 2025)
LAT_MAX, DEP0, DEP1, MINCOV = (float(a) for a in sys.argv[6:10]) if len(sys.argv) > 9 else (-60.0, 1000.0, 3000.0, 0.5)
files = [mdir/f"{y}"/f"{m:02d}.nc" for y in range(y0, y1+1) for m in range(1, 13)]
files = [f for f in files if f.exists()]
with xr.open_dataset(files[0]) as d0:
    lat, lon = d0.latitude.values, d0.longitude.values
dep = -xr.open_dataset(bat_f).elevation.interp(latitude=xr.DataArray(lat, dims="latitude"),
                                                longitude=xr.DataArray(lon, dims="longitude")).values
dep = np.nan_to_num(dep, nan=-1.0)
lband = (np.mod(lon, 360)//30).astype(int)[None, :]*np.ones(dep.shape, int)
w = np.cos(np.radians(lat))[:, None]*np.ones(dep.shape)
reg = (lat[:, None] <= LAT_MAX) & (dep >= DEP0) & (dep <= DEP1)
import os
if os.environ.get("H0_CONTIGUOUS", "1") != "0":
    from scipy import ndimage
    shallow = (dep < DEP1) & (lat[:, None] <= LAT_MAX + 5)        # land, shelf and slope (dep < 0 on land)
    lab, n = ndimage.label(shallow, structure=np.ones((3, 3), int))
    # periodic longitude: merge labels that touch across the dateline seam
    parent = np.arange(n + 1)
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for i in range(lab.shape[0]):
        for di in (-1, 0, 1):
            k = i + di
            if 0 <= k < lab.shape[0] and lab[i, 0] and lab[k, -1]:
                parent[find(lab[i, 0])] = find(lab[k, -1])
    roots = np.array([find(a) for a in range(n + 1)])
    lab = roots[lab]
    srow = 0 if lat[0] < lat[-1] else -1                            # southernmost row
    keep = np.unique(lab[srow][lab[srow] > 0])
    cont = np.isin(lab, keep) & (lab > 0)
    print(f"H0 region: {reg.sum()} points in depth range; {(reg & cont).sum()} contiguous with Antarctica kept, "
          f"{(reg & ~cont).sum()} isolated (rises/plateaus/ridges) removed")
    reg = reg & cont
dbin = np.clip(((dep - DEP0)//100).astype(int), 0, int((DEP1-DEP0)//100) - 1)
ii, jj = np.nonzero(reg)                    # keep only region points (memory-light)
A = np.empty((len(files), ii.size)); ICE = np.empty_like(A)
for k, f in enumerate(files):
    with xr.open_dataset(f) as d:
        A[k] = d.adt_mean.values[ii, jj]; ICE[k] = np.nan_to_num(d.ice_frac.values[ii, jj])
ok = np.isfinite(A) & (ICE <= float(os.environ.get("H0_ICEMAX", "0.10")))
anom = A - np.nanmean(np.where(ok, A, np.nan), 0)[None]
w = w[ii, jj]; dbin = dbin[ii, jj]; lband = lband[ii, jj]
rows = []
for k, f in enumerate(files):
    r = dict(year=int(f.parent.name), month=int(f.stem), coverage=float(ok[k].mean()))
    vals = []
    for b in np.unique(dbin):
        mk = ok[k] & (dbin == b)
        if mk.sum() > 0:
            vals.append(np.sum(anom[k][mk]*w[mk])/np.sum(w[mk]))
    r["h_cm"] = 100*float(np.mean(vals)) if vals else np.nan
    for L in range(12):
        mk = ok[k] & (lband == L)
        r[f"h_band{L*30:03d}_cm"] = 100*float(np.sum(anom[k][mk]*w[mk])/np.sum(w[mk])) if mk.sum() > 20 else np.nan
    rows.append(r)
t = pd.DataFrame(rows); t["x"] = t.year + (t.month-0.5)/12
t.loc[t.coverage < MINCOV, [c for c in t.columns if c.startswith("h")]] = np.nan
# remove bias and annual + semiannual harmonics, then long-period fit (10, 20, 40 yr) -> h_f
x = t.x.values; y = t.h_cm.values; g = np.isfinite(y)
def design(x, periods):
    cols = [np.ones_like(x)]
    for P in periods:
        cols += [np.cos(2*np.pi*x/P), np.sin(2*np.pi*x/P)]
    return np.column_stack(cols)
Xs = design(x, [1.0, 0.5]); cs = np.linalg.lstsq(Xs[g], y[g], rcond=None)[0]
h = y - Xs @ cs
Xl = design(x, [10.0, 20.0, 40.0]); cl = np.linalg.lstsq(Xl[g], h[g], rcond=None)[0]
hf = Xl @ cl
t["h_noseas_cm"] = h; t["h_f_cm"] = hf
t["T_hogg_model_calibrated_Sv"] = -3.69*(h - hf) - 1.11*hf
t.attrs = {}
out.parent.mkdir(parents=True, exist_ok=True); t.to_csv(out, index=False)
print(f"region points {int(reg.sum())}; months {len(t)}; mean coverage {t.coverage.mean():.2f}; "
      f"months below MINCOV {int((t.coverage < MINCOV).sum())}")
