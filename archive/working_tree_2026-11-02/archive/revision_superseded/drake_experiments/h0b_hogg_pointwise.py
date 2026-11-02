"""H0b: Antarctic-slope sea-level transport index, Hogg et al. (2015, Sect. 2.3) treatment of sea-ice gaps.
Same region, depth bins and conversion as h0_antarctic_slope_index.py (unchanged); differences, following Hogg et al.:
  "A bias ... and a seasonal cycle were removed from the time series at each point as part of the averaging process
   ... to reduce artifacts due to the incomplete sampling in regions of seasonal sea ice cover."
  -> per point: anomaly relative to that point's calendar-month mean (bias + seasonal cycle), computed from ice-free
     samples only (ice_frac <= 0.10; DUACS vNov2024 fills ADT under ice, so ice points are masked, as in older AVISO);
     a point-month climatology needs >= MINYRS valid years, otherwise that point is unused in that calendar month.
  -> no 50 % coverage threshold (Hogg et al. average whatever is sampled); a month is kept if it has >= MINPTS valid
     points and >= half of the 100-m depth bins populated. Coverage and number of populated 30-deg bands are written.
Then, as in h0: equal-weight mean of depth bins; annual + semiannual harmonics removed from the index; h_f from a fit of
10/20/40-yr sines/cosines; T = -3.69 (h - h_f) - 1.11 h_f  [Sv, h in cm]  (signs checked against the Hogg et al. 2015
PDF; factors model-derived, Hughes et al. 2014, "preliminary").
Optional (20 Oct 2026): env H0_GMSL=<csv with year,month,gmsl_noGIA_mm> subtracts global-mean sea level (no GIA
correction, as DUACS ADT) from ADT at every point and month before the anomalies, so h is slope sea level relative to
the global mean (the Hughes et al. 2014 factors come from Boussinesq models without GMSL rise).
Usage: python h0b_hogg_pointwise.py <monthly dir> <bathymetry nc> <out csv> [y0 y1 LAT_MAX DEP0 DEP1 MINYRS MINPTS]"""
import sys, os
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy import ndimage
mdir, bat_f, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
y0, y1 = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (2001, 2025)
LAT_MAX, DEP0, DEP1, MINYRS, MINPTS = ([float(a) for a in sys.argv[6:11]] if len(sys.argv) > 10
                                       else [-60.0, 1000.0, 3000.0, 5, 100])
files = [f for f in (mdir/f"{y}"/f"{m:02d}.nc" for y in range(y0, y1+1) for m in range(1, 13)) if f.exists()]
with xr.open_dataset(files[0]) as d0:
    lat, lon = d0.latitude.values, d0.longitude.values
dep = -xr.open_dataset(bat_f).elevation.interp(latitude=xr.DataArray(lat, dims="latitude"),
                                                longitude=xr.DataArray(lon, dims="longitude")).values
dep = np.nan_to_num(dep, nan=-1.0)
lband = (np.mod(lon, 360)//30).astype(int)[None, :]*np.ones(dep.shape, int)
w = np.cos(np.radians(lat))[:, None]*np.ones(dep.shape)
reg = (lat[:, None] <= LAT_MAX) & (dep >= DEP0) & (dep <= DEP1)
if os.environ.get("H0_CONTIGUOUS", "1") != "0":          # identical to h0
    shallow = (dep < DEP1) & (lat[:, None] <= LAT_MAX + 5)
    lab, n = ndimage.label(shallow, structure=np.ones((3, 3), int))
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
    lab = np.array([find(a) for a in range(n + 1)])[lab]
    srow = 0 if lat[0] < lat[-1] else -1
    keep = np.unique(lab[srow][lab[srow] > 0]); reg = reg & np.isin(lab, keep) & (lab > 0)
nb = int((DEP1-DEP0)//100)
dbin = np.clip(((dep - DEP0)//100).astype(int), 0, nb - 1)
ii, jj = np.nonzero(reg)
A = np.empty((len(files), ii.size), np.float32); ICE = np.empty_like(A)
for k, f in enumerate(files):
    with xr.open_dataset(f) as d:
        A[k] = d.adt_mean.values[ii, jj]; ICE[k] = np.nan_to_num(d.ice_frac.values[ii, jj], nan=1.0)
if os.environ.get("H0_GMSL"):
    G = pd.read_csv(os.environ["H0_GMSL"], comment="#").set_index(["year", "month"]).gmsl_noGIA_mm
    A -= np.array([G.loc[(int(f.parent.name), int(f.stem))] for f in files], np.float32)[:, None]/1000.0
    print("GMSL removed:", os.environ["H0_GMSL"])
ok = np.isfinite(A) & (ICE <= 0.10)
Av = np.where(ok, A, np.nan)
mon = np.array([int(f.stem) for f in files])
anom = np.full_like(Av, np.nan)
for m in range(1, 13):                                   # per-point bias + seasonal cycle (calendar-month mean)
    s = mon == m
    nv = np.isfinite(Av[s]).sum(0)
    clim = np.where(nv >= MINYRS, np.nanmean(Av[s], 0), np.nan)
    anom[s] = Av[s] - clim[None]
use = np.isfinite(anom)
w = w[ii, jj]; dbin = dbin[ii, jj]; lband = lband[ii, jj]
rows = []
for k, f in enumerate(files):
    r = dict(year=int(f.parent.name), month=int(f.stem), coverage_icefree=float(ok[k].mean()),
             coverage_used=float(use[k].mean()), n_points=int(use[k].sum()),
             n_bands_30deg=int(len(np.unique(lband[use[k]]))))
    vals = [np.sum(anom[k][mk]*w[mk])/np.sum(w[mk]) for b in range(nb) for mk in [use[k] & (dbin == b)] if mk.sum() > 0]
    r["n_depth_bins"] = len(vals)
    r["h_cm"] = 100*float(np.mean(vals)) if (len(vals) >= nb/2 and use[k].sum() >= MINPTS) else np.nan
    for L in range(12):
        mk = use[k] & (lband == L)
        r[f"h_band{L*30:03d}_cm"] = 100*float(np.sum(anom[k][mk]*w[mk])/np.sum(w[mk])) if mk.sum() > 20 else np.nan
    rows.append(r)
t = pd.DataFrame(rows); t["x"] = t.year + (t.month-0.5)/12
x = t.x.values; y = t.h_cm.values; g = np.isfinite(y)
def design(x, P):
    return np.column_stack([np.ones_like(x)] + [f(2*np.pi*x/p) for p in P for f in (np.cos, np.sin)])
Xs = design(x, [1.0, 0.5]); h = y - Xs @ np.linalg.lstsq(Xs[g], y[g], rcond=None)[0]
Xl = design(x, [10.0, 20.0, 40.0]); hf = Xl @ np.linalg.lstsq(Xl[g], h[g], rcond=None)[0]
t["h_noseas_cm"] = h; t["h_f_cm"] = np.where(g, hf, np.nan)
t["T_hogg_model_calibrated_Sv"] = -3.69*(h - hf) - 1.11*hf
t["coverage"] = t.coverage_used                          # column name expected by h2_hogg_vs_insitu.py
out.parent.mkdir(parents=True, exist_ok=True); t.to_csv(out, index=False)
print(f"region points {ii.size}; months {len(t)}, with index {int(g.sum())}; coverage used: mean {t.coverage_used.mean():.2f}, "
      f"Jan-Mar {t[t.month<=3].coverage_used.mean():.2f}, Jul-Sep {t[(t.month>=7)&(t.month<=9)].coverage_used.mean():.2f}; "
      f"median 30-deg bands sampled Jul-Sep {t[(t.month>=7)&(t.month<=9)].n_bands_30deg.median():.0f}/12")
