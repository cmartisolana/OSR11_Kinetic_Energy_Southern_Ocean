"""dk9: half-power positions around Drake Passage (110-20W) in three periods (2001-2008, 2009-2016, 2017-2025).
Baseline half-power (200 cm2 s-2, 25 %, 100 km) on monthly KE of the monthly-mean velocity (Product 1), monthly envelopes
matched to the 2001-2025 mean envelope at each longitude (as a2_run.py / dk6); envelopes in the ACC band (mean ADT -0.6..0.2 m)
present in >= 50 % of months; Malvinas-Current envelopes (62-50W north of 52S) excluded. Per envelope: mean latitude in each
period (>= 50 % of its months valid), 2001-2025 mean and monthly standard deviation; trend over 2001-2025 (OLS, km/decade).

Purpose : half-power positions around Drake Passage (110-20 W) in three periods, per envelope (SI S3 table of envelope positions).
Inputs  : $OSR11_DATA/monthly_allsat (00a), outputs/base/mean_adt.nc
Outputs : outputs/drake/dk9/dk9_hp_periods.csv (distributed)
Run     : python scripts/03_drake/03p_halfpower_periods.py      (import of dk6_common fixed as in 03o)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, glob, numpy as np, pandas as pd, xarray as xr
from osr11.drake_halfpower import F, FI, P, MAXE, detect_all, match  # noqa: F401
from osr11 import paths as PA
OUT = str(PA.dk(9))
BOX = dict(longitude=slice(-112, -18), latitude=slice(-66, -35))
fs = [str(f) for f in PA.monthly_files("allsat")]
KE = []; 
for f in fs:
    with xr.open_dataset(f) as d:
        s = d.sel(**BOX); KE.append((0.5 * (s.u_mean.values ** 2 + s.v_mean.values ** 2)).astype(np.float32)); lat, lon = s.latitude.values.astype(float), s.longitude.values.astype(float)
KE = np.array(KE)
C, TB = detect_all(KE, lat)
cc, ctb = F.detect_field(np.ascontiguousarray(np.nanmean(KE, 0).astype(float)), lat, F.meridional_dists(lat), P["thr"], P["min_w"], P["frac"], MAXE)
meta, arr = match(C, TB, cc, ctb, lat)
adt = xr.open_dataset(PA.MEAN_ADT).mean_adt.sel(**BOX).interp(latitude=lat, longitude=lon).values
M = pd.DataFrame(meta, columns=["j", "c", "lo", "hi", "hpc", "ls", "ln"]); M["j"] = M.j.astype(int); M["lon"] = lon[M.j]
ia = np.clip(np.searchsorted(lat, M.hpc), 0, lat.size - 1); M["adt"] = adt[ia, M.j]
hp = arr[:, :, FI["lat_hp"]]; M["cov"] = np.isfinite(hp).mean(0)
yr = np.repeat(np.arange(2001, 2026), 12); t = 2001 + (np.arange(300) + 0.5) / 12
PER = {"2001-2008": (2001, 2008), "2009-2016": (2009, 2016), "2017-2025": (2017, 2025)}
for k, (a, b) in PER.items():
    m = (yr >= a) & (yr <= b); v = hp[m]; ok = np.isfinite(v).mean(0) >= 0.5
    M[k] = np.where(ok, np.nanmean(v, 0), np.nan)
M["mean"] = np.nanmean(hp, 0); M["sd_km"] = np.nanstd(hp, 0) * 111.195
tr = []
for q in range(hp.shape[1]):
    y = hp[:, q]; g = np.isfinite(y)
    tr.append(np.polyfit(t[g], y[g], 1)[0] * 1111.95 if g.sum() > 60 else np.nan)
M["trend_km_dec"] = tr
keep = (M.adt >= -0.6) & (M.adt <= 0.2) & (M["cov"] >= 0.5) & ~((M.lon >= -62) & (M.lon < -50) & (M.hpc > -52))
M[keep].to_csv(f"{OUT}/dk9_hp_periods.csv", index=False)
print(keep.sum(), "envelopes kept;", M[keep][["2001-2008", "2009-2016", "2017-2025"]].describe().round(2).to_string())
d = M[keep]; print("mean period differences (km): P2-P1", round(111.195 * np.nanmean(d["2009-2016"] - d["2001-2008"]), 1), "P3-P1", round(111.195 * np.nanmean(d["2017-2025"] - d["2001-2008"]), 1), "median monthly SD km", round(d.sd_km.median(), 1))
