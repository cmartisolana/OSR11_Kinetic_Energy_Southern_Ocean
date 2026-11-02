"""dk14: circumpolar occurrence of monthly half-power points (Fig. 2b around the ACC; same rules as dk10):
all envelopes detected each month on the KE of the monthly-mean velocity (baseline D1: 200 cm2 s-2, 100 km, 25 %),
points kept where the 2001-2025 mean ADT is in the ACC band (-0.6..0.2 m); Malvinas shelf-break branch (62-56 W north of
52 S) excluded. Density = % of months with a point in each 1 deg lon x 0.5 deg lat cell, per grid longitude of the column.
Saved for 2001-2025, 2001-2015 and 2016-2025. Product 1.
Usage: python dk14_hp_density_circ.py <revision dir> <out dir>"""
import sys, numpy as np, xarray as xr
from pathlib import Path
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import fronts as F
FI = {k: i for i, k in enumerate(F.FIELDS)}; P = dict(thr=0.02, min_w=100.0, frac=0.25); MAXE = 8
fs = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(fs[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.interp(latitude=lat, longitude=lon).values
dist = F.meridional_dists(lat)
le = np.arange(-180, 180.01, 1.0); ae = np.arange(-65, -35.01, 0.5)
H = {k: np.zeros((ae.size - 1, le.size - 1)) for k in ("all", "early", "late")}
NM = {"all": 300, "early": 180, "late": 120}
for k, f in enumerate(fs):
    with xr.open_dataset(f) as d:
        ke = (0.5 * (d.u_mean.values.astype(float) ** 2 + d.v_mean.values.astype(float) ** 2))
    cnt, tab = F.detect_field(ke, lat, dist, P["thr"], P["min_w"], P["frac"], MAXE)
    jj, cc = np.nonzero(np.arange(MAXE)[None, :] < cnt[:, None]); la = tab[jj, cc, FI["lat_hp"]]; lo = lon[jj]
    ok = np.isfinite(la)
    i = np.clip(np.searchsorted(lat, la[ok]), 0, lat.size - 1); a_ok = (adt[i, jj[ok]] >= -0.6) & (adt[i, jj[ok]] <= 0.2)
    la, lo = la[ok][a_ok], lo[ok][a_ok]
    mal = (lo >= -62) & (lo < -56) & (la > -52); la, lo = la[~mal], lo[~mal]
    h, _, _ = np.histogram2d(la, lo, bins=[ae, le])
    H["all"] += h; H["early" if k < 180 else "late"] += h
ncol = np.array([((lon >= a) & (lon < a + 1.0)).sum() for a in le[:-1]])
np.savez(OUT / "dk14_hp_density_circ.npz", lon_edges=le, lat_edges=ae,
         **{k: 100 * v / (NM[k] * np.maximum(ncol, 1))[None, :] for k, v in H.items()})
print("ok")
