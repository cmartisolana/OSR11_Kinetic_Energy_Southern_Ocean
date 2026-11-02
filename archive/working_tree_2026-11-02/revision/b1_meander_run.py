"""B1: monthly pathway and energy series for further standing meanders (same metrics as b1_pilot.py).

Sections were fixed on 2026-10-08 from time-mean fields only (b1_regions.py: mean speed,
mean EKE, mean ADT, GEBCO 3000 m isobath), before any time series was computed:

SEIR150  Southeast Indian Ridge crossing, the meander TNG14 analyse at 150-155E in OFES
         upstream    125-138E  (zonal flow over the deep basin, before the ridge crossing)
         meander     140-150E  (ridge crossing and southward bend)
         downstream  150-158E  (eastern half of the EKE maximum, west of Macquarie Ridge)
         latitude band 60-47S; meander box 60-46S
PAR      Pacific-Antarctic Ridge, identical to b1_pilot.py (Liu et al. 2026 sections), rerun here
         so that all meanders share one code path and the ice rule below
KERG     Fawn Trough branch of the Kerguelen Plateau
         upstream    70-77E, 58.5-54S (approach to Fawn Trough)
         meander     76-84E  (trough exit and northward turn along the plateau's eastern flank)
         downstream  84-94E  (EKE maximum east of the plateau)
         latitude band 57-48S; meander box 58-48S (excludes the SAF branch north of 47.5S)

EKE boxes exclude pixels with monthly ice fraction > 0.10 in any month (same rule as a5_energy.py).
Pathway: relative-ADT contour at the climatological jet-core level of the meander box
(b1_meander.core_level), spanning >= 90 % of the box width (gate: >= 90 % of months).
Usage: python b1_meander_run.py <revision/data dir> <PAR|SEIR150|KERG> <out_dir>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b1_meander as B  # noqa: E402

CONF = {
    "SEIR150": dict(reg=(120.0, 162.0, -62.0, -44.0), up=(125.0, 138.0), mea=(140.0, 150.0), down=(150.0, 158.0),
                    band=(-60.0, -47.0), up_band=(-60.0, -47.0), mbox_lat=(-60.0, -46.0)),
    "PAR": dict(reg=(-152.0, -78.0, -62.0, -45.0), up=(-150.0, -135.13), mea=(-135.13, -109.38), down=(-109.38, -80.0),
                band=(-60.0, -48.0), up_band=(-60.0, -48.0), mbox_lat=(-62.0, -45.0)),
    "KERG": dict(reg=(66.0, 98.0, -62.0, -44.0), up=(70.0, 77.0), mea=(76.0, 84.0), down=(84.0, 94.0),
                 band=(-57.0, -48.0), up_band=(-58.5, -54.0), mbox_lat=(-58.0, -48.0)),
}
data, name, outdir = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
c = CONF[name]; outdir.mkdir(parents=True, exist_ok=True)
files = [data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
t["x"] = t.year + (t.month - 0.5) / 12
lo0, lo1, la0, la1 = c["reg"]
stack = {v: [] for v in ("u_mean", "v_mean", "adt_mean", "var_intra", "ice_frac")}
for f in files:
    with xr.open_dataset(f) as ds:
        s = ds.sel(longitude=slice(lo0, lo1), latitude=slice(la0, la1))
        for v in stack:
            stack[v].append(s[v].values.astype(np.float64))
lat = s.latitude.values.astype(float); lon = s.longitude.values.astype(float)
U, V, A, S2, ICE = (np.array(stack[v]) for v in ("u_mean", "v_mean", "adt_mean", "var_intra", "ice_frac"))
ub, vb, Ab = np.nanmean(U, 0), np.nanmean(V, 0), np.nanmean(A, 0)
mbox = dict(lon0=c["mea"][0], lon1=c["mea"][1], lat0=c["mbox_lat"][0], lat1=c["mbox_lat"][1])
core, core_spd = B.core_level(Ab, ub, vb, lat, lon, mbox)
print(name, "core level (relative ADT, m):", core, "speed", round(core_spd, 3), flush=True)
w = np.cos(np.radians(lat))[:, None] * np.ones((1, lon.size))
bm = (lat >= c["band"][0]) & (lat <= c["band"][1]); ubm = (lat >= c["up_band"][0]) & (lat <= c["up_band"][1])
lm = lambda a: (lon >= a[0]) & (lon <= a[1])
dm = bm[:, None] & lm(c["down"])[None, :]; mm = bm[:, None] & lm(c["mea"])[None, :]
ice_any = np.nanmax(np.nan_to_num(ICE), 0) > 0.10
print("EKE-box pixels dropped for ice:", int((ice_any & (dm | mm)).sum()), "of", int((dm | mm).sum()), flush=True)
dm, mm = dm & ~ice_any, mm & ~ice_any
rows = []
for k in range(len(files)):
    rec = {"i": k}
    rel, la, lo = B.box_relative(A[k], lat, lon, mbox)
    for lname, lev in (("core", core), ("core-0.10", core - 0.10), ("core+0.10", core + 0.10)):
        line = B.spanning_contour(rel, la, lo, lev)
        for sm in (100.0, 150.0):
            m = None if line is None else B.path_metrics(line, smooth=sm)
            for key in ("length_km", "sinuosity", "amp_km", "rms_kappa", "int_abs_kappa"):
                rec[f"{lname}_s{int(sm)}_{key}"] = np.nan if m is None else m[key]
    spd = np.hypot(U[k], V[k])
    rec["up_core_speed"] = float(np.nanmean(np.nanmax(spd[np.ix_(ubm, lm(c["up"]))], axis=0)))
    rec["mea_core_speed"] = float(np.nanmean(np.nanmax(spd[np.ix_(bm, lm(c["mea"]))], axis=0)))
    eke = 0.5 * ((U[k] - ub) ** 2 + (V[k] - vb) ** 2) + 0.5 * S2[k]
    for nm, msk in (("down_eke", dm), ("mea_eke", mm)):
        ok = msk & np.isfinite(eke)
        rec[nm] = float(1e4 * np.sum(eke[ok] * w[ok]) / np.sum(w[ok]))
    rows.append(rec)
P = pd.concat([t, pd.DataFrame(rows)], axis=1)
P.to_csv(outdir / f"b1_series_{name}.csv", index=False)
cov = {cc: round(float(np.isfinite(P[cc]).mean()), 3) for cc in P.columns if cc.endswith("_s150_length_km")}
print(name, "pathway coverage (s150):", cov)
