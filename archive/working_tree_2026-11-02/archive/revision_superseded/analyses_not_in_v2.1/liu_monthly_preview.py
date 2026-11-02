#!/usr/bin/env python3
"""Monthly-data analogue of Liu et al. (2026) PAR metrics, 2001-2025 (preview, NOT a reproduction).

From monthly-mean fields (extract_monthly.py), for each month m, using months m-1..m+2:
  gradient  G = sum over the 4 months of |grad ADT_month| (m/100 km)
  axis      latitude (60S-48S) of max G at each longitude (150W-80W)   [Liu: modal latitude of
            daily threshold exceedances]
  speed     mean over the 4 months of |U_month| at the axis              [Liu: mean daily speed]
  band      meridional extent of the contiguous band around the axis where G exceeds 20 % of the
            4-month domain maximum (domain 150W-60W, 65S-35S; Liu: 70S-30S)
            [NOT Liu's zero-to-zero occurrence width, which measures day-to-day wandering]
  position, peaks/troughs (Liu Table 1 windows), Tr2-Pk3 amplitude as in liu_reproduce.py
Trends: Liu-style (Sen, MK, their modified MK), Hamed-Rao, and our calibrated test.
Usage: python liu_monthly_preview.py <revision/data dir> <allsat|twosat> <out_dir> [--quick]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import liu_reproduce as LR  # noqa: E402

data, product, outdir = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]); quick = "--quick" in sys.argv
outdir.mkdir(parents=True, exist_ok=True)
keys = [(y, m) for y in range(2001, 2026) for m in range(1, 13)]
Gs, Ss, gmax = [], [], []
for y, m in keys:
    with xr.open_dataset(data / f"monthly_{product}" / f"{y}" / f"{m:02d}.nc") as ds:
        d = ds.sel(longitude=slice(-150, -60), latitude=slice(-65, -35))
        lat = d.latitude.values.astype(float); lon = d.longitude.values.astype(float)
        g = LR_grad = None
        A = d.adt_mean.values.astype(float)
        phi = np.radians(lat)[:, None]
        dy = 6.371e6 * np.radians(np.gradient(lat))[:, None]
        dx = 6.371e6 * np.cos(phi) * np.radians(np.gradient(lon))[None, :]
        g = np.hypot(np.gradient(A, axis=1) / dx, np.gradient(A, axis=0) / dy) * 1e5
        Gs.append(g); Ss.append(np.hypot(d.u_mean.values, d.v_mean.values))
Gs, Ss = np.array(Gs), np.array(Ss)
core_lat = (lat >= -60) & (lat <= -48)
jl = np.nonzero((lon >= -150) & (lon <= -80))[0]
nt = len(keys)
axis = np.full((nt, jl.size), np.nan); band = axis.copy(); speed = axis.copy()
for k in range(1, nt - 2):
    G = np.nansum(Gs[k - 1:k + 3], 0); Sp = np.nanmean(Ss[k - 1:k + 3], 0)
    thr = 0.2 * np.nanmax(G)
    for q, j in enumerate(jl):
        col = np.where(core_lat, G[:, j], -np.inf)
        if not np.isfinite(col).any() or np.nanmax(col) <= 0:
            continue
        i0 = int(np.nanargmax(col)); axis[k, q] = lat[i0]; speed[k, q] = Sp[i0, j]
        above = G[:, j] > thr
        if above[i0]:
            n = i0
            while n + 1 < lat.size and above[n + 1]:
                n += 1
            s = i0
            while s - 1 >= 0 and above[s - 1]:
                s -= 1
            band[k, q] = (lat[n] - lat[s]) * LR.KM_PER_DEG
t = pd.DataFrame(keys, columns=["year", "month"]); M = t.copy(); L = lon[jl]
for s, (a, b) in LR.SECT.items():
    msk = (L >= a) & (L <= b)
    M[f"{s}_position"] = np.nanmean(axis[:, msk], 1)
    M[f"{s}_band_width"] = np.nanmean(band[:, msk], 1)
    M[f"{s}_speed"] = np.nanmean(speed[:, msk], 1)
for nme, (a, b) in LR.PEAKS.items():
    M[nme] = np.nanmin(axis[:, (L >= a) & (L <= b)], 1)
for nme, (a, b) in LR.TROUGHS.items():
    M[nme] = np.nanmax(axis[:, (L >= a) & (L <= b)], 1)
M["amplitude_Tr2_Pk3"] = 0.5 * (M.Tr2 - M.Pk3)
M.to_csv(outdir / f"liu_monthly_metrics_{product}.csv", index=False)
xr.Dataset({"axis_lat": (("time", "lon"), axis), "band_km": (("time", "lon"), band), "speed": (("time", "lon"), speed)},
           coords=dict(time=np.arange(nt), lon=L)).to_netcdf(outdir / f"liu_monthly_axis_{product}.nc")
print(product, "means:", M.drop(columns=["year", "month"]).mean().round(3).to_dict(), flush=True)
T = LR.trends(M, quick, only=["2001-2025"])
T.to_csv(outdir / f"liu_monthly_trends_{product}.csv", index=False)
print("done", product)
