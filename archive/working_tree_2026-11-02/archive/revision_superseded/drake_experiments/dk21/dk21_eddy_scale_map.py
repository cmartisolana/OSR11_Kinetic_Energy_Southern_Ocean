"""dk21: scale at which the surface kinetic energy becomes eddy dominated, mapped (Product 1, monthly means 2001-2025).
Coarse-graining as dk20 (Gaussian with the second moment of a top-hat of width l, land as zero velocity). Because the
filter is linear, u_l = ubar_l + u'_l, and the time mean of the filtered KE splits exactly into a mean-flow part
0.5|ubar_l|^2 and an eddy part <0.5|u'_l|^2> (u' = monthly anomaly from the 2001-2025 mean). For each band
l1 < l2 (0-100, 100-150, 150-200, 200-300, 300-500, 500-1000, >1000 km) the eddy fraction is
f = EKE_band / KE_band, with EKE_band = EKE_l1 - EKE_l2. The eddy-dominance scale l_e is where f falls through 0.5,
interpolated in log l between band centres (geometric means; 50 and 1400 km for the open bands). Maps are smoothed
with a 1 deg Gaussian for display only. Sub-monthly variance is not included (no scale information).
Usage: python dk21_eddy_scale_map.py <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, xarray as xr
from scipy.ndimage import gaussian_filter1d
from multiprocessing import Pool
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
L = [0, 100, 150, 200, 300, 500, 1000]
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.interp(latitude=lat, longitude=lon).values
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
RE = 6371.0; dy = RE * np.radians(np.abs(np.diff(lat).mean())); dx = RE * np.cos(np.radians(lat)) * np.radians(np.abs(np.diff(lon).mean()))
k = np.fft.rfftfreq(lon.size)


def filt(a, l):
    if l == 0: return a
    s = l / np.sqrt(12.0)
    A = np.fft.rfft(a, axis=1) * np.exp(-0.5 * (2 * np.pi * k[None, :] * (s / dx)[:, None]) ** 2)
    return gaussian_filter1d(np.fft.irfft(A, n=lon.size, axis=1), s / dy, axis=0, mode="nearest")


nv = np.zeros(madt.shape); su = np.zeros(madt.shape); sv = np.zeros(madt.shape)
for f in files:
    with xr.open_dataset(f) as d: u, v = d.u_mean.values.astype(float), d.v_mean.values.astype(float)
    ok = np.isfinite(u) & np.isfinite(v); nv += ok; su += np.where(ok, u, 0); sv += np.where(ok, v, 0)
valid = (nv == len(files)) & ~ice & np.isfinite(madt)
ub = np.where(valid, su / np.maximum(nv, 1), 0.0); vb = np.where(valid, sv / np.maximum(nv, 1), 0.0)


def month(f):
    with xr.open_dataset(f) as d: u, v = d.u_mean.values.astype(float), d.v_mean.values.astype(float)
    up = np.where(valid, u, 0.0) - ub; vp = np.where(valid, v, 0.0) - vb
    return np.stack([0.5 * (filt(up, l) ** 2 + filt(vp, l) ** 2) for l in L]).astype(np.float32)


if __name__ == "__main__":
    E = np.zeros((len(L),) + madt.shape)
    with Pool(2) as P:
        for i, e in enumerate(P.imap_unordered(month, files, chunksize=4)):
            E += e
            if i % 50 == 0: print("month", i, flush=True)
    E = 1e4 * E / len(files)
    M = 1e4 * np.stack([0.5 * (filt(ub, l) ** 2 + filt(vb, l) ** 2) for l in L])
    K = M + E
    edges = L + [None]; cent = [50, np.sqrt(100 * 150), np.sqrt(150 * 200), np.sqrt(200 * 300), np.sqrt(300 * 500), np.sqrt(500 * 1000), 1400]
    FB = []
    for i in range(len(L)):
        eb = E[i] - (E[i + 1] if i + 1 < len(L) else 0); kb = K[i] - (K[i + 1] if i + 1 < len(L) else 0)
        FB.append(np.where(kb > 1e-3, eb / np.where(kb > 1e-3, kb, 1), np.nan))
    FB = np.array(FB)
    le = np.full(madt.shape, np.nan); lc = np.log(cent)
    for i in range(len(L) - 1):
        a, b = FB[i], FB[i + 1]; cross = np.isnan(le) & (a >= 0.5) & (b < 0.5)
        t = (a - 0.5) / np.where(a - b != 0, a - b, np.nan)
        le = np.where(cross, np.exp(lc[i] + t * (lc[i + 1] - lc[i])), le)
    le = np.where(valid & np.isnan(le) & (FB[0] < 0.5), 0.0, le)       # mean-flow dominated at all scales
    le[~valid] = np.nan
    np.savez(OUT / "dk21_eddy_scale.npz", lat=lat, lon=lon, L=np.array(L), cent=np.array(cent), E=E.astype(np.float32), M=M.astype(np.float32),
             FB=FB.astype(np.float32), le=le.astype(np.float32), valid=valid)
    acc = valid & (madt >= -0.6) & (madt <= 0.2); w = np.cos(np.radians(lat))[:, None] * acc
    for i in range(len(L)):
        eb = E[i] - (E[i + 1] if i + 1 < len(L) else 0); kb = K[i] - (K[i + 1] if i + 1 < len(L) else 0)
        print(f"band from {L[i]} km: ACC-band eddy fraction {(w * eb).sum() / (w * kb).sum():.3f}")
    q = le[acc]; print("l_e ACC band: median", np.nanmedian(q), "p25", np.nanpercentile(q, 25), "p75", np.nanpercentile(q, 75))
