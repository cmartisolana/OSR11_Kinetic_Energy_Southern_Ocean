"""dk20: scale-aware KE of the surface geostrophic flow by coarse-graining (Aluie et al. 2018; applied to altimetry by
Storer et al. 2022), Products 1 and 4, monthly-mean fields 2001-2025, ACC band.
Filter: Gaussian kernel with the second moment of a top-hat of width l (sigma = l / sqrt(12)), applied separately to u and v
(Cartesian components in longitude-latitude, an approximation adequate for l <= 1000 km at these latitudes), zonally by FFT
on each latitude row (periodic; sigma in grid points from the local zonal spacing) and meridionally with scipy (mode
'nearest' at 65 and 35 S). Land, ice-masked and missing pixels are set to zero velocity before filtering (the usual
convention: land as zero velocity), so the filtered KE near coasts is reduced.
For each month and l in L: KE_l = 0.5 |u_l|^2 of (a) the full monthly-mean velocity and (b) its anomaly from the 2001-2025
mean (u' = u - ubar: the monthly-mean part of EKE), averaged over the ACC band (mean ADT -0.6..0.2 m, ice-free, as a5).
Band KE between l1 < l2 = KE_l1 - KE_l2 (filtering spectrum); KE at l = 0 is the unfiltered KE.
Trends: trend_bcz on monthly series (300 null series), % of the band mean per decade; absolute in cm2 s-2 per decade.
Usage: python dk20_coarse_grain.py <revision dir> <out dir> <product: allsat|twosat>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy.ndimage import gaussian_filter1d
from multiprocessing import Pool
REV, OUT, PROD = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]; OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
L = [0, 100, 150, 200, 300, 500, 1000]                     # km; 0 = unfiltered
files = [REV / "data" / f"monthly_{PROD}" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.interp(latitude=lat, longitude=lon).values
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.interp(latitude=lat, longitude=lon, method="nearest").values > 0.5
RE = 6371.0; dlat = np.abs(np.diff(lat).mean()); dlon = np.abs(np.diff(lon).mean())
dy = RE * np.radians(dlat); dx = RE * np.cos(np.radians(lat)) * np.radians(dlon)
k = np.fft.rfftfreq(lon.size)                               # cycles per grid point
w = np.cos(np.radians(lat))[:, None]


def filt(a, l):
    if l == 0: return a
    s = l / np.sqrt(12.0)
    A = np.fft.rfft(a, axis=1) * np.exp(-0.5 * (2 * np.pi * k[None, :] * (s / dx)[:, None]) ** 2)
    a = np.fft.irfft(A, n=lon.size, axis=1)
    return gaussian_filter1d(a, s / dy, axis=0, mode="nearest")


sf = OUT / f"dk20_series_{PROD}.npz"
if not sf.exists():
    nv = np.zeros(madt.shape); su = np.zeros(madt.shape); sv = np.zeros(madt.shape)
    for f in files:
        with xr.open_dataset(f) as d: u, v = d.u_mean.values.astype(float), d.v_mean.values.astype(float)
        ok = np.isfinite(u) & np.isfinite(v); nv += ok; su += np.where(ok, u, 0); sv += np.where(ok, v, 0)
    valid = (nv == len(files)) & ~ice & np.isfinite(madt)
    ub, vb = su / np.maximum(nv, 1), sv / np.maximum(nv, 1)
    acc = valid & (madt >= -0.6) & (madt <= 0.2); wa = w * acc; W = wa.sum()

    def month(f):
        with xr.open_dataset(f) as d: u, v = d.u_mean.values.astype(float), d.v_mean.values.astype(float)
        u = np.where(valid, u, 0.0); v = np.where(valid, v, 0.0); up = np.where(valid, u - ub, 0.0); vp = np.where(valid, v - vb, 0.0)
        tot = [float((wa * 0.5 * (filt(u, l) ** 2 + filt(v, l) ** 2)).sum() / W) * 1e4 for l in L]
        ano = [float((wa * 0.5 * (filt(up, l) ** 2 + filt(vp, l) ** 2)).sum() / W) * 1e4 for l in L]
        return tot, ano
    with Pool(2) as P: res = P.map(month, files)
    np.savez(sf, KE=np.array([r[0] for r in res]), EKEm=np.array([r[1] for r in res]), L=np.array(L))
    print("series saved", flush=True)
Z = np.load(sf); KE, EK = Z["KE"], Z["EKEm"]
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
BANDS = [(0, 100), (100, 200), (200, 300), (300, 500), (500, 1000), (1000, None)]
jobs = []
for nm, X in (("KE", KE), ("EKEm", EK)):
    for a, b in BANDS:
        ia = L.index(a); y = X[:, ia] - (X[:, L.index(b)] if b else 0.0)
        jobs.append((nm, f"{a}-{b}" if b else f">{a}", y))
    jobs.append((nm, "all", X[:, 0]))


def tr(j):
    nm, band, y = j; r = R.trend_bcz(T, y, month=MON, n_outer=300); mu = float(y.mean())
    return dict(product=PROD, field=nm, band_km=band, mean=mu, abs_dec=10 * r["slope"], abs_lo=10 * r["ci_lo"], abs_hi=10 * r["ci_hi"],
                pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"])


if __name__ == "__main__":
    with Pool(2) as P: rows = P.map(tr, jobs)
    D = pd.DataFrame(rows); D.to_csv(OUT / f"dk20_trends_{PROD}.csv", index=False); print(D.round(3).to_string(index=False))
