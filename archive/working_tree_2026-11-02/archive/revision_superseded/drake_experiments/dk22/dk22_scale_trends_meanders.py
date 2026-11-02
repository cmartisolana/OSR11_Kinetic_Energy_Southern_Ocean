"""dk22: KE and the monthly eddy part of KE (u' = monthly anomaly from the 2001-2025 mean) by length-scale band in the
upstream, ridge and downstream sectors of the five standing meanders (dk13 regions, ACC band) and away from them,
Product 1, monthly means 2001-2025. Coarse-graining as dk20/dk21. Bands: <100 (below effective resolution), 100-300
(eddy dominated, dk21), 300-1000 and >1000 km (mean-flow dominated). Trends: trend_bcz (300 null series), % of the
band mean per decade, FDR 10 % over all tests of each band.
Usage: python dk22_scale_trends_meanders.py <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy.ndimage import gaussian_filter1d
from multiprocessing import Pool
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
L = [0, 100, 300, 1000]
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.interp(latitude=lat, longitude=lon).values
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
RE = 6371.0; dy = RE * np.radians(np.abs(np.diff(lat).mean())); dx = RE * np.cos(np.radians(lat)) * np.radians(np.abs(np.diff(lon).mean()))
k = np.fft.rfftfreq(lon.size)
MEANDERS = {"Kerguelen Plateau": ((55, 65), (65, 80), (80, 95)), "Southeast Indian Ridge": ((120, 130), (130, 140), (140, 150)),
            "Macquarie Ridge": ((150, 155), (155, 160), (160, 175)), "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
            "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35))}


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
acc = valid & (madt >= -0.6) & (madt <= 0.2); w = np.cos(np.radians(lat))[:, None]
REG = {}; anym = np.zeros(lon.size, bool)
for nm, (u_, r_, d_) in MEANDERS.items():
    anym |= (lon > u_[0]) & (lon < d_[1])
    for part, (a, b) in (("upstream", u_), ("ridge", r_), ("downstream", d_)):
        m = (lon > a) & (lon < b)
        if nm.startswith("Drake") and part == "downstream": m &= ~((lon > -60) & (lon < -50))
        REG[f"{nm}|{part}"] = w * (acc & m[None, :])
REG["Away from the five meanders|-"] = w * (acc & ~anym[None, :])
KEYS = list(REG); WS = np.stack([REG[q] for q in KEYS]); WN = WS.sum(axis=(1, 2))


def month(f):
    with xr.open_dataset(f) as d: u, v = d.u_mean.values.astype(float), d.v_mean.values.astype(float)
    u = np.where(valid, u, 0.0); v = np.where(valid, v, 0.0); up, vp = u - ub, v - vb
    out = np.zeros((2, len(L), len(KEYS)))
    for j, l in enumerate(L):
        for i, (a, b) in enumerate(((u, v), (up, vp))):
            e = 0.5 * (filt(a, l) ** 2 + filt(b, l) ** 2) * 1e4
            out[i, j] = np.tensordot(WS, e, axes=([1, 2], [0, 1])) / WN
    return out


if __name__ == "__main__":
    sf = OUT / "dk22_series.npz"
    if not sf.exists():
        with Pool(2) as P: S = np.array(P.map(month, files, chunksize=4))
        np.savez(sf, S=S, keys=np.array(KEYS), L=np.array(L))
    S = np.load(sf)["S"]
    T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
    BANDS = [("<100", 0, 1), ("100-300", 1, 2), ("300-1000", 2, 3), (">1000", 3, None)]
    jobs = []
    for i, fld in enumerate(("KE", "EKEm")):
        for bn, a, b in BANDS:
            for r, key in enumerate(KEYS):
                y = S[:, i, a, r] - (S[:, i, b, r] if b is not None else 0.0); jobs.append((fld, bn, key, y))

    def tr(j):
        fld, bn, key, y = j; r = R.trend_bcz(T, y, month=MON, n_outer=300); mu = float(y.mean())
        return dict(field=fld, band_km=bn, region=key.split("|")[0], sector=key.split("|")[1], mean=mu, abs_dec=10 * r["slope"],
                    pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"])
    with Pool(2) as P: rows = P.map(tr, jobs)
    D = pd.DataFrame(rows); D["fdr10"] = False
    for (f_, b_), g in D.groupby(["field", "band_km"]): D.loc[g.index, "fdr10"] = R.bh_fdr(g.p.values, 0.10)
    D.to_csv(OUT / "dk22_trends.csv", index=False); pd.set_option("display.width", 200); print(D.round(3).to_string(index=False))
