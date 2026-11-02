"""A8b: longitude-resolved surface geostrophic transport proxy across the ACC band (as a8_surface_transport.py:
T(lon, t) = integral over the fixed ACC band of the monthly-mean eastward geostrophic velocity u dy, m2 s-1; band = mean ADT
-0.6..+0.2 m, ice-free in every month), averaged in 5 deg longitude bins; trend_bcz (300 null series), % of the mean per
decade, FDR 10 %. Also the five meander composites (ridge + downstream, dk13 regions) and the circumpolar mean (check
against a8_trends_allsat.csv). Surface proxy per unit depth only; not a volume transport.
Usage: python a8b_transport_bins.py <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from multiprocessing import Pool
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True); sys.path.insert(0, str(REV)); import rstats as R
src = open(Path(__file__).parent / "b5_meander_shape.py").read(); exec(src[src.index("MEANDERS = {"):src.index("files = [")])
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
band = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice
dy = 6371e3 * np.radians(np.abs(np.gradient(lat)))[:, None]
sf = OUT / "a8b_series.npz"
if not sf.exists():
    TL = np.full((300, lon.size), np.nan)
    for k, f in enumerate(files):
        with xr.open_dataset(f) as d: u = d.u_mean.values.astype(float)
        TL[k] = np.nansum(np.where(band, u, 0.0) * dy, 0)
    TL[:, band.sum(0) == 0] = np.nan
    np.savez(sf, TL=TL, lon=lon); print("series saved", flush=True)
TL = np.load(sf)["TL"]
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
EDGES = np.arange(-180, 181, 5.0)
jobs = [("ACC", np.nanmean(TL, 1))]
for a in EDGES[:-1]:
    m = (lon >= a) & (lon < a + 5); jobs.append((f"{a + 2.5:g}", np.nanmean(TL[:, m], 1) if np.isfinite(TL[:, m]).any() else np.full(300, np.nan)))
for nm, (u_, r_, d_) in MEANDERS.items():
    m = ((lon > r_[0]) & (lon < r_[1])) | ((lon > d_[0]) & (lon < d_[1]))
    if nm.startswith("Drake"): m &= ~((lon > -60) & (lon < -50))
    jobs.append((nm, np.nanmean(TL[:, m], 1)))


def tr(args):
    key, y = args
    if not np.isfinite(y).all(): return dict(key=key)
    r = R.trend_bcz(T, y, month=MON, n_outer=300); mu = float(np.mean(y))
    return dict(key=key, mean=mu, pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"])


if __name__ == "__main__":
    with Pool(2) as P: rows = P.map(tr, jobs)
    D = pd.DataFrame(rows); D["lon"] = pd.to_numeric(D.key, errors="coerce")
    mk = D.lon.notna() & D.p.notna(); pv = D.p[mk].values; o = np.argsort(pv); q = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; f = np.zeros(len(pv), bool); f[o] = q <= 0.10; D["fdr10"] = False; D.loc[mk, "fdr10"] = f
    D.to_csv(OUT / "a8b_trends_allsat.csv", index=False); b = D[mk]
    print("bins", len(b), "positive", int((b.pct > 0).sum()), "p<0.05", int((b.p.round(3) < 0.05).sum()), "FDR", int(b.fdr10.sum()))
    print(D[D.lon.isna()].round(3).to_string(index=False)); print("done")
