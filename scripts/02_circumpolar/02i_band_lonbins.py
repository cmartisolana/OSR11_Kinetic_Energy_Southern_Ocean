"""dk12: circumpolar extension of dk8 (Fig. 1c-d) for the circumpolar figure.
Band = between the mean SAF and PF streamlines of Fig. 1a (central streamlines of the GV23 frontal zones, psi = -0.5 and
-1.125, mapped to the 2001-2025 mean ADT via dk5_stream_edges.csv), pixels ice-free (ice fraction <= 0.10) in every month.
Monthly EKE (Eq. 1: departures of the monthly-mean velocity from the 2001-2025 mean plus the intra-monthly variance),
cos-lat weighted means in 2 deg longitude bins around the ACC (Product 1). Outputs:
  dk12_band_series_allsat.npz  monthly series (300 x 180) and bin centres
  dk12_band_trends_allsat.csv  trend per bin (Theil-Sen on the deseasonalised series, trend_bcz, 150 null series),
                               mean and shallowest depth in the band per bin
  dk12_band_annual_allsat.csv  annual means as anomaly in % of each bin's 2001-2025 mean (Hovmoller)

Purpose : EKE between the mean SAF and PF in 2-deg longitude bins around the ACC, trends per bin (Fig. 3d; 91 % positive bins).
Inputs  : $OSR11_DATA/monthly_allsat (00a), outputs/base/clim_allsat.nc, outputs/base/bathymetry_so.nc,
          outputs/drake/dk5/dk5_stream_edges.csv (03b; distributed, so this step can run before stage 03)
Outputs : outputs/drake/dk12/dk12_band_series_allsat.npz, dk12_band_trends_allsat.csv, dk12_band_annual_allsat.csv (distributed)
Run     : python scripts/02_circumpolar/02i_band_lonbins.py
Seeds   : rstats.trend_bcz default seed (0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, numpy as np, pandas as pd, xarray as xr
from pathlib import Path
from multiprocessing import Pool
from osr11 import rstats, paths as PA
D5, OUT = PA.dk(5), PA.dk(12)
fs = PA.monthly_files("allsat")
assert all(f.exists() for f in fs), "missing monthly files"
with xr.open_dataset(fs[0]) as d0:
    lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
su = np.zeros((lat.size, lon.size)); sv = np.zeros_like(su); n = np.zeros_like(su); icemax = np.zeros_like(su)
for f in fs:                                                   # pass 1: 2001-2025 mean velocity, maximum ice fraction
    with xr.open_dataset(f) as d:
        u, v = d.u_mean.values.astype(float), d.v_mean.values.astype(float); ok = np.isfinite(u) & np.isfinite(v)
        su[ok] += u[ok]; sv[ok] += v[ok]; n[ok] += 1; icemax = np.maximum(icemax, np.nan_to_num(d.ice_frac.values))
ub, vb = su / np.where(n > 0, n, np.nan), sv / np.where(n > 0, n, np.nan)
e = pd.read_csv(D5 / "dk5_stream_edges.csv").set_index("psi").adt_edge
LEV_SAF, LEV_PF = (float(np.interp(p, e.index[::-1], e.values[::-1])) for p in (-0.5, -1.125))
adt = xr.open_dataset(PA.CLIM).adt.interp(latitude=lat, longitude=lon).values
band = np.isfinite(adt) & (adt >= LEV_PF) & (adt <= LEV_SAF) & (icemax <= 0.10) & (n == len(fs))
w = np.cos(np.radians(lat))[:, None] * band
BINS = np.arange(-180, 180, 2.0); jb = np.clip(((lon + 180) // 2).astype(int), 0, BINS.size - 1)
Y = np.full((len(fs), BINS.size), np.nan)
for k, f in enumerate(fs):                                     # pass 2: monthly band-mean EKE per bin
    with xr.open_dataset(f) as d:
        E = (0.5 * ((d.u_mean.values - ub) ** 2 + (d.v_mean.values - vb) ** 2) + 0.5 * d.var_intra.values) * 1e4
    ww = np.where(np.isfinite(E), w, 0.0); num = np.nansum(E * ww, 0); den = ww.sum(0)
    Y[k] = np.bincount(jb, num, BINS.size) / np.where(np.bincount(jb, den, BINS.size) > 0, np.bincount(jb, den, BINS.size), np.nan)
bat = xr.open_dataset(PA.BATHY).elevation.interp(latitude=lat, longitude=lon).values
np.savez(OUT / "dk12_band_series_allsat.npz", eke=Y, lon=BINS + 1.0, lev_saf=LEV_SAF, lev_pf=LEV_PF)
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)


def one(b):
    y = Y[:, b]
    if not np.isfinite(y).all(): return dict(lon=BINS[b] + 1.0)
    r = rstats.trend_bcz(T, y, month=MON, n_outer=150); mu = float(np.mean(y)); m = (jb == b)
    dep = -bat[:, m][band[:, m]]
    return dict(lon=BINS[b] + 1.0, mean=mu, slope_dec=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"],
                pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu,
                depth_mean=float(np.nanmean(dep)) if dep.size else np.nan, depth_min=float(np.nanmin(dep)) if dep.size else np.nan,
                n_pix=int(band[:, m].sum()))


if __name__ == "__main__":
    yr = np.repeat(np.arange(2001, 2026), 12)
    A = pd.DataFrame(Y, index=yr).groupby(level=0).mean()
    (100 * (A / A.mean() - 1)).set_axis(BINS + 1.0, axis=1).to_csv(OUT / "dk12_band_annual_allsat.csv")
    print("levels", LEV_SAF, LEV_PF, "band pixels", int(band.sum()), flush=True)
    with Pool(2) as P:
        rows = P.map(one, range(BINS.size))
    pd.DataFrame(rows).to_csv(OUT / "dk12_band_trends_allsat.csv", index=False); print("done")
