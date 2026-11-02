"""dk18: surface transport proxy between mean streamlines (replaces the along-stream speed in Fig. 3b, d), Product 1.
(1) For each 0.05 m bin of the 2001-2025 mean ADT (-0.9..+0.5 m; pixels valid in all months and ice-free, as dk17/s2):
    T_bin(t) = longitude mean of the meridional integral of the monthly-mean eastward geostrophic velocity over the
    pixels of the bin, sum(u dy) (m2 s-1), as the ACC-band proxy of a8 (Sect. 2.4b). Fixed bins only: between streamlines
    that follow the flow, the surface geostrophic transport is set by the ADT difference (g/f dADT) and is nearly constant
    by construction, so moving bins carry no information on transport.
(2) Away from the five meanders: a8 proxy (ACC band) averaged over the longitudes outside the upstream-to-downstream range
    of all five meanders (dk13 regions), to complete a8b_trends_allsat.csv (meanders: ridge + downstream).
Trends: trend_bcz on monthly series (300 null series), % of the mean per decade.
Usage: python dk18_streamline_transport.py <revision dir> <a8b out dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from multiprocessing import Pool
REV, A8B, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
MEANDERS = {"Kerguelen Plateau": ((55, 65), (65, 80), (80, 95)), "Southeast Indian Ridge": ((120, 130), (130, 140), (140, 150)),
            "Macquarie Ridge": ((150, 155), (155, 160), (160, 175)), "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
            "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35))}
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
edges = np.round(np.arange(-0.9, 0.5001, 0.05), 3); NB = edges.size - 1
dy = 6371e3 * np.radians(np.abs(np.gradient(lat)))[:, None] * np.ones(lon.size)[None, :]
sf = OUT / "dk18_series.npz"
if not sf.exists():
    nv = np.zeros(madt.shape)
    for f in files:
        with xr.open_dataset(f) as d: nv += np.isfinite(d.u_mean.values)
    valid = (nv == len(files)) & ~ice & np.isfinite(madt)
    lab = np.where(valid, np.digitize(madt, edges) - 1, -1); lab[(lab < 0) | (lab >= NB)] = -1
    TB = np.zeros((len(files), NB))
    for k, f in enumerate(files):
        with xr.open_dataset(f) as d: u = d.u_mean.values.astype(float)
        x = np.where(lab >= 0, np.nan_to_num(u) * dy, 0.0)
        for b in range(NB):
            TB[k, b] = np.where(lab == b, x, 0.0).sum() / lon.size      # longitude mean of sum(u dy) (zero where the bin is absent)
        if k % 60 == 0: print("month", k, flush=True)
    np.savez(sf, TB=TB, edges=edges)
TB = np.load(sf)["TB"]
TL = np.load(A8B / "a8b_series.npz")["TL"]
anym = np.zeros(lon.size, bool)
for nm, (u_, r_, d_) in MEANDERS.items(): anym |= (lon > u_[0]) & (lon < d_[1])
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
jobs = [(f"bin|{0.5 * (edges[b] + edges[b + 1]):.3f}", TB[:, b]) for b in range(NB)] + [("Away from the five meanders", np.nanmean(TL[:, ~anym], 1))]


def tr(a):
    key, y = a; r = R.trend_bcz(T, y, month=MON, n_outer=300); mu = float(np.mean(y))
    return dict(key=key, mean=mu, pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"])


if __name__ == "__main__":
    with Pool(2) as P: rows = P.map(tr, jobs)
    D = pd.DataFrame(rows); D.to_csv(OUT / "dk18_trends.csv", index=False); print(D.round(3).to_string(index=False))
