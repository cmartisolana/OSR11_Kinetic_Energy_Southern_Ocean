"""dk21: streamline (mean-ADT) profiles of the 2001-2025 trends at the five standing meanders and away from them, Product 1.
Requested 21 Oct 2026: the circumpolar mean-ADT profile (dk17/s5, now SI) at specific locations, with the surface
transport between streamlines (dk18) in place of the along-stream speed.
Regions (longitudes as dk17/s2 and dk18, defined from the bathymetry before any composite): ridge + downstream ranges of
the Kerguelen Plateau, Southeast Indian Ridge, Macquarie Ridge, Pacific-Antarctic Ridge and Drake Passage-Scotia Sea (60-50 W,
the Malvinas Current, excluded), and 'away' = longitudes outside the upstream-to-downstream range of all five.
Pixels: finite in all 300 months, outside the sea-ice mask, mean ADT in [-0.9, +0.5) m (dk17 convention).
Per region, month and 0.1 m bin of the 2001-2025 mean ADT (fixed streamlines):
  EKE = 0.5|U - Ubar|^2 + 0.5 var_intra, KE = 0.5|U|^2 (cm2 s-2; cos-lat area means, as s2);
  T = surface transport between the bin's mean streamlines, longitude mean over the region of the meridional integral of
      the monthly-mean eastward geostrophic velocity over the bin's pixels, sum(u dy) (m2 s-1; as dk18).
Moving bins (monthly ADT minus its domain-mean anomaly, as s2 'iadt') for EKE and KE only: between moving streamlines the
surface transport is fixed by their ADT difference, so T is not computed for them (dk18).
Band totals: T and area-mean EKE/KE over the bins inside the ACC band (mean ADT -0.6..+0.2 m).
Trends: rstats.trend_bcz on the monthly series (deseasonalised Theil-Sen, calibrated AR noise, n_outer=300); % of the mean
per decade. Bins with < 2 % of the region's valid area are dropped.
Usage: python dk21_meander_profiles.py <revision dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from multiprocessing import Pool
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
MEANDERS = {"Kerguelen Plateau": ((55, 65), (65, 80), (80, 95)), "Southeast Indian Ridge": ((120, 130), (130, 140), (140, 150)),
            "Macquarie Ridge": ((150, 155), (155, 160), (160, 175)), "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
            "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35))}
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
madt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values.astype(float)
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
edges = np.round(np.arange(-0.9, 0.5001, 0.1), 3); NB = edges.size - 1
shape = madt.shape
dy = 6371e3 * np.radians(np.abs(np.gradient(lat)))[:, None] * np.ones(lon.size)[None, :]
w = np.cos(np.radians(lat))[:, None] * np.ones(lon.size)[None, :]
# regions by longitude
regmask = {}; anym = np.zeros(lon.size, bool)
for nm, (u_, r_, d_) in MEANDERS.items():
    m = ((lon > r_[0]) & (lon < r_[1])) | ((lon > d_[0]) & (lon < d_[1]))
    if nm.startswith("Drake"): m &= ~((lon > -60) & (lon < -50))
    regmask[nm] = m; anym |= (lon > u_[0]) & (lon < d_[1])
regmask["Away from the five meanders"] = ~anym
REG = list(regmask); NR = len(REG)
sf = OUT / "dk21_series.npz"
if not sf.exists():
    nv = np.zeros(shape); su = np.zeros(shape); sv = np.zeros(shape); sa = np.zeros(shape)
    for f in files:
        with xr.open_dataset(f) as d:
            u = d.u_mean.values.astype(float); v = d.v_mean.values.astype(float); a = d.adt_mean.values.astype(float)
        ok = np.isfinite(u) & np.isfinite(v) & np.isfinite(a); nv += ok
        su += np.where(ok, u, 0); sv += np.where(ok, v, 0); sa += np.where(ok, a, 0)
    valid = (nv == len(files)) & ~ice & np.isfinite(madt) & (madt >= edges[0]) & (madt < edges[-1])
    ub, vb, abar = su / np.maximum(nv, 1), sv / np.maximum(nv, 1), sa / np.maximum(nv, 1)
    V = valid.ravel(); wv = w.ravel()[V]; dyv = dy.ravel()[V]
    LONI = np.broadcast_to(np.arange(lon.size)[None, :], shape).ravel()[V]
    rid = np.full(V.sum(), -1)
    for i, nm in enumerate(REG): rid[regmask[nm][LONI]] = i
    nlon = np.array([regmask[nm].sum() for nm in REG], float)
    labF = np.digitize(madt.ravel()[V], edges) - 1
    mean_dom = np.sum(abar.ravel()[V] * wv) / wv.sum()
    key = lambda r, b: r * NB + b
    nt = len(files); NK = NR * NB
    AREA = np.zeros((nt, NK)); EKE = np.zeros((nt, NK)); KE = np.zeros((nt, NK)); T = np.zeros((nt, NK))
    AREAi = np.zeros((nt, NK)); EKEi = np.zeros((nt, NK)); KEi = np.zeros((nt, NK))
    okF = (rid >= 0) & (labF >= 0) & (labF < NB); kF = key(rid, labF)[okF]
    AF = np.bincount(kF, weights=wv[okF], minlength=NK)
    for t, f in enumerate(files):
        with xr.open_dataset(f) as d:
            u = d.u_mean.values.astype(float).ravel()[V]; v = d.v_mean.values.astype(float).ravel()[V]
            s2 = d.var_intra.values.astype(float).ravel()[V]; a = d.adt_mean.values.astype(float).ravel()[V]
        eke = (0.5 * ((u - ub.ravel()[V]) ** 2 + (v - vb.ravel()[V]) ** 2) + 0.5 * s2) * 1e4
        ke = 0.5 * (u ** 2 + v ** 2) * 1e4
        AREA[t] = AF
        EKE[t] = np.bincount(kF, weights=(wv * eke)[okF], minlength=NK) / np.where(AF > 0, AF, np.nan)
        KE[t] = np.bincount(kF, weights=(wv * ke)[okF], minlength=NK) / np.where(AF > 0, AF, np.nan)
        T[t] = np.bincount(kF, weights=(u * dyv)[okF], minlength=NK) / np.repeat(nlon, NB)
        an = np.sum(a * wv) / wv.sum() - mean_dom
        labI = np.digitize(a - an, edges) - 1; okI = (rid >= 0) & (labI >= 0) & (labI < NB); kI = key(rid, labI)[okI]
        AI = np.bincount(kI, weights=wv[okI], minlength=NK); AREAi[t] = AI
        EKEi[t] = np.bincount(kI, weights=(wv * eke)[okI], minlength=NK) / np.where(AI > 0, AI, np.nan)
        KEi[t] = np.bincount(kI, weights=(wv * ke)[okI], minlength=NK) / np.where(AI > 0, AI, np.nan)
        if t % 60 == 0: print("month", t, flush=True)
    np.savez(sf, edges=edges, reg=np.array(REG), AREA=AREA, EKE=EKE, KE=KE, T=T, AREAi=AREAi, EKEi=EKEi, KEi=KEi)
Z = np.load(sf)
X = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
mid = 0.5 * (edges[:-1] + edges[1:]); inband = (edges[:-1] >= -0.6 - 1e-9) & (edges[1:] <= 0.2 + 1e-9)
jobs = []
for r, nm in enumerate(REG):
    sl = slice(r * NB, (r + 1) * NB); A = Z["AREA"][:, sl].mean(0); frac = A / A.sum()
    for b in range(NB):
        if frac[b] < 0.02: continue
        for v in ("EKE", "KE", "T"):
            jobs.append(dict(region=nm, kind="fixed", var=v, mid=mid[b], area_frac=frac[b], y=Z[v][:, sl][:, b]))
        Ai = Z["AREAi"][:, sl][:, b]
        if np.all(Ai > 0):
            for v in ("EKE", "KE"):
                jobs.append(dict(region=nm, kind="moving", var=v, mid=mid[b], area_frac=frac[b], y=Z[v + "i"][:, sl][:, b]))
    Ab = Z["AREA"][:, sl][:, inband]
    jobs.append(dict(region=nm, kind="band", var="T", mid=np.nan, area_frac=np.nan, y=Z["T"][:, sl][:, inband].sum(1)))
    for v in ("EKE", "KE"):
        jobs.append(dict(region=nm, kind="band", var=v, mid=np.nan, area_frac=np.nan, y=np.nansum(Ab * Z[v][:, sl][:, inband], 1) / Ab.sum(1)))


def run(j):
    y = j.pop("y"); o = R.trend_bcz(X, y, month=MON, n_outer=300); m = float(np.mean(y))
    return dict(**j, mean=m, slope=o["slope"], pct=1000 * o["slope"] / m, pct_lo=1000 * o["ci_lo"] / m, pct_hi=1000 * o["ci_hi"] / m, p=o["p"])


if __name__ == "__main__":
    print(len(jobs), "trend calls", flush=True)
    with Pool(2) as P: rows = P.map(run, jobs, chunksize=4)
    D = pd.DataFrame(rows); D.to_csv(OUT / "dk21_profiles.csv", index=False)
    pd.set_option("display.width", 200)
    print(D[D.kind == "band"].round(3).to_string(index=False))
    print(D[D.kind == "fixed"].pivot_table(index=["region", "mid"], columns="var", values="pct").round(1).to_string())
