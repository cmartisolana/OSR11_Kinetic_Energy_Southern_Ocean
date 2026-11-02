"""B5d: north-south flow inside the ACC from MONTHLY velocity anomalies (replaces the annual-mean version B5c), with the
paper's trend test (deseasonalised Theil-Sen, calibrated, trend_bcz; 300 null series).
ACC band: mean ADT -0.6..+0.2 m, ice-free (Sect. 2.2), all months valid. Monthly-mean geostrophic velocity (u, v), its
2001-2025 mean (ub, vb) and anomalies u' = u - ub, v' = v - vb.
(1) Standing north-south flow s = v' sign(vb) on cells with |vb| >= 2 cm/s (standing branches): s > 0 strengthens the
    standing northward/southward branches. Signed, so transient eddies average out. cm/s.
(2) North-south share of the month-to-month eddy energy: q = sum(w v'^2) / sum(w (u'^2 + v'^2)) over all ACC-band cells
    (0.5 = isotropic); and the north-south eddy variance sum(w v'^2)/sum(w) (cm2/s2). From monthly means only: the
    intra-monthly variance in the files is not split into u and v.
Outputs: 1 deg map of the OLS slope of deseasonalised monthly s (slope only, as Fig. 2a); 5 deg longitude profiles and the
five meander composites (ridge + downstream, dk13 regions) of s, q and v'^2 with trend_bcz; FDR 10 % over the bins.
Usage: python b5d_meridional_monthly.py <revision dir> <b5 out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from multiprocessing import Pool
REV, B5 = Path(sys.argv[1]), Path(sys.argv[2]); sys.path.insert(0, str(REV)); import rstats as R
src = open(Path(__file__).parent / "b5_meander_shape.py").read(); exec(src[src.index("MEANDERS = {"):src.index("files = [")])
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
Zf = np.load(B5 / "b5_annual_fields.npz"); ub = np.nanmean(Zf["U"], 0).astype(float) * 100; vb = np.nanmean(Zf["V"], 0).astype(float) * 100
allok = np.all(np.isfinite(Zf["V"]), 0); del Zf
mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
acc = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice & allok
stand = acc & (np.abs(vb) >= 2.0); sg = np.nan_to_num(np.sign(vb))
cw = np.cos(np.radians(lat))[:, None]; wS = cw * stand; wA = cw * acc
EDGES = np.arange(-180, 181, 5.0); NB = 72; jb = np.clip(((lon + 180) // 5).astype(int), 0, NB - 1)
ny = (lat.size // 8) * 8
def blk(a): return a[:ny].reshape(ny // 8, 8, lon.size // 8, 8).sum(axis=(1, 3))
REG = {}
for nm, (u_, r_, d_) in MEANDERS.items():
    m = ((lon > r_[0]) & (lon < r_[1])) | ((lon > d_[0]) & (lon < d_[1]))
    if nm.startswith("Drake"): m &= ~((lon > -60) & (lon < -50))
    REG[nm] = m
NT = len(files)
S_bin = np.full((NT, NB), np.nan); Q_bin = np.full((NT, NB), np.nan); V2_bin = np.full((NT, NB), np.nan)
S_reg = {nm: np.full(NT, np.nan) for nm in MEANDERS}; Q_reg = {nm: np.full(NT, np.nan) for nm in MEANDERS}; V2_reg = {nm: np.full(NT, np.nan) for nm in MEANDERS}
den1 = blk(wS); S_map = np.full((NT,) + den1.shape, np.nan)
dS = np.bincount(jb, wS.sum(0), NB); dA = np.bincount(jb, wA.sum(0), NB)
for k, f in enumerate(files):
    with xr.open_dataset(f) as d: u = d.u_mean.values.astype(float) * 100; v = d.v_mean.values.astype(float) * 100
    up, vp = np.nan_to_num(u - ub), np.nan_to_num(v - vb); s = vp * sg
    S_map[k] = blk(s * wS) / np.where(den1 > 0, den1, np.nan)
    nS = np.bincount(jb, (s * wS).sum(0), NB); nV = np.bincount(jb, (vp ** 2 * wA).sum(0), NB); nU = np.bincount(jb, (up ** 2 * wA).sum(0), NB)
    S_bin[k] = nS / np.where(dS > 0, dS, np.nan); V2_bin[k] = nV / np.where(dA > 0, dA, np.nan); Q_bin[k] = nV / np.where(nU + nV > 0, nU + nV, np.nan)
    for nm, m in REG.items():
        S_reg[nm][k] = (s * wS)[:, m].sum() / wS[:, m].sum(); a = (vp ** 2 * wA)[:, m].sum(); b = (up ** 2 * wA)[:, m].sum()
        V2_reg[nm][k] = a / wA[:, m].sum(); Q_reg[nm][k] = a / (a + b)
    if k % 50 == 0: print("month", k, flush=True)
# map: OLS slope of deseasonalised monthly s (cm/s per decade)
T = 2001 + (np.arange(NT) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
Sd = S_map - np.stack([np.nanmean(S_map[MON == mm], 0) for mm in range(1, 13)])[MON - 1]
ta = T - T.mean(); slope = 10 * np.nansum(ta[:, None, None] * Sd, 0) / (ta ** 2).sum(); slope[den1 <= 0.25 * 64 * np.cos(np.radians(55))] = np.nan
np.savez(B5 / "b5d_map.npz", lat=lat[:ny].reshape(-1, 8).mean(1), lon=lon.reshape(-1, 8).mean(1), slope=slope,
         vmean=blk(vb * wS) / np.where(den1 > 0, den1, np.nan))
np.savez(B5 / "b5d_series.npz", S_bin=S_bin, Q_bin=Q_bin, V2_bin=V2_bin, lon=EDGES[:-1] + 2.5,
         **{f"S|{nm}": S_reg[nm] for nm in MEANDERS}, **{f"Q|{nm}": Q_reg[nm] for nm in MEANDERS}, **{f"V2|{nm}": V2_reg[nm] for nm in MEANDERS})


def tr(args):
    key, y, absolute = args
    if not np.isfinite(y).all(): return dict(key=key)
    r = R.trend_bcz(T, y, month=MON, n_outer=300); mu = float(np.mean(y))
    o = dict(key=key, mean=mu, slope_dec=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"])
    if not absolute: o.update(pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu)
    return o


if __name__ == "__main__":
    jobs = [(f"S|{b}", S_bin[:, b], True) for b in range(NB)] + [(f"Q|{b}", Q_bin[:, b], False) for b in range(NB)] + [(f"V2|{b}", V2_bin[:, b], False) for b in range(NB)]
    jobs += [(f"S|{nm}", S_reg[nm], True) for nm in MEANDERS] + [(f"Q|{nm}", Q_reg[nm], False) for nm in MEANDERS] + [(f"V2|{nm}", V2_reg[nm], False) for nm in MEANDERS]
    import os
    only = os.environ.get("B5D_ONLY")                      # e.g. "S": recompute only that variable, keep the others from the csv
    if only: jobs = [j for j in jobs if j[0].split("|")[0] == only]
    with Pool(2) as P: rows = P.map(tr, jobs)
    if only and (B5 / "b5d_trends.csv").exists():
        old = pd.read_csv(B5 / "b5d_trends.csv"); old = old[old.key.str.split("|").str[0] != only]
        rows = old.drop(columns=[c for c in ("var", "where", "lon", "fdr10") if c in old]).to_dict("records") + rows
    D = pd.DataFrame(rows); D["var"] = D.key.str.split("|").str[0]; D["where"] = D.key.str.split("|").str[1]
    isbin = D["where"].str.fullmatch(r"\d+"); D["lon"] = np.where(isbin, EDGES[:-1][pd.to_numeric(D["where"], errors="coerce").fillna(0).astype(int)] + 2.5, np.nan)
    D["fdr10"] = False
    for v in ("S", "Q", "V2"):
        mk = isbin & (D["var"] == v) & D.p.notna(); pv = D.p[mk].values; o = np.argsort(pv); q = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]; f = np.zeros(len(pv), bool); f[o] = q <= 0.10; D.loc[mk, "fdr10"] = f
        b = D[mk]; print(v, "bins", len(b), "positive", int((b.slope_dec > 0).sum()), "p<0.05", int((b.p.round(3) < 0.05).sum()), "FDR", int(b.fdr10.sum()), flush=True)
    D.to_csv(B5 / "b5d_trends.csv", index=False)
    print(D[~isbin].drop(columns=["key", "lon", "fdr10"]).round(4).to_string(index=False)); print("done")
