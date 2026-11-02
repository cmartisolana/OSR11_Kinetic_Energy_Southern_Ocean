"""B5e: RMS of the month-to-month velocity anomalies inside the ACC band (north-south v' and east-west u'), from the
monthly-mean geostrophic velocity minus its 2001-2025 mean (as B5d). ACC band: mean ADT -0.6..+0.2 m, ice-free, all
months valid. rms_v = sqrt(cos-lat weighted mean of v'^2) over a region and month (cm/s).
Map: 1 deg cells, 2001-2025 mean rms_v and OLS slope of the deseasonalised monthly rms_v (cm/s per decade; slope only).
Profile: 5 deg longitude bins and the five meander composites (ridge + downstream, dk13 regions), from the B5d series
(V2 = mean v'^2, Q = v'^2/(u'^2+v'^2)): rms_v = sqrt(V2), rms_u = sqrt(V2 (1 - Q) / Q); trend_bcz (300 null series), FDR 10 %.
Intra-monthly variability is not split into u and v in the monthly files, so only month-to-month anomalies are used.
Usage: python b5e_rms.py <revision dir> <b5 out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from multiprocessing import Pool
REV, B5 = Path(sys.argv[1]), Path(sys.argv[2]); sys.path.insert(0, str(REV)); import rstats as R
src = open(Path(__file__).parent / "b5_meander_shape.py").read(); exec(src[src.index("MEANDERS = {"):src.index("files = [")])
files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
with xr.open_dataset(files[0]) as d0: lat, lon = d0.latitude.values.astype(float), d0.longitude.values.astype(float)
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
mf = B5 / "b5e_map.npz"
if not mf.exists():
    Zf = np.load(B5 / "b5_annual_fields.npz"); vb = np.nanmean(Zf["V"], 0).astype(float) * 100; allok = np.all(np.isfinite(Zf["V"]), 0); del Zf
    mean_adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc").mean_adt.values
    ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
    acc = np.isfinite(mean_adt) & (mean_adt >= -0.6) & (mean_adt <= 0.2) & ~ice & allok
    w = np.cos(np.radians(lat))[:, None] * acc; ny = (lat.size // 8) * 8
    blk = lambda a: a[:ny].reshape(ny // 8, 8, lon.size // 8, 8).sum(axis=(1, 3))
    den = blk(w); RMS = np.full((300,) + den.shape, np.nan)
    for k, f in enumerate(files):
        with xr.open_dataset(f) as d: vp = np.nan_to_num(d.v_mean.values.astype(float) * 100 - vb)
        RMS[k] = np.sqrt(blk(vp ** 2 * w) / np.where(den > 0, den, np.nan))
    Rd = RMS - np.stack([np.nanmean(RMS[MON == mm], 0) for mm in range(1, 13)])[MON - 1]; ta = T - T.mean()
    slope = 10 * np.nansum(ta[:, None, None] * Rd, 0) / (ta ** 2).sum(); mean = np.sqrt(np.nanmean(RMS ** 2, 0))
    bad = den <= 0.25 * 64 * np.cos(np.radians(55)); slope[bad] = np.nan; mean[bad] = np.nan
    np.savez(mf, lat=lat[:ny].reshape(-1, 8).mean(1), lon=lon.reshape(-1, 8).mean(1), mean=mean, slope=slope); print("map saved", flush=True)
Z = np.load(B5 / "b5d_series.npz"); BIN = Z["lon"]
def series(v2, q): return np.sqrt(v2), np.sqrt(v2 * (1 - q) / q)


def tr(args):
    key, y = args
    if not np.isfinite(y).all(): return dict(key=key)
    r = R.trend_bcz(T, y, month=MON, n_outer=300); mu = float(np.mean(y))
    return dict(key=key, mean=mu, slope_dec=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"],
                pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu)


if __name__ == "__main__":
    jobs = []
    for b in range(BIN.size):
        rv, ru = series(Z["V2_bin"][:, b], Z["Q_bin"][:, b]); jobs += [(f"rms_v|{b}", rv), (f"rms_u|{b}", ru)]
    for nm in MEANDERS:
        rv, ru = series(Z[f"V2|{nm}"], Z[f"Q|{nm}"]); jobs += [(f"rms_v|{nm}", rv), (f"rms_u|{nm}", ru)]
    with Pool(2) as P: rows = P.map(tr, jobs)
    D = pd.DataFrame(rows); D["var"] = D.key.str.split("|").str[0]; D["where"] = D.key.str.split("|").str[1]
    isbin = D["where"].str.fullmatch(r"\d+"); D["lon"] = np.nan; D.loc[isbin, "lon"] = BIN[D.loc[isbin, "where"].astype(int).values]
    D["fdr10"] = False
    for v in ("rms_v", "rms_u"):
        mk = isbin & (D["var"] == v) & D.p.notna(); pv = D.p[mk].values; o = np.argsort(pv); q = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]; f = np.zeros(len(pv), bool); f[o] = q <= 0.10; D.loc[mk, "fdr10"] = f
        b = D[mk]; print(v, "bins", len(b), "positive", int((b.slope_dec > 0).sum()), "p<0.05", int((b.p.round(3) < 0.05).sum()), "FDR", int(b.fdr10.sum()),
                         "median %/dec", round(b.pct.median(), 2), flush=True)
    D.to_csv(B5 / "b5e_trends.csv", index=False)
    print(D[~isbin][["var", "where", "mean", "slope_dec", "lo", "hi", "p", "pct"]].round(3).to_string(index=False)); print("done")
