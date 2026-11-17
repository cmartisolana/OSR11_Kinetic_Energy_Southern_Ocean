"""Step 2: EKE between the mean SAF and PF in 2-degree longitude bins, Product 4 (counterpart of main Figs. 2c and 3d; P1 source:
outputs/drake/dk12 from scripts/02_circumpolar/02i_band_lonbins.py), and the Drake-domain late-versus-early increases (SI S2.3).

Purpose : apply 02i_band_lonbins.py unchanged to the Product 4 monthly fields.
Inputs  : Product 4 monthly fields; outputs/base/clim_allsat.nc (Product 1 mean ADT), outputs/base/bathymetry_so.nc,
          outputs/drake/dk5/dk5_stream_edges.csv (SAF/PF levels), outputs/drake/dk12/* (Product 1, for the check)
Outputs : outputs/p4/p4_bins_series[<tag>].npz  monthly bin series (300 x 180), bin centres, levels
          outputs/p4/p4_bins_trends[<tag>].csv  lon, slope, lo, hi, p, mean (+ pct, pct_lo, pct_hi, depth_mean, depth_min, n_pix, p_aic, p_ar12);
                                        slope/lo/hi in cm2 s-2 per decade
          outputs/p4/p4_bins_annual[<tag>].csv  annual bin means as anomaly in % of each bin's 2001-2025 mean (as dk12_band_annual_allsat.csv)
          outputs/p4/p4_drake_domains[<tag>].csv  upstream 86-68 W, passage 68-58 W, downstream 58-40 W: 2018-2025 over 2001-2017 (%), both products
          outputs/p4/p4_bins_p1check.csv        (with --check-p1) trend test re-run on Product 1 bins from dk12_band_series_allsat.npz vs dk12 csv
Method  : exactly 02i: levels SAF = +0.1048 m, PF = -0.5424 m (psi = -0.5 and -1.125 of dk5_stream_edges.csv); ubar = unweighted mean
          of the monthly velocities; band = pixels whose mean ADT lies between the levels, maximum monthly ice fraction <= 0.10 and
          finite velocity in all 300 months; monthly EKE = 0.5|U_m - ubar|^2 + 0.5 var_intra (cm2 s-2), cos(lat)-weighted mean per
          2-degree bin (bin = floor((lon + 180) / 2)); trend = rstats.trend_bcz(n_outer=150, seed 0) on the monthly bin series.
          Mean ADT used for the levels: default tag '' = Product 1 clim_allsat.nc linearly interpolated to the 1/4 deg grid (what the
          02i code does when it is given another grid: `xr.open_dataset(PA.CLIM).adt.interp(latitude=lat, longitude=lon)`; same
          geographic band as Product 1); tag '_p4adt' (--p4adt) = Product 4's own 2001-2025 mean ADT (sensitivity; as 02a does for
          the ACC band). Drake domains as figures_old/figS2_timing.py: from the annual anomalies f (fraction of the bin mean), the
          domain series is the mean of f over the bins with centres strictly inside the domain, increase = 100 ((1 + mean_late) /
          (1 + mean_early) - 1), early 2001-2017, late 2018-2025.
Run     : python scripts/06_product4/p4_02_lonbins.py [--p4adt] [--check-p1] [--no-trends]
Seeds   : rstats.trend_bcz default seed (0).
"""
import argparse
from multiprocessing import Pool
import p4_common as C
import numpy as np, pandas as pd, xarray as xr

ap = argparse.ArgumentParser(); ap.add_argument("--p4adt", action="store_true"); ap.add_argument("--check-p1", action="store_true")
ap.add_argument("--no-trends", action="store_true")
args = ap.parse_args() if __name__ == "__main__" else ap.parse_args([])
TAG = "_p4adt" if args.p4adt else ""
lat, lon = C.grid()
U, V, S2, A, I = (C.stack(v) for v in ("u_mean", "v_mean", "var_intra", "adt_mean", "ice_frac"))
NT = U.shape[0]
su = np.zeros((lat.size, lon.size)); sv = np.zeros_like(su); n = np.zeros_like(su); icemax = np.zeros_like(su); sa = np.zeros_like(su)
for k in range(NT):
    u, v = np.asarray(U[k], float), np.asarray(V[k], float); ok = np.isfinite(u) & np.isfinite(v)
    su[ok] += u[ok]; sv[ok] += v[ok]; n[ok] += 1; icemax = np.maximum(icemax, np.nan_to_num(np.asarray(I[k], float)))
    sa += np.nan_to_num(np.asarray(A[k], float))
ub, vb = su / np.where(n > 0, n, np.nan), sv / np.where(n > 0, n, np.nan)
e = pd.read_csv(C.PA.dk(5, mkdir=False) / "dk5_stream_edges.csv").set_index("psi").adt_edge
LEV_SAF, LEV_PF = (float(np.interp(p, e.index[::-1], e.values[::-1])) for p in (-0.5, -1.125))
if args.p4adt:
    adt = sa / NT
else:
    adt = xr.open_dataset(C.PA.CLIM).adt.interp(latitude=lat, longitude=lon).values
band = np.isfinite(adt) & (adt >= LEV_PF) & (adt <= LEV_SAF) & (icemax <= 0.10) & (n == NT)
w = np.cos(np.radians(lat))[:, None] * band
BINS = np.arange(-180, 180, 2.0); jb = np.clip(((lon + 180) // 2).astype(int), 0, BINS.size - 1)
Y = np.full((NT, BINS.size), np.nan)
for k in range(NT):
    E = (0.5 * ((np.asarray(U[k], float) - ub) ** 2 + (np.asarray(V[k], float) - vb) ** 2) + 0.5 * np.asarray(S2[k], float)) * 1e4
    ww = np.where(np.isfinite(E), w, 0.0); num = np.nansum(E * ww, 0); den = ww.sum(0)
    dd = np.bincount(jb, den, BINS.size)
    Y[k] = np.bincount(jb, num, BINS.size) / np.where(dd > 0, dd, np.nan)
bat = xr.open_dataset(C.PA.BATHY).elevation.interp(latitude=lat, longitude=lon).values
TT = 2001 + (np.arange(NT) + 0.5) / 12; MO = np.tile(np.arange(1, 13), 25)


def one(b, Yser=None):
    Yser = Y if Yser is None else Yser
    y = Yser[:, b]
    if not np.isfinite(y).all(): return dict(lon=BINS[b] + 1.0)
    r = C.R.trend_bcz(TT, y, month=MO, n_outer=150); mu = float(np.mean(y)); m = (jb == b)
    dep = -bat[:, m][band[:, m]]
    return dict(lon=BINS[b] + 1.0, mean=mu, slope=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"],
                pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu,
                depth_mean=float(np.nanmean(dep)) if dep.size else np.nan, depth_min=float(np.nanmin(dep)) if dep.size else np.nan,
                n_pix=int(band[:, m].sum()), p_aic=r["p_aic"], p_ar12=r["p_ar12"])


def p1_one(b):
    Z = np.load(C.PA.dk(12, mkdir=False) / "dk12_band_series_allsat.npz")
    y = Z["eke"][:, b]; r = C.R.trend_bcz(TT, y, month=MO, n_outer=150)
    return dict(lon=float(Z["lon"][b]), slope_dec=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"])


def domains(Aan, tag):
    lo_ = Aan.columns.astype(float).values; f = Aan.values / 100; yrs = Aan.index.values.astype(int)
    e, l = yrs < 2018, yrs >= 2018; out = []
    for nm, (a, b) in (("upstream 86-68W", (-86, -68)), ("passage 68-58W", (-68, -58)), ("downstream 58-40W", (-58, -40))):
        m = (lo_ > a) & (lo_ < b); s_ = f[:, m].mean(1)
        out.append(dict(product=tag, domain=nm, n_bins=int(m.sum()), pct_increase_2018_2025_vs_2001_2017=100 * ((1 + s_[l].mean()) / (1 + s_[e].mean()) - 1)))
    return out


if __name__ == "__main__":
    yr = np.repeat(np.arange(2001, 2026), 12)
    np.savez(C.HERE / f"p4_bins_series{TAG}.npz", eke=Y, lon=BINS + 1.0, lev_saf=LEV_SAF, lev_pf=LEV_PF, band=band)
    Aan = pd.DataFrame(Y, index=yr).groupby(level=0).mean()
    Aan = (100 * (Aan / Aan.mean() - 1)).set_axis(BINS + 1.0, axis=1); Aan.to_csv(C.HERE / f"p4_bins_annual{TAG}.csv")
    A1 = pd.read_csv(C.PA.dk(12, mkdir=False) / "dk12_band_annual_allsat.csv", index_col=0); A1.index = A1.index.astype(int)
    D = pd.DataFrame(domains(A1, "P1 (dk12_band_annual_allsat)") + domains(Aan, "P4 (levels on P4 mean ADT)" if args.p4adt else "P4 (levels on P1 mean ADT)"))
    D.to_csv(C.HERE / f"p4_drake_domains{TAG}.csv", index=False); print(D.round(1).to_string(index=False))
    print("levels", LEV_SAF, LEV_PF, "band pixels", int(band.sum()), "bins with data", int(np.isfinite(Y).all(0).sum()), flush=True)
    if args.check_p1:
        chk = [47, 30, 167]                      # bins centred at -85, -119 and +155 deg
        with Pool(2) as P:
            got = P.map(p1_one, chk)
        ref = pd.read_csv(C.PA.dk(12, mkdir=False) / "dk12_band_trends_allsat.csv")
        rows = []
        for g in got:
            r = ref[np.isclose(ref.lon, g["lon"])].iloc[0]
            rows.append(dict(lon=g["lon"], slope_here=g["slope_dec"], slope_dk12=r.slope_dec, lo_here=g["lo"], lo_dk12=r.lo,
                             hi_here=g["hi"], hi_dk12=r.hi, p_here=g["p"], p_dk12=r.p))
        pd.DataFrame(rows).to_csv(C.HERE / "p4_bins_p1check.csv", index=False); print(pd.DataFrame(rows).to_string(index=False), flush=True)
    if not args.no_trends:
        with Pool(2) as P:
            rows = P.map(one, range(BINS.size))
        T = pd.DataFrame(rows)
        cols = ["lon", "slope", "lo", "hi", "p", "mean", "pct", "pct_lo", "pct_hi", "depth_mean", "depth_min", "n_pix", "p_aic", "p_ar12"]
        T[cols].to_csv(C.HERE / f"p4_bins_trends{TAG}.csv", index=False)
        ok = T.slope.notna()
        print("bins", int(ok.sum()), "positive", int((T.slope[ok] > 0).sum()), "p<0.05", int((T.p[ok] < 0.05).sum()), "done", flush=True)
