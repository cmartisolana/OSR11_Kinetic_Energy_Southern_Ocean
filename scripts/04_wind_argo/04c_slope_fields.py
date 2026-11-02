"""F4: trends of the isopycnal tilt in and around the ACC from the Roemmich-Gilson Argo fields (candidate Figure 4, 28 Oct 2026).

Input: dk32_iso.npz (plan1019/drake/dk32): annual-mean pressure of sigma0 = 27.3, 27.5, 27.7 kg m-3 on the RG 1 deg grid,
64.5-35.5 S, 2004-2025 (TEOS-10; NaN where the surface outcrops or lies below 1975 dbar), and the 2001-2025 mean ADT
(altimetry, Product 1) averaged onto the same grid.
Domain: ACC band (mean ADT -0.6 to +0.2 m, as in the main text) widened by a 100 km buffer: every 1 deg cell whose centre lies
within 100 km (great circle) of a band pixel of the 1/8 deg mean ADT (so the buffer is measured from the band edge). Cells with any sea-ice-masked 1/8 deg altimetry pixel are removed
(same rule as the altimetric analyses), and a cell enters a surface only if that surface is defined in all 22 years, so each
regional mean uses a fixed set of cells.
Tilt: signed cross-stream slope s = grad(p) . n, with n the unit vector along grad(mean ADT) (pointing to higher ADT, i.e.
equatorward across the streamlines), in dbar per 100 km; positive = the isopycnal deepens equatorward. The unsigned
magnitude |grad p| (the Fig. 2c diagnostic) is computed for comparison. Centred differences on the 1 deg grid (periodic in
longitude).
Trends: Theil-Sen, % of the mean per decade. Regional series (buffered band, sectors, 5 deg longitude bins): rstats.trend_bcz
(annual, AR(1) with iterated bias correction, 300 null series). Maps: the same z statistic per cell, calibrated per
persistence class as in a1_maps_bc.py (10 classes of the bias-corrected AR(1) coefficient, 25 cells x 40 null series each,
kappa = max(q95/1.96, q99/2.576)); BH-FDR at 10 % over the domain cells. Sectors: Atlantic 70 W-20 E, Indian 20-147 E,
Pacific 147 E-70 W.
       ../outputs/sea_ice_mask.nc)

Purpose : signed cross-stream slope of the isopycnals on the ACC band + 100 km domain (annual); domain mask used by Fig. 4.
Inputs  : outputs/drake/dk32/dk32_iso.npz (04a), outputs/base/{sea_ice_mask,mean_adt}.nc
Outputs : outputs/f4/f4_tilt_maps.nc (distributed), f4_regional_trends.csv, f4_regional_series.npz, f4_lon5_trends.csv
Run     : python scripts/04_wind_argo/04c_slope_fields.py [buffer_km=100]     (no raw data needed beyond dk32_iso.npz)
Seeds   : default_rng(0)/(1).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy.stats import norm

from osr11 import rstats as R, paths as PA  # noqa: E402
OUT = PA.out("f4")

BUF_KM = float(sys.argv[1]) if len(sys.argv) > 1 else 100.0
SURF = ("27.3", "27.5", "27.7")
I = np.load(PA.dk(32) / "dk32_iso.npz")
yrs = I["years"].astype(float) + 0.5; lat = I["lat"]; lon0 = I["lon"]
lon = np.where(lon0 > 180, lon0 - 360, lon0); o = np.argsort(lon); lon = lon[o]
adt = I["adt"][:, o]
band = np.isfinite(adt) & (adt >= -0.6) & (adt <= 0.2)
LON, LAT = np.meshgrid(lon, lat)

# sea ice: any masked 1/8 deg pixel inside the 1 deg cell
ice8 = xr.open_dataset(PA.ICE_MASK).ice_mask
icec = np.zeros(band.shape, bool)
v = ice8.values > 0
ii, jj = np.nonzero(v)
ilat = np.searchsorted(lat - 0.5, ice8.latitude.values[ii], side="right") - 1
jlon = np.searchsorted(lon - 0.5, ice8.longitude.values[jj], side="right") - 1
ok = (ilat >= 0) & (ilat < lat.size) & (jlon >= 0) & (jlon < lon.size)
icec[ilat[ok], jlon[ok]] = True


from scipy.spatial import cKDTree
def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


# band at full altimetric resolution (1/8 deg mean ADT, ../outputs/mean_adt.nc): the buffer is measured from its edge
madt = xr.open_dataset(PA.MEAN_ADT).mean_adt
A8 = madt.values; b8 = np.isfinite(A8) & (A8 >= -0.6) & (A8 <= 0.2)
L8, B8 = np.meshgrid(madt.longitude.values, madt.latitude.values)
tree = cKDTree(xyz(B8[b8], L8[b8]))
dch, _ = tree.query(xyz(LAT.ravel(), LON.ravel()))
dmin = (2 * 6371.0 * np.arcsin(np.clip(dch / 2, 0, 1))).reshape(band.shape)
buf = dmin <= BUF_KM + 1e-6
print(f"band cells (1 deg ADT) {band.sum()}, buffered {buf.sum()}, ice-removed {int((buf & icec).sum())}, buffered minus band {int((buf & ~band).sum())}", flush=True)

# gradients (dbar or m per 100 km); lon periodic
KM = 111.195
def grad(F):
    dy = np.gradient(F, lat, axis=-2) / KM * 100
    dx = (np.roll(F, -1, axis=-1) - np.roll(F, 1, axis=-1)) / 2.0 / (KM * np.cos(np.radians(lat))[:, None]) * 100
    return dx, dy


ax_, ay_ = grad(adt); an = np.hypot(ax_, ay_); nx, ny = ax_ / an, ay_ / an
fields = {}
for s in SURF:
    P = I[f"p{s}"][:, :, o]
    dx, dy = grad(P)
    fields[(s, "signed")] = dx * nx + dy * ny
    fields[(s, "mag")] = np.hypot(dx, dy)
    fields[(s, "depth")] = P

w = np.cos(np.radians(LAT))
SECT = {"Atlantic": (LON >= -70) & (LON < 20), "Indian": (LON >= 20) & (LON < 147), "Pacific": (LON >= 147) | (LON < -70)}


def domain(s, base):
    F = fields[(s, "signed")]
    return base & ~icec & np.isfinite(F).all(0) & np.isfinite(fields[(s, "mag")]).all(0)


def regmean(F, m):
    return np.array([np.sum(F[t][m] * w[m]) / np.sum(w[m]) for t in range(F.shape[0])])


rows, series = [], {}
def tr(name, y, x=yrs, meta=None):
    r = R.trend_bcz(x, y, monthly=False, n_outer=300)
    mu = float(np.mean(y))
    rows.append(dict(series=name, mean=mu, slope_dec=10 * r["slope"], lo_dec=10 * r["ci_lo"], hi_dec=10 * r["ci_hi"],
                     pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"],
                     n=r["n"], **(meta or {})))
    print(f"{name}: {rows[-1]['pct']:+.2f} % per decade [{rows[-1]['pct_lo']:+.2f}, {rows[-1]['pct_hi']:+.2f}] p = {r['p']:.3f}", flush=True)


for s in SURF:
    for kind in ("signed", "mag", "depth"):
        F = fields[(s, kind)]
        for dname, base in (("buffer100", buf), ("band", band)):
            m = domain(s, base)
            y = regmean(F, m); series[f"{s}|{kind}|{dname}|ACC"] = y
            tr(f"{s}|{kind}|{dname}|ACC", y, meta=dict(surface=s, kind=kind, domain=dname, region="ACC", ncell=int(m.sum())))
            if dname == "buffer100" and kind != "depth":
                for sec, ms in SECT.items():
                    y = regmean(F, m & ms); series[f"{s}|{kind}|{dname}|{sec}"] = y
                    tr(f"{s}|{kind}|{dname}|{sec}", y, meta=dict(surface=s, kind=kind, domain=dname, region=sec, ncell=int((m & ms).sum())))
        if kind == "signed":   # later start (sparser early Argo sampling)
            m = domain(s, buf); y = regmean(F, m)[3:]
            tr(f"{s}|signed|buffer100|ACC|2007-2025", y, x=yrs[3:], meta=dict(surface=s, kind=kind, domain="buffer100", region="ACC 2007-2025", ncell=int(m.sum())))
pd.DataFrame(rows).to_csv(OUT / "f4_regional_trends.csv", index=False)
np.savez(OUT / "f4_regional_series.npz", years=yrs, **{k.replace("|", "__"): v for k, v in series.items()})

# 5 deg longitude bins, signed slope, all surfaces
lrows = []
edges = np.arange(-180, 181, 5)
for s in SURF:
    F = fields[(s, "signed")]; m0 = domain(s, buf)
    for k in range(edges.size - 1):
        m = m0 & (LON >= edges[k]) & (LON < edges[k + 1])
        if m.sum() < 2:
            continue
        y = regmean(F, m); r = R.trend_bcz(yrs, y, monthly=False, n_outer=300); mu = y.mean()
        lrows.append(dict(surface=s, lon=edges[k] + 2.5, ncell=int(m.sum()), mean=mu, pct=1000 * r["slope"] / mu,
                          lo=1000 * r["ci_lo"] / mu, hi=1000 * r["ci_hi"] / mu, p=r["p"]))
    print(s, "longitude bins done", flush=True)
L5 = pd.DataFrame(lrows)
for s in SURF:
    q = L5.surface == s; L5.loc[q, "fdr10"] = R.bh_fdr(L5.loc[q, "p"].values, 0.10)
L5.to_csv(OUT / "f4_lon5_trends.csv", index=False)

# maps: per-cell z, class-calibrated p, FDR within the domain
def cellmap(F, m):
    idx = np.argwhere(m); n = len(idx)
    b = np.full(n, np.nan); z = np.full(n, np.nan); ph = np.full(n, np.nan); sg = np.full(n, np.nan); mu = np.full(n, np.nan); models = []
    rng = np.random.default_rng(0)
    for q, (i, j) in enumerate(idx):
        y = F[:, i, j]
        bb, a, sd, phis, phi_bc, sig = R._z_stat(yrs, y, False, 200, rng)
        b[q], z[q], ph[q], sg[q], mu[q] = bb, bb / sd, float(phi_bc.sum()), sig, y.mean(); models.append(phi_bc)
    e = np.quantile(ph, np.linspace(0, 1, 11)); e[-1] += 1e-9
    cls = np.clip(np.digitize(ph, e) - 1, 0, 9); kappa = np.zeros(10); rng = np.random.default_rng(1)
    for k in range(10):
        pool = np.nonzero(cls == k)[0]
        pick = rng.choice(pool, min(25, pool.size), replace=False); zz = []
        for q in pick:
            for s_ in R._simulate_ar(models[q], sg[q], yrs.size, 40, rng, burn=500):
                bi, _, sdi, *_ = R._z_stat(yrs, s_, False, 200, rng); zz.append(abs(bi / sdi))
        zz = np.array(zz); kappa[k] = max(np.quantile(zz, 0.95) / 1.96, np.quantile(zz, 0.99) / 2.576)
    p = 2 * norm.sf(np.abs(z) / kappa[cls])
    fdr = R.bh_fdr(p, 0.10).astype(bool)
    out = {k_: np.full(m.shape, np.nan) for k_ in ("pct", "p", "mean", "slope_dec")}
    out["fdr"] = np.zeros(m.shape, bool)
    for q, (i, j) in enumerate(idx):
        out["pct"][i, j] = 1000 * b[q] / mu[q]; out["p"][i, j] = p[q]; out["mean"][i, j] = mu[q]; out["slope_dec"][i, j] = 10 * b[q]
        out["fdr"][i, j] = fdr[q]
    print(f"  cells {n}, p<0.05 {int((p < 0.05).sum())}, FDR {int(fdr.sum())}, positive {int((b > 0).sum())}, kappa {kappa.round(2)}", flush=True)
    return out


DS = {}
for s in SURF:
    for kind in ("signed", "mag"):
        print("map", s, kind, flush=True)
        r = cellmap(fields[(s, kind)], domain(s, buf))
        for k_, v_ in r.items():
            DS[f"{kind}_{s.replace('.', 'p')}_{k_}"] = (("lat", "lon"), v_.astype(float))
for s in SURF:
    DS[f"meanslope_{s.replace('.', 'p')}"] = (("lat", "lon"), np.nanmean(fields[(s, "signed")], 0))
    DS[f"meandepth_{s.replace('.', 'p')}"] = (("lat", "lon"), np.nanmean(fields[(s, "depth")], 0))
DS["band"] = (("lat", "lon"), band.astype(np.int8)); DS["buffer100"] = (("lat", "lon"), buf.astype(np.int8))
DS["ice"] = (("lat", "lon"), icec.astype(np.int8)); DS["adt"] = (("lat", "lon"), adt)
xr.Dataset(DS, coords={"lat": lat, "lon": lon}).to_netcdf(OUT / "f4_tilt_maps.nc")
print("done", flush=True)
