"""Half-power baseline (200 cm2 s-2, 100 km, 25 %) on the two-satellite product (1/4 deg), same matching and
trend test as a2_run / retest; compared with the all-satellite tracks.

Purpose : half-power baseline on the two-satellite product and comparison with the all-satellite tracks (Fig. S6).
Inputs  : $OSR11_DATA/monthly_half-power_points_two_sat/YYYY/MM.nc (01c), outputs/base/mean_adt.nc,
          outputs/a2/tracks_baseline.npz (02c), $OSR11_CACHE/lon.npy (02b/02c)
Outputs : outputs/a2/tracks_baseline_twosat.npz, outputs/a2/a2h_twosat_vs_allsat.csv
Run     : python scripts/02_circumpolar/02e_halfpower_twosat.py
Seeds   : rstats.trend_bcz default seed (0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd, xarray as xr
from osr11 import halfpower as A, fronts as F, rstats as R, paths as PA  # noqa: F401
base = str(PA.HP_POINTS_TWOSAT); A.OUT.mkdir(parents=True, exist_ok=True)
ke = np.empty((300, 120, 1440), np.float32); k = 0
for y in range(2001, 2026):
    for m in range(1, 13):
        with xr.open_dataset(f"{base}/{y}/{m:02d}.nc") as ds:
            ke[k] = ds.ke.values
            if k == 0: lat = ds.latitude.values.astype(float); lon = ds.longitude.values.astype(float)
        k += 1
adt = xr.open_dataset(PA.MEAN_ADT)["mean_adt"].interp(latitude=lat, longitude=lon, method="nearest").values
acc = np.isfinite(adt) & (adt >= -0.6) & (adt <= 0.2)
kemean = np.nanmean(ke.astype(np.float64), 0); clim_max = np.nanmax(kemean, 0); P = A.BASE
counts, tab = A.run_set(ke, lat, lon, P, clim_max)
cc, ct = A.clim_envelopes(kemean, lat, lon, P, clim_max)
arr, meta = A.matched_tracks(counts, tab, cc, ct, lat, lon, acc)
keep = meta.keep.values; print("two-sat matched tracks kept", int(keep.sum()), "of", len(meta), f"({keep.sum()/lon.size:.0%} of longitudes)")
np.savez_compressed(A.OUT / "tracks_baseline_twosat.npz", **arr, j=meta.j.values, keep=keep, lon=lon)
x = 2001 + (np.arange(300) + .5) / 12; month = np.tile(np.arange(1, 13), 25); K = 111.195
def series(L, PK, lo, sel):
    y = np.nanmean(L[:, sel] - np.nanmean(L[:, sel], 0, keepdims=True), 1); pk = np.nanmean(PK[:, sel], 1); return y, pk
def report(tag, L, PK, lo):
    sec = {"ACC": np.ones(lo.size, bool), "Atlantic": (lo >= -70) & (lo < 20), "Indian": (lo >= 20) & (lo < 147), "Pacific": (lo >= 147) | (lo < -70)}
    rows = []
    for nm, sel in sec.items():
        y, pk = series(L, PK, lo, sel)
        o = R.trend_bcz(x, y, month=month, n_outer=300); q = R.trend_bcz(x, pk, month=month, n_outer=300); mp = np.nanmean(pk)
        rows.append(dict(product=tag, region=nm, n_lon=int(sel.sum()), pos_km_dec=o["slope"]*K*10, pos_lo=o["ci_lo"]*K*10, pos_hi=o["ci_hi"]*K*10, pos_p=o["p"],
                         peak_mean=mp*1e4, peak_pct_dec=100*q["slope"]*10/mp, peak_lo=100*q["ci_lo"]*10/mp, peak_hi=100*q["ci_hi"]*10/mp, peak_p=q["p"]))
        r = rows[-1]; print(f"{tag:8s} {nm:9s} n={r['n_lon']:4d}  pos {r['pos_km_dec']:+5.1f} [{r['pos_lo']:+5.1f},{r['pos_hi']:+5.1f}] p={r['pos_p']:.2f} | peak KE {r['peak_mean']:5.0f}  {r['peak_pct_dec']:+5.1f} %/dec [{r['peak_lo']:+5.1f},{r['peak_hi']:+5.1f}] p={r['peak_p']:.3f}", flush=True)
    return rows
lo2 = lon[meta.j.values[keep]]
rows = report("two-sat", arr["lat_hp"][:, keep].astype(float), arr["peak_val"][:, keep].astype(float), lo2)
d = np.load(A.OUT / "tracks_baseline.npz", allow_pickle=True); lon1 = np.load(PA.CACHE / "lon.npy"); k1 = d["keep"] & d["acc"]
rows += report("all-sat", d["lat_hp"][:, k1].astype(float), d["peak_val"][:, k1].astype(float), lon1[d["j"][k1]])
pd.DataFrame(rows).to_csv(A.OUT / "a2h_twosat_vs_allsat.csv", index=False)
