"""x6b (5 Nov 2026): 2-degree longitude bins of the relative-threshold half-power tracks (outputs/x6/x6_relthr_tracks.npz, ACC
envelopes retained in >= 90 % of months): bin-mean anomaly series of the half-power latitude and of the fixed ADT contour through
its mean position, Theil-Sen trends with the calibrated monthly test; meander amplitude by period and trend of the annual amplitude;
also the circumpolar-mean series and their trends. Outputs: outputs/x6/x6_relthr_bins.csv, x6_relthr_circumpolar.csv

Origin  : x6b_relthr_bins.py (5 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : position trends and meander amplitude by 2-degree longitude bin (Sect. 3.2; Fig. S7; SI S4.3).
Inputs  : outputs/x6/x6_relthr_tracks.npz (02n)
Outputs : outputs/x6/x6_relthr_bins.csv, x6_relthr_circumpolar.csv (distributed)
Run     : python scripts/02_circumpolar/02o_pathway_bins.py   (no raw data; about 1 h on 2 cores)
Seeds   : rstats.trend_bcz defaults.
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import os
import numpy as np, pandas as pd
from multiprocessing import Pool
from osr11 import rstats as R, paths as PA
OUT = PA.out("x6")
KM = 111.195; Z = np.load(OUT / "x6_relthr_tracks.npz"); keep = Z["keep"]; lon = Z["lon"][keep]
t = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25); yr = np.repeat(np.arange(2001, 2026), 12); YRS = np.arange(2001, 2026)
hp = (Z["lat_hp"][:, keep] - np.nanmean(Z["lat_hp"][:, keep], 0)) * KM
lcraw = Z["lat_c"][:, keep]; okc = np.isfinite(lcraw).sum(0) >= 200; lc = (lcraw - np.nanmean(lcraw, 0)) * KM; lc[:, ~okc] = np.nan
edges = np.arange(-180, 181, 2)
def fit(y, monthly=True):
    if np.isfinite(y).sum() < (200 if monthly else 20): return dict(slope=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan)
    return R.trend_bcz(t if monthly else YRS.astype(float), y, month=month if monthly else None, monthly=monthly)
def one(b):
    a, bb = edges[b], edges[b + 1]; m = (lon >= a) & (lon < bb)
    if m.sum() == 0: return None
    row = dict(lon=a + 1, n_tracks=int(m.sum()), lat_mean=float(np.nanmean(Z["lat_hp_clim"][keep][m])))
    for key, arr in (("hp", hp), ("contour", lc)):
        r = fit(np.nanmean(arr[:, m], 1)); row.update({f"{key}_slope": 10 * r["slope"], f"{key}_lo": 10 * r["ci_lo"], f"{key}_hi": 10 * r["ci_hi"], f"{key}_p": r["p"]})
    row.update(sd_all=np.nanmean(np.nanstd(hp[:, m], 0)), sd_early=np.nanmean(np.nanstd(hp[yr <= 2017][:, m], 0)), sd_late=np.nanmean(np.nanstd(hp[yr >= 2018][:, m], 0)))
    ann = np.array([np.nanmean(np.nanstd(hp[yr == y][:, m], 0)) for y in YRS]); ra = fit(ann, monthly=False)
    row.update(sd_pct_change=100 * (row["sd_late"] / row["sd_early"] - 1), sd_slope_pct=1000 * ra["slope"] / np.nanmean(ann), sd_lo_pct=1000 * ra["ci_lo"] / np.nanmean(ann), sd_hi_pct=1000 * ra["ci_hi"] / np.nanmean(ann), sd_p=ra["p"])
    print(row["lon"], row["n_tracks"], round(row["hp_slope"], 1), round(row["contour_slope"], 1), flush=True); return row
if __name__ == "__main__":
    rows = []
    for key, arr in (("lat_hp", hp), ("lat_contour", lc)):
        y = np.nanmean(arr, 1); r = fit(y); rows.append(dict(series=key, unit="km per decade", slope=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"], n_tracks=int(np.isfinite(arr).any(0).sum())))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUT / "x6_relthr_circumpolar.csv", index=False)
    with Pool(int(os.environ.get("NPROC", "2"))) as pool: out = [r for r in pool.map(one, range(180)) if r is not None]
    pd.DataFrame(out).to_csv(OUT / "x6_relthr_bins.csv", index=False); print("done", len(out))
