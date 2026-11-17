"""Review check C8 (4 Nov 2026): does the sea-level-driven drift of the fixed ADT bins affect the Fig. 2d cross-track velocity trends?
Inputs : dk24_section_fields_allsat.npz (monthly ADT, u, v, ice at the 121 track points), dk5_stream_edges.csv (the 17 bin edges),
         dk5_trends.csv (reported values). Reads only; writes c04_results.csv, c04_log.txt.
Method : reproduce the streamwise binning of 03b (points whose monthly ADT falls between two fixed edges; month NaN if any point of the
         band is ice-covered; trend from the complete monthly series or from annual means of >= 6 months), then repeat with the bin
         edges following the secular rise of ADT: edges + linear trend of the section-mean ADT (fitted on the 2001-2025 monthly
         section means, deseasonalised), and, as a second variant, edges + the full monthly section-mean ADT anomaly (12-month running
         mean). Reports the two significant interfrontal bins and the whole profile.

Origin  : review check c04_adt_bin_drift.py (4 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : Fig. S1c and SI S2.2 (drift of the ADT bins of Fig. 2d).
Outputs : outputs/c04/c04_results.csv (distributed), c04_var{0,1,2}.csv, c04_log.txt
Run     : python scripts/03_drake/03w_adt_bin_drift.py 0; ... 1; ... 2; then without argument (no raw data; about 2 min per variant)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
import numpy as np, pandas as pd
from osr11 import rstats as R, paths as PA
HERE = PA.out("c04")
Z = np.load(PA.dk(24, mkdir=False) / "dk24_section_fields_allsat.npz"); lo, la, dsx = Z["lon"], Z["lat"], Z["dsx"]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
ice = Z["ice"] > 0.10; un = np.where(ice, np.nan, Z["u"] * nx + Z["v"] * ny) * 100; adt = Z["adt"]
E = pd.read_csv(PA.dk(5, mkdir=False) / "dk5_stream_edges.csv"); edges = E.adt_edge.values; psic = 0.5 * (E.psi.values[1:] + E.psi.values[:-1])
REP = pd.read_csv(PA.dk(5, mkdir=False) / "dk5_trends.csv"); REP = REP[(REP.coord == "stream") & (REP.src == "allsat") & (REP["var"] == "u") & (REP.period == "2001-2025")]
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - .5) / 12
x = T.x.values; mo = T.month.values
# section-mean ADT: secular rise
smean = np.nanmean(np.where(ice, np.nan, adt), 1)            # may be NaN-free? ice months: use the always-ice-free points instead
perm = ~ice.any(0); smean_n = adt[:, perm].mean(1)
ds = R.deseasonalise(smean_n, mo); b, a = R.theil_sen(x, ds)
rise = b * (x - x.mean())                                      # m, zero-mean secular rise
anom = pd.Series(ds - ds.mean()).rolling(12, center=True, min_periods=6).mean().values
log = [f"section-mean ADT (always ice-free points, {perm.sum()} of 121): Theil-Sen rise {100*b*10:.2f} cm per decade; total over the record {100*b*25:.1f} cm; bin widths 5-15 cm"]
print(log[-1])
def annual(years, vals, minn):
    df = pd.DataFrame(dict(year=years, v=vals)).dropna(); g = df.groupby("year").v; aa = g.mean()[g.count() >= minn]
    return aa.index.values.astype(float), aa.values
def bin_trends(shift, label, n_outer=300):
    out = []
    for j in range(len(edges) - 1):
        lo_e, hi_e = sorted((edges[j], edges[j + 1]))
        inb = (adt - shift[:, None] >= lo_e) & (adt - shift[:, None] < hi_e)
        icy = (inb & ice).any(1)
        y = np.array([np.nanmean(un[k, inb[k]]) if inb[k].any() else np.nan for k in range(300)]); y[icy] = np.nan
        if np.isfinite(y).all(): r = R.trend_bcz(x, y, month=mo, n_outer=n_outer); kind = "monthly"
        else:
            xa, ya = annual(T.year.values, y, 6); r = R.trend_bcz(xa, ya, monthly=False, n_outer=n_outer); kind = f"annual n={len(xa)}"
        out.append(dict(variant=label, bin_psic=round(psic[j], 3), adt_lo=lo_e, adt_hi=hi_e, mean_u=np.nanmean(y), slope_dec=10 * r["slope"], lo=10 * r["ci_lo"], hi=10 * r["ci_hi"], p=r["p"], kind=kind))
    return pd.DataFrame(out)
VAR = sys.argv[1] if len(sys.argv) > 1 else "all"      # run one variant per call (each ~2 min): 0 fixed, 1 rise-following, 2 anomaly-following
if VAR == "0": bin_trends(np.zeros(300), "fixed edges (as 03b)", 150).to_csv(HERE / "c04_var0.csv", index=False); sys.exit()
if VAR == "1": bin_trends(rise, "edges follow the secular ADT rise", 150).to_csv(HERE / "c04_var1.csv", index=False); sys.exit()
if VAR == "2": bin_trends(np.nan_to_num(anom), "edges follow the 12-month section-mean ADT anomaly", 150).to_csv(HERE / "c04_var2.csv", index=False); sys.exit()
D0, D1, D2 = (pd.read_csv(HERE / f"c04_var{k}.csv") for k in "012")
D = pd.concat([D0, D1, D2]); D.to_csv(HERE / "c04_results.csv", index=False)
pd.set_option("display.width", 250)
chk = D0.merge(REP[["bin", "slope", "p"]].rename(columns={"bin": "bin_psic", "slope": "rep_slope_per_yr", "p": "rep_p"}), on="bin_psic", how="left")
log.append("reproduction of dk5_trends (slope per decade here vs 10 x reported per-year slope):")
log.append(chk.assign(rep_dec=10 * chk.rep_slope_per_yr)[["bin_psic", "slope_dec", "rep_dec", "p", "rep_p", "kind"]].round(3).to_string(index=False))
for j in range(len(D0)):
    log.append(f"bin psi {psic[j]:+.3f} (ADT {D0.adt_lo.iloc[j]:+.2f}..{D0.adt_hi.iloc[j]:+.2f}): fixed {D0.slope_dec.iloc[j]:+.2f} [{D0.lo.iloc[j]:+.2f},{D0.hi.iloc[j]:+.2f}] p={D0.p.iloc[j]:.3f} | "
               f"rise-following {D1.slope_dec.iloc[j]:+.2f} [{D1.lo.iloc[j]:+.2f},{D1.hi.iloc[j]:+.2f}] p={D1.p.iloc[j]:.3f} | anomaly-following {D2.slope_dec.iloc[j]:+.2f} p={D2.p.iloc[j]:.3f}")
(HERE / "c04_log.txt").write_text("\n".join(log)); print("\n".join(log[1:]))
