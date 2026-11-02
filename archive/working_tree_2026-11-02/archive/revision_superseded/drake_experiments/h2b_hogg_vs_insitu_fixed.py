"""h2 with one row per transect (fixes the year-month many-to-many merge with d1b_insitu_split.csv in h2_hogg_vs_insitu.py,
which gave 519 rows from 248 transects). Split rows are aligned to transects by order, checked on U_tot.
Annual means: (a) index over matched transect months only, in situ over transects; (b) index over all months of the year.
Usage: python h2b_fixed.py <h0b csv> <d1 dir> <GV repo dir> <out csv>"""
import sys, numpy as np, pandas as pd, scipy.io as sio
from pathlib import Path
from scipy import stats
h0, D1, GV, OUTF = (Path(a) for a in sys.argv[1:5])
m = sio.loadmat(GV/"Datasets"/"TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
tr = np.array(m.transp, float); dist = np.array(m.dist, float); tr[dist > 840, :] = np.nan
t = pd.to_datetime(np.array(m.time, float)-719529, unit="D")
ins = pd.DataFrame(dict(year=t.year, month=t.month, U_tot=np.nansum(tr, 0)/1e6, n=np.isfinite(tr).sum(0)))
sp = pd.read_csv(D1/"d1b_insitu_split.csv")
assert len(sp) == len(ins) and np.allclose(sp.U_tot.values, ins.U_tot.values, atol=0.05), "split not aligned"
ins["U_north"], ins["U_south"] = sp.U_north.values, sp.U_south.values
ins = ins[ins.n >= 30]
H = pd.read_csv(h0); out = {}
for idx in ("h_noseas_cm", "T_hogg_model_calibrated_Sv"):
    M = ins.merge(H[["year", "month", idx]], on=["year", "month"], how="left").dropna(subset=[idx])
    out[f"{idx}|n_transects"] = len(M); out[f"{idx}|period"] = f"{M.year.min()}-{M.month.iloc[0]:02d} to {M.year.max()}-{M.month.iloc[-1]:02d}"
    A = M.groupby("year")[[idx, "U_tot", "U_north", "U_south"]].mean()
    Afull = H.groupby("year")[idx].mean().loc[A.index]
    for c in ("U_tot", "U_north", "U_south"):
        r, p = stats.pearsonr(M[idx], M[c]); out[f"{idx}|r_transects_{c}"] = f"{r:.2f} (p={p:.3f}, n={len(M)})"
        r, p = stats.pearsonr(A[idx], A[c]); out[f"{idx}|r_annual_matched_{c}"] = f"{r:.2f} (p={p:.2f}, n={len(A)})"
        r, p = stats.pearsonr(Afull, A[c]); out[f"{idx}|r_annual_fullyear_{c}"] = f"{r:.2f} (p={p:.2f}, n={len(A)})"
pd.Series(out).to_csv(OUTF); print(pd.Series(out).to_string())
