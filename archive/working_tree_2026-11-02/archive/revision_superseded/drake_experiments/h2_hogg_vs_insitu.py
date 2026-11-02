"""H2: Antarctic-slope sea-level transport index (Hogg et al. 2015 / Hughes et al. 2014) against the Drake Passage records.
Inputs: the output CSV of h0_antarctic_slope_index.py (h_noseas_cm, h_f_cm, T_hogg_model_calibrated_Sv, coverage),
the GV2023 total 0-760 m transport (per transect) and the altimetric surface proxies (d1b_series.csv).
Compares, over 2005-2019: (i) monthly matched and annual means of the index (z-scores) and of T_hogg (Sv) with the in situ
total and north/south parts; (ii) the same with the surface proxy T_perm on the GV track; trends (trend_bcz) of T_hogg over
the full index record. The Sv conversion uses model-derived factors (Hughes et al. 2014) and is labelled as such; the in
situ series covers 0-760 m only, so only the variability, not the magnitude, is comparable.
Usage: python h2_hogg_vs_insitu.py <h0 csv> <d1 output dir> <GV repo dir> <revision dir> <out csv>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
h0, D1, GV, REV, OUTF = (Path(a) for a in sys.argv[1:6])
sys.path.insert(0, str(REV)); import rstats as R
H = pd.read_csv(h0)
m = sio.loadmat(GV/"Datasets"/"TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
tr = np.array(m.transp, float); dist = np.array(m.dist, float); tr[dist > 840, :] = np.nan
time = pd.to_datetime(np.array(m.time, float)-719529, unit="D")
ins = pd.DataFrame(dict(year=time.year, month=time.month, U_tot=np.nansum(tr, 0)/1e6, n=np.isfinite(tr).sum(0)))
ins = ins[ins.n >= 30]
sp = pd.read_csv(D1/"d1b_series.csv")[["year", "month", "Tperm_gv", "Tall_gv"]]
split = pd.read_csv(D1/"d1b_insitu_split.csv")[["year", "month", "U_north", "U_south"]]
ins = ins.merge(split, on=["year", "month"], how="left")
z = lambda s: (s - s.mean())/s.std()
out = {}
for idx in ("h_noseas_cm", "T_hogg_model_calibrated_Sv"):
    M = ins.merge(H[["year", "month", idx]], on=["year", "month"]).dropna(subset=[idx])
    A = M.groupby("year")[[idx, "U_tot", "U_north", "U_south"]].mean().loc[2005:2019]
    for c in ("U_tot", "U_north", "U_south"):
        out[f"{idx}|r_transects_{c}"] = np.corrcoef(M[idx], M[c])[0, 1] if len(M) > 10 else np.nan
        out[f"{idx}|r_annual_{c}"] = np.corrcoef(A[idx], A[c])[0, 1] if len(A) > 8 else np.nan
    S = sp.merge(H[["year", "month", idx]], on=["year", "month"]).dropna()
    out[f"{idx}|r_monthly_Tperm_gv"] = np.corrcoef(S[idx], S.Tperm_gv)[0, 1] if len(S) > 10 else np.nan
    out[f"{idx}|n_transects_matched"] = len(M)
G = H.dropna(subset=["T_hogg_model_calibrated_Sv"])
if len(G) > 120:
    o = R.trend_bcz(G.x.values, G.T_hogg_model_calibrated_Sv.values, month=G.month.values)
    out["T_hogg|trend_Sv_per_decade"] = 10*o["slope"]; out["T_hogg|ci"] = (10*o["ci_lo"], 10*o["ci_hi"]); out["T_hogg|p"] = o["p"]
    out["T_hogg|record"] = f"{int(G.year.min())}-{int(G.year.max())}"
out["index|mean_coverage"] = H.coverage.mean()
pd.Series(out).to_csv(OUTF); print(pd.Series(out).to_string())
