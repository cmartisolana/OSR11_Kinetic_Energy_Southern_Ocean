"""DK5b: in situ cross-track velocity (70-118 m) per 25-km bin in distance FROM NORTH, as in GV23 (sshbinning_transp.m):
northbound transects are flipped (dist -> max(dist) - dist + 12.5) and the ADCP offset uoffs is applied
(u - uoffs for northbound, u + uoffs for southbound). Replaces the in situ 'u' rows of dk5_trends.csv.

Purpose : recompute the in situ cross-track velocity rows of dk5_trends.csv with GV23's north-referenced distance and ADCP offset.
Inputs  : $GV23_DIR (Gutierrez-Villanueva et al. 2023 archive, doi:10.5281/zenodo.10044261); outputs/drake/dk5/dk5_trends.csv (03b)
Outputs : outputs/drake/dk5/dk5_trends.csv (in situ 'u' rows replaced; distributed)
Run     : python scripts/03_drake/03c_insitu_velocity_fix.py     (always after 03b)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, numpy as np, pandas as pd, scipy.io as sio
from pathlib import Path
from osr11 import rstats as R, paths as PA
GV, OUT = PA.GV23_DIR, PA.dk(5)
LV = sio.loadmat(GV/"Datasets"/"lmgvelos38nb_50_300.mat", squeeze_me=True, struct_as_record=False)["lmgdata"]
TH = sio.loadmat(GV/"Datasets"/"trandsmiss_theta.mat", squeeze_me=True)
uoffs = np.mean(np.abs(np.r_[np.atleast_1d(TH["uoffs_north"]), np.atleast_1d(TH["uoffs_south"])]))
U = np.array(LV.u, float); z = np.array(LV.z, float); d = np.array(LV.dist, float); nb = np.array(LV.northbound) == 1
U[:, :, nb] -= uoffs; U[:, :, ~nb] += uoffs
for k in np.where(nb)[0]:
    o = np.argsort(np.abs(np.nanmax(d[:, k]) - d[:, k]) + 12.5); d[:, k] = (np.abs(np.nanmax(d[:, k]) - d[:, k]) + 12.5)[o]; U[:, :, k] = U[:, o, k]
zs = (z >= 70) & (z <= 118)
t = pd.to_datetime(np.array(LV.time, float) - 719529, unit="D")
D = pd.read_csv(OUT/"dk5_trends.csv"); D = D[~((D.src == "insitu") & (D["var"] == "u"))]
rows = []
for c in np.arange(12.5, 850, 25):
    uu = np.full(U.shape[2], np.nan)
    for k in range(U.shape[2]):
        j = np.where(np.abs(d[:, k] - c) < 12.6)[0]
        if j.size: uu[k] = np.nanmean(U[zs][:, j[0], k]) * 100
    df = pd.DataFrame(dict(year=t.year, v=uu)).dropna(); df = df[(df.year >= 2005) & (df.year <= 2019)]
    g = df.groupby("year").v; a = g.mean()[g.count() >= 3]
    o = R.trend_bcz(a.index.values.astype(float), a.values, monthly=False)
    rows.append(dict(coord="passage", src="insitu", var="u", bin=c, period="2005-2019", mean=np.nanmean(uu), slope=o["slope"], ci_lo=o["ci_lo"], ci_hi=o["ci_hi"], p=o["p"], n=len(a)))
D = pd.concat([D, pd.DataFrame(rows)], ignore_index=True); D.to_csv(OUT/"dk5_trends.csv", index=False)
print("uoffs", uoffs); print(pd.DataFrame(rows)[["bin", "mean", "slope", "p"]].round(3).to_string())
