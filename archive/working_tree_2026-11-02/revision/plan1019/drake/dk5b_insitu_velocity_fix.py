"""DK5b: in situ cross-track velocity (70-118 m) per 25-km bin in distance FROM NORTH, as in GV23 (sshbinning_transp.m):
northbound transects are flipped (dist -> max(dist) - dist + 12.5) and the ADCP offset uoffs is applied
(u - uoffs for northbound, u + uoffs for southbound). Replaces the in situ 'u' rows of dk5_trends.csv.
Usage: python dk5b_insitu_velocity_fix.py <revision dir> <GV repo dir> <dk5 out dir>"""
import sys, numpy as np, pandas as pd, scipy.io as sio
from pathlib import Path
REV, GV, OUT = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV)); import rstats as R
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
