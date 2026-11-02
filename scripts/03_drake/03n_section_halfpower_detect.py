"""dk6 T1 detectability (Product 1, no ice mask, coverage >= 50 %): climatological KE shifted northward linearly to S km over
2001-2025, plus observed monthly KE anomalies resampled by calendar year; detection, matching and passage-mean series as in
dk6_run.py; calibrated test (trend_bcz). Usage: python dk6_detect.py <S_km> <n_real> <seed>

Purpose : detectability of a synthetic northward shift of S km in Drake Passage (one third of 25 km shifts detected; SI S3).
Inputs  : outputs/drake/dk6/fields_allsat.npz (03l), outputs/drake/dk6/dk6_tracks_allsat_ice0.csv (03m)
Outputs : outputs/drake/dk6/dk6_detect_S<S>_seed<seed>.csv (distributed)
Run     : for S in 0 10 25; do python scripts/03_drake/03n_section_halfpower_detect.py $S 30 $S; done
Seeds   : the third argument (the shipped files used seed = S).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, time, numpy as np, pandas as pd
from osr11.drake_halfpower import *  # noqa: F401,F403
Skm, NR, seed = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); KM = 111.195
D = load("allsat", False); lat, lon = D["lat"], D["lon"]; dists = F.meridional_dists(lat)
ke = D["ke"]; kec = np.nanmean(ke, 0); A = ke - kec[None]
cc, ctb = F.detect_field(np.ascontiguousarray(kec), lat, dists, P["thr"], P["min_w"], P["frac"], MAXE)
trk = pd.read_csv(f"{OUT}/dk6_tracks_allsat_ice0.csv"); trk = trk[trk["cov"] >= 0.5]
rng = np.random.default_rng(seed); out = []
s_t = Skm * (T - T[0]) / (T[-1] - T[0])          # km, linear ramp
# shifted climatology per month (interpolate along latitude, column by column)
for r in range(NR):
    t0 = time.time()
    yrs = rng.integers(0, 25, 25)
    syn = np.empty_like(ke)
    for k in range(300):
        y, m = divmod(k, 12)
        dl = s_t[k] / KM
        sh = np.array([np.interp(lat - dl, lat, kec[:, j], left=np.nan, right=np.nan) for j in range(lon.size)]).T
        syn[k] = sh + A[yrs[y] * 12 + m]
    C, TB = detect_all(syn, lat); meta, arr = match(C, TB, cc, ctb, lat)
    q = [int(np.where((meta[:, 0] == rw.j) & (meta[:, 1] == rw.c))[0][0]) for rw in trk.itertuples()]
    hp = arr[:, q, FI["lat_hp"]]; yser = np.nanmean((hp - np.nanmean(hp, 0, keepdims=True)) * KM, 1)
    rr = tr(yser)
    out.append(dict(S_km=Skm, real=r, trend_25yr_km=rr["slope"] * 25, p=rr["p"], detected=rr["p"] < 0.05, secs=time.time() - t0))
    print(out[-1], flush=True)
pd.DataFrame(out).to_csv(f"{OUT}/dk6_detect_S{int(Skm)}_seed{seed}.csv", index=False)
