"""dk7: in situ EKE (0-780 m depth-integrated, Gutierrez-Villanueva et al. 2023, eke_uvprime_all_os38.mat) in 25 km bins of
distance from the north along the LM Gould transect; annual means (>= 3 transects per year), Oct 2005-Apr 2019; trend with the
calibrated annual test (trend_bcz, AR(1)); relative trend in % of the bin mean per decade (units differ from surface EKE)."""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[3])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, numpy as np, pandas as pd, scipy.io as sio
REV = _REV; sys.path.insert(0, REV); import rstats
m = sio.loadmat((_GV + "/Datasets/eke_uvprime_all_os38.mat"), squeeze_me=True, struct_as_record=False)["lmgdata"]
eke = np.array(m.eke, float); dist = np.array(m.dist, float); t = pd.to_datetime(np.array(m.time, float) - 719529, unit="D")
ok = np.array(m.flag) == 1
rows = []
for b0 in np.arange(0, 850, 25):
    c = b0 + 12.5; y = np.full(eke.shape[1], np.nan)
    for k in range(eke.shape[1]):
        j = np.where(np.abs(dist[:, k] - c) < 12.6)[0]
        if j.size and ok[k]: y[k] = eke[j[0], k]
    s = pd.Series(y, index=t.year); s = s[(s.index >= 2005) & (s.index <= 2019)]
    g = s.groupby(level=0).agg(["mean", "count"]); g = g[g["count"] >= 3]
    if len(g) < 8: continue
    r = rstats.trend_bcz(g.index.values.astype(float), g["mean"].values, monthly=False)
    mu = np.nanmean(s.values)
    rows.append(dict(bin=c, mean=mu, n_years=len(g), pct=1000 * r["slope"] / mu, pct_lo=1000 * r["ci_lo"] / mu, pct_hi=1000 * r["ci_hi"] / mu, p=r["p"]))
R = pd.DataFrame(rows); R.to_csv((_DK_OUT + "/dk7/dk7_insitu_eke_passage.csv"), index=False); print(R.round(2).to_string())
