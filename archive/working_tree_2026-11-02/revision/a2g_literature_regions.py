"""A2g: half-power latitude trends 2001-2025 for the regions used in earlier studies (for the literature comparison).
Same matched baseline tracks and preparation as a2_run (mean of per-longitude anomalies within the region), revised
trend test. Regions (longitudes, deg E in -180..180): Chambers (2018) 90-170E and 160-90W (200-270E); Kim and Orsi (2014)
southeast Indian 75-150E and southeast Pacific 150-90W; Liu et al. (2026) Pacific-Antarctic Ridge meander 150-80W;
circumpolar. Fixed before computing (2026-10-15).
Usage: python a2g_literature_regions.py"""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[0])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import rstats as R
OUT = HERE / "outputs"
d = np.load(OUT / "a2" / "tracks_baseline.npz", allow_pickle=True)
lon = np.load((_CACHE + "/lon.npy"))
keep = d["keep"] & d["acc"]; lo = lon[d["j"][keep]]; L = d["lat_hp"][:, keep].astype(float)
x = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25)
REG = {"Circumpolar": (-180, 180), "90-170E (Chambers 2018)": (90, 170), "160-90W (Chambers 2018)": (-160, -90),
       "75-150E SE Indian (Kim and Orsi 2014)": (75, 150), "150-90W SE Pacific (Kim and Orsi 2014)": (-150, -90),
       "150-80W PAR meander (Liu et al. 2026)": (-150, -80)}
rows = []
for name, (a, b) in REG.items():
    sel = (lo >= a) & (lo < b)
    y = np.nanmean(L[:, sel] - np.nanmean(L[:, sel], 0, keepdims=True), 1)
    o = R.trend_bcz(x, y, month=month, n_outer=300); K = 111.195 * 10
    rows.append(dict(region=name, n_lon=int(sel.sum()), km_dec=o["slope"] * K, lo=o["ci_lo"] * K, hi=o["ci_hi"] * K, p=o["p"]))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "a2" / "a2g_literature_regions.csv", index=False)
