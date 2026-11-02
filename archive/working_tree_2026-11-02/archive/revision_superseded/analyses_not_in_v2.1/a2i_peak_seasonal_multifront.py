import sys, numpy as np, pandas as pd
sys.path.insert(0, ".")
import rstats as R, a2_run as A
d = np.load("outputs/a2/tracks_baseline.npz", allow_pickle=True)
lon = np.load("/home/claude/w/cache/lon.npy")
keep = d["keep"] & d["acc"]
j = d["j"][keep]; clim = d["lat_hp_clim"][keep]
L = d["lat_hp"][:, keep].astype(float); PK = d["peak_val"][:, keep].astype(float) * 1e4
x = 2001 + (np.arange(300) + .5) / 12; month = np.tile(np.arange(1, 13), 25)
K = 111.195
def tr(y, scale=1.0):
    o = R.trend_bcz(x, y, month=month, n_outer=300); return o["slope"]*scale*10, o["ci_lo"]*scale*10, o["ci_hi"]*scale*10, o["p"]
# 1 peak KE within envelope, % per decade
pk = np.nanmean(PK, 1); m = np.nanmean(pk)
s, lo, hi, p = tr(pk); print(f"peak KE ACC: {100*s/m:+.1f} %/dec [{100*lo/m:+.1f},{100*hi/m:+.1f}] p={p:.3f}  mean {m:.0f}")
# 2 seasonal cycle of position (km, + north)
y = np.nanmean(L - np.nanmean(L, 0, keepdims=True), 1) * K
sc = pd.Series(y).groupby(month).mean(); print("seasonal cycle km:", np.round(sc.values, 1), "range", round(sc.max()-sc.min(), 1))
an = pd.Series(y).groupby(np.repeat(np.arange(25), 12)).mean(); print("sd annual means km", round(an.std(), 1), "sd monthly deseason", round((y - sc.values[month-1]).std(), 1))
# 3 northern vs southern envelopes where >=2 matched envelopes at a longitude
df = pd.DataFrame(dict(j=j, clim=clim, i=np.arange(len(j))))
multi = df.groupby("j").filter(lambda g: len(g) >= 2)
N = multi.loc[multi.groupby("j").clim.idxmax()].sort_values("j"); S = multi.loc[multi.groupby("j").clim.idxmin()].sort_values("j")
print("longitudes with >=2 envelopes:", len(N))
for name, sel in (("northern", N.i.values), ("southern", S.i.values)):
    yy = np.nanmean(L[:, sel] - np.nanmean(L[:, sel], 0, keepdims=True), 1)
    s, lo, hi, p = tr(yy, K); print(f"{name}: mean lat {np.nanmean(L[:, sel]):.2f}  {s:+.1f} km/dec [{lo:+.1f},{hi:+.1f}] p={p:.3f}")
sep = L[:, N.i.values] - L[:, S.i.values]
yy = np.nanmean(sep - np.nanmean(sep, 0, keepdims=True), 1)
s, lo, hi, p = tr(yy, K); print(f"separation N-S: mean {np.nanmean(sep)*K:.0f} km  {s:+.1f} km/dec [{lo:+.1f},{hi:+.1f}] p={p:.3f}")
for name, sel in (("northern", N.i.values), ("southern", S.i.values)):
    pk = np.nanmean(PK[:, sel], 1); m = np.nanmean(pk); s, lo, hi, p = tr(pk)
    print(f"peak KE {name}: {100*s/m:+.1f} %/dec [{100*lo/m:+.1f},{100*hi/m:+.1f}] p={p:.3f}")
