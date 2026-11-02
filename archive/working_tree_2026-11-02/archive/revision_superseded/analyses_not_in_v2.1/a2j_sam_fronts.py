import sys, numpy as np, pandas as pd
sys.path.insert(0, "."); import rstats as R
rng = np.random.default_rng(20261016)
sam = pd.read_csv("data/sam_marshall_monthly.txt", sep=r"\s+", comment="#", header=None, na_values="NA")
sam.columns = ["year"] + list(range(1, 13)); sam = sam.melt(id_vars="year", var_name="month", value_name="sam"); sam["month"] = sam.month.astype(int)
sam = sam[(sam.year >= 2001) & (sam.year <= 2025)].sort_values(["year", "month"]).sam.values.astype(float)
d = np.load("outputs/a2/tracks_baseline.npz", allow_pickle=True); lon = np.load("/home/claude/w/cache/lon.npy")
keep = d["keep"] & d["acc"]; lo = lon[d["j"][keep]]; L = d["lat_hp"][:, keep].astype(float) * 111.195; PK = d["peak_val"][:, keep].astype(float)
x = 2001 + (np.arange(300) + .5) / 12; month = np.tile(np.arange(1, 13), 25)
def prep(y):
    y = R.deseasonalise(y, month); b, a = R.theil_sen(x, y); return y - (a + b * x)
def lagc(s, y, Lg): return np.corrcoef(s[:300 - Lg], y[Lg:])[0, 1]
def surr(s):
    f = np.fft.rfft(s); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0; ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=s.size)
fs = prep(sam); SU = [surr(fs) for _ in range(2000)]
sec = {"ACC": np.ones(lo.size, bool), "Atlantic": (lo >= -70) & (lo < 20), "Indian": (lo >= 20) & (lo < 147), "Pacific": (lo >= 147) | (lo < -70)}
for name, sel in sec.items():
    for var, A in (("position", L), ("peak KE", PK)):
        y = np.nanmean(A[:, sel] - np.nanmean(A[:, sel], 0, keepdims=True), 1) if var == "position" else np.nanmean(A[:, sel], 1)
        fy = prep(y); rr = [lagc(fs, fy, g) for g in range(13)]; b = int(np.argmax(np.abs(rr)))
        smax = np.array([max(abs(lagc(s, fy, g)) for g in range(13)) for s in SU]); p = (1 + (smax >= abs(rr[b])).sum()) / 2001
        ya = pd.Series(fy).groupby(np.repeat(np.arange(25), 12)).mean(); sa = pd.Series(fs).groupby(np.repeat(np.arange(25), 12)).mean()
        print(f"{name:8s} {var:9s} r(lag0)={rr[0]:+.2f}  best lag {b:2d} r={rr[b]:+.2f} p_max={p:.3f}  annual r={np.corrcoef(ya, sa)[0,1]:+.2f}")
