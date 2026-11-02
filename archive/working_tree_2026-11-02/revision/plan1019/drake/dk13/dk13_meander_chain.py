"""dk13: the standing-meander chain of Thompson and Naveira Garabato (2014) at five meanders, from the along-stream
monthly series of b3 (jet-core path, 5 deg bins; b3_bin_series.npz). Requested by Referee Gutierrez-Villanueva (major
comment 3): test the links at specific meanders instead of sector averages.
Regions per meander, defined from bathymetry (the 5 deg bin trends of b3 had already been seen; the composites had not) (ridge = where the path or the plateau
constrains the flow; upstream = 10 deg west of it; downstream = 15 deg east, as the 3 downstream bins of b3_relations):
  links: (1) upstream jet speed, (2) sinuosity and (2b) RMS curvature of the path (B4) over the ridge and downstream (flexing),
         (3) EKE downstream, (4) jet speed downstream (final link: slowing / compensation).
Composite series = mean of the bin series in the region (relative to each bin's mean, so bins weigh equally);
trend: trend_bcz (deseasonalised Theil-Sen, calibrated, 300 null series), % of the mean per decade; FDR (BH, 10 %)
over all tests (25 with the curvature link). Malvinas bins (60-50 W) excluded as in b3. Product 1.
Usage: python dk13_meander_chain.py <revision dir> <out dir> [b4_bin_series.npz (adds the curvature link, B4)]"""
import sys, numpy as np, pandas as pd
from pathlib import Path
REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
Z = dict(np.load(REV / "outputs" / "b3" / "b3_bin_series.npz")); Z["rms_kappa"] = np.load(Path(sys.argv[3]))["rms_kappa"] if len(sys.argv) > 3 else None; C = pd.read_csv(REV / "outputs" / "b3" / "b3_climatology.csv")
lon = C.lon.values
MEANDERS = {   # name: (upstream, ridge, downstream) longitude ranges (deg E, -180..180)
    "Kerguelen Plateau":       ((55, 65), (65, 80), (80, 95)),
    "Southeast Indian Ridge":  ((120, 130), (130, 140), (140, 150)),
    "Macquarie Ridge":         ((150, 155), (155, 160), (160, 175)),
    "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
    "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35)),
}
LINKS = [("upstream speed", "speed", 0), ("sinuosity, ridge and downstream", "sinu", (1, 2)),
         ("curvature, ridge and downstream", "rms_kappa", (1, 2)), ("downstream EKE", "eke", 2),
         ("downstream speed", "speed", 2)]
T = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
excl = (lon > -60) & (lon < -50)
rows = []
for name, rngs in MEANDERS.items():
    for lab, var, which in LINKS:
        if Z.get(var) is None: continue
        idx = which if isinstance(which, tuple) else (which,)
        m = np.zeros(lon.size, bool)
        for i in idx:
            a, b = rngs[i]; m |= (lon > a) & (lon < b)
        m &= ~excl
        X = Z[var][:, m].astype(float); X = X / np.nanmean(X, 0, keepdims=True)
        y = 100 * np.nanmean(X, 1)
        r = R.trend_bcz(T, y, month=MON, n_outer=300); mu = np.nanmean(y)
        rows.append(dict(meander=name, link=lab, var=var, bins=",".join(f"{v:g}" for v in lon[m]), n_bins=int(m.sum()),
                         pct=1000 * r["slope"] / mu, lo=1000 * r["ci_lo"] / mu, hi=1000 * r["ci_hi"] / mu, p=r["p"]))
        print(name, lab, round(rows[-1]["pct"], 2), round(r["p"], 3), flush=True)
D = pd.DataFrame(rows)
p = D.p.values; o = np.argsort(p); q = p[o] * len(p) / (np.arange(len(p)) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
fdr = np.zeros(len(p), bool); fdr[o] = q <= 0.10; D["fdr10"] = fdr
D.to_csv(OUT / "dk13_meander_chain_allsat.csv", index=False); print(D.to_string(index=False))
