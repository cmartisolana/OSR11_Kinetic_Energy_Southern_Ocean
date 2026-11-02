"""Quick look at B2 (box-free meander -> EKE): (a) best-lag r by window centre for the three links;
(b) along-stream response of band EKE to the path length of the three meander windows (alt. 1)."""
import numpy as np, pandas as pd
from figstyle import *  # noqa
T = pd.read_csv(OUT / "b2" / "b2_window_tests_v2.csv"); A = pd.read_csv(OUT / "b2" / "b2_alt1_lonbins.csv")
fig, (a, b) = plt.subplots(2, 1, figsize=(10, 7.2), gridspec_kw=dict(hspace=0.45))
for link, lab, c in (("E", "EKE in the same window (lag 0–3, max |r|)", "#1F2933"), ("D", "EKE 15–35° downstream", "#C0392B"), ("B", "eddy births per 1000 km of path", "#2A9D8F")):
    r = T[f"rmax_{link}"].values; p = T[f"p_{link}"].values
    sgn = np.sign(T[[f"r_{link}_lag{k}" for k in range(4)]].values[np.arange(len(T)), np.abs(T[[f"r_{link}_lag{k}" for k in range(4)]].values).argmax(1)])
    a.plot(T.centre, r * sgn, "-", color=c, lw=1.4, label=lab)
    ok = p < 0.05
    a.plot(T.centre[ok], (r * sgn)[ok], "o", color=c, ms=5); a.plot(T.centre[~ok], (r * sgn)[~ok], "o", mfc="white", color=c, ms=5)
for c0 in T.centre[T.coverage < 0.8]:
    a.axvspan(c0 - 5, c0 + 5, color="#E9ECEF", lw=0, zorder=0)
a.axhline(0, color="#495057", lw=0.8); a.set_xlim(-185, 180); a.set_xlabel("Window centre (°E)"); a.set_ylabel("r (signed, best lag)")
a.legend(fontsize=7.5, frameon=False, loc="upper left", ncol=1); base_timeseries_style(a)
a.set_title("a) Path length → EKE / births, 36 windows of 30° (grey: path not tracked in ≥ 20 % of months)", loc="left", fontsize=10, fontweight="bold")
for m, c in (("KERG", "#457B9D"), ("SEIR150", "#2A9D8F"), ("PAR", "#E07A5F")):
    s = A[A.meander == m].sort_values("rel_lon")
    for k, ls in ((0, "-"), (2, "--")):
        y = s[f"r_lag{k}"].rolling(5, center=True, min_periods=3).mean()
        b.plot(s.rel_lon, y, ls, color=c, lw=1.5, label=f"{m} lag {k}")
b.axvspan(-15, 15, color="#E9ECEF", lw=0, zorder=0); b.axhline(0, color="#495057", lw=0.8)
b.set_xlabel("Longitude relative to the window centre (°; shaded = the window whose path length is used)"); b.set_ylabel("r (path length, band EKE)")
b.legend(fontsize=7, frameon=False, ncol=3); base_timeseries_style(b)
b.set_title("b) Along-stream response (1° bins, 5° running mean)", loc="left", fontsize=10, fontweight="bold")
fig.savefig(OUT / "b2" / "b2_quicklook.png", dpi=150); print("saved")
