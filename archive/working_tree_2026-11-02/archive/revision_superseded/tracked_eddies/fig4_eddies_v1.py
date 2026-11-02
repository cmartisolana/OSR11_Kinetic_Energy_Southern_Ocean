"""Candidate Figure 4 (v2.1, draft 1): tracked-eddy distributions and their trends by longitude in the ACC band.
a-d  Distributions (density, log scale) of eddy amplitude, speed radius, rotational speed (all band observations) and
     track lifetime (tracks born in the band) in 2001-2011 (grey) and 2012-2023 (green); dashed lines: medians.
     Below each: ratio of the two densities (late / early).
e-g  Trends by 10 deg longitude bin, % of the bin mean per decade with 95 % intervals (filled: p < 0.05; ring: FDR 10 %):
     e mean amplitude and mean rotational speed; f sum(amplitude^2) per day (eddy-energy proxy) and eddies present per day;
     g mean speed radius. Horizontal lines: ACC-band values. Meander ranges (dk13 ridge + downstream) shaded.
Usage: python fig4_eddies_v1.py <e6 out dir> <revision dir> <figure dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

out, rev, figdir = (Path(a) for a in sys.argv[1:4]); figdir.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(rev))
import style_nc as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

D = pd.read_csv(out / "e6_distributions.csv"); Tr = pd.read_csv(out / "e6_lonbin_trends.csv")
MEANDERS = {"KP": (65, 95), "SEIR": (130, 150), "MR": (155, 175), "PAR": (-145, -100), "DP–SS": (-70, -35)}
C_E, C_L = "#9A9A9A", st.C["te"]

fig = plt.figure(figsize=(st.W2, 150 * st.MM))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.35], hspace=0.45)
gtop = gs[0].subgridspec(2, 4, height_ratios=[3, 1], hspace=0.08, wspace=0.45)
panels = [("amplitude", "amplitude (cm)", 100, "a"), ("speed_radius", "speed radius (km)", 1e-3, "b"),
          ("speed_average", "rotational speed (cm s$^{-1}$)", 100, "c"), ("lifetime_d", "lifetime (d)", 1, "d")]
for k, (v, lab, sc, let) in enumerate(panels):
    ax = fig.add_subplot(gtop[0, k]); axr = fig.add_subplot(gtop[1, k], sharex=ax)
    e = D[(D["var"] == v) & (D.period == "early") & D.left.notna()]; l = D[(D["var"] == v) & (D.period == "late") & D.left.notna()]
    x = (e.left.values + e.right.values) / 2 * sc
    ax.plot(x, e.density.values / sc, color=C_E, lw=0.9); ax.plot(x, l.density.values / sc, color=C_L, lw=0.9)
    for per, c in (("early", C_E), ("late", C_L)):
        med = D[(D["var"] == v) & (D.period == per) & D.left.isna()].iloc[0]["median"] * sc
        ax.axvline(med, color=c, lw=0.6, ls="--")
    ax.set_yscale("log"); ax.set_ylabel("density" if k == 0 else ""); ax.tick_params(labelbottom=False)
    st.plabel(ax, let, lab)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = l.density.values / e.density.values
    r[(e.density.values * len(e) < 1e-6)] = np.nan
    axr.plot(x, r, color=st.C["ink"], lw=0.8); axr.axhline(1, color="#9A9A9A", lw=0.5)
    axr.set_ylim(0.6, 1.6); axr.set_yticks([0.8, 1.0, 1.2, 1.4]); axr.set_xlabel(lab)
    if k == 0: axr.set_ylabel("late / early")
    if v == "lifetime_d": ax.set_xlim(0, 250)
    if v == "amplitude": ax.set_xlim(0, 35)
    if v == "speed_radius": ax.set_xlim(20, 120)
    if v == "speed_average": ax.set_xlim(0, 70)
ax0 = fig.axes[0]
ax0.legend(handles=[Line2D([], [], color=C_E, label="2001–2011"), Line2D([], [], color=C_L, label="2012–2023"),
                    Line2D([], [], color=st.C["ink"], ls="--", lw=0.6, label="median")], loc="upper right")

gbot = gs[1].subgridspec(3, 1, hspace=0.35)
spec = [("e", [("amp", st.C["te"], "o", "mean amplitude"), ("spd", st.C["ke"], "s", "mean rotational speed")]),
        ("f", [("amp2_per_day", st.C["te"], "D", "Σ amplitude² per day"), ("present_per_day", st.C["tr"], "^", "eddies present per day")]),
        ("g", [("rad", st.C["pos"], "v", "mean speed radius")])]
lonc = np.arange(-175, 180, 10)
for k, (let, items) in enumerate(spec):
    ax = fig.add_subplot(gbot[k])
    for lo, hi in MEANDERS.values():
        ax.axvspan(lo, hi, color=st.C["light"], lw=0, zorder=0)
    st.zero(ax)
    ns = len(items)
    for i, (v, c, m, lab) in enumerate(items):
        d = Tr[(Tr.region == "bin") & (Tr["var"] == v)].sort_values("bin")
        x = d.lon_c.values + (i - (ns - 1) / 2) * 2.0
        ax.vlines(x, d.pct_lo, d.pct_hi, color=c, lw=0.6)
        sig = d.p.values < 0.05
        ax.plot(x[~sig], d.pct.values[~sig], m, mfc="white", mec=c, mew=0.7, ms=3, lw=0)
        ax.plot(x[sig], d.pct.values[sig], m, color=c, ms=3, lw=0)
        f = d.fdr10.values.astype(bool)
        ax.plot(x[f], d.pct.values[f], "o", mfc="none", mec=st.C["ink"], mew=0.6, ms=6, lw=0)
        a = Tr[(Tr.region == "ACC") & (Tr["var"] == v)].iloc[0]
        ax.axhline(a.pct, color=c, lw=0.6, ls=":")
    ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 30))
    ax.set_xticklabels([f"{abs(t)}°{'W' if t < 0 else 'E'}" if t not in (0, 180, -180) else ("0°" if t == 0 else "180°") for t in np.arange(-180, 181, 30)])
    if k < 2: ax.tick_params(labelbottom=False)
    ax.set_ylabel("% per decade")
    ax.legend(handles=[Line2D([], [], color=c, marker=m, lw=0, ms=3, label=lab) for v, c, m, lab in items], loc="upper left", ncol=2)
    st.plabel(ax, let)
    if k == 0:
        for name, (lo, hi) in MEANDERS.items():
            ax.annotate(name, xy=((lo + hi) / 2, 1), xycoords=("data", "axes fraction"), xytext=(0, 1), textcoords="offset points",
                        ha="center", va="bottom", fontsize=5.5, color=st.C["grey"])
fig.axes[-1].set_xlabel("longitude")
fig.savefig(figdir / "Figure4_eddies_v1.png"); fig.savefig(figdir / "Figure4_eddies_v1.pdf")
print("saved")
