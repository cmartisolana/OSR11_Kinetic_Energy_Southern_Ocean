"""v2.3 Figure 3: transport proxies against eddy energy, all in percent of the mean.
(a) Drake Passage: transect (solid: ice-free points used) and the domain for wind stress, KE and EKE.
(b) Transport proxies, 12-month running means in % of the mean: Drake altimetric surface proxy, circumpolar surface
    proxy (fixed ACC band, a8), in situ 0-760 m transport (annual means, 2005-2019).
(c) Drake-domain EKE and KE, 12-month running means in % of the mean.
(d) Trends in % of the mean per decade with 95 % CI (revised test; filled p < 0.05): transport proxies (Drake altimetric,
    Drake in situ, circumpolar and sectors), KE and EKE (Drake domain, ACC band), eastward wind stress (Drake box,
    circumpolar 45-60S). All-satellite product.
Usage: python v23_fig3.py <repository outputs/trends dir>
"""
import sys
import numpy as np, pandas as pd
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

trends_dir = sys.argv[1]
FIG = OUT / "figures_v23"; FIG.mkdir(exist_ok=True)
RT = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")
Ttr = pd.read_csv(OUT / "a6" / "a6c_transport_fixed_series.csv")
PTS = pd.read_csv(OUT / "a6" / "a6c_section_points.csv").set_index("line")
BX = pd.read_csv(OUT / "a6" / "a6b_drake_box_series.csv")
WS = pd.read_csv(OUT / "a6" / "a6_wind_stress_series.csv"); WS = WS[(WS.year >= 2001) & (WS.year <= 2025)].reset_index(drop=True)
A8 = pd.read_csv(OUT / "a8" / "a8_series_allsat.csv"); A8T = pd.read_csv(OUT / "a8" / "a8_trends_allsat.csv").set_index("series")
obs = pd.read_csv(f"{trends_dir}/drake_passage_total_transport_observations.csv", parse_dates=["date"]).dropna(subset=["total_transport_sv"])
ann = obs.groupby(obs.date.dt.year).total_transport_sv.mean(); ann = ann[(ann.index >= 2005) & (ann.index <= 2019)]
x = (Ttr.year + (Ttr.month - 0.5) / 12).values
C_DP, C_ACC, C_OBS, C_EKE, C_KE, C_W = "#1D6FA3", "#1F2933", "#D9822B", STYLE["eke"], STYLE["ke"], "#7A8288"


def title(ax, s):
    ax.set_title(s, loc="left", fontsize=11.5, fontweight="bold")


def pct_row(key, mean):
    r = RT.loc[key]; f = 1000 / mean
    return dict(pct=r.slope * f, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p_new)


def a8_row(s):
    r = A8T.loc[s]; return dict(pct=r.pct_per_decade, pct_lo=r.pct_lo, pct_hi=r.pct_hi, p=r.p)


pa = lambda y: 100 * (y / np.nanmean(y) - 1)
fig = plt.figure(figsize=(11.0, 10.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], width_ratios=[0.85, 1.6], hspace=0.28, wspace=0.22)
gr = gs[0, 1].subgridspec(2, 1, hspace=0.42)

ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo(central_longitude=-65))
ax.set_extent([-84, -46, -66, -51], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor="#111111", edgecolor="#111111", zorder=3)
ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.4, color="#7A8288", alpha=0.6, linestyle="--")
bx = np.r_[np.c_[np.linspace(-80, -52, 50), np.full(50, -63)], np.c_[np.full(20, -52), np.linspace(-63, -54, 20)],
           np.c_[np.linspace(-52, -80, 50), np.full(50, -54)], np.c_[np.full(20, -80), np.linspace(-54, -63, 20)]]
ax.plot(bx[:, 0], bx[:, 1], color=STYLE["box"], lw=1.8, transform=ccrs.PlateCarree(), label="Domain (wind, KE, EKE)", zorder=4)
la_s = PTS.loc["v2", "perm_southern_lat"]; lo_s = -65 + (la_s + 55) / (-63 + 55) * 5
ax.plot([-65, lo_s], [-55, la_s], color=C_DP, lw=3, transform=ccrs.Geodetic(), zorder=5, label="Transect (ice-free part)")
ax.plot([lo_s, -60], [la_s, -63], color=C_DP, lw=1.5, ls=(0, (2, 1.5)), transform=ccrs.Geodetic(), zorder=5, label="Seasonal ice (excluded)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=7.5, frameon=False, handlelength=2.2)
title(ax, "a) Drake Passage")

ax = fig.add_subplot(gr[0])
ax.plot(x, running(pa(Ttr.Tperm_v2.values)), color=C_DP, lw=2.0, label="Drake, altimetric surface proxy")
xa = (A8.year + (A8.month - 0.5) / 12).values
ax.plot(xa, running(pa(A8.ACC_T.values)), color=C_ACC, lw=2.0, label="ACC band, circumpolar surface proxy")
ax.plot(ann.index + 0.5, pa(ann.values), "o-", color=C_OBS, ms=4, lw=1.4, label="Drake, in situ 0–760 m (annual)")
ax.axhline(0, color="#6b7280", lw=0.8, ls="--"); ax.set_xlim(2000.5, 2026); ax.set_ylabel("% of mean")
base_timeseries_style(ax); ax.legend(loc="upper left", fontsize=7.5, frameon=False, ncol=1)
yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.6 * (yl[1] - yl[0]))
title(ax, "b) Transport proxies")

ax = fig.add_subplot(gr[1])
xb = (BX.year + (BX.month - 0.5) / 12).values
ax.plot(xb, running(pa(BX.ice_any_EKE.values)), color=C_EKE, lw=2.0, label="EKE")
ax.plot(xb, running(pa(BX.ice_any_KE.values)), color=C_KE, lw=2.0, label="KE")
ax.axhline(0, color="#6b7280", lw=0.8, ls="--"); ax.set_xlim(2000.5, 2026); ax.set_ylabel("% of mean"); ax.set_xlabel("Year")
base_timeseries_style(ax); ax.legend(loc="upper left", fontsize=7.5, frameon=False, ncol=2)
yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.3 * (yl[1] - yl[0]))
title(ax, "c) Drake Passage KE and EKE")

ax = fig.add_subplot(gs[1, :])
mean_dp = float(np.nanmean(Ttr.Tperm_v2.values)); mean_obs = float(ann.mean())
mBX = {v: float(np.nanmean(BX[f"ice_any_{v}"])) for v in ("EKE", "KE")}
mW = {c: float(WS[c].mean()) for c in ("Drake_box", "Circumpolar_45_60S")}
T, E, K, W_ = "transport proxy", "eddy kinetic energy", "kinetic energy", "eastward wind stress"
rows = [("Drake: altimetric surface proxy (2001–2025)", {T: pct_row("Tperm_v2|2001-2025", mean_dp)}),
        ("Drake: in situ 0–760 m (2005–2019)", {T: pct_row("insitu_annual|2005-2019", mean_obs)}),
        ("ACC band: circumpolar surface proxy", {T: a8_row("ACC_T")}),
        ("   Pacific sector", {T: a8_row("Pacific_T")}),
        ("   Atlantic sector", {T: a8_row("Atlantic_T")}),
        ("   Indian sector", {T: a8_row("Indian_T")}),
        ("Drake: EKE", {E: pct_row("Drake_EKE_ice_any|2001-2025", mBX["EKE"])}),
        ("Drake: KE", {K: pct_row("Drake_KE_ice_any|2001-2025", mBX["KE"])}),
        ("ACC band: EKE", {E: retest_row("allsat|ACC_EKE|2001-2025")}),
        ("ACC band: KE", {K: retest_row("allsat|ACC_KE|2001-2025")}),
        ("Wind stress: Drake domain", {W_: pct_row("tau_Drake_box|2001-2025", mW["Drake_box"])}),
        ("Wind stress: 45–60°S circumpolar", {W_: pct_row("tau_Circumpolar_45_60S|2001-2025", mW["Circumpolar_45_60S"])})]
forest_panel(ax, rows, series_colors={T: C_DP, E: C_EKE, K: C_KE, W_: C_W}, series_markers={T: "s", E: "o", K: "o", W_: "D"},
             xlim=(-18, 18), label_size=9, dodge=False)
for yb in (-5.5, -9.5):
    ax.axhline(yb, color="#CED4DA", lw=0.8)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=8, frameon=False, ncol=5, handletextpad=0.3, columnspacing=1.0)
title(ax, "d) Trends 2001–2025, % of the mean per decade (in situ: 2005–2019)")
fig.savefig(FIG / "Figure3.png"); fig.savefig(FIG / "Figure3.pdf")
print("saved")
for gl, d in rows:
    for s, r in d.items(): print(f"{gl:45s} {r['pct']:+6.1f} [{r['pct_lo']:+6.1f}, {r['pct_hi']:+6.1f}] p={r['p']:.3f}")
