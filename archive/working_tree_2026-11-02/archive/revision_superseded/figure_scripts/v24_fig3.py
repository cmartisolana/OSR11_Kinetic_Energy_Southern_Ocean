"""v2.4 Figure 3. Changes from v23_fig3.py (author comments, 2026-10-15): the trend panel is split into
(d) transport proxies and (e) KE, EKE and wind stress; new (f, g) comparison of the altimetric and in situ Drake proxies:
standardised anomalies at the 244 transects (altimetric proxy of the transect month), scatter with r and RMSE,
and histogram of the differences. Same data and statistics as v2.3.
--- v2.3 description ---
v2.3 Figure 3: transport proxies against eddy energy, all in percent of the mean.
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
FIG = OUT / "figures_v24"; FIG.mkdir(exist_ok=True)
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
fig = plt.figure(figsize=(12.0, 14.2))
outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.72, 0.62], hspace=0.48)
gs = outer[0].subgridspec(1, 2, width_ratios=[0.85, 1.6], wspace=0.22)
gr = gs[0, 1].subgridspec(2, 1, hspace=0.42)
gmid = outer[1].subgridspec(1, 2, wspace=1.05)
gbot = outer[2].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.3)

ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo(central_longitude=-65))
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

mean_dp = float(np.nanmean(Ttr.Tperm_v2.values)); mean_obs = float(ann.mean())
mBX = {v: float(np.nanmean(BX[f"ice_any_{v}"])) for v in ("EKE", "KE")}
mW = {c: float(WS[c].mean()) for c in ("Drake_box", "Circumpolar_45_60S")}
T, E, K, W_ = "transport proxy", "eddy kinetic energy", "kinetic energy", "eastward wind stress"
rows_t = [("Drake: altimetric surface\nproxy (2001–2025)", {T: pct_row("Tperm_v2|2001-2025", mean_dp)}),
          ("Drake: in situ 0–760 m\n(2005–2019)", {T: pct_row("insitu_annual|2005-2019", mean_obs)}),
          ("ACC band: circumpolar\nsurface proxy", {T: a8_row("ACC_T")}),
          ("Pacific sector", {T: a8_row("Pacific_T")}),
          ("Atlantic sector", {T: a8_row("Atlantic_T")}),
          ("Indian sector", {T: a8_row("Indian_T")})]
rows_e = [("Drake: EKE", {E: pct_row("Drake_EKE_ice_any|2001-2025", mBX["EKE"])}),
          ("Drake: KE", {K: pct_row("Drake_KE_ice_any|2001-2025", mBX["KE"])}),
          ("ACC band: EKE", {E: retest_row("allsat|ACC_EKE|2001-2025")}),
          ("ACC band: KE", {K: retest_row("allsat|ACC_KE|2001-2025")}),
          ("Wind stress:\nDrake domain", {W_: pct_row("tau_Drake_box|2001-2025", mW["Drake_box"])}),
          ("Wind stress:\n45–60°S", {W_: pct_row("tau_Circumpolar_45_60S|2001-2025", mW["Circumpolar_45_60S"])})]
ax = fig.add_subplot(gmid[0, 0])
forest_panel(ax, rows_t, series_colors={T: C_DP}, series_markers={T: "s"}, xlim=(-18, 12), label_size=8.5, value_size=7.5)
ax.axhline(-2.5, color="#CED4DA", lw=0.8)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8, frameon=False, ncol=2, handletextpad=0.3, columnspacing=1.0)
title(ax, "d) Transport proxy trends, % dec$^{-1}$")
ax = fig.add_subplot(gmid[0, 1])
forest_panel(ax, rows_e, series_colors={E: C_EKE, K: C_KE, W_: C_W}, series_markers={E: "o", K: "o", W_: "D"},
             xlim=(-4, 18), label_size=8.5, value_size=7.5, dodge=False)
ax.axhline(-3.5, color="#CED4DA", lw=0.8)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8, frameon=False, ncol=2, handletextpad=0.3, columnspacing=1.0)
title(ax, "e) KE, EKE and wind stress trends, % dec$^{-1}$")
rows = rows_t + rows_e

# (f) altimetric vs in situ at the transects
obs["year"] = obs.date.dt.year; obs["month"] = obs.date.dt.month
mm = obs.merge(Ttr[["year", "month", "Tperm_v2"]], on=["year", "month"])
zs = lambda v: (v - v.mean()) / v.std()
za, zi = zs(mm.Tperm_v2.values), zs(mm.total_transport_sv.values)
r_tr = float(np.corrcoef(za, zi)[0, 1]); rmse = float(np.sqrt(np.mean((za - zi) ** 2)))
ax = fig.add_subplot(gbot[0, 0])
ax.scatter(za, zi, s=14, color=C_OBS, alpha=0.65, lw=0)
lim = 1.05 * max(np.abs(za).max(), np.abs(zi).max())
ax.plot([-lim, lim], [-lim, lim], color="#6b7280", lw=0.8, ls="--")
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xlabel("Altimetric surface proxy (z)"); ax.set_ylabel("In situ 0–760 m (z)"); base_timeseries_style(ax)
ax.text(0.03, 0.97, f"n = {len(mm)} transects\nr = {r_tr:.2f}\nRMSE = {rmse:.2f}", transform=ax.transAxes, va="top", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))
title(ax, "f) Drake transport at the transects")
ax = fig.add_subplot(gbot[0, 1])
dz = za - zi
ax.hist(dz, bins=np.arange(-5, 5.01, 0.5), color=C_OBS, alpha=0.75, edgecolor="white")
ax.axvline(0, color="#6b7280", lw=0.8, ls="--")
ax.set_xlabel("Altimetric − in situ (z)"); ax.set_ylabel("Transects"); base_timeseries_style(ax)
ax.text(0.97, 0.97, f"sd = {dz.std():.2f}\n(√2 = 1.41 for\nuncorrelated series)", transform=ax.transAxes, ha="right", va="top", fontsize=8.5)
title(ax, "g) Differences at the transects")
print(f"transects: n={len(mm)} r={r_tr:.4f} rmse={rmse:.3f} sd_diff={dz.std():.3f}")
fig.savefig(FIG / "Figure3.png"); fig.savefig(FIG / "Figure3.pdf")
print("saved")
for gl, d in rows:
    for s, r in d.items(): print(f"{gl:45s} {r['pct']:+6.1f} [{r['pct_lo']:+6.1f}, {r['pct_hi']:+6.1f}] p={r['p']:.3f}")
