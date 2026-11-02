"""Main Fig. 4 (v1 style): Drake Passage, 2001-2025.
(a) Drake box (80-52W, 63-54S) and altimetric transport-proxy section (v2 line; thin: variants;
    solid part: points ice-free in all months, used for the proxy)
(b) standardised trends (sigma per yr, 95 % CI, calibrated p): altimetric proxy (v2 line, fixed
    ice-free points, 2001-2025; diamond: all points incl. under ice),
    in situ ADCP 0-760 m (annual means 2005-2019), zonal wind stress 45-60S (2001-2025)
(c) Drake-box EKE and KE trends (ice mask; diamonds: no mask)
(d) standardised anomalies: altimetric proxy (fixed ice-free points; monthly, 12-month running mean, band = line variants),
    in situ (transects, 12-transect running mean), wind stress 45-60S (12-month running mean)
(e) Drake-box KE and EKE (ice mask: raw monthly, 12-month running mean, trend; dotted: no mask)
Usage: python fig4_main.py <revision/data dir> <repository outputs/trends dir>
"""
import sys
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

data, trends_dir = sys.argv[1], sys.argv[2]
TR = {"v2": (-65.0, -55.0, -60.0, -63.0)}
for dl in (-1.0, -0.5, 0.5, 1.0):
    TR[f"v2_lon{dl:+.1f}"] = (-65.0 + dl, -55.0, -60.0 + dl, -63.0)
TR["v2_south+0.5"] = (-65.0, -55.0, -60.0, -62.5); TR["v2_north-0.5"] = (-65.0, -55.5, -60.0, -63.0)
TR["gould_approx"] = (-65.0, -55.0, -58.5, -62.0)
# transport proxy: fixed geometry (a6c_transport_fixed.py); a6_transport_series (monthly ice exclusion) is not used
Ttr = pd.read_csv(OUT / "a6" / "a6c_transport_fixed_series.csv")
A6c = pd.read_csv(OUT / "a6" / "a6c_trends.csv"); A6c.index = A6c.series + "_" + A6c.period
A6c = retest_override(A6c, lambda i, r: f"{r.series}|{r.period}")
PTS = pd.read_csv(OUT / "a6" / "a6c_section_points.csv").set_index("line")
A6 = pd.read_csv(OUT / "a6" / "a6_trends.csv").set_index("series")
A6 = retest_override(A6, lambda i, r: {"insitu_annual_2005-2019": "insitu_annual|2005-2019"}.get(i) or
                     (f"{i[:-10]}|{i[-9:]}" if i.startswith(("tau_", "dADT")) else None))
A6b = pd.read_csv(OUT / "a6" / "a6b_drake_box_trends.csv").query("period == '2001-2025'").set_index("series")
A6b = retest_override(A6b, lambda i, r: f"{i}|2001-2025")
BX = pd.read_csv(OUT / "a6" / "a6b_drake_box_series.csv")
WS = pd.read_csv(OUT / "a6" / "a6_wind_stress_series.csv"); WS = WS[(WS.year >= 2001) & (WS.year <= 2025)].reset_index(drop=True)
AG = pd.read_csv(OUT / "a6" / "a6c_agreement.csv", index_col=0)
obs = pd.read_csv(f"{trends_dir}/drake_passage_total_transport_observations.csv", parse_dates=["date"]).dropna(subset=["total_transport_sv"])
x = (Ttr.year + (Ttr.month - 0.5) / 12).values

fig = plt.figure(figsize=(10.0, 10.2))
gs = fig.add_gridspec(3, 3, height_ratios=[0.95, 1.0, 1.0], hspace=0.38, wspace=0.42)
# (a) map
ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo(central_longitude=-65))
ax.set_extent([-84, -46, -66, -51], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor="#111111", edgecolor="#111111", zorder=3)
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.4, color="#7A8288", alpha=0.6, linestyle="--")
bx = np.r_[np.c_[np.linspace(-80, -52, 50), np.full(50, -63)], np.c_[np.full(20, -52), np.linspace(-63, -54, 20)],
           np.c_[np.linspace(-52, -80, 50), np.full(50, -54)], np.c_[np.full(20, -80), np.linspace(-54, -63, 20)]]
ax.plot(bx[:, 0], bx[:, 1], color=STYLE["box"], lw=1.8, transform=ccrs.PlateCarree(), label="Analysis box", zorder=4)
for k, (a, b, c_, d_) in TR.items():
    if k != "v2":
        ax.plot([a, c_], [b, d_], color=STYLE["transport"], lw=0.6, alpha=0.5, transform=ccrs.Geodetic(), zorder=4)
la_s = PTS.loc["v2", "perm_southern_lat"]; lo_s = -65 + (la_s + 55) / (-63 + 55) * 5
ax.plot([-65, lo_s], [-55, la_s], color=STYLE["transport"], lw=3, transform=ccrs.Geodetic(), zorder=5,
        label="Section (ice-free points)")
ax.plot([lo_s, -60], [la_s, -63], color=STYLE["transport"], lw=2, ls=(0, (2, 1.5)), transform=ccrs.Geodetic(), zorder=5,
        label="Section (seasonal ice)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=7.5, frameon=False, ncol=1, handlelength=2.2)
panel_label(ax, "a", x=-0.14, y=1.08)

# (b) standardised trends
ax = fig.add_subplot(gs[0, 1])
def std_trend(key, y):
    r = A6c.loc[key] if key.startswith("T") else A6.loc[key]; s = r.slope / r.sigma_per_yr
    return r.sigma_per_yr, r.ci_cal_lo / s, r.ci_cal_hi / s, r.p_cal
items = [("Altimetry", "Tperm_v2_2001-2025"), ("In situ", "insitu_annual_2005-2019"), ("Wind 45–60°S", "tau_Circumpolar_45_60S_2001-2025")]
for i, (lab, key) in enumerate(items):
    v, lo, hi, p = std_trend(key, None)
    col = (STYLE["transport"] if i == 0 else STYLE["transport_obs"] if i == 1 else STYLE["wind"]) if p < 0.05 else STYLE["ns"]
    ax.bar(i, v, width=0.62, color=col, edgecolor="white", zorder=3)
    ax.errorbar(i, v, yerr=[[v - lo], [hi - v]], color="#111111", lw=1.0, capsize=2.8, capthick=1.0, zorder=4)
    ax.text(i, hi + 0.004, stars(p), ha="center", va="bottom", fontsize=10)
ax.plot(0.22, A6c.loc["Tall_v2_2001-2025"].sigma_per_yr, "D", ms=4.5, mfc="white", mec="#1F2933", mew=1.0, zorder=5,
        label="all points\n(incl. under ice)")
ax.legend(loc="lower right", fontsize=7, frameon=False, markerfirst=False, bbox_to_anchor=(1.04, 0.0))
ax.axhline(0, color="#7f8a96", lw=0.8, ls="--"); ax.set_xticks(range(3)); ax.set_xticklabels([i[0] for i in items], fontsize=8.5)
ax.set_ylabel(r"Standardized slope ($\sigma$ yr$^{-1}$)", fontsize=9); base_timeseries_style(ax); panel_label(ax, "b", x=-0.32, y=1.08)
ax.spines["bottom"].set_visible(False)

# (c) Drake-box EKE / KE trends
ax = fig.add_subplot(gs[0, 2])
for i, v in enumerate(("EKE", "KE")):
    r = A6b.loc[f"Drake_{v}_ice_any"]; rn = A6b.loc[f"Drake_{v}_nomask"]
    col = (STYLE["eke"] if v == "EKE" else STYLE["ke"]) if r.p_cal < 0.05 else STYLE["ns"]
    ax.bar(i, r.slope, width=0.62, color=col, edgecolor="white", zorder=3)
    ax.errorbar(i, r.slope, yerr=[[r.slope - r.ci_cal_lo], [r.ci_cal_hi - r.slope]], color="#111111", lw=1.0, capsize=2.8, capthick=1.0, zorder=4)
    ax.plot(i + 0.22, rn.slope, "D", ms=4.5, mfc="white", mec="#1F2933", mew=1.0, zorder=5, label="no ice mask" if i == 0 else None)
    ax.text(i, r.ci_cal_hi + 0.08, stars(r.p_cal), ha="center", va="bottom", fontsize=10)
ax.axhline(0, color="#7f8a96", lw=0.8, ls="--"); ax.set_xticks([0, 1]); ax.set_xticklabels(["EKE", "KE"], fontsize=9)
ax.set_ylabel(r"cm$^2$ s$^{-2}$ yr$^{-1}$", fontsize=9); base_timeseries_style(ax); panel_label(ax, "c", x=-0.32, y=1.08)
ax.legend(loc="upper right", fontsize=7.5, frameon=False); ax.spines["bottom"].set_visible(False)

# (d) standardised anomalies
ax = fig.add_subplot(gs[1, :])
z = lambda s_: (s_ - np.nanmean(s_)) / np.nanstd(s_)
Z = np.array([running(z(Ttr[f"Tperm_{k}"].values)) for k in TR])
ax.fill_between(x, np.nanmin(Z, 0), np.nanmax(Z, 0), color=STYLE["transport"], alpha=0.15, lw=0)
ax.plot(x, z(Ttr["Tperm_v2"].values), color=STYLE["transport"], alpha=STYLE["raw_alpha"], lw=STYLE["raw_lw"])
ax.plot(x, Z[0], color=STYLE["transport"], lw=3.0, label="Altimetry (surface proxy)")
ob = obs.sort_values("date"); xo = (ob.date.dt.year + (ob.date.dt.dayofyear - 0.5) / 365.25).values
ax.plot(xo, pd.Series(z(ob.total_transport_sv.values)).rolling(12, center=True, min_periods=6).mean(),
        color=STYLE["transport_obs"], lw=2.4, label="In situ (0–760 m ADCP)")
xw = (WS.year + (WS.month - 0.5) / 12).values
ax.plot(xw, running(z(WS["Circumpolar_45_60S"].values)), color=STYLE["wind"], lw=1.4, ls="--", label="Zonal wind stress 45–60°S")
ax.axhline(0, color="#6b7280", lw=1.0, ls="--")
ra = A6c.loc["Tperm_v2_2001-2025"]; ri = A6.loc["insitu_annual_2005-2019"]; rw = A6.loc["tau_Circumpolar_45_60S_2001-2025"]
ax.text(0.01, 0.04, f"Altimetry: {ra.pct_per_decade:+.1f} % dec$^{{-1}}$ (p={ra.p_cal:.2f});  in situ: {ri.pct_per_decade:+.1f} % dec$^{{-1}}$ "
        f"(p={ri.p_cal:.2f});  wind: {rw.pct_per_decade:+.1f} % dec$^{{-1}}$ (p={rw.p_cal:.2f});  "
        f"altimetry vs in situ: r = {AG.loc['Tperm_v2', 'r_transects']:.2f} (transects), {AG.loc['Tperm_v2', 'r_annual']:.2f} (annual)",
        transform=ax.transAxes, fontsize=7.6, color="#1f2937")
ax.set_ylabel(r"Standardized anomaly ($\sigma$)"); base_timeseries_style(ax)
ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=3); ax.set_xlim(2000.5, 2026.5)
yl = ax.get_ylim(); ax.set_ylim(yl[0] - 0.4, yl[1] + 0.6); panel_label(ax, "d", x=-0.06, y=1.08)

# (e) Drake-box energy
ax = fig.add_subplot(gs[2, :])
xb = BX.year + (BX.month - 0.5) / 12
for v, col in (("KE", STYLE["ke"]), ("EKE", STYLE["eke"])):
    r = A6b.loc[f"Drake_{v}_ice_any"]; y = BX[f"ice_any_{v}"].values
    ts_panel(ax, xb.values, y, col, label=v, slope=r.slope, intercept=theil_line(xb.values, deseason(BX, f"ice_any_{v}"), r.slope))
    ax.plot(xb, running(BX[f"nomask_{v}"].values), color=col, lw=1.0, ls=(0, (1.5, 1.2)))
ax.plot([], [], color="0.4", lw=1.0, ls=(0, (1.5, 1.2)), label="no ice mask")
rk, re = A6b.loc["Drake_KE_ice_any"], A6b.loc["Drake_EKE_ice_any"]
ax.text(0.01, 0.04, f"KE: {rk.slope:+.2f} cm$^2$ s$^{{-2}}$ yr$^{{-1}}$ (p={rk.p_cal:.3f})\nEKE: {re.slope:+.2f} cm$^2$ s$^{{-2}}$ yr$^{{-1}}$ "
        f"(p={re.p_cal:.3f}){stars(re.p_cal)}", transform=ax.transAxes, fontsize=8.5, color="#1f2937")
ax.set_xlabel("Year"); ax.set_ylabel(r"Energy (cm$^2$ s$^{-2}$)"); base_timeseries_style(ax)
ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=3); ax.set_xlim(2000.5, 2026.5)
yl = ax.get_ylim(); ax.set_ylim(yl[0] - 60, yl[1] + 40); panel_label(ax, "e", x=-0.06, y=1.08)
fig.savefig(FIG / "Fig4_Drake.png"); fig.savefig(FIG / "Fig4_Drake.pdf")
print("saved")
