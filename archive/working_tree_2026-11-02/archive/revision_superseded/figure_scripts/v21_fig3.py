"""v2.1 Figure 3: v2 Fig. 3 layout (a-e) with corrected inputs.
(a) Drake box (80-52W, 63-54S; KE, EKE and wind stress) and transport section (65W 55S - 60W 63S);
    solid: points ice-free in every month (used for the proxy); dashed: seasonally ice-covered part.
(b) standardised trends (sigma per yr, 95 % CI, revised p): altimetric surface transport proxy
    (fixed ice-free points, 2001-2025) and in-situ total ADCP transport 0-760 m (annual means, 2005-2019).
(c) Drake-box EKE and KE trends 2001-2025 (95 % CI, revised p).
(d) standardised transport anomalies (altimetry monthly and 12-month running mean; in situ 12-transect running
    mean) and eastward wind stress over the Drake box (12-month running mean, right axis).
(e) Drake-box EKE (monthly, 12-month running mean, trend).
Usage: python v21_fig3.py <repository outputs/trends dir>
"""
import sys
import numpy as np, pandas as pd
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

trends_dir = sys.argv[1]
FIGV3 = OUT / "figures_v21"; FIGV3.mkdir(exist_ok=True)
Ttr = pd.read_csv(OUT / "a6" / "a6c_transport_fixed_series.csv")
A6c = pd.read_csv(OUT / "a6" / "a6c_trends.csv"); A6c.index = A6c.series + "_" + A6c.period
A6c = retest_override(A6c, lambda i, r: f"{r.series}|{r.period}")
PTS = pd.read_csv(OUT / "a6" / "a6c_section_points.csv").set_index("line")
A6 = pd.read_csv(OUT / "a6" / "a6_trends.csv").set_index("series")
A6 = retest_override(A6, lambda i, r: {"insitu_annual_2005-2019": "insitu_annual|2005-2019"}.get(i) or
                     (f"{i[:-10]}|{i[-9:]}" if i.startswith("tau_") else None))
A6b = pd.read_csv(OUT / "a6" / "a6b_drake_box_trends.csv").query("period == '2001-2025'").set_index("series")
A6b = retest_override(A6b, lambda i, r: f"{i}|2001-2025")
BX = pd.read_csv(OUT / "a6" / "a6b_drake_box_series.csv")
WS = pd.read_csv(OUT / "a6" / "a6_wind_stress_series.csv"); WS = WS[(WS.year >= 2001) & (WS.year <= 2025)].reset_index(drop=True)
AG = pd.read_csv(OUT / "a6" / "a6c_agreement.csv", index_col=0)
obs = pd.read_csv(f"{trends_dir}/drake_passage_total_transport_observations.csv", parse_dates=["date"]).dropna(subset=["total_transport_sv"])
x = (Ttr.year + (Ttr.month - 0.5) / 12).values


def title(ax, s):
    ax.set_title(s, loc="left", fontsize=12, fontweight="bold")


fig = plt.figure(figsize=(10.0, 10.6))
gs = fig.add_gridspec(3, 3, height_ratios=[0.95, 1.0, 0.9], hspace=0.42, wspace=0.75)
ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo(central_longitude=-65))
ax.set_extent([-84, -46, -66, -51], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor="#111111", edgecolor="#111111", zorder=3)
ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.4, color="#7A8288", alpha=0.6, linestyle="--")
bx = np.r_[np.c_[np.linspace(-80, -52, 50), np.full(50, -63)], np.c_[np.full(20, -52), np.linspace(-63, -54, 20)],
           np.c_[np.linspace(-52, -80, 50), np.full(50, -54)], np.c_[np.full(20, -80), np.linspace(-54, -63, 20)]]
ax.plot(bx[:, 0], bx[:, 1], color=STYLE["box"], lw=1.8, transform=ccrs.PlateCarree(), label="Analysis box", zorder=4)
la_s = PTS.loc["v2", "perm_southern_lat"]; lo_s = -65 + (la_s + 55) / (-63 + 55) * 5
ax.plot([-65, lo_s], [-55, la_s], color=STYLE["transport"], lw=3, transform=ccrs.Geodetic(), zorder=5, label="Transport section")
ax.plot([lo_s, -60], [la_s, -63], color=STYLE["transport"], lw=1.5, ls=(0, (2, 1.5)), transform=ccrs.Geodetic(), zorder=5,
        label="(seasonal ice, excluded)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=7.5, frameon=False, handlelength=2.2)
title(ax, "a) Drake Passage")

ax = fig.add_subplot(gs[0, 1])
def sig_row(r):
    s = r.slope / r.sigma_per_yr
    return dict(pct=r.sigma_per_yr * 10, pct_lo=r.ci_cal_lo / s * 10, pct_hi=r.ci_cal_hi / s * 10, p=r.p_cal)
rows = [("Altimetry\n(surface, 2001–2025)", {"transport": sig_row(A6c.loc["Tperm_v2_2001-2025"])}),
        ("In-situ\n(0–760 m, 2005–2019)", {"transport": sig_row(A6.loc["insitu_annual_2005-2019"])})]
forest_panel(ax, rows, unit=r"$\sigma$ dec$^{-1}$", series_colors={"transport": STYLE["transport"]}, xlim=(-6, 4), label_size=8)
title(ax, "b) Transport Trends")

ax = fig.add_subplot(gs[0, 2])
rows = [(v, {"energy": retest_row(f"Drake_{v}_ice_any|2001-2025")}) for v in ("EKE", "KE")]
rows[0][1]["energy"]["p"] = A6b.loc["Drake_EKE_ice_any"].p_cal
forest_panel(ax, rows, series_colors={"energy": STYLE["eke"]}, xlim=(-6, 18))
title(ax, "c) KE and EKE Trends")

ax = fig.add_subplot(gs[1, :])
z = lambda s_: (s_ - np.nanmean(s_)) / np.nanstd(s_)
ax.plot(x, z(Ttr["Tperm_v2"].values), color=STYLE["transport"], alpha=STYLE["raw_alpha"], lw=STYLE["raw_lw"])
ax.plot(x, running(z(Ttr["Tperm_v2"].values)), color=STYLE["transport"], lw=3.0, label="Altimetry")
ob = obs.sort_values("date"); xo = (ob.date.dt.year + (ob.date.dt.dayofyear - 0.5) / 365.25).values
ax.plot(xo, pd.Series(z(ob.total_transport_sv.values)).rolling(12, center=True, min_periods=6).mean(),
        color=STYLE["transport_obs"], lw=2.4, label="In-situ")
ax.axhline(0, color="#6b7280", lw=1.0, ls="--")
ax.set_ylabel(r"Standardized transport anomaly ($\sigma$)"); base_timeseries_style(ax); ax.set_xlim(2000.5, 2026.5)
a2 = ax.twinx()
xw = (WS.year + (WS.month - 0.5) / 12).values
a2.plot(xw, running(WS["Drake_box"].values), color=STYLE["wind"], lw=1.4, ls="--", label="Eastward wind stress (Drake box)")
a2.set_ylabel(r"Eastward wind stress (N m$^{-2}$)", color=STYLE["wind"]); a2.tick_params(axis="y", colors=STYLE["wind"])
a2.spines["right"].set_visible(True); a2.spines["right"].set_color(STYLE["wind"]); a2.grid(False)
ra = A6c.loc["Tperm_v2_2001-2025"]; ri = A6.loc["insitu_annual_2005-2019"]; rw = A6.loc["tau_Drake_box_2001-2025"]
ax.text(0.01, 0.04, f"Altimetry: {ra.sigma_per_yr:+.3f} $\\sigma$ yr$^{{-1}}$ (p={ra.p_cal:.2f});  In-situ: {ri.sigma_per_yr:+.3f} $\\sigma$ yr$^{{-1}}$ "
        f"(p={ri.p_cal:.2f});  wind stress: {rw.pct_per_decade:+.1f} % dec$^{{-1}}$ (p={rw.p_cal:.2f})",
        transform=ax.transAxes, fontsize=8, color="#1f2937")
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8.5, frameon=False, ncol=3)
yl = ax.get_ylim(); ax.set_ylim(yl[0] - 0.4, yl[1] + 0.6)
title(ax, "d) Transport and Wind Stress Time Series (2001–2025)")

ax = fig.add_subplot(gs[2, :])
xb = (BX.year + (BX.month - 0.5) / 12).values
r = A6b.loc["Drake_EKE_ice_any"]; y = BX["ice_any_EKE"].values
ts_panel(ax, xb, y, STYLE["eke"], label="EKE", slope=r.slope, intercept=theil_line(xb, deseason(BX, "ice_any_EKE"), r.slope))
ax.text(0.01, 0.05, f"EKE: {r.slope:+.2f} cm$^2$ s$^{{-2}}$ yr$^{{-1}}$ (p={r.p_cal:.4f}){stars(r.p_cal)}", transform=ax.transAxes,
        fontsize=8.5, color="#1f2937")
ax.set_xlabel("Year"); ax.set_ylabel(r"Energy (cm$^2$ s$^{-2}$)"); base_timeseries_style(ax); ax.set_xlim(2000.5, 2026.5)
ax.legend(loc="upper left", fontsize=8.5, frameon=False)
title(ax, "e) EKE Time Series (2001–2025)")
fig.savefig(FIGV3 / "Figure3.png"); fig.savefig(FIGV3 / "Figure3.pdf")
print("saved")
