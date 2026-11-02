"""Figure 3 (v2.6, style_nc): Drake Passage and transport proxies, Product 1.
(a) Drake Passage: mean KE 2001-2025 (shading), transect with the mean cross-track surface geostrophic velocity
    (arrows; solid part: points free of ice in every month, used for the proxy), domain of the wind stress / KE / EKE series;
(b) mean KE and EKE along the transect vs distance from its northern end, 2001-2025;
(c) trends of KE and EKE along the transect (50 km bins, % per decade, 2001-2025);
(d) trends of the cross-track velocity along the transect (cm s-1 per decade), 2001-2025 and October 2005-April 2019,
    with the zones where Gutierrez-Villanueva et al. (2023) found acceleration / deceleration of the 0-760 m transport;
(e) Drake transport proxy, circumpolar proxy, Drake-domain KE and EKE, 12-month running means in % of the mean;
(f) transport proxy trends (in situ: trend of the 0-760 m transport only); (g) Drake KE, EKE and wind stress trends.
Usage: python fig3.py <repository outputs/trends dir>   (OSR11_CACHE for the KE map)
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.colors import LogNorm
from style_nc import *  # noqa
from matplotlib.lines import Line2D

trends_dir = sys.argv[1]
FIG = OUT / "figures_final"
RT = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")
Ttr = pd.read_csv(OUT / "a6" / "a6c_transport_fixed_series.csv")
BX = pd.read_csv(OUT / "a6" / "a6b_drake_box_series.csv")
WS = pd.read_csv(OUT / "a6" / "a6_wind_stress_series.csv"); WS = WS[(WS.year >= 2001) & (WS.year <= 2025)].reset_index(drop=True)
A8 = pd.read_csv(OUT / "a8" / "a8_series_allsat.csv"); A8T = pd.read_csv(OUT / "a8" / "a8_trends_allsat.csv").set_index("series")
SM = pd.read_csv(OUT / "a6" / "a6d_section_mean.csv")
ST = pd.read_csv(OUT / "a6" / "a6d_section_trends.csv").set_index("key")
obs = pd.read_csv(f"{trends_dir}/drake_passage_total_transport_observations.csv", parse_dates=["date"]).dropna(subset=["total_transport_sv"])
ann = obs.groupby(obs.date.dt.year).total_transport_sv.mean(); ann = ann[(ann.index >= 2005) & (ann.index <= 2019)]
cache = Path(os.environ.get("OSR11_CACHE", HERE / "cache"))
x = (Ttr.year + (Ttr.month - 0.5) / 12).values
C_DP, C_ACC, C_EKE, C_KE, C_W = C["tr"], C["ink"], C["eke"], C["ke"], C["wind"]
pa = lambda y: 100 * (y / np.nanmean(y) - 1)
PC = ccrs.PlateCarree()


def pct_row(key, mean):
    r = RT.loc[key]; f = 1000 / mean
    return dict(pct=r.slope * f, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p_new)


def a8_row(s):
    r = A8T.loc[s]; return dict(pct=r.pct_per_decade, pct_lo=r.pct_lo, pct_hi=r.pct_hi, p=r.p)


def prof(v, period):
    k = [i for i in ST.index if i.startswith(f"{v}|") and i.endswith(f"|{period}")]
    d = ST.loc[k].copy(); d["dist"] = [float(i.split("|")[1]) for i in k]
    return d.dropna(subset=["tr"]).sort_values("dist")


def front_marks(ax, ymax=None):
    for dkm, lab in FRONTS:
        ax.axvline(dkm, color="#9A9A9A", lw=0.5, ls=":", zorder=0)
        ax.text(dkm + 8, 0.98, lab, transform=ax.get_xaxis_transform(), ha="left", va="top", fontsize=5.5, color=C["grey"])


# front positions on the transect: local maxima of the mean cross-track velocity (SAF, PF, SACCF)
vn = SM.vn.values; dist = SM.dist_km.values
FRONTS = []
for lab, (a, b) in (("SAF", (0, 120)), ("PF", (300, 500)), ("SACCF", (700, 850))):
    s = (dist >= a) & (dist <= b); FRONTS.append((float(dist[s][np.argmax(vn[s])]), lab))

fig = plt.figure(figsize=(W2, 7.3))
g = fig.add_gridspec(4, 2, height_ratios=[1.05, 0.8, 0.62, 1.0], hspace=0.62, wspace=0.32, width_ratios=[0.9, 1.1])
gb = g[3, :].subgridspec(1, 2, wspace=0.75)
perm = SM.ice_months.values == 0
DX_ = "Distance from northern end (km)"

# a map
ax = fig.add_subplot(g[0, 0], projection=ccrs.LambertConformal(central_longitude=-63, standard_parallels=(-55, -62)))
ax.set_extent([-73, -52, -65.5, -52.5], crs=PC)
KM_ = xr.open_dataset(OUT / "a6" / "a6d_drake_mean_ke.nc").ke_mean
im = ax.pcolormesh(KM_.longitude, KM_.latitude, KM_.values, transform=PC, cmap="YlGnBu", norm=LogNorm(30, 3000), shading="auto",
                   rasterized=True, zorder=0)
ax.add_feature(cfeature.LAND, facecolor="#3A3A3A", edgecolor="none", zorder=3)
gl = ax.gridlines(crs=PC, draw_labels=False, linewidth=0.25, color="#9A9A9A", alpha=0.6, linestyle=":")
gl.xlocator = plt.FixedLocator([-70, -65, -60, -55]); gl.ylocator = plt.FixedLocator([-64, -60, -56])
bx = np.r_[np.c_[np.linspace(-80, -52, 50), np.full(50, -63)], np.c_[np.full(20, -52), np.linspace(-63, -54, 20)],
           np.c_[np.linspace(-52, -80, 50), np.full(50, -54)], np.c_[np.full(20, -80), np.linspace(-54, -63, 20)]]
ax.plot(bx[:, 0], bx[:, 1], color=C["ink"], lw=0.6, ls=(0, (3, 1.5)), transform=PC, zorder=4, label="Domain (wind, KE, EKE)")
ax.plot(SM.lon[perm], SM.lat[perm], color=C_DP, lw=1.3, transform=PC, zorder=5, label="Transect, ice-free part")
ax.plot(SM.lon[~perm], SM.lat[~perm], color=C_DP, lw=0.8, ls=(0, (1.5, 1)), transform=PC, zorder=5, label="Transect under seasonal ice")
R_E = 6.371e6; latr = np.radians(-59.0)
DXm = R_E * np.cos(latr) * np.radians(5.0); DYm = R_E * np.radians(-8.0); L = np.hypot(DXm, DYm); nx, ny = DYm / L, -DXm / L
if nx < 0: nx, ny = -nx, -ny
k = np.arange(0, len(SM), 4)
q = ax.quiver(SM.lon.values[k], SM.lat.values[k], SM.vn.values[k] * nx, SM.vn.values[k] * ny, transform=PC, color=C["ink"],
              scale=2.2, width=0.005, headwidth=3.5, headlength=4, zorder=6)
ax.quiverkey(q, 0.84, 0.32, 0.25, "25 cm s$^{-1}$", labelpos="N", fontproperties=dict(size=5.5), coordinates="axes")
cb = fig.colorbar(im, ax=ax, orientation="vertical", shrink=0.8, pad=0.03, aspect=22)
cb.set_label("Mean KE (cm$^2$ s$^{-2}$)", fontsize=6); cb.ax.tick_params(labelsize=5.5, width=0.4, length=2); cb.outline.set_linewidth(0.4)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=2, fontsize=5.5)
plabel(ax, "a", "Drake Passage")

# b mean distribution along the transect
ax = fig.add_subplot(g[0, 1])
for v, col in (("KE", C_KE), ("EKE", C_EKE)):
    ax.plot(dist, SM[v], color=col, lw=1.0, label=v)
ax.axvspan(dist[perm].max(), dist.max(), color=C["ice"], lw=0, zorder=0)
front_marks(ax)
ax.set_xlim(0, dist.max()); ax.set_ylim(0, 1300); ax.set_xlabel(DX_); ax.set_ylabel("cm$^2$ s$^{-2}$")
ax.legend(loc="center right")
plabel(ax, "b", "Mean KE and EKE along the transect, 2001–2025", dx=-30)


def trend_profile(ax, rows):
    zero(ax)
    for d, col, mk, off, lab in rows:
        sig = (d.p < 0.05).values; xx = (d.dist + off).values
        ax.vlines(xx, d.lo, d.hi, color=col, lw=0.7)
        ax.plot(xx[sig], d.tr.values[sig], mk, color=col, ms=3.0, mew=0.7, label=lab)
        ax.plot(xx[~sig], d.tr.values[~sig], mk, mfc="white", color=col, ms=3.0, mew=0.7)
    ax.axvspan(dist[perm].max(), dist.max(), color=C["ice"], lw=0, zorder=0)
    front_marks(ax); ax.set_xlim(0, dist.max()); ax.set_xlabel(DX_)


ax = fig.add_subplot(g[1, 0])
trend_profile(ax, [(prof("KE", "2001-2025"), C_KE, "o", -6, "KE"), (prof("EKE", "2001-2025"), C_EKE, "o", 6, "EKE")])
ax.set_ylabel("% per decade"); ax.legend(loc="lower left", ncol=2)
plabel(ax, "c", "KE and EKE trends along the transect, 2001–2025", dx=-30)
ax = fig.add_subplot(g[1, 1])
for (a, b), col in (((150, 300), C["ke"]), ((350, 550), C["tr"])):
    ax.axvspan(a, b, color=col, alpha=0.12, lw=0, zorder=0.5)
trend_profile(ax, [(prof("vn", "2001-2025"), C_DP, "o", -6, "2001–2025"), (prof("vn", "2005-2019"), C["obs"], "s", 6, "Oct 2005–Apr 2019")])
ax.set_ylabel("cm s$^{-1}$ per decade"); ax.legend(loc="lower right", ncol=1)
ax.text(225, 0.03, "G-V23\nacceleration", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.3, color=C["ke"])
ax.text(450, 0.03, "G-V23\ndeceleration", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.3, color=C["tr"])
plabel(ax, "d", "Cross-track velocity trends along the transect", dx=-30)

# e series
ax = fig.add_subplot(g[2, :])
xa = (A8.year + (A8.month - 0.5) / 12).values; xb = (BX.year + (BX.month - 0.5) / 12).values
ax.plot(x, running(pa(Ttr.Tperm_v2.values)), color=C_DP, lw=1.0, label="Drake Passage transport proxy")
ax.plot(xa, running(pa(A8.ACC_T.values)), color=C_DP, lw=1.0, ls=(0, (3, 1.5)), label="Circumpolar transport proxy")
ax.plot(xb, running(pa(BX.ice_any_EKE.values)), color=C_EKE, lw=1.0, label="Drake Passage EKE")
ax.plot(xb, running(pa(BX.ice_any_KE.values)), color=C_KE, lw=1.0, label="Drake Passage KE")
zero(ax); ts_style(ax); ax.set_ylabel("% of mean")
ax.legend(loc="upper left", ncol=4)
yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.25 * (yl[1] - yl[0]))
plabel(ax, "e", "Transport proxies and Drake Passage KE and EKE, 12-month running means", dx=-24)

# f, g trends
mean_dp = float(np.nanmean(Ttr.Tperm_v2.values)); mean_obs = float(ann.mean())
mBX = {v: float(np.nanmean(BX[f"ice_any_{v}"])) for v in ("EKE", "KE")}
mW = {c: float(WS[c].mean()) for c in ("Drake_box", "Circumpolar_45_60S")}
T, E, K, W_ = "transport proxy", "EKE", "KE", "eastward wind stress"
rows_t = [("Drake altimetric, 2001–2025", {T: pct_row("Tperm_v2|2001-2025", mean_dp)}),
          ("Drake altimetric, 2005–2019", {T: pct_row("Tperm_v2|2005-2019", mean_dp)}),
          ("Drake in situ, 2005–2019", {T: pct_row("insitu_annual|2005-2019", mean_obs)}),
          ("Circumpolar", {T: a8_row("ACC_T")}),
          ("Pacific sector", {T: a8_row("Pacific_T")}), ("Atlantic sector", {T: a8_row("Atlantic_T")}),
          ("Indian sector", {T: a8_row("Indian_T")})]
rows_e = [("Drake Passage EKE", {E: pct_row("Drake_EKE_ice_any|2001-2025", mBX["EKE"])}),
          ("Drake Passage KE", {K: pct_row("Drake_KE_ice_any|2001-2025", mBX["KE"])}),
          ("Wind stress, Drake Passage", {W_: pct_row("tau_Drake_box|2001-2025", mW["Drake_box"])}),
          ("Wind stress, 45–60°S", {W_: pct_row("tau_Circumpolar_45_60S|2001-2025", mW["Circumpolar_45_60S"])})]
ax = fig.add_subplot(gb[0, 0])
forest(ax, rows_t, {T: C_DP}, {T: "s"}, xlim=(-24, 34), xlabel="Trend (% per decade)")
ax.axhline(-2.5, color=C["light"], lw=0.5)
ax.legend(handles=sig_handles({"Transport proxy": C_DP}, {"Transport proxy": "s"}), loc="upper center", bbox_to_anchor=(0.4, -0.3), ncol=2)
plabel(ax, "f", "Transport proxy trends", dx=-80)
ax = fig.add_subplot(gb[0, 1])
forest(ax, rows_e, {E: C_EKE, K: C_KE, W_: C_W}, {W_: "D"}, xlim=(-4, 18), dodge=False, xlabel="Trend 2001–2025 (% per decade)")
ax.axhline(-1.5, color=C["light"], lw=0.5)
ax.legend(handles=sig_handles({"EKE": C_EKE, "KE": C_KE, "Wind stress": C_W}, {"Wind stress": "D"}), loc="upper center",
          bbox_to_anchor=(0.4, -0.3), ncol=4)
plabel(ax, "g", "Drake Passage KE, EKE and wind stress", dx=-80)
save(fig, FIG, "Figure3")
print("fronts (km):", FRONTS, "ice-free to km", dist[perm].max())
