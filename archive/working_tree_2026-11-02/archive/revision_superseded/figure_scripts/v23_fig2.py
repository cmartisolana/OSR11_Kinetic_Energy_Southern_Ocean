"""v2.3 Figure 2: (c) now compares this study with earlier estimates over the same regions (a2g_literature_regions.py).
v2.2 Figure 2 (three panels: a detections, b circumpolar position, c sector trends; the 90E profile moves to the SI), recomputed with the parameters stated
in the Methods (200 cm2 s-2, 100 km, 25 %; a2_run.py baseline) and the matched fixed sample.
(a) half-power points detected in January 2001 (ACC band) over the monthly KE field;
(b) KE profile at 90°E, January 2001: envelope, peak, half-power point, threshold;
(c) ACC median-free matched half-power latitude (mean of per-longitude anomalies + climatological mean),
    raw monthly, 12-month running mean, deseasonalised Theil-Sen trend, revised p;
(d) trends by sector with 95 % CI and revised p.
Usage: OSR11_CACHE=<cache> OSR11_OUTPUTS=<repo outputs> python v21_fig2.py
"""
import os
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403
import a2_run as A
import fronts as F

FIGV3 = OUT / "figures_v23"; FIGV3.mkdir(exist_ok=True)
ke, lat, lon, acc, t = A.load()
P = A.BASE
k0 = 0  # January 2001
fld = np.asarray(ke[k0], float)
d = F.meridional_dists(lat)
cnt, tab = F.detect_field(fld, lat, d, P["thr"], P["min_w"], P["frac"], A.MAX_ENV)
FI = A.FI
jj, kk = np.nonzero(np.arange(A.MAX_ENV)[None, :] < cnt[:, None])
lhp = tab[jj, kk, FI["lat_hp"]]
inacc = A.in_acc(acc, lat, lon, lhp, jj)
ice = xr.open_dataset(A.REPO_OUT / "sea_ice_mask.nc").ice_mask.values > 0.5

dtr = np.load(OUT / "a2" / "tracks_baseline.npz", allow_pickle=True)
keep = dtr["keep"] & dtr["acc"]
sec = A.sector(lon[dtr["j"][keep]])
lat_hp = dtr["lat_hp"][:, keep].astype(float)
S = pd.read_csv(OUT / "a2" / "a2_summary.csv"); S = S[S.set == "baseline"].set_index("region")
S = retest_override(S, lambda g, r: f"baseline|{g}|lat_hp", {"m_lat_hp_p": "p_new", "m_lat_hp_lo": "ci_lo", "m_lat_hp_hi": "ci_hi"})
x = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25)


def title(ax, s):
    ax.set_title(s, loc="left", fontsize=12, fontweight="bold")


fig = plt.figure(figsize=(12.0, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.25, 1.15], wspace=0.32)
ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo())
polar_ax_style(ax, ccrs, cfeature)
LON, LAT = np.meshgrid(lon, lat)
ax.pcolormesh(LON, LAT, np.where(ice, np.nan, fld * 1e4), transform=ccrs.PlateCarree(), cmap="Greys", vmin=0, vmax=1500,
              shading="auto", rasterized=True, zorder=0)
add_ice_mask(ax, ccrs, lon, lat, ice)
ax.scatter(lon[jj[inacc]], lhp[inacc], s=2.2, c="#C0392B", lw=0, transform=ccrs.PlateCarree(), zorder=5)
J90 = int(np.argmin(np.abs(lon - 90.0)))
title(ax, "a) Half-power points, Jan 2001")

# (c) matched position series
ax = fig.add_subplot(gs[0, 1])
y = np.nanmean(lat_hp - np.nanmean(lat_hp, 0, keepdims=True), 1) + np.nanmean(lat_hp)
r = S.loc["ACC"]
yd = pd.Series(y); yd = (yd - yd.groupby(month).transform("mean") + yd.mean()).values
ts_panel(ax, x, y, STYLE["lat"], slope=r.m_lat_hp_slope, intercept=theil_line(x, yd, r.m_lat_hp_slope))
ax.text(0.98, 0.04, f"{r.m_lat_hp_slope:+.4f} ° yr$^{{-1}}$ [{r.m_lat_hp_lo:+.4f}, {r.m_lat_hp_hi:+.4f}], p={r.m_lat_hp_p:.2f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=STYLE["neutral"])
ax.set_xlabel("Year"); ax.set_ylabel("Latitude (°)"); base_timeseries_style(ax)
title(ax, "b) Circumpolar-mean position")

# (c) this study vs earlier estimates over the same regions (km per decade, + north)
ax = fig.add_subplot(gs[0, 2])
KM = 111.195
G = pd.read_csv(OUT / "a2" / "a2g_literature_regions.csv").set_index("region")
def g(name):
    r = G.loc[name]; return dict(pct=r.km_dec, pct_lo=r.lo, pct_hi=r.hi, p=r.p)
def sec(reg):
    q = S.loc[reg]; return dict(pct=q.m_lat_hp_slope * KM * 10, pct_lo=q.m_lat_hp_lo * KM * 10, pct_hi=q.m_lat_hp_hi * KM * 10, p=q.m_lat_hp_p)
TS = "this study, 2001–2025"
rows = [("Circumpolar", {TS: g("Circumpolar")}), ("Atlantic", {TS: sec("Atlantic")}), ("Indian", {TS: sec("Indian")}),
        ("Pacific", {TS: sec("Pacific")}), ("90–170°E", {TS: g("90-170E (Chambers 2018)")}),
        ("160–90°W", {TS: g("160-90W (Chambers 2018)")}), ("PAR meander\n150–80°W", {TS: g("150-80W PAR meander (Liu et al. 2026)")})]
forest_panel(ax, rows, unit=r"km dec$^{-1}$, + north", series_colors={TS: STYLE["lat"]}, xlim=(-45, 50), label_size=9, value_size=7.5)
ax.get_legend().remove() if ax.get_legend() else None
LIT = [("Sokolov and Rintoul (2009b), SSH contours, 1992–2007", 0, -40, "v", "#7B2D26"),
       ("Gille (2014), SSH contours, ~1993–2012", 0, -30, "^", "#B5651D"),
       ("Chambers (2018), along-track half-power, 1993–2016", 4, -29, "P", "#5B5F97"),
       ("Chambers (2018)", 5, 46, "P", "#5B5F97"),
       ("Liu et al. (2026), PAR meander, 1993–2023", 6, 3.3, "X", "#2E7D32")]
seen = set()
for lab, row, val, mk, col in LIT:
    lab2 = lab if lab.split(",")[0] not in seen else None; seen.add(lab.split(",")[0])
    ax.plot(val, -row, mk, ms=8, color=col, mec="white", mew=0.6, zorder=6, label=lab2)
ax.legend(loc="upper center", bbox_to_anchor=(0.2, -0.13), fontsize=7, frameon=False, ncol=1, handletextpad=0.3)
title(ax, "c) Trends compared with earlier estimates")
fig.savefig(FIGV3 / "Figure2.png"); fig.savefig(FIGV3 / "Figure2.pdf")
print("saved; Jan 2001 ACC detections:", int(inacc.sum()))
