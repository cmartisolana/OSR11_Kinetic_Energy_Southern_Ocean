"""v2.1 Figure 2: v2 Fig. 2 layout without the width panels (v2 d, f), recomputed with the parameters stated
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

FIGV3 = OUT / "figures_v21"; FIGV3.mkdir(exist_ok=True)
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


fig = plt.figure(figsize=(10.6, 10.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.8], width_ratios=[1.5, 1], hspace=0.32, wspace=0.22)
ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo())
polar_ax_style(ax, ccrs, cfeature)
LON, LAT = np.meshgrid(lon, lat)
ax.pcolormesh(LON, LAT, np.where(ice, np.nan, fld * 1e4), transform=ccrs.PlateCarree(), cmap="Greys", vmin=0, vmax=1500,
              shading="auto", rasterized=True, zorder=0)
add_ice_mask(ax, ccrs, lon, lat, ice)
ax.scatter(lon[jj[inacc]], lhp[inacc], s=2.2, c="#C0392B", lw=0, transform=ccrs.PlateCarree(), zorder=5)
J90 = int(np.argmin(np.abs(lon - 90.0)))
ax.plot(np.full(50, lon[J90]), np.linspace(-65, -35, 50), color="#3B3B98", lw=1.4, ls=(0, (4, 2)), transform=ccrs.PlateCarree(), zorder=6)
title(ax, "a) Half-power Point Monthly Detection (Jan 2001)")

# (b) profile at 90E
ax = fig.add_subplot(gs[0, 1])
prof = fld[:, J90] * 1e4
sel = (lat >= -60) & (lat <= -35)
ax.plot(prof[sel], lat[sel], color="#C0392B", lw=1.4, label="KE profile")
rows = [tab[J90, c] for c in range(cnt[J90]) if A.in_acc(acc, lat, lon, np.array([tab[J90, c, FI["lat_hp"]]]), np.array([J90]))[0]]
if rows:
    r = max(rows, key=lambda q: q[FI["peak_val"]])
    ax.axhspan(r[FI["lat_s"]], r[FI["lat_n"]], color="#E8A0A0", alpha=0.35, lw=0, label="Envelope")
    ax.plot(r[FI["peak_val"]] * 1e4, r[FI["peak_lat"]], "^", ms=9, color="#B22222", mec="#111111", label="Envelope peak")
    php = np.interp(r[FI["lat_hp"]], lat, fld[:, J90]) * 1e4
    ax.plot(php, r[FI["lat_hp"]], "*", ms=13, color="#F4B400", mec="#3B3B98", label="Half-power point")
ax.axvline(P["thr"] * 1e4, color="#555555", ls="--", lw=1.0, label="Threshold (200 cm$^2$ s$^{-2}$)")
ax.set_xlabel(r"KE (cm$^2$ s$^{-2}$)"); ax.set_ylabel("Latitude (°)"); base_timeseries_style(ax)
ax.legend(loc="lower right", fontsize=7.5, frameon=False)
title(ax, "b) Latitudinal Section at 90°E")

# (c) matched position series
ax = fig.add_subplot(gs[1, 0])
y = np.nanmean(lat_hp - np.nanmean(lat_hp, 0, keepdims=True), 1) + np.nanmean(lat_hp)
r = S.loc["ACC"]
yd = pd.Series(y); yd = (yd - yd.groupby(month).transform("mean") + yd.mean()).values
ts_panel(ax, x, y, STYLE["lat"], slope=r.m_lat_hp_slope, intercept=theil_line(x, yd, r.m_lat_hp_slope))
ax.text(0.98, 0.04, f"{r.m_lat_hp_slope:+.4f} ° yr$^{{-1}}$ [{r.m_lat_hp_lo:+.4f}, {r.m_lat_hp_hi:+.4f}], p={r.m_lat_hp_p:.2f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=STYLE["neutral"])
ax.set_xlabel("Year"); ax.set_ylabel("Latitude (°)"); base_timeseries_style(ax)
title(ax, "c) Mean Half-power Point Position (2001–2025)")

# (d) sector trends (forest, km per decade)
ax = fig.add_subplot(gs[1, 1])
KM = 111.195
rows = []
for reg in ("ACC", "Pacific", "Atlantic", "Indian"):
    q = S.loc[reg]
    rows.append((reg, {"position": dict(pct=q.m_lat_hp_slope * KM * 10, pct_lo=q.m_lat_hp_lo * KM * 10, pct_hi=q.m_lat_hp_hi * KM * 10, p=q.m_lat_hp_p)}))
forest_panel(ax, rows, unit=r"km dec$^{-1}$, + north", series_colors={"position": STYLE["lat"]}, xlim=(-12, 12))
title(ax, "d) Regional Half-power Point Trends")
fig.savefig(FIGV3 / "Figure2.png"); fig.savefig(FIGV3 / "Figure2.pdf")
print("saved; Jan 2001 ACC detections:", int(inacc.sum()))
