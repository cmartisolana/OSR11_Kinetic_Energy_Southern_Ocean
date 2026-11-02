"""v2.4 Figure 2 (five panels). Changes from v23_fig2.py (author comments, 2026-10-15):
layout redistributed (two rows); new (b) example half-power profile (moved back from SI S3); this study's
trends (d) separated from the comparison with earlier estimates (e).
(a) half-power points detected in January 2001 (ACC band) over the monthly KE field; star = longitude of (b);
(b) meridional KE profile at the (b) longitude, January 2001: 200 cm2 s-2 threshold, segment, peak, 25 % level,
    envelope bounds, half-power latitude; dashed: 2001-2025 mean profile at the same longitude;
(c) circumpolar-mean matched half-power latitude, monthly, 12-month running mean, trend, revised p;
(d) trends (km per decade, + north) for the circumpolar mean and the sectors;
(e) this study over the regions of earlier studies vs their published estimates (a2g_literature_regions.csv).
Usage: OSR11_CACHE=<cache> OSR11_OUTPUTS=<repo outputs> python v24_fig2.py
"""
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403
import a2_run as A
import fronts as F

FIGV = OUT / "figures_v24"; FIGV.mkdir(exist_ok=True)
LON_EX = 147.9          # example longitude for (b): single, wide, matched envelope south of Tasmania in Jan 2001
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
lat_hp = dtr["lat_hp"][:, keep].astype(float)
S = pd.read_csv(OUT / "a2" / "a2_summary.csv"); S = S[S.set == "baseline"].set_index("region")
S = retest_override(S, lambda g, r: f"baseline|{g}|lat_hp", {"m_lat_hp_p": "p_new", "m_lat_hp_lo": "ci_lo", "m_lat_hp_hi": "ci_hi"})
x = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25)
KM = 111.195
G = pd.read_csv(OUT / "a2" / "a2g_literature_regions.csv").set_index("region")


def title(ax, s):
    ax.set_title(s, loc="left", fontsize=12, fontweight="bold")


def g(name):
    r = G.loc[name]; return dict(pct=r.km_dec, pct_lo=r.lo, pct_hi=r.hi, p=r.p)


def sec(reg):
    q = S.loc[reg]; return dict(pct=q.m_lat_hp_slope * KM * 10, pct_lo=q.m_lat_hp_lo * KM * 10, pct_hi=q.m_lat_hp_hi * KM * 10, p=q.m_lat_hp_p)


fig = plt.figure(figsize=(13.0, 10.8))
outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.78], hspace=0.50)
top = outer[0].subgridspec(1, 3, width_ratios=[1.0, 1.05, 1.25], wspace=0.30)
bot = outer[1].subgridspec(1, 2, width_ratios=[0.8, 1.2], wspace=0.55)

# (a) detections, Jan 2001
ax = fig.add_subplot(top[0, 0], projection=ccrs.SouthPolarStereo())
polar_ax_style(ax, ccrs, cfeature)
LON, LAT = np.meshgrid(lon, lat)
ax.pcolormesh(LON, LAT, np.where(ice, np.nan, fld * 1e4), transform=ccrs.PlateCarree(), cmap="Greys", vmin=0, vmax=1500,
              shading="auto", rasterized=True, zorder=0)
add_ice_mask(ax, ccrs, lon, lat, ice)
ax.scatter(lon[jj[inacc]], lhp[inacc], s=2.2, c="#C0392B", lw=0, transform=ccrs.PlateCarree(), zorder=5)
JEX = int(np.argmin(np.abs(lon - LON_EX)))
ex = tab[JEX, :cnt[JEX]]
ex = ex[A.in_acc(acc, lat, lon, ex[:, FI["lat_hp"]], np.full(len(ex), JEX))][0]
ax.plot(lon[JEX], ex[FI["lat_hp"]], "*", ms=13, color="#F4A261", mec="#1F2933", mew=0.8, transform=ccrs.PlateCarree(), zorder=7)
title(ax, "a) Half-power points, Jan 2001")

# (b) example profile
ax = fig.add_subplot(top[0, 1])
prof = fld[:, JEX] * 1e4
clim = np.nanmean(np.asarray(ke[:, :, JEX], float), 0) * 1e4
s, e = int(ex[FI["i_start"]]), int(ex[FI["i_end"]])
ls_, ln_, pk_lat, pv, hp = ex[FI["lat_s"]], ex[FI["lat_n"]], ex[FI["peak_lat"]], ex[FI["peak_val"]] * 1e4, ex[FI["lat_hp"]]
ax.plot(lat, clim, color="#7A8288", lw=1.3, ls="--", label="2001–2025 mean")
ax.plot(lat, prof, color="#1F2933", lw=1.8, label="January 2001")
sel = (lat >= ls_) & (lat <= ln_)
ax.fill_between(lat[sel], 0, prof[sel], color=STYLE["lat"], alpha=0.25, lw=0, label="envelope (integrated)")
ax.axhline(P["thr"] * 1e4, color="#C0392B", lw=1.0, ls=":", label="threshold, 200 cm$^2$ s$^{-2}$")
ax.plot([lat[s], lat[e]], [P["thr"] * 1e4] * 2, color="#C0392B", lw=3.0, solid_capstyle="butt", alpha=0.6)
ax.plot([lat[s] - 0.6, lat[e] + 0.6], [P["frac"] * pv] * 2, color="#E07A5F", lw=1.0, ls="-.", label="25 % of peak")
ax.plot(pk_lat, pv, "o", ms=6, color="#1F2933")
ax.axvline(hp, color=STYLE["lat"], lw=2.0, label=r"half-power latitude $\varphi_{hp}$")
for v, lab in ((ls_, r"$\varphi_S$"), (ln_, r"$\varphi_N$")):
    ax.axvline(v, color=STYLE["lat"], lw=0.9, ls="--")
    ax.text(v + (-0.15 if lab.endswith("S$") else 0.15), pv * 0.72, lab, ha="right" if lab.endswith("S$") else "left",
            va="center", fontsize=11, color=STYLE["lat"])
ax.set_xlim(hp - 7, hp + 7); ax.set_ylim(0, pv * 1.12)
ax.set_xlabel("Latitude (°)"); ax.set_ylabel(r"KE (cm$^2$ s$^{-2}$)"); base_timeseries_style(ax)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=7, frameon=False, handlelength=1.8, columnspacing=0.8)
title(ax, f"b) KE profile at {lon[JEX]:.1f}°E")

# (c) matched position series
ax = fig.add_subplot(top[0, 2])
y = np.nanmean(lat_hp - np.nanmean(lat_hp, 0, keepdims=True), 1) + np.nanmean(lat_hp)
r = S.loc["ACC"]
yd = pd.Series(y); yd = (yd - yd.groupby(month).transform("mean") + yd.mean()).values
ts_panel(ax, x, y, STYLE["lat"], slope=r.m_lat_hp_slope, intercept=theil_line(x, yd, r.m_lat_hp_slope))
ax.text(0.98, 0.03, f"{r.m_lat_hp_slope * KM * 10:+.1f} km dec$^{{-1}}$ [{r.m_lat_hp_lo * KM * 10:+.1f}, {r.m_lat_hp_hi * KM * 10:+.1f}], p={r.m_lat_hp_p:.2f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=STYLE["neutral"])
ax.set_xlabel("Year"); ax.set_ylabel("Latitude (°)"); base_timeseries_style(ax)
title(ax, "c) Circumpolar-mean position")

# (d) this study: circumpolar and sectors
ax = fig.add_subplot(bot[0, 0])
TS = "this study, 2001–2025"
rows = [("Circumpolar", {TS: g("Circumpolar")}), ("Atlantic", {TS: sec("Atlantic")}), ("Indian", {TS: sec("Indian")}),
        ("Pacific", {TS: sec("Pacific")})]
forest_panel(ax, rows, unit=r"km dec$^{-1}$, + north", series_colors={TS: STYLE["lat"]}, xlim=(-15, 15), label_size=9.5, value_size=8)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=8, frameon=False, ncol=1, handletextpad=0.3)
title(ax, "d) Position trends, 2001–2025")

# (e) comparison with earlier estimates over the same regions
ax = fig.add_subplot(bot[0, 1])
rows = [("Circumpolar", {TS: g("Circumpolar")}), ("90–170°E", {TS: g("90-170E (Chambers 2018)")}),
        ("160–90°W", {TS: g("160-90W (Chambers 2018)")}), ("PAR meander\n150–80°W", {TS: g("150-80W PAR meander (Liu et al. 2026)")})]
forest_panel(ax, rows, unit=r"km dec$^{-1}$, + north", series_colors={TS: STYLE["lat"]}, xlim=(-45, 50), label_size=9.5, value_size=8)
ax.get_legend().remove() if ax.get_legend() else None
for ln in list(ax.lines):
    if ln.get_label() in (TS, "not significant (p ≥ 0.05)"):
        ln.remove()
LIT = [("Sokolov and Rintoul (2009b), SSH contours, 1992–2007", 0, -40, "v", "#7B2D26"),
       ("Gille (2014), SSH contours, ~1993–2012", 0, -30, "^", "#B5651D"),
       ("Chambers (2018), along-track half-power, 1993–2016", 1, -29, "P", "#5B5F97"),
       ("Chambers (2018)", 2, 46, "P", "#5B5F97"),
       ("Liu et al. (2026), PAR meander, 1993–2023", 3, 3.3, "X", "#2E7D32")]
ax.plot([], [], "o", color=STYLE["lat"], mfc="white", mew=1.4, ms=6.5, label="this study, 2001–2025 (open: p ≥ 0.05)")
seen = set()
for lab, row, val, mk, col in LIT:
    lab2 = lab if lab.split(",")[0] not in seen else None; seen.add(lab.split(",")[0])
    ax.plot(val, -row - 0.28, mk, ms=9, color=col, mec="white", mew=0.6, zorder=6, label=lab2)
ax.legend(loc="upper center", bbox_to_anchor=(0.4, -0.13), fontsize=7.5, frameon=False, ncol=2, handletextpad=0.3, columnspacing=1.0)
title(ax, "e) Compared with earlier estimates")
fig.savefig(FIGV / "Figure2.png"); fig.savefig(FIGV / "Figure2.pdf")
print("saved; example lon", lon[JEX], "phi_hp", hp, "phi_S", ls_, "phi_N", ln_, "peak", pv)
