"""Figure 2 (v2.6): half-power positions of the ACC KE envelope (style_nc).
a half-power points detected in January 2001 (ACC band) over the monthly KE field; star = longitude of b;
b meridional KE profile at that longitude, January 2001, with threshold, envelope, 25 % level and half-power latitude;
c circumpolar-mean matched half-power latitude, monthly and 12-month running mean, with the Theil-Sen trend;
d trends (km per decade, + north) for the circumpolar mean, the sectors and the regions of earlier studies, with the
  published estimates (a2g_literature_regions.csv).
Usage: OSR11_CACHE=<cache> OSR11_OUTPUTS=<repo outputs> python fig2.py
"""
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from style_nc import *  # noqa
import a2_run as A
import fronts as F

FIGD = OUT / "figures_final"
LON_EX = 147.9
ke, lat, lon, acc, t = A.load()
P = A.BASE
fld = np.asarray(ke[0], float)
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
PC = ccrs.PlateCarree()


def g(name):
    r = G.loc[name]; return dict(pct=r.km_dec, pct_lo=r.lo, pct_hi=r.hi, p=r.p)


def sec(reg):
    q = S.loc[reg]; return dict(pct=q.m_lat_hp_slope * KM * 10, pct_lo=q.m_lat_hp_lo * KM * 10, pct_hi=q.m_lat_hp_hi * KM * 10, p=q.m_lat_hp_p)


fig = plt.figure(figsize=(W2, 4.7))
gs = fig.add_gridspec(2, 2, width_ratios=[0.8, 1.2], height_ratios=[1, 1], hspace=0.45, wspace=0.28)

# a
ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo())
polar(ax, ccrs, cfeature)
LON, LAT = np.meshgrid(lon, lat)
ax.pcolormesh(LON, LAT, np.where(ice, np.nan, fld * 1e4), transform=PC, cmap="Greys", vmin=0, vmax=1500, shading="auto",
              rasterized=True, zorder=0)
add_ice_mask(ax, ccrs, lon, lat, ice)
ax.scatter(lon[jj[inacc]], lhp[inacc], s=0.5, c=C["pos"], lw=0, transform=PC, zorder=5)
JEX = int(np.argmin(np.abs(lon - LON_EX)))
ex = tab[JEX, :cnt[JEX]]
ex = ex[A.in_acc(acc, lat, lon, ex[:, FI["lat_hp"]], np.full(len(ex), JEX))][0]
ax.plot(lon[JEX], ex[FI["lat_hp"]], "*", ms=7, color=C["obs"], mec="k", mew=0.5, transform=PC, zorder=7)
plabel(ax, "a", "Half-power points, January 2001")

# b
ax = fig.add_subplot(gs[0, 1])
prof = fld[:, JEX] * 1e4
clim = np.nanmean(np.asarray(ke[:, :, JEX], float), 0) * 1e4
s, e = int(ex[FI["i_start"]]), int(ex[FI["i_end"]])
ls_, ln_, pk_lat, pv, hp = ex[FI["lat_s"]], ex[FI["lat_n"]], ex[FI["peak_lat"]], ex[FI["peak_val"]] * 1e4, ex[FI["lat_hp"]]
ax.plot(lat, clim, color=C["grey"], lw=0.8, ls=(0, (3, 1.5)), label="2001–2025 mean")
ax.plot(lat, prof, color=C["ink"], lw=1.0, label="January 2001")
sel = (lat >= ls_) & (lat <= ln_)
ax.fill_between(lat[sel], 0, prof[sel], color=C["pos"], alpha=0.18, lw=0, label="Envelope")
ax.axhline(P["thr"] * 1e4, color=C["ke"], lw=0.7, ls=":", label="Threshold (200 cm$^2$ s$^{-2}$)")
ax.plot([lat[s] - 0.5, lat[e] + 0.5], [P["frac"] * pv] * 2, color=C["ke"], lw=0.7, ls="-.", label="25 % of peak")
ax.axvline(hp, color=C["pos"], lw=1.0, label="Half-power latitude")
for v_, lab in ((ls_, r"$\varphi_S$"), (ln_, r"$\varphi_N$")):
    ax.axvline(v_, color=C["pos"], lw=0.5, ls="--")
    ax.text(v_, pv * 1.02, lab, ha="center", va="bottom", fontsize=6, color=C["pos"])
ax.set_xlim(hp - 5.5, hp + 11); ax.set_ylim(0, pv * 1.12)
ax.set_xlabel("Latitude (°)"); ax.set_ylabel("KE (cm$^2$ s$^{-2}$)")
ax.legend(loc="upper right", fontsize=5.6)
plabel(ax, "b", f"KE profile at {lon[JEX]:.1f}°E", dx=-30)

# c
ax = fig.add_subplot(gs[1, 0])
y = np.nanmean(lat_hp - np.nanmean(lat_hp, 0, keepdims=True), 1) + np.nanmean(lat_hp)
r = S.loc["ACC"]
yd = pd.Series(y); yd = (yd - yd.groupby(month).transform("mean") + yd.mean()).values
ax.plot(x, y, color=C["pos"], lw=0.35, alpha=0.4)
ax.plot(x, running(y), color=C["pos"], lw=1.0)
b0 = theil_line(x, yd, r.m_lat_hp_slope); xf = np.array([x[0], x[-1]])
ax.plot(xf, b0 + r.m_lat_hp_slope * xf, color=C["ink"], lw=0.7, ls="--")
ts_style(ax); ax.set_ylabel("Latitude (°)")
plabel(ax, "c", "Circumpolar-mean half-power latitude", dx=-30)

# d
ax = fig.add_subplot(gs[1, 1])
TS = "This study, 2001–2025"
rows = [("Circumpolar", {TS: g("Circumpolar")}), ("Atlantic sector", {TS: sec("Atlantic")}), ("Indian sector", {TS: sec("Indian")}),
        ("Pacific sector", {TS: sec("Pacific")}), ("90–170°E", {TS: g("90-170E (Chambers 2018)")}),
        ("75–150°E", {TS: g("75-150E SE Indian (Kim and Orsi 2014)")}), ("160–90°W", {TS: g("160-90W (Chambers 2018)")}),
        ("150–80°W (PAR meander)", {TS: g("150-80W PAR meander (Liu et al. 2026)")})]
ys = forest(ax, rows, {TS: C["pos"]}, xlim=(-45, 50), xlabel="Trend (km per decade, positive northward)")
ax.axhline(-3.5, color=C["light"], lw=0.5)
LIT = [("Sokolov and Rintoul (2009b), 1992–2007", 0, -40, "v", "#7B3F00"), ("Gille (2014), ~1993–2012", 0, -30, "^", "#B07D2B"),
       ("Chambers (2018), 1993–2016", 4, -29, "P", "#444444"), (None, 6, 46, "P", "#444444"),
       ("Liu et al. (2026), 1993–2023", 7, 3.3, "X", "#56B4E9")]
for lab, row, val, mk, col in LIT:
    ax.plot(val, ys[row] - 0.3, mk, ms=4.2, color=col, mec="white", mew=0.4, zorder=6)
h = [Line2D([], [], ls="-", lw=0.8, marker="o", color=C["pos"], mfc="white", ms=3.6, label="This study, 2001–2025 (open: p ≥ 0.05)")]
h += [Line2D([], [], ls="none", marker=mk, color=col, mec="white", mew=0.4, ms=4.2, label=lab) for lab, _, _, mk, col in LIT if lab]
ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.35, -0.25), ncol=2, fontsize=5.6)
plabel(ax, "d", "Half-power latitude trends", dx=-80)
save(fig, FIGD, "Figure2")
