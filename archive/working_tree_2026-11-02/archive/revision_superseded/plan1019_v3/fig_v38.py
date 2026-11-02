"""Main-text figures for manuscript v3.3, Product 1 (all-satellite) unless stated; style_nc house style.
Redesigned toward the physical displays of the reference studies (EKE-trend maps with fronts, annual anomaly
time series of forcing and response as in Meredith and Hogg 2006 / Hogg et al. 2015, trend profiles across the
passage and across streamlines as in Gutierrez-Villanueva et al. 2023, front-position maps and longitude profiles
with topography as in Chambers 2018 / Sokolov and Rintoul 2009). No regions or values of earlier studies are drawn.
Figure 1  Drake Passage: where the flow changed (map, forcing/response series, transect, streamlines).
Figure 2  Circumpolar: EKE trend map, forcing and response series, sector partition, tracked eddies.
Figure 3  Pathway and standing meanders: half-power positions (map, circumpolar series, longitude profile)
          and along-stream EKE, jet speed, sinuosity, with the depth under the path.
Trend lines in the time-series panels are the reported trends (Theil-Sen on deseasonalised monthly series),
drawn through the anomaly mean: solid if p < 0.05 (p rounded to three decimals), dashed otherwise.
Usage: python fig_v33.py <revision dir> <drake out dir> <plan1019 out dir> <figure dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, DO, PO, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import retest_override, add_ice_mask, add_acc_contours, add_sector_boundaries
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import LogNorm, TwoSlopeNorm
plt = S.plt; C = S.C
OUT = REV / "outputs"; ROOT_OUT = REV.parent / "outputs"
PC = ccrs.PlateCarree(); KM = 111.195
FEAT = [(75, -47.5, "Kerguelen\nPlateau"), (122, -46.0, "Southeast\nIndian Ridge"), (166, -49.0, "Campbell\nPlateau"),
        (-140, -60.5, "Pacific–Antarctic\nRidge"), (-66, -60.0, "Drake\nPassage")]


def sigp(p):
    return np.isfinite(p) and round(float(p), 3) < 0.05


def ann(df, col, mean=None, y0=2001, y1=2025):
    s = df[(df.year >= y0) & (df.year <= y1)].groupby("year")[col].mean()
    mu = df[col].mean() if mean is None else mean
    return 100 * (s / mu - 1)


def trend_line(ax, pct, p, col, x0=2001, x1=2025.99, xc=None, lw=0.9):
    xc = (x0 + x1) / 2 if xc is None else xc
    xx = np.array([x0, x1]); yy = pct / 10 * (xx - xc)
    ax.plot(xx, yy, color=col, lw=lw, ls="-" if sigp(p) else (0, (3, 1.6)), zorder=3)


def series(ax, s, col, lab, ms=1.8, lw=0.6, mk="o", ls="-"):
    ax.plot(s.index + 0.5, s.values, ls, marker=mk, color=col, lw=lw, ms=ms, label=lab, alpha=0.9, zorder=4)


def bar_trends(ax, groups, series_def, width=0.17, ylim=None, ylabel=None):
    """groups: [(label, {key: dict(pct, lo, hi, p)})]; series_def: [(key, label, colour)]. Filled bar if p < 0.05,
    white with coloured edge otherwise; whiskers are 95 % intervals."""
    for gi, (gl, d) in enumerate(groups):
        keys = [k for k, _, _ in series_def if k in d]
        for k_i, k in enumerate(keys):
            col = dict((a, c) for a, _, c in series_def)[k]; r = d[k]
            x = gi + (k_i - (len(keys) - 1) / 2) * width
            ax.bar(x, r["pct"], width * 0.9, color=col if sigp(r["p"]) else "white", edgecolor=col, lw=0.8, zorder=3)
            ax.plot([x, x], [r["lo"], r["hi"]], color="#333333", lw=0.6, zorder=4)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g for g, _ in groups]); ax.tick_params(axis="x", length=0)
    S.zero(ax)
    if ylim: ax.set_ylim(*ylim)
    if ylabel: ax.set_ylabel(ylabel)
    h = [Patch(facecolor=c, edgecolor=c, label=l) for _, l, c in series_def]
    h.append(Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05"))
    return h


def lonfmt(ticks):
    return [("180°" if abs(t) == 180 else "0°" if t == 0 else f"{abs(t)}°{'W' if t < 0 else 'E'}") for t in ticks]


bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation
clim = xr.open_dataset(OUT / "clim_allsat.nc")
ice_ds = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc"); ice_g = ice_ds.ice_mask.values > 0.5
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
RT = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")


def rt(key):
    r = RT.loc[key]; f = 1000 / (r.slope * 1000 / r.pct_per_decade)
    return dict(pct=r.pct_per_decade, lo=r.ci_lo * f, hi=r.ci_hi * f, p=r.p_new)


# ================================================================== Figure 1: Drake Passage
tr1 = pd.read_csv(DO / "dk1" / "dk1_trends.csv"); tr1 = tr1[tr1.period == "2001-2025"].set_index("name")
ser1 = pd.read_csv(DO / "dk1" / "dk1_series.csv")
d1s = pd.read_csv(PO / "d1" / "d1b_series.csv"); d1t = pd.read_csv(PO / "d1" / "d1b_trends.csv").set_index(["series", "period"])
mf = sorted((REV / "data" / "monthly_allsat").glob("*/*.nc")); mf = [f for f in mf if 2001 <= int(f.parent.name) <= 2025]
box = dict(longitude=slice(-86, -44), latitude=slice(-66, -50))
with xr.open_dataset(mf[0]) as d0:
    sub = d0.sel(**box); dlat, dlon = sub.latitude.values, sub.longitude.values
us, vs, vis, yrs = [], [], [], []; icemax = np.zeros((len(dlat), len(dlon)))
for f in mf:
    with xr.open_dataset(f) as d:
        s = d.sel(**box)
        us.append(s.u_mean.values.astype(np.float32)); vs.append(s.v_mean.values.astype(np.float32))
        vis.append(s.var_intra.values.astype(np.float32)); icemax = np.maximum(icemax, np.nan_to_num(s.ice_frac.values))
        yrs.append(int(f.parent.name))
us, vs, vis, yrs = np.array(us), np.array(vs), np.array(vis), np.array(yrs)
ub, vb = np.nanmean(us, 0), np.nanmean(vs, 0)
ekem = (0.5 * ((us - ub) ** 2 + (vs - vb) ** 2) + 0.5 * vis) * 1e4          # monthly EKE, cm2 s-2
del us, vs, vis
YY = np.arange(2001, 2026); ea = np.array([np.nanmean(ekem[yrs == y], 0) for y in YY]); del ekem
xa = YY - YY.mean(); slope_map = np.tensordot(xa, ea - ea.mean(0), axes=(0, 0)) / (xa ** 2).sum() * 10   # per decade
slope_map[icemax > 0.10] = np.nan
cs = 2; ny, nx = (len(dlat) // cs) * cs, (len(dlon) // cs) * cs
sm = np.nanmean(slope_map[:ny, :nx].reshape(ny // cs, cs, nx // cs, cs), axis=(1, 3))
slat = dlat[:ny].reshape(-1, cs).mean(1); slon = dlon[:nx].reshape(-1, cs).mean(1)
cl = clim.sel(**box); adt_d = cl.adt.values
trk = pd.read_csv(PO / "d1" / "d1b_track.csv")
edges = pd.read_csv(DO / "dk5" / "dk5_stream_edges.csv").set_index("psi").adt_edge

fig = plt.figure(figsize=(S.W2, 1.0 * S.W2))
gsF = fig.add_gridspec(2, 1, height_ratios=[1.2, 1], hspace=0.3)
gtop = gsF[0].subgridspec(1, 2, width_ratios=[1, 1.3], wspace=0.22)
gbot = gsF[1].subgridspec(1, 2, width_ratios=[1.2, 1], wspace=0.32)
ax = fig.add_subplot(gtop[0])
im = ax.pcolormesh(slon, slat, sm, cmap="RdBu_r", norm=TwoSlopeNorm(0, -80, 80), shading="auto", rasterized=True)
bs = bat.sel(**box)
ax.contour(bs.longitude, bs.latitude, -bs.values, levels=[3000], colors="#8A8A8A", linewidths=0.35)
ax.contourf(bs.longitude, bs.latitude, bs.values, levels=[0, 1e5], colors=["#3A3A3A"])
# mean positions of the SAF and PF cores: central streamlines of the GV23 zones, mapped to ADT
for psi, lab, ls in ((-0.5, "SAF", "-"), (-1.125, "PF", (0, (4, 1.5)))):
    lev = float(np.interp(psi, edges.index[::-1], edges.values[::-1]))
    ax.contour(cl.longitude, cl.latitude, adt_d, levels=[lev], colors="#1A1A1A", linewidths=1.0, linestyles="solid", zorder=5)
ax.add_patch(Rectangle((-80, -63), 28, 9, fill=False, ec="#444444", lw=0.7, ls=(0, (1, 1)), zorder=6))
ax.text(-56.5, -62.3, "Drake box", ha="center", va="center", fontsize=5.5, zorder=7,
        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))
ax.plot(trk.lon[trk.perm], trk.lat[trk.perm], color=C["obs"], lw=1.6, zorder=8, solid_capstyle="butt")
ax.plot(trk.lon[~trk.perm], trk.lat[~trk.perm], color=C["obs"], lw=1.2, ls=(0, (2, 1.2)), zorder=8)
tk6 = pd.read_csv(DO / "dk6" / "dk6_tracks_allsat_ice0.csv"); tk6 = tk6[tk6["cov"] >= 0.5]
ax.plot(tk6.lon, tk6.hpc, ".", color=C["pos"], ms=1.6, zorder=7)
ax.text(-85, -55.4, "SAF", fontsize=5.5, fontweight="bold", zorder=8); ax.text(-85, -60.4, "PF", fontsize=5.5, fontweight="bold", zorder=8)
ax.set_xlim(-86, -44); ax.set_ylim(-66, -50); ax.set_aspect(1 / np.cos(np.radians(58)))
ax.set_xticks(np.arange(-80, -44, 10)); ax.set_xticklabels([f"{-v}°W" for v in np.arange(-80, -44, 10)])
ax.set_yticks([-65, -60, -55, -50]); ax.set_yticklabels([f"{-v}°S" for v in (-65, -60, -55, -50)])
ax.tick_params(labelsize=5.5)
S.hcbar(fig, im, ax, "EKE trend, 2001–2025 (cm$^2$ s$^{-2}$ per decade)", shrink=0.7, pad=0.12, extend="both")
S.plabel(ax, "a", "Where the eddy energy changed")

# b, c: monthly deseasonalised anomalies (thin) with annual means (markers) and, to the right, the trends with 95 % intervals
def mon_anom(df, col, pct=True):
    y = df[col].astype(float); m = df["month"] if "month" in df else ((df.x - df.x.astype(int)) * 12 + 0.5).astype(int)
    a = y - y.groupby(m.values).transform("mean")
    if pct: a = 100 * a / y.mean()
    yr = df["year"] if "year" in df else df.x.astype(int)
    return pd.Series(a.values, index=(yr + (m - 0.5) / 12).values), a.groupby(yr.values).mean()
DK6 = DO / "dk6"
s6 = pd.read_csv(DK6 / "dk6_series_allsat_ice0.csv"); t6 = pd.read_csv(DK6 / "dk6_trends_allsat_ice0.csv")
t6 = t6[t6.coverage_min == 0.5].set_index("series")
def T1(k): r = tr1.loc[k]; return (r.pct, r.pct_lo, r.pct_hi, r.p)
rd = d1t.loc[("Tperm_gv", "2001-2025 monthly")]
H1 = "T1 half-power latitude (km, +N)"; H2 = "T2 latitude of the climatological psi_hp streamline (km, +N)"
ROWS = [  # (entries, ylabel, row title, trend unit)
    ([(ser1, "allsat_box_KE", True, C["ke"], "KE", "o", T1("allsat_box_KE")), (ser1, "allsat_box_EKE", True, C["eke"], "EKE", "o", T1("allsat_box_EKE"))],
     "%", "Kinetic energy", "% dec$^{-1}$"),
    ([(s6, f"0.5|{H1}", False, C["pos"], "energetic band (half-power latitude, 70–55°W)", "o", tuple(t6.loc[H1, ["trend_dec", "lo", "hi", "p"]])),
      (s6, f"0.5|{H2}", False, C["ink"], "sea-surface-height streamline through its mean position", "D", tuple(t6.loc[H2, ["trend_dec", "lo", "hi", "p"]]))],
     "km, + north", "Position", "km dec$^{-1}$")]
gts0 = gtop[1].subgridspec(2, 1, height_ratios=[1, 1], hspace=0.42)
gtsA = gts0[0].subgridspec(1, 2, width_ratios=[1, 0.13], wspace=0.05)
gtsB = gts0[1].subgridspec(1, 2, width_ratios=[1, 0.13], wspace=0.05)
cells = [(gtsA[0, 0], gtsA[0, 1]), (gtsB[0, 0], gtsB[0, 1])]
for r_, (entries, ylab, ttl, tunit) in enumerate(ROWS):
    ax = fig.add_subplot(cells[r_][0]); axt = fig.add_subplot(cells[r_][1])
    lim = 0
    for k, (df, col, pct, colr, lab, mk, trd) in enumerate(entries):
        ms_, an = mon_anom(df, col, pct)
        ax.plot(ms_.index, ms_.values, color=colr, lw=0.35, alpha=0.45, zorder=2)
        ax.plot(an.index + 0.5, an.values, "-", marker=mk, color=colr, lw=0.9, ms=2.2, mec="white", mew=0.3, zorder=4)
        lim = max(lim, np.nanpercentile(np.abs(ms_.values), 96))
        ax.text(0.01 + 0.0 * k, 0.97 - 0.15 * k, lab, transform=ax.transAxes, fontsize=5.0, color=colr, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.05", fc="white", ec="none", alpha=0.7), zorder=6)
        v, lo_, hi_, pp = trd; sg = sigp(pp); x0 = k - (len(entries) - 1) / 2
        axt.plot([x0, x0], [lo_, hi_], color=colr, lw=0.9, solid_capstyle="butt")
        axt.plot(x0, v, mk, color=colr, mfc=colr if sg else "white", mec=colr, mew=0.8, ms=3.2, zorder=3)
    lim = np.ceil(lim / 10) * 10 if lim > 10 else np.ceil(lim)
    S.zero(ax); S.ts_style(ax); ax.set_ylim(-lim, lim); ax.set_ylabel(ylab, labelpad=1)
    axt.axhline(0, color="#9A9A9A", lw=0.5, zorder=0); axt.set_xlim(-1, 1); axt.set_xticks([])
    axt.yaxis.tick_right(); axt.yaxis.set_label_position("right"); axt.spines["left"].set_visible(False); axt.spines["right"].set_visible(True)
    axt.spines["bottom"].set_visible(False); axt.tick_params(axis="y", labelsize=5.0, length=2)
    axt.set_ylabel(tunit, fontsize=5.0, labelpad=1)
    if r_ == 0:
        axt.set_title("trend", fontsize=5.5, pad=2)
        S.plabel(ax, "b", "Eddy and mean-flow energy in the Drake box", dy=10)
    if r_ == 1:
        S.plabel(ax, "c", "Position of the energetic band in the passage", dy=4)
# c: trends along the LM Gould transect, passage coordinates (as GV23)
K5 = pd.read_csv(DO / "dk5" / "dk5_trends.csv")
def prof(coord, src, var, per):
    return K5[(K5.coord == coord) & (K5.src == src) & (K5["var"] == var) & (K5.period == per)].sort_values("bin")
REC = ((("insitu", "2005-2019"), C["obs"], "in situ ADCP, Oct 2005–Apr 2019", 3.0),
       (("allsat", "2005-2019"), C["eke"], "altimetry, same period", 2.4),
       (("allsat", "2001-2025"), C["ink"], "altimetry, 2001–2025", 2.4))
# d: EKE trends along the LM Gould transect (relative, % of the bin mean per decade): altimetry (surface) and in situ (0-780 m, GV23)
ax = fig.add_subplot(gbot[0]); ax2 = ax.twinx()
me = prof("passage", "allsat", "eke", "2001-2025")
ax2.fill_between(me.bin, 0, me["mean"], color="#EDEDED", lw=0, zorder=0); ax2.set_ylim(0, 700)
ax2.set_ylabel("Mean surface EKE (cm$^2$ s$^{-2}$)", color="#8A8A8A", labelpad=1); ax2.spines["right"].set_visible(True); ax2.tick_params(colors="#8A8A8A")
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.axvspan(408, 850, facecolor="none", hatch="....", edgecolor="#C9DDF0", lw=0, zorder=0)
q = prof("passage", "allsat", "eke", "2001-2025")
ax.fill_between(q.bin, 1000 * q.ci_lo / q["mean"], 1000 * q.ci_hi / q["mean"], color=C["eke"], alpha=0.15, lw=0, zorder=1)
INS = pd.read_csv(DO / "dk7" / "dk7_insitu_eke_passage.csv")
curves = [(INS.bin, INS.pct, INS.p, C["obs"], "in situ ADCP, 0–780 m, Oct 2005–Apr 2019", "-", 3.0)]
for per, ls, lab in (("2005-2019", (0, (3, 1.5)), "altimetry, Oct 2005–Apr 2019"), ("2001-2025", "-", "altimetry, 2001–2025 (shading: 95 % interval)")):
    q = prof("passage", "allsat", "eke", per); curves.append((q.bin, 1000 * q.slope / q["mean"], q.p, C["eke"], lab, ls, 2.4))
for x, y, pp, col, lab, ls, ms in curves:
    x, y, pp = np.asarray(x), np.asarray(y), np.asarray(pp); sg = np.round(pp, 3) < 0.05
    ax.plot(x, y, color=col, lw=0.9, ls=ls, label=lab, zorder=3)
    ax.plot(x[sg], y[sg], "o", color=col, ms=ms, zorder=4); ax.plot(x[~sg], y[~sg], "o", mfc="white", color=col, ms=ms, mew=0.6, zorder=4)
S.zero(ax); ax.set_xlim(0, 850); ax.set_ylim(-60, 60); ax.set_xlabel("Distance across the passage from the north (km)")
ax.set_ylabel("EKE trend (% of mean per decade)")
for x_, lab in ((40, "SAF"), (400, "PF"), (790, "SACCF")):
    ax.text(x_, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.5, color="#555555", fontweight="bold")
ax.legend(loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "d", "EKE trend along the LM Gould transect", dx=-24, dy=10)

# e: trends between pairs of SSH streamlines in the frontal zones of GV23 (streamwise coordinates)
ZA = pd.read_csv(DO / "dk5" / "dk5_zones.csv"); ZA = ZA[ZA.coord == "stream"]
ZONES = [("SAF (-0.40 to -0.60 m)", "SAF"), ("SAF-PF zone (-0.60 to -1.00 m)", "SAF–PF"), ("PF (-1.00 to -1.25 m)", "PF")]
gse = gbot[1].subgridspec(1, 2, wspace=0.35)
for k, (var, ttl, ylim) in enumerate((("transport", "Transport", (-6, 6)), ("eke", "EKE", (-30, 45)))):
    ax = fig.add_subplot(gse[k])
    for zi, (zk, zl) in enumerate(ZONES):
        for ri, ((src, per), col, lab, ms) in enumerate(REC[2:]):
            z = ZA[(ZA.zone == zk) & (ZA.src == src) & (ZA.period == per) & (ZA["var"].str.lower().str.startswith(var))]
            if z.empty: continue
            z = z.iloc[0]; x = zi
            ax.bar(x, z.pct, 0.5, color=col if sigp(z.p) else "white", edgecolor=col, lw=0.8, zorder=3)
            ax.plot([x, x], [z.pct_lo, z.pct_hi], color="#333333", lw=0.6, zorder=4)
    S.zero(ax); ax.set_xticks(range(3)); ax.set_xticklabels([z[1] for z in ZONES]); ax.tick_params(axis="x", length=0)
    ax.set_ylim(*ylim); ax.set_ylabel("% per decade" if k == 0 else "")
    ax.text(0.5, 1.01, ttl, transform=ax.transAxes, ha="center", va="bottom", fontsize=6)
    if k == 0: S.plabel(ax, "e", "Trends between streamlines, 2001–2025", dx=-24, dy=12)
S.save(fig, FIGD, "Figure1"); plt.close(fig)

# ================================================================== Figure 2: circumpolar
M = xr.open_dataset(OUT / "a1_maps_bc.nc")
A5 = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv"); A5t = pd.read_csv(OUT / "a5" / "a5_series_twosat.csv")
A8 = pd.read_csv(OUT / "a8" / "a8_series_allsat.csv"); A8t = pd.read_csv(OUT / "a8" / "a8_trends_allsat.csv").set_index("series")
W = pd.read_csv(OUT / "a6" / "a6_wind_stress_series.csv")
E2 = pd.read_csv(OUT / "eddies" / "e2_series.csv"); TE = pd.read_csv(OUT / "eddies" / "e2_e3_trends.csv").set_index("series")
for df in (A5, A5t):
    df["year"] = 2001 + df.i // 12
lat, lon = ice_ds.latitude.values, ice_ds.longitude.values

fig = plt.figure(figsize=(S.W2, 0.88 * S.W2))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.1], height_ratios=[1.15, 0.85], wspace=0.3, hspace=0.3)
ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo())
S.polar(ax, ccrs, cfeature)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
im = ax.pcolormesh(L1, A1, 10 * M["EKE_slope"].values, transform=PC, cmap=tc, vmin=-60, vmax=60, shading="nearest", rasterized=True, zorder=0)
bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
add_ice_mask(ax, ccrs, lon, lat, ice_g); add_sector_boundaries(ax, ccrs)
jj, ii = np.nonzero(M["EKE_fdr"].values == 1)
ax.scatter(M.lon.values[ii], M.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
for lo_, la_, s_ in FEAT:
    ax.text(lo_, la_, s_, transform=PC, fontsize=4.6, ha="center", va="center", color="#222222", zorder=9, linespacing=0.9,
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
for lo_, s_ in ((-110, "Pacific"), (65, "Indian"), (-15, "Atlantic")):
    ax.text(lo_, -38.5, s_, transform=PC, fontsize=5.2, ha="center", va="center", color="#555555", style="italic", zorder=9)
S.plabel(ax, "a", "Where the eddy energy changed, 2001–2025")
S.hcbar(fig, im, ax, "EKE trend (cm$^2$ s$^{-2}$ per decade)", shrink=0.7, extend="both")

gsb = gs[0, 1].subgridspec(2, 1, hspace=0.2)
ax = fig.add_subplot(gsb[0])
W["year"] = W["year"].astype(int)
series(ax, ann(A5, "ACC_SPD"), C["tr"], "along-stream speed"); r = rt("allsat|ACC_SPD|2001-2025"); trend_line(ax, r["pct"], r["p"], C["tr"])
series(ax, ann(A8, "ACC_T"), "#56B4E9", "surface transport proxy", mk="s", ms=1.6); r = A8t.loc["ACC_T"]; trend_line(ax, r.pct_per_decade, r.p, "#56B4E9")
S.zero(ax); S.ts_style(ax); ax.set_ylim(-6, 6); ax.set_xticklabels([]); ax.set_ylabel("Anomaly (%)")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "b", "Mean flow, ACC band", dy=20)
ax = fig.add_subplot(gsb[1])
series(ax, ann(A5, "ACC_EKE"), C["eke"], "EKE"); r = rt("allsat|ACC_EKE|2001-2025"); trend_line(ax, r["pct"], r["p"], C["eke"])
series(ax, ann(A5, "ACC_KE"), C["ke"], "KE"); r = rt("allsat|ACC_KE|2001-2025"); trend_line(ax, r["pct"], r["p"], C["ke"])
series(ax, ann(A5t, "ACC_EKE"), C["two"], "EKE, two-satellite product", mk="^", ms=1.8); r = rt("twosat|ACC_EKE|2001-2025"); trend_line(ax, r["pct"], r["p"], C["two"])
S.zero(ax); S.ts_style(ax); ax.set_ylim(-20, 20); ax.set_ylabel("Anomaly (%)")
ax.legend(loc="upper left", ncol=3, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "c", "Kinetic energy, ACC band", dy=2)

ax = fig.add_subplot(gs[1, 0])
groups = []
for reg, gl in (("ACC", "ACC band"), ("Pacific", "Pacific"), ("Indian", "Indian"), ("Atlantic", "Atlantic")):
    t8 = A8t.loc[f"{reg}_T"]
    d = {"EKE": rt(f"allsat|{reg}_EKE|2001-2025"), "KE": rt(f"allsat|{reg}_KE|2001-2025"),
         "SPD": rt(f"allsat|{reg}_SPD|2001-2025"), "T": dict(pct=t8.pct_per_decade, lo=t8.pct_lo, hi=t8.pct_hi, p=t8.p)}
    groups.append((gl, d))
h = bar_trends(ax, groups, [("EKE", "EKE", C["eke"]), ("KE", "KE", C["ke"]), ("SPD", "along-stream speed", C["tr"]),
                            ("T", "surface transport proxy", "#56B4E9")], ylim=(-3, 17), ylabel="Trend 2001–2025 (% per decade)")
ax.legend(handles=h, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "d", "Eddies and mean flow by sector", dy=20)

ax = fig.add_subplot(gs[1, 1])
E2["year"] = E2.year.astype(int)
for key, col, lab, mk in (("ACC_amp", C["te"], "amplitude", "D"), ("ACC_spd", "#CC79A7", "rotational speed", "v"),
                          ("ACC_present_per_day", C["ink"], "number present", "o"), ("ACC_rad", C["two"], "speed radius", "s")):
    s = ann(E2[E2.year <= 2022], key, mean=E2[E2.year <= 2022][key].mean(), y1=2022)
    series(ax, s, col, lab, mk=mk, ms=1.8)
    r = TE.loc[key]; trend_line(ax, r.pct_per_decade, r.p, col, x0=2001, x1=2023.67)
S.zero(ax); S.ts_style(ax, (2000.8, 2024)); ax.set_ylim(-12, 12); ax.set_ylabel("Annual anomaly (%)")
ax.legend(loc="upper left", ncol=2, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "e", "Tracked eddies, ACC band: stronger, not more numerous", dy=20)
S.save(fig, FIGD, "Figure2"); plt.close(fig)

# ================================================================== Figure 3: pathway and standing meanders
TB = np.load(OUT / "a2" / "tracks_baseline.npz")
keep = TB["keep"]; jlon = clim.longitude.values[TB["j"]]; latc = TB["lat_hp_clim"]
hp = TB["lat_hp"][:, keep].astype(float); hp_series = np.nanmean(hp - np.nanmean(hp, 0, keepdims=True), axis=1)
t_m = 2001 + (np.arange(hp.shape[0]) + 0.5) / 12
LB = RT[RT.index.str.startswith("lonbin|")].copy(); LB["lon0"] = LB.index.str.split("|").str[1].astype(float); LB["c"] = LB.lon0 + 2.5
LB["tr"] = LB.slope * KM * 10; LB["lo"] = LB.ci_lo * KM * 10; LB["hi"] = LB.ci_hi * KM * 10
CL = pd.read_csv(OUT / "b3" / "b3_climatology.csv")
T = {v: pd.read_csv(OUT / "b3" / f"b3_trends_{v}.csv") for v in ("eke", "speed", "sinu", "lat")}
SHALLOW = CL.lon[CL.depth_min < 3000].values
_dk2 = pd.read_csv(DO / "dk2" / "dk2_trends.csv"); _r = _dk2[(_dk2["product"] == "allsat") & (_dk2.period == "2001-2025")].iloc[0]
rc = RT.loc["baseline|ACC|lat_hp"]

fig = plt.figure(figsize=(S.W2, 1.22 * S.W2))
gs = fig.add_gridspec(6, 1, height_ratios=[2.25, 1, 1, 1, 1, 0.55], hspace=0.5)
gtop = gs[0].subgridspec(1, 2, width_ratios=[1.25, 1], wspace=0.22)
ax = fig.add_subplot(gtop[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
spd = np.hypot(clim.u.values, clim.v.values) * 100
c4 = clim.coarsen(latitude=2, longitude=4, boundary="trim").mean(); sp4 = np.hypot(c4.u, c4.v) * 100
im = ax.pcolormesh(c4.longitude, c4.latitude, sp4, transform=PC, cmap="YlGnBu", vmin=0, vmax=40, shading="nearest", rasterized=True, zorder=0)
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
# climatological half-power positions: grey = in the ACC band but not tracked in >= 90 % of months; coloured = tracked,
# coloured by the 5-degree bin trend
bin_tr = dict(zip(LB.lon0, LB.tr))
cc = np.array([bin_tr.get(np.floor(l / 5) * 5, np.nan) for l in jlon])
nk = TB["acc"] & ~keep
ax.scatter(jlon[nk][::2], latc[nk][::2], s=1.0, c="#6F6F6F", lw=0, transform=PC, zorder=3)
ax.scatter(jlon[keep], latc[keep], s=3.0, c="#7B2D8E", lw=0, transform=PC, zorder=4)
ax.plot(np.r_[CL.lon, CL.lon.iloc[0] + 360], np.r_[CL.lat_mean, CL.lat_mean.iloc[0]], color="#D7191C", lw=0.8, transform=PC, zorder=5)
ax.plot(trk.lon, trk.lat, color=C["obs"], lw=1.2, transform=PC, zorder=6)
for lo_, la_, s_ in FEAT:
    ax.text(lo_, la_ + (3 if s_.startswith("Drake") else 0), s_, transform=PC, fontsize=4.6, ha="center", va="center", zorder=9, linespacing=0.9,
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
S.hcbar(fig, im, ax, "Mean surface speed, 2001–2025 (cm s$^{-1}$)", shrink=0.7, extend="max")
S.plabel(ax, "a", "Energetic pathway and jet-core path")
ax = fig.add_subplot(gtop[1])
ax.plot(t_m, hp_series * KM, color="#BDBDBD", lw=0.5)
rm = pd.Series(hp_series * KM).rolling(12, center=True).mean()
ax.plot(t_m, rm, color=C["pos"], lw=1.1, label="12-month running mean")
tl = rc.slope * KM * (t_m - t_m.mean()); ax.plot(t_m, tl, color="k", lw=0.9, ls=(0, (3, 1.6)))
ax.fill_between(t_m, rc.ci_lo * KM * (t_m - t_m.mean()), rc.ci_hi * KM * (t_m - t_m.mean()), color="k", alpha=0.08, lw=0)
ax.fill_between(t_m, -12.5, 12.5, color=C["pos"], alpha=0.0)
S.zero(ax); S.ts_style(ax); ax.set_ylim(-60, 60); ax.set_ylabel("Half-power latitude anomaly (km, + north)")
ax.text(0.02, 0.04, f"trend {rc.slope * KM * 10:+.1f} km per decade ({rc.ci_lo * KM * 10:+.1f} to {rc.ci_hi * KM * 10:+.1f})\n"
        f"a 25 km shift over 2001–2025 would be detected", transform=ax.transAxes, fontsize=5.3, va="bottom")
ax.legend(loc="upper right", fontsize=5.0)
S.plabel(ax, "b", "Circumpolar mean over the tracked longitudes", dx=-26)

def bands(ax):
    for l in SHALLOW: ax.axvspan(l - 2.5, l + 2.5, color="#EDEDED", lw=0, zorder=0)
    ax.axvspan(-60, -50, facecolor="none", edgecolor="#BDBDBD", hatch="////", lw=0, zorder=0)
def pts(ax, x, y, lo, hi, p, col, ms=2.6, mk="o", fdr=None, dx=0.0):
    sg = np.round(p, 3) < 0.05
    ax.vlines(x + dx, lo, hi, color=col, lw=0.6)
    ax.plot((x + dx)[sg], y[sg], mk, color=col, ms=ms, zorder=4)
    ax.plot((x + dx)[~sg], y[~sg], mk, mfc="white", color=col, ms=ms, mew=0.6, zorder=4)
    if fdr is not None: ax.plot((x + dx)[fdr], y[fdr], "o", mfc="none", mec="#111111", ms=5.5, mew=0.7, zorder=5)
TICK = np.arange(-180, 181, 30)
axs = [fig.add_subplot(gs[k]) for k in range(1, 6)]
ax = axs[0]; bands(ax); S.zero(ax)
pts(ax, LB.c.values, LB.tr.values, LB.lo.values, LB.hi.values, LB.p_new.values, C["pos"], dx=-0.8)
d = T["lat"].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]
pts(ax, d.lon.values, d.tr.values, d.lo.values, d.hi.values, d.p.values, "#555555", mk="D", ms=2.2, dx=0.8,
    fdr=d.fdr10.astype(bool).values if "fdr10" in d else None)
ax.errorbar([-64], [_r.trend_km_dec], yerr=[[_r.trend_km_dec - _r.ci_lo], [_r.ci_hi - _r.trend_km_dec]], fmt="*", color=C["obs"],
            mfc="white" if not sigp(_r.p) else C["obs"], ms=5.5, lw=0.7, zorder=6)
ax.set_ylim(-60, 60); ax.set_ylabel("km per decade")
ax.legend(handles=[Line2D([], [], ls="-", marker="o", color=C["pos"], ms=2.6, lw=0.6, label="half-power latitude (5° bins)"),
                   Line2D([], [], ls="-", marker="D", color="#555555", ms=2.2, lw=0.6, label="jet-core path latitude"),
                   Line2D([], [], ls="-", marker="*", color=C["obs"], mfc="white", ms=5, lw=0.6, label="half-power, LM Gould section")],
          loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "c", "Meridional displacement (positive northward)", dx=-28, dy=12)
spec = [("eke", "d", "EKE within 150 km of the jet-core path", C["eke"]), ("speed", "e", "Jet speed along the path", C["tr"]),
        ("sinu", "f", "Path sinuosity (meander amplitude)", C["sin"])]
for ax, (v, l, ttl, col) in zip(axs[1:4], spec):
    bands(ax); S.zero(ax)
    d = T[v].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]
    pts(ax, d.lon.values, d.tr.values, d.lo.values, d.hi.values, d.p.values, col, fdr=d.fdr10.astype(bool).values if "fdr10" in d else None)
    lo_, hi_ = np.nanmin(d.lo), np.nanmax(d.hi); pad = 0.06 * (hi_ - lo_); ax.set_ylim(lo_ - pad, hi_ + pad)
    ax.set_ylabel("% per decade"); S.plabel(ax, l, ttl, dx=-28)
for ax in axs[:4]:
    ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels([])
ax = axs[4]
dep = CL.sort_values("lon")
ax.fill_between(dep.lon, -dep.depth_min / 1000, -6.5, color="#C8C8C8", step="mid", lw=0, zorder=0)
ax.fill_between(dep.lon, -dep.depth / 1000, -6.5, color="#6A6A6A", step="mid", lw=0, zorder=1)
bands(ax); ax.set_ylim(-6.0, 0); ax.set_yticks([-5, -3, -1]); ax.set_yticklabels(["5", "3", "1"]); ax.set_ylabel("Depth (km)")
ax.set_xlim(-180, 180); ax.set_xticks(TICK); ax.set_xticklabels(lonfmt(TICK)); ax.set_xlabel("Longitude (flow from left to right)")
for x_, s_ in [(-142, "PAR"), (-65, "Drake"), (37, "SWIR"), (78, "Kerguelen"), (138, "SEIR"), (160, "Macquarie")]:
    ax.text(x_, -0.3, s_, ha="center", va="top", fontsize=5.0, color="#222222")
S.plabel(ax, "g", "Depth under the jet-core path (dark: bin mean; light: shallowest point)", dx=-28)
S.save(fig, FIGD, "Figure3"); plt.close(fig)
print("done")
