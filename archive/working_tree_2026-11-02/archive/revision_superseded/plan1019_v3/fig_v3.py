"""Main-text figures for manuscript v3 (Drake-first structure), Product 1 (all-satellite) only; style_nc house style.
Figure 1  Drake Passage setting and trends (box, upstream, downstream), 2001-2025.
Figure 2  Altimetry compared with the LM Gould ADCP record (2005-2019) and the extended record (2001-2025);
          panel e = transport trends (Fig. 3c of the manuscript while the Antarctic-slope index is pending).
Figure 4  Circumpolar EKE trends and comparison with earlier estimates; ACC/sector trends; tracked eddies.
Figure 5  Pathway (half-power latitude trends) and along-stream structure (EKE, jet speed, sinuosity, latitude).
(Figure 3, the Antarctic-slope index, is pending the southern extraction.)
Significance symbols: filled if p < 0.05 after rounding p to three decimals (p = 0.050 -> open, "at the threshold").
Usage: python fig_v3.py <revision dir> <drake out dir> <plan1019 out dir> <GV repo dir> <m1_trends.csv> <figure dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, scipy.io as sio
REV, DO, PO, GV, M1F, FIGD = (Path(a) for a in sys.argv[1:7])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import retest_override, running, add_ice_mask, add_acc_contours, add_sector_boundaries
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
plt = S.plt; C = S.C
OUT = REV / "outputs"; ROOT_OUT = REV.parent / "outputs"
PC = ccrs.PlateCarree()


def rp(d):
    """round p to 3 decimals so that p = 0.0498 (reported 0.050, at the threshold) plots as open."""
    d = dict(d); d["p"] = round(float(d["p"]), 3) if np.isfinite(d["p"]) else np.nan; return d


def pt(ax, r, yy, c, m="o", ms=3.4):
    r = rp(r)
    ax.plot([r["pct_lo"], r["pct_hi"]], [yy, yy], color=c, lw=0.8, solid_capstyle="butt", zorder=2)
    sig = np.isfinite(r["p"]) and r["p"] < 0.05
    ax.plot(r["pct"], yy, m, color=c, mfc=c if sig else "white", mec=c, mew=0.8, ms=ms, zorder=3)


def grouped_forest(ax, groups, xlim, xlabel, label_x=-0.01, head_size=6.0, tick_size=5.8):
    """groups: [(header, [(label, [(row dict, colour, marker), ...]), ...])]; several entries per label are dodged."""
    y = 0.0; yt, yl = [], []
    for head, items in groups:
        if head:
            ax.text(label_x, y + 0.75, head, transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=head_size, fontweight="bold")
        for lab, entries in items:
            n = len(entries)
            for k, (r, c, m) in enumerate(entries):
                pt(ax, r, y + (0.16 * ((n - 1) / 2 - k) if n > 1 else 0), c, m)
            yt.append(y); yl.append(lab); y -= 1
        y -= 0.6
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=tick_size); ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False); ax.set_ylim(y + 1.0, 1.3)
    ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0); ax.set_xlim(*xlim); ax.set_xlabel(xlabel)
    return yt


def sig_legend(ax, cols, markers=None, **kw):
    ax.legend(handles=S.sig_handles(cols, markers, ns_label="open: p ≥ 0.05"), **kw)


# ------------------------------------------------------------------ shared inputs
tr1 = pd.read_csv(DO / "dk1" / "dk1_trends.csv"); ser1 = pd.read_csv(DO / "dk1" / "dk1_series.csv")
def dk(name):
    r = tr1[(tr1.name == name) & (tr1.period == "2001-2025")].iloc[0]
    return dict(pct=r.pct, pct_lo=r.pct_lo, pct_hi=r.pct_hi, p=r.p)
bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation
clim = xr.open_dataset(OUT / "clim_allsat.nc")
ice_ds = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc"); ice_g = ice_ds.ice_mask.values > 0.5
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt

# ================================================================== Figure 1: Drake Passage
mf = sorted((REV / "data" / "monthly_allsat").glob("*/*.nc"))
mf = [f for f in mf if 2001 <= int(f.parent.name) <= 2025]
box = dict(longitude=slice(-106, -28), latitude=slice(-66, -49))
U = V = None; n = 0; VI = None; ICE = None; UU = []
with xr.open_dataset(mf[0]) as d0:
    sub = d0.sel(**box); dlat, dlon = sub.latitude.values, sub.longitude.values
um = np.zeros((len(dlat), len(dlon))); vm = np.zeros_like(um); vi = np.zeros_like(um); icemax = np.zeros_like(um)
us, vs = [], []
for f in mf:
    with xr.open_dataset(f) as d:
        s = d.sel(**box)
        u, v = s.u_mean.values, s.v_mean.values
        us.append(u.astype(np.float32)); vs.append(v.astype(np.float32))
        vi += np.nan_to_num(s.var_intra.values); icemax = np.maximum(icemax, np.nan_to_num(s.ice_frac.values))
us, vs = np.array(us), np.array(vs); ub, vb = np.nanmean(us, 0), np.nanmean(vs, 0)
eke_map = (0.5 * np.nanmean((us - ub) ** 2 + (vs - vb) ** 2, 0) + 0.5 * vi / len(mf)) * 1e4   # cm2 s-2
eke_map[icemax > 0.10] = np.nan
del us, vs
cl = clim.sel(**box); adt_d = cl.adt.values
trk = pd.read_csv(PO / "d1" / "d1b_track.csv")

fig = plt.figure(figsize=(S.W2, 1.0 * S.W2))
gs = fig.add_gridspec(3, 2, width_ratios=[1.25, 1], height_ratios=[1.35, 1.0, 0.95], wspace=0.95, hspace=0.75)
ax = fig.add_subplot(gs[0, 0])
from matplotlib.colors import LogNorm
im = ax.pcolormesh(dlon, dlat, eke_map, cmap="YlGnBu", norm=LogNorm(30, 2000), shading="auto", rasterized=True)
bs = bat.sel(**box)
ax.contour(bs.longitude, bs.latitude, -bs.values, levels=[3000], colors="#8A8A8A", linewidths=0.35)
ax.contourf(bs.longitude, bs.latitude, bs.values, levels=[0, 1e5], colors=["#3A3A3A"])
ax.contour(cl.longitude, cl.latitude, adt_d, levels=[-0.6, 0.2], colors=["#5f5f5f", "#2f2f2f"], linewidths=[0.6, 0.9])
ax.add_patch(Rectangle((-80, -63), 28, 9, fill=False, ec="k", lw=0.9, zorder=6))
acc_d = np.isfinite(adt_d) & (adt_d >= -0.6) & (adt_d <= 0.2)
for (l0, l1), lab in (((-100, -80), "upstream"), ((-55, -35), "downstream")):
    msk = np.where(acc_d & (cl.longitude.values[None, :] >= l0) & (cl.longitude.values[None, :] < l1), 1.0, np.nan)
    ax.contourf(cl.longitude, cl.latitude, msk, levels=[0.5, 1.5], colors="none", hatches=["////"], zorder=4)
    ax.contour(cl.longitude, cl.latitude, np.nan_to_num(msk), levels=[0.5], colors="k", linewidths=0.6, zorder=5)
    ax.text((l0 + l1) / 2, -50.3, lab, ha="center", va="center", fontsize=5.5,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8), zorder=7)
ax.text(-66, -62.4, "Drake box", ha="center", va="center", fontsize=5.5, zorder=7,
        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
ax.plot(trk.lon[trk.perm], trk.lat[trk.perm], color=C["obs"], lw=1.6, zorder=8, solid_capstyle="butt")
ax.plot(trk.lon[~trk.perm], trk.lat[~trk.perm], color=C["obs"], lw=1.2, ls=(0, (2, 1.2)), zorder=8)
ax.set_xlim(-106, -28); ax.set_ylim(-66, -49); ax.set_aspect(1 / np.cos(np.radians(57.5)))
ax.set_xticks(np.arange(-100, -29, 15)); ax.set_xticklabels([f"{-v}°W" for v in np.arange(-100, -29, 15)])
ax.set_yticks([-65, -60, -55, -50]); ax.set_yticklabels([f"{-v}°S" for v in (-65, -60, -55, -50)])
ax.tick_params(labelsize=5.5)
S.hcbar(fig, im, ax, "Mean EKE (cm$^2$ s$^{-2}$)", shrink=0.7, pad=0.12)
S.plabel(ax, "a", "Drake Passage regions and the LM Gould track")

ax = fig.add_subplot(gs[0:2, 1])
R_ = {"box": "Drake box (80–52°W, 63–54°S)", "up": "Upstream (100–80°W, ACC band)", "down": "Downstream (55–35°W, ACC band)"}
groups = []
for reg, rl in R_.items():
    items = [(vl, [(dk(f"allsat_{reg}_{var}"), C[col], "o")]) for var, vl, col in
             (("KE", "KE", "ke"), ("EKE", "EKE", "eke"), ("SPD", "along-stream speed", "tr"))]
    items.append(("eastward wind stress", [(dk(f"wind_{reg}_TAUX"), C["wind"], "o")]))
    groups.append((rl, items))
_d1 = pd.read_csv(PO / "d1" / "d1b_trends.csv").set_index(["series", "period"])
def _t(s_, p_):
    r = _d1.loc[(s_, p_)]; return dict(pct=r.pct_per_decade, pct_lo=r.ci_lo_pct, pct_hi=r.ci_hi_pct, p=r.p)
groups.append(("Drake Passage transport", [("surface proxy, north, 2001–2025", [(_t("Tperm_gv", "2001-2025 monthly"), C["tr"], "s")]),
                                           ("in situ 0–760 m, 2005–2019", [(_t("U_tot", "2005-2019 annual"), C["obs"], "s")])]))
grouped_forest(ax, groups, (-15, 28), "Trend (% of mean per decade); 2001–2025 unless stated")
sig_legend(ax, {"KE": C["ke"], "EKE": C["eke"], "speed": C["tr"], "wind": C["wind"], "in situ": C["obs"]}, loc="lower left",
           bbox_to_anchor=(-0.55, 1.0), ncol=5, fontsize=5.4)
S.plabel(ax, "b", x=-0.62, dy=16)

# c: trends along the LM Gould transect, passage coordinates (as Gutierrez-Villanueva et al., 2023, Figs. 3-4)
K5 = pd.read_csv(DO / "dk5" / "dk5_trends.csv"); K5 = K5[K5.coord == "passage"]
def prof(src, var, per):
    q = K5[(K5.src == src) & (K5["var"] == var) & (K5.period == per)].sort_values("bin"); return q
ax = fig.add_subplot(gs[1, 0]); ax2 = ax.twinx()
mu = prof("allsat", "u", "2001-2025")
ax2.fill_between(mu.bin, 0, mu["mean"], color="#EDEDED", lw=0, zorder=0); ax2.set_ylim(0, 45)
ax2.set_ylabel("Mean (cm s$^{-1}$)", color="#8A8A8A", labelpad=1); ax2.spines["right"].set_visible(True); ax2.tick_params(colors="#8A8A8A")
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.axvspan(408, 850, facecolor="none", hatch="....", edgecolor="#C9DDF0", lw=0, zorder=0)
for (src, per), col, lab, ms in ((("insitu", "2005-2019"), C["obs"], "in situ ADCP (70–118 m), 2005–2019", 3.0),
                                (("allsat", "2005-2019"), C["eke"], "altimetry, 2005–2019", 2.6),
                                (("allsat", "2001-2025"), C["ink"], "altimetry, 2001–2025", 2.6)):
    q = prof(src, "u", per); sig = q.p.round(3) < 0.05
    ax.plot(q.bin, 10 * q.slope, "-", color=col, lw=0.9, label=lab)
    ax.plot(q.bin[sig], 10 * q.slope[sig], "o", color=col, ms=ms); ax.plot(q.bin[~sig], 10 * q.slope[~sig], "o", mfc="white", color=col, ms=ms, mew=0.6)
S.zero(ax); ax.set_xlim(0, 850); ax.set_xlabel("Distance from the north (km)"); ax.set_ylabel("Trend (cm s$^{-1}$ per decade)")
for x_, lab in ((40, "SAF"), (400, "PF"), (790, "SACCF")):
    ax.text(x_, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.5, color="#555555")
ax.legend(loc="upper left", bbox_to_anchor=(0, -0.28), ncol=3, fontsize=5.0, columnspacing=0.8)
S.plabel(ax, "c", "Cross-track velocity trends along the transect", dx=-24, dy=10)

# d, e: streamwise frontal zones (pairs of SSH streamlines, GV23 Figs. 5 and 7)
ZN = pd.read_csv(DO / "dk5" / "dk5_zones.csv"); ZN = ZN[ZN.coord == "stream"]
def zr(zone, src, var, per):
    q = ZN[(ZN.zone == zone) & (ZN.src == src) & (ZN["var"].str.startswith(var)) & (ZN.period == per)].iloc[0]
    return dict(pct=q.pct, pct_lo=q.pct_lo, pct_hi=q.pct_hi, p=q.p)
zones = [("SAF (-0.40 to -0.60 m)", "SAF"), ("SAF-PF zone (-0.60 to -1.00 m)", "SAF–PF zone"), ("PF (-1.00 to -1.25 m)", "PF")]
gsd = gs[2, :].subgridspec(1, 2, wspace=0.9)
for k, (var, letter, ttl, xl) in enumerate((("transport", "d", "Transport between streamlines", (-25, 25)), ("EKE", "e", "EKE between streamlines", (-80, 80)))):
    ax = fig.add_subplot(gsd[0, k])
    groups = [(None, [(lab, [(zr(z, "insitu", var, "2005-2019"), C["obs"], "o"), (zr(z, "allsat", var, "2005-2019"), C["eke"], "o"),
                             (zr(z, "allsat", var, "2001-2025"), C["ink"], "o")])]) for z, lab in zones]
    groups = [(None, [it for g in groups for it in g[1]])]
    grouped_forest(ax, groups, xl, "Trend (% of mean per decade)", tick_size=5.8)
    S.plabel(ax, letter, ttl, x=-0.3, dy=4)
S.save(fig, FIGD, "Figure1"); plt.close(fig)

# ================================================================== Figure 2: comparison with in situ
m = sio.loadmat(GV / "Datasets" / "TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
trs = np.array(m.transp, float); blat = np.array(m.lat, float); dist = np.array(m.dist, float); trs[dist > 840, :] = np.nan
time = pd.to_datetime(np.array(m.time, float) - 719529, unit="D"); nval = np.isfinite(trs).sum(0)
def ins_pos(a, b):
    out = np.full(trs.shape[1], np.nan)
    for k in range(trs.shape[1]):
        if nval[k] < 30: continue
        bb, qq = blat[:, k], trs[:, k]; s = np.isfinite(bb) & np.isfinite(qq) & (bb >= min(a, b)) & (bb <= max(a, b))
        if s.sum() < 3: continue
        o = np.argsort(-bb[s]); xx, q = bb[s][o], np.clip(qq[s][o], 0, None); c = np.cumsum(q)
        if c[-1] > 0: out[k] = np.interp(0.5 * c[-1], c, xx)
    return out
ins = pd.DataFrame(dict(year=time.year, pos=ins_pos(-55.1, -59.33), n=nval)); ins = ins[ins.n >= 30]
split = pd.read_csv(PO / "d1" / "d1b_insitu_split.csv")
ekei = pd.read_csv(DO / "dk3" / "dk3_eke_transects_allsat.csv"); ekei = ekei[ekei.band == "north"]
ser = pd.read_csv(DO / "dk4" / "dk4_section_series.csv"); pos = pd.read_csv(DO / "dk2" / "dk2_positions.csv")
d1t = pd.read_csv(PO / "d1" / "d1b_trends.csv").set_index(["series", "period"])
OV = (2005, 2019)
def rel(s):
    mu = s[(s.index >= OV[0]) & (s.index <= OV[1])].mean(); return 100 * (s / mu - 1)
def z(s):
    s2 = s[(s.index >= OV[0]) & (s.index <= OV[1])]; return (s - s2.mean()) / s2.std()

fig = plt.figure(figsize=(S.W2, 0.74 * S.W2))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], wspace=0.55, hspace=0.6)
ax = fig.add_subplot(gs[0, 0])
pr = pd.read_csv(DO / "dk3" / "dk3_velocity_profile_allsat.csv"); pr = pr[pr.uin_count >= 30]
ax.plot(pr.latb, pr.uin_mean, color=C["obs"], lw=1.2, label="in situ ADCP, 70–118 m")
ax.plot(pr.latb, pr.ualt_mean, color=C["ink"], lw=1.0, label="altimetry")
ax.set_xlim(-55, -63); S.zero(ax); ax.set_xlabel("Latitude (°)"); ax.set_ylabel("Mean cross-track velocity (cm s$^{-1}$)")
ax.legend(loc="upper right", fontsize=5.5); S.plabel(ax, "a", "Mean velocity profile", dx=-26)

def series_panel(ax, alt, insitu, col, ylab, letter, title, relative=True):
    y = rel(alt) if relative else alt
    ax.plot(alt.index, y, "-o", color=col, lw=1.0, ms=1.8, label="altimetry, 2001–2025")
    ii = insitu.loc[OV[0]:OV[1]]; yi = rel(ii) if relative else ii
    ax.plot(ii.index, yi, "-o", color=C["obs"], lw=1.2, ms=2.4, zorder=4, label="in situ, 2005–2019")
    ax.axvspan(OV[0] - 0.5, OV[1] + 0.5, color=C["light"], lw=0, zorder=0)
    if relative: S.zero(ax)
    S.ts_style(ax, (2000.5, 2025.5)); ax.set_ylabel(ylab); S.plabel(ax, letter, title, dx=-26)

ax = fig.add_subplot(gs[0, 1])
series_panel(ax, ser.groupby("year").allsat_EKE_north.mean(), ekei.groupby("year").eke_in.mean(), C["eke"],
             "EKE (% of 2005–2019 mean)", "b", "EKE north of 58.5°S")
ax.legend(loc="upper left", fontsize=5.3)
ax = fig.add_subplot(gs[0, 2])
series_panel(ax, pos.groupby("year").allsat_front1_lat.mean(), ins.groupby("year").pos.mean(), C["pos"],
             "Latitude (°)", "c", "Position of the energetic band", relative=False)

ax = fig.add_subplot(gs[1, 0:2])
pa = z(ser.groupby("year").allsat_Tperm.mean())
sp = split[(split.year >= OV[0]) & (split.year <= OV[1])].groupby("year")[["U_tot", "U_north", "U_south"]].mean()
ax.axvspan(OV[0] - 0.5, OV[1] + 0.5, color=C["light"], lw=0, zorder=0)
ax.plot(pa.index, pa, "-o", color=C["tr"], lw=1.1, ms=1.8, label="surface transport proxy, north of 58.55°S (altimetry)")
for c_, lab, col, ls in (("U_north", "in situ 0–760 m, north of 58.55°S", C["obs"], "-"),
                         ("U_south", "in situ 0–760 m, south of 58.55°S", "#8C5A00", (0, (3, 1.2))),
                         ("U_tot", "in situ 0–760 m, full section", "#222222", (0, (1, 1)))):
    ax.plot(sp.index, z(sp[c_]), marker="o", color=col, lw=1.0, ls=ls, ms=2.0, label=lab, zorder=4)
S.zero(ax); S.ts_style(ax, (2000.5, 2025.5)); ax.set_ylabel("Standardised annual anomaly")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.3)
S.plabel(ax, "d", x=-0.08, dy=26)

ax = fig.add_subplot(gs[1, 2])
def d1(s, p):
    r = d1t.loc[(s, p)]; return dict(pct=r.pct_per_decade, pct_lo=r.ci_lo_pct, pct_hi=r.ci_hi_pct, p=r.p)
groups = [("Altimetry, 2001–2025", [("proxy, north", [(d1("Tperm_gv", "2001-2025 monthly"), C["tr"], "o")])]),
          ("In situ 0–760 m, 2005–2019", [("full section", [(d1("U_tot", "2005-2019 annual"), "#222222", "o")]),
                                           ("north of 58.55°S", [(d1("U_north", "2005-2019 annual"), C["obs"], "o")]),
                                           ("south of 58.55°S", [(d1("U_south", "2005-2019 annual"), "#8C5A00", "o")])])]
grouped_forest(ax, groups, (-60, 60), "Trend (% of mean per decade)", head_size=5.6, tick_size=5.6)
S.plabel(ax, "e", "Transport trends", x=-0.55, dy=10)
S.save(fig, FIGD, "FigureS_comparison"); plt.close(fig)

# ================================================================== Figure 4: circumpolar energy, earlier estimates
M = xr.open_dataset(OUT / "a1_maps_bc.nc")
m1 = pd.read_csv(M1F).set_index("series")
TE = pd.read_csv(OUT / "eddies" / "e2_e3_trends.csv").set_index("series")
A8 = pd.read_csv(OUT / "a8" / "a8_trends_allsat.csv").set_index("series")
RT = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")
lat, lon = ice_ds.latitude.values, ice_ds.longitude.values

fig = plt.figure(figsize=(S.W2, 0.92 * S.W2))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.05], height_ratios=[1.15, 1], wspace=0.75, hspace=0.42)
ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo())
S.polar(ax, ccrs, cfeature)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
im = ax.pcolormesh(L1, A1, M["EKE_slope"].values, transform=PC, cmap=tc, vmin=-6, vmax=6, shading="nearest", rasterized=True, zorder=0)
bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
add_ice_mask(ax, ccrs, lon, lat, ice_g); add_sector_boundaries(ax, ccrs)
jj, ii = np.nonzero(M["EKE_fdr"].values == 1)
ax.scatter(M.lon.values[ii], M.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
def outline(l0, l1, a0, a1, **kw):
    lo = np.r_[np.linspace(l0, l1, 60), np.full(10, l1), np.linspace(l1, l0, 60), np.full(10, l0)]
    la = np.r_[np.full(60, a0), np.linspace(a0, a1, 10), np.full(60, a1), np.linspace(a1, a0, 10)]
    ax.plot(lo, la, transform=PC, **kw)
for l0, l1, a0, a1 in ((40, 150, -57, -44), (150, 288, -62, -48), (-35, 10, -56, -46)):
    outline(l0, l1, a0, a1, color="#6A3D9A", lw=0.9, zorder=8)
for l0, l1 in ((75, 95), (145, 175), (-150, -125), (-55, -35)):     # wake boxes (downstream), within the ACC band
    outline(l0, l1, -62, -45, color="#E69F00", lw=0.8, ls=(0, (2, 1)), zorder=8)
S.plabel(ax, "a", "EKE trend, 2001–2025")
S.hcbar(fig, im, ax, "EKE trend (cm$^2$ s$^{-2}$ yr$^{-1}$)", shrink=0.7, extend="both")
ax.legend(handles=[Line2D([], [], color="#6A3D9A", lw=0.9, label="Hogg et al. (2015) boxes"),
                   Line2D([], [], color="#E69F00", lw=0.8, ls=(0, (2, 1)), label="wake boxes (ACC band only)"),
                   Line2D([], [], ls="none", marker="o", ms=2, color="k", label="significant after FDR")],
          loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=1, fontsize=5.3)

def m1r(key):
    r = m1.loc[key]; return dict(pct=r.pct_per_decade, pct_lo=r.pct_ci_lo, pct_hi=r.pct_ci_hi, p=r.p)
ax = fig.add_subplot(gs[0, 1])
TS20, TS25, LIT = C["two"], C["eke"], "#222222"
groups = []
for box_, head, lit in (("Hogg_Indian", "Indian sector box", ("Hogg et al. (2015), from 1993", 6.0, 1.6)),
                         ("Hogg_Pacific", "Pacific sector box", ("Hogg et al. (2015), from 1993", 6.4, 1.8)),
                         ("Hogg_Atlantic", "Atlantic sector box", None),
                         ("box35_65", "35–65°S", ("Martínez-Moreno et al. (2021), 1993–2020", 5.1, np.nan)),
                         ("ACC", "ACC band", None)):
    items = [("2001–2020", [(m1r(f"allsat|{box_}_EKE|2001-2020"), TS25, "o")]),
             ("2001–2025", [(m1r(f"allsat|{box_}_EKE|2001-2025"), TS25, "o")])]
    if lit:
        v, e = lit[1], lit[2]
        items.append((lit[0], [(dict(pct=v, pct_lo=v - e if np.isfinite(e) else v, pct_hi=v + e if np.isfinite(e) else v, p=0.0), LIT, "D")]))
    groups.append((head, items))
grouped_forest(ax, groups, (-12, 22), "EKE trend (% per decade)", head_size=5.8, tick_size=5.5)
ax.legend(handles=[Line2D([], [], ls="-", lw=0.8, marker="o", color=TS25, ms=3.4, label="this study (open: p ≥ 0.05)"),
                   Line2D([], [], ls="-", lw=0.8, marker="D", color=LIT, ms=3.4, label="earlier estimate, as published")],
          loc="lower left", bbox_to_anchor=(-0.9, 1.0), ncol=2, fontsize=5.3)
S.plabel(ax, "b", x=-0.95, dy=16)

def rt(key):
    r = RT.loc[key]; mean = r.slope * 1000 / r.pct_per_decade; f = 1000 / mean
    return dict(pct=r.pct_per_decade, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p_new)
def a8(reg):
    r = A8.loc[f"{reg}_T"]; return dict(pct=r.pct_per_decade, pct_lo=r.pct_lo, pct_hi=r.pct_hi, p=r.p)
ax = fig.add_subplot(gs[1, 0])
groups = []
for reg, head in (("ACC", "ACC band"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian")):
    items = [("EKE", [(rt(f"allsat|{reg}_EKE|2001-2025"), C["eke"], "o")]), ("KE", [(rt(f"allsat|{reg}_KE|2001-2025"), C["ke"], "o")]),
             ("along-stream speed", [(rt(f"allsat|{reg}_SPD|2001-2025"), C["tr"], "o")]),
             ("surface transport proxy", [(a8(reg), C["tr"], "s")])]
    if reg == "ACC":
        items.append(("wind stress, 45–60°S", [(rt("tau_Circumpolar_45_60S|2001-2025"), C["wind"], "o")]))
    groups.append((head, items))
grouped_forest(ax, groups, (-4, 19), "Trend 2001–2025 (% per decade)", head_size=5.8, tick_size=5.5)
S.plabel(ax, "c", x=-0.75, dy=6)

def erow(k):
    r = TE.loc[k]; f = 1000 / r["mean"]; return dict(pct=r.pct_per_decade, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p)
ax = fig.add_subplot(gs[1, 1])
groups = [("Tracked eddies, ACC band", [("mean amplitude", [(erow("ACC_amp"), C["te"], "D")]),
                                         ("rotational speed", [(erow("ACC_spd"), C["te"], "D")]),
                                         ("number present", [(erow("ACC_present_per_day"), C["te"], "D")]),
                                         ("births", [(erow("ACC_births"), C["te"], "D")]),
                                         ("speed radius", [(erow("ACC_rad"), C["te"], "D")])]),
          ("Gridded, same months", [("EKE", [(erow("ACC_EKE_allsat"), C["eke"], "o")])])]
groups = [(h, [(l, e) for l, e in it if e]) for h, it in groups]
grouped_forest(ax, groups, (-4, 16), "Trend Jan 2001–Aug 2023 (% per decade)", head_size=5.8, tick_size=5.5)
S.plabel(ax, "d", x=-0.95, dy=6)
S.save(fig, FIGD, "Figure2"); plt.close(fig)

# ================================================================== Figure 5: pathway and meanders
SM = pd.read_csv(OUT / "a2" / "a2_summary.csv"); SM = SM[SM.set == "baseline"].set_index("region")
SM = retest_override(SM, lambda g, r: f"baseline|{g}|lat_hp", {"m_lat_hp_p": "p_new", "m_lat_hp_lo": "ci_lo", "m_lat_hp_hi": "ci_hi"})
G = pd.read_csv(OUT / "a2" / "a2g_literature_regions.csv").set_index("region"); KM = 111.195
def gr(name):
    r = G.loc[name]; return dict(pct=r.km_dec, pct_lo=r.lo, pct_hi=r.hi, p=r.p)
def sec(reg):
    q = SM.loc[reg]; return dict(pct=q.m_lat_hp_slope * KM * 10, pct_lo=q.m_lat_hp_lo * KM * 10, pct_hi=q.m_lat_hp_hi * KM * 10, p=q.m_lat_hp_p)
CL = pd.read_csv(OUT / "b3" / "b3_climatology.csv")
T = {v: pd.read_csv(OUT / "b3" / f"b3_trends_{v}.csv") for v in ("eke", "speed", "sinu", "lat")}
SHALLOW = CL.lon[CL.depth_min < 3000].values; MALV = CL.lon[(CL.lon > -60) & (CL.lon < -50)].values

fig = plt.figure(figsize=(S.W2, 1.18 * S.W2))
gs = fig.add_gridspec(6, 1, height_ratios=[1.75, 1.0, 1, 1, 1, 1], hspace=0.75)
ax = fig.add_subplot(gs[0])
_dk2 = pd.read_csv(DO / "dk2" / "dk2_trends.csv"); _r = _dk2[(_dk2["product"] == "allsat") & (_dk2.period == "2001-2025")].iloc[0]
drake_sec = dict(pct=_r.trend_km_dec, pct_lo=_r.ci_lo, pct_hi=_r.ci_hi, p=_r.p)
rows = [("Persistently tracked longitudes", {"t": gr("Circumpolar")}), ("Atlantic sector", {"t": sec("Atlantic")}),
        ("Indian sector", {"t": sec("Indian")}), ("Pacific sector", {"t": sec("Pacific")}),
        ("90–170°E", {"t": gr("90-170E (Chambers 2018)")}), ("75–150°E", {"t": gr("75-150E SE Indian (Kim and Orsi 2014)")}),
        ("160–90°W", {"t": gr("160-90W (Chambers 2018)")}), ("150–80°W (PAR meander)", {"t": gr("150-80W PAR meander (Liu et al. 2026)")}),
        ("Drake Passage, LM Gould section (merged SAF–PF band)", {"t": drake_sec})]
rows = [(l, {k: rp(v) for k, v in d.items()}) for l, d in rows]
ys = S.forest(ax, rows, {"t": C["pos"]}, xlim=(-45, 50), xlabel="Half-power latitude trend (km per decade, positive northward)")
LITM = [("Sokolov and Rintoul (2009b), 1992–2007", 0, -40, "v", "#7B3F00"), ("Gille (2014), SSH contours", 0, -30, "^", "#B07D2B"),
        ("Chambers (2018), 1993–2016", 4, -29, "P", "#444444"), (None, 6, 46, "P", "#444444"),
        ("Liu et al. (2026), 1993–2023", 7, 3.3, "X", "#56B4E9")]
for lab, row, val, mk, col in LITM:
    ax.plot(val, ys[row] - 0.3, mk, ms=4.2, color=col, mec="white", mew=0.4, zorder=6)
ax.text(44, ys[5], "Kim and Orsi (2014):\nsouthward drift, no rate", fontsize=5.0, va="center", ha="right", color="#555555")
h = [Line2D([], [], ls="-", lw=0.8, marker="o", color=C["pos"], mfc="white", ms=3.6, label="this study, 2001–2025 (open: p ≥ 0.05)")]
h += [Line2D([], [], ls="none", marker=mk, color=col, mec="white", mew=0.4, ms=4.2, label=lab) for lab, _, _, mk, col in LITM if lab]
ax.legend(handles=h, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=1, fontsize=5.3)
S.plabel(ax, "a", "Position of the energetic band (half-power latitude)", x=-0.28)
# b: 5-degree longitude bins (revised test) with the detection frequency of the matched envelope
LB = pd.read_csv(OUT / "retest" / "retest_summary.csv"); LB = LB[LB.key.str.startswith("lonbin|")].copy()
LB["lon0"] = LB.key.str.split("|").str[1].astype(float); LB["c"] = LB.lon0 + 2.5
DF = pd.read_csv(OUT / "a2" / "detection_frequency_lon.csv")
ax = fig.add_subplot(gs[1]); ax2 = ax.twinx()
ax2.fill_between(DF.lon, 0, 100 * DF.freq, color="#E6E6E6", lw=0, zorder=0, step="mid"); ax2.axhline(90, color="#9A9A9A", lw=0.5, ls="--")
ax2.set_ylim(0, 100); ax2.set_ylabel("Detection (% of months)", color="#7F7F7F"); ax2.spines["right"].set_visible(True)
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
sigb = LB.p_new.round(3) < 0.05
ax.vlines(LB.c, LB.ci_lo * KM * 10, LB.ci_hi * KM * 10, color=C["pos"], lw=0.6)
ax.plot(LB.c[~sigb], LB.slope[~sigb] * KM * 10, "o", mfc="white", color=C["pos"], ms=2.6, mew=0.6)
ax.plot(LB.c[sigb], LB.slope[sigb] * KM * 10, "o", color=C["pos"], ms=2.6)
S.zero(ax); ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 30)); ax.set_xticklabels([])
ax.set_ylabel("km per decade")
S.plabel(ax, "b", "Half-power latitude trend in 5° longitude bins (shading: detection frequency)", dx=-28)

def bands(ax):
    for l in SHALLOW: ax.axvspan(l - 2.5, l + 2.5, color="#EDEDED", lw=0, zorder=0)
    for l in MALV: ax.axvspan(l - 2.5, l + 2.5, facecolor="none", edgecolor="#BDBDBD", hatch="////", lw=0, zorder=0)
spec = [("eke", "c", "EKE within 150 km of the jet-core path", "% per decade", C["eke"]),
        ("speed", "d", "Jet speed along the path", "% per decade", C["tr"]),
        ("sinu", "e", "Path sinuosity", "% per decade", C["sin"]),
        ("lat", "f", "Path latitude (positive northward)", "km per decade", C["pos"])]
axs = [fig.add_subplot(gs[k + 2]) for k in range(4)]
for ax, (v, l, ttl, unit, col) in zip(axs, spec):
    bands(ax); S.zero(ax)
    d = T[v].dropna(subset=["tr"]); d = d[~d.lon.between(-60, -50)]; sig = d.p.round(3) < 0.05
    ax.vlines(d.lon, d.lo, d.hi, color=col, lw=0.6)
    ax.plot(d.lon[sig], d.tr[sig], "o", color=col, ms=2.6, zorder=4)
    ax.plot(d.lon[~sig], d.tr[~sig], "o", mfc="white", color=col, ms=2.6, mew=0.6, zorder=4)
    f = d.fdr10.astype(bool) if "fdr10" in d else np.zeros(len(d), bool)
    ax.plot(d.lon[f], d.tr[f], "o", mfc="none", mec="#111111", ms=5.5, mew=0.7, zorder=5)
    lo_, hi_ = np.nanmin(d.lo), np.nanmax(d.hi); pad = 0.06 * (hi_ - lo_); ax.set_ylim(lo_ - pad, hi_ + pad)
    ax.set_ylabel(unit); ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 30))
    ax.set_xticklabels([] if v != "lat" else [f"{abs(t)}°{'W' if t < 0 else 'E' if t > 0 else ''}" if abs(t) != 180 else "180°" for t in np.arange(-180, 181, 30)])
    S.plabel(ax, l, ttl, dx=-28)
for x_, s in [(-135, "PAR"), (-65, "Drake"), (40, "SWIR"), (78, "Kerguelen"), (133, "SEIR"), (163, "Macquarie")]:
    axs[0].text(x_, axs[0].get_ylim()[1], s, ha="center", va="bottom", fontsize=5.0, color="#555555")
axs[-1].set_xlabel("Longitude (flow from left to right)")
S.save(fig, FIGD, "Figure3"); plt.close(fig)
print("done")
