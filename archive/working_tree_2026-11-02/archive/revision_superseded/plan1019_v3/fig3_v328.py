"""Figure 3 (v3.28): where the surface eddy field intensified, in streamline coordinates. Product 1 (all-satellite).
a  EKE trend on 1 deg cells (a1 maps; dots: FDR), mean SAF and PF streamlines, 3000 m isobath, sea-ice mask, the five
   meander regions (ridge + downstream, dk13).
b  Trends of EKE and KE (fixed and moving bins) and of the surface transport proxy (fixed bins, dk18) against the 2001-2025 mean ADT (0.05 m streamline bins, circumpolar), fixed
   bins (solid, 95 % interval shaded for EKE) and bins following the monthly ADT minus its domain-mean anomaly (dashed);
   above: mean EKE per bin. After the streamwise analysis of Gutierrez-Villanueva et al.
   (2023), applied around the ACC (dk17/s2, s5; paper test on monthly series).
c  Circumpolar mean anomaly of the half-power latitude over the persistently tracked longitudes (a2 baseline; km,
   positive northward): monthly, 12-month running mean, trend with 95 % interval (retest).
d  Trends of EKE, KE and the surface transport proxy (a8b, dk18) in the ACC band over the ridge and downstream of the five standing meanders
   and away from them (after Thompson and Naveira Garabato 2014, Fig. 1b).
Usage: python fig3_v328.py <revision dir> <drake out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, DO, FIGD = (Path(a) for a in sys.argv[1:4])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C
OUT = REV / "outputs"; ROOT_OUT = REV.parent / "outputs"; PC = ccrs.PlateCarree(); KM = 111.195
CSPD = "#0072B2"


def sigp(p):
    return np.isfinite(p) and round(float(p), 3) < 0.05


M = xr.open_dataset(OUT / "a1_maps_bc.nc")
clim = xr.open_dataset(OUT / "clim_allsat.nc")
Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z12["lev_saf"]), float(Z12["lev_pf"])
ice_ds = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc"); ice_g = ice_ds.ice_mask.values > 0.5
lat, lon = ice_ds.latitude.values, ice_ds.longitude.values
bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation; bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
P = pd.read_csv(DO / "dk17" / "s5_profiles.csv")
T18 = pd.read_csv(DO / "dk18" / "dk18_trends.csv"); TA8 = pd.read_csv(DO / "a8b" / "a8b_trends_allsat.csv").set_index("key")
TBIN = T18[T18.key.str.startswith("bin|")].copy(); TBIN["mid"] = TBIN.key.str.split("|").str[1].astype(float)
TBIN = TBIN.rename(columns={"pct": "pct_dec", "pct_lo": "pct_dec_lo", "pct_hi": "pct_dec_hi"}).sort_values("mid")
RT = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")
TB = np.load(OUT / "a2" / "tracks_baseline.npz")
MEANDERS = {"Kerguelen Plateau": ((55, 65), (65, 80), (80, 95)), "Southeast Indian Ridge": ((120, 130), (130, 140), (140, 150)),
            "Macquarie Ridge": ((150, 155), (155, 160), (160, 175)), "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
            "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35))}
SHORT = {"Kerguelen Plateau": "Kerguelen\nPlateau", "Southeast Indian Ridge": "Southeast\nIndian Ridge", "Macquarie Ridge": "Macquarie\nRidge",
         "Pacific–Antarctic Ridge": "Pacific–Antarctic\nRidge", "Drake Passage–Scotia Sea": "Drake Passage–\nScotia Sea",
         "Away from the five meanders": "Away from\nthe meanders"}

fig = plt.figure(figsize=(S.W2, 0.92 * S.W2))
gR = fig.add_gridspec(2, 1, height_ratios=[1.25, 0.8], hspace=0.30)
g1 = gR[0].subgridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.16)
# ---------------------------------------------------------------- a map
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
im = ax.pcolormesh(L1, A1, 10 * M["EKE_slope"].values, transform=PC, cmap=tc, vmin=-60, vmax=60, shading="nearest", rasterized=True, zorder=0)
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[LEV_PF, LEV_SAF], colors="#1A1A1A", linewidths=0.6, linestyles="solid", transform=PC, zorder=6)
add_ice_mask(ax, ccrs, lon, lat, ice_g)
jj, ii = np.nonzero(M["EKE_fdr"].values == 1)
ax.scatter(M.lon.values[ii], M.lat.values[jj], s=2.5, c="k", lw=0, transform=PC, zorder=7)
for nm, (u_, r_, d_) in MEANDERS.items():          # ridge + downstream arcs just north of the band
    for a_, b_ in (r_, d_):
        if nm.startswith("Drake") and a_ == -50: a_ = -50
        xs = np.linspace(a_, b_, 30); ax.plot(xs, np.full_like(xs, -38.2), color="#A2508A", lw=1.6, transform=PC, zorder=8, solid_capstyle="butt")
for lo_, la_, s_ in ((75, -44.0, "KP"), (135, -42.5, "SEIR"), (160, -44.5, "MR"), (-130, -44.0, "PAR"), (-55, -44.0, "DP–SS")):
    ax.text(lo_, la_, s_, transform=PC, fontsize=4.8, ha="center", va="center", color="#A2508A", zorder=9,
            bbox=dict(boxstyle="round,pad=0.06", fc="white", ec="none", alpha=0.8))
S.plabel(ax, "a", "EKE trend, 2001–2025")
S.hcbar(fig, im, ax, "EKE trend (cm$^2$ s$^{-2}$ per decade)", shrink=0.7, extend="both")
# ---------------------------------------------------------------- b streamline profile
gb = g1[1].subgridspec(2, 1, height_ratios=[0.32, 1.0], hspace=0.06)
axm = fig.add_subplot(gb[0]); ax = fig.add_subplot(gb[1], sharex=axm)
q = lambda b, v: P[(P["bin"] == b) & (P.sector == "ACC") & (P["var"] == v)].sort_values("mid")
for a_ in (axm, ax):
    a_.axvspan(-0.6, 0.2, color="#F1F1F1", zorder=0)
    for lv, lab in ((LEV_PF, "PF"), (LEV_SAF, "SAF")):
        a_.axvline(lv, color="#555555", lw=0.5, ls=(0, (2, 2)), zorder=1)
e = q("madt", "EKE")
axm.fill_between(e.mid, 0, e["mean"], color=C["eke"], alpha=0.25, lw=0); axm.plot(e.mid, e["mean"], color=C["eke"], lw=0.8)
axm.set_ylim(0, 340); axm.set_ylabel("EKE\n(cm$^2$ s$^{-2}$)", color=C["eke"], fontsize=5.5); axm.tick_params(axis="y", labelsize=5, colors=C["eke"])
plt.setp(axm.get_xticklabels(), visible=False)
for lv, lab in ((LEV_PF, "PF"), (LEV_SAF, "SAF")):
    axm.text(lv, 345, lab, ha="center", va="bottom", fontsize=5.5, color="#333333")
axm.text(-0.2, 345, "ACC band", ha="center", va="bottom", fontsize=5.0, color="#777777")
S.plabel(axm, "b", "Trends along streamlines, circumpolar", dy=12)
S.zero(ax)
for v, col in (("EKE", C["eke"]), ("KE", C["ke"]), ("T", CSPD)):
    f = q("madt", v) if v != "T" else TBIN; m = q("iadt", v) if v != "T" else None
    if v == "EKE":
        ax.fill_between(f.mid, f.pct_dec_lo, f.pct_dec_hi, color=col, alpha=0.12, lw=0, zorder=2)
    ax.plot(f.mid, f.pct_dec, color=col, lw=0.9, zorder=4)
    if m is not None: ax.plot(m.mid, m.pct_dec, color=col, lw=0.7, ls=(0, (3, 1.6)), zorder=4)
    sg = f.p.round(3) < 0.05
    ax.scatter(f.mid[sg], f.pct_dec[sg], s=7, color=col, zorder=5, lw=0)
    ax.scatter(f.mid[~sg], f.pct_dec[~sg], s=7, facecolor="white", edgecolor=col, lw=0.5, zorder=5)
ax.set_xlim(-0.9, 0.5); ax.set_ylim(-6, 20)
ax.set_xlabel("2001–2025 mean ADT (m)   south ← → north"); ax.set_ylabel("trend (% per decade)")
ax.legend(handles=[Line2D([], [], color=C["eke"], lw=0.9, marker="o", ms=2.5, label="EKE"), Line2D([], [], color=C["ke"], lw=0.9, marker="o", ms=2.5, label="KE"),
                   Line2D([], [], color=CSPD, lw=0.9, marker="o", ms=2.5, label="surface transport proxy"),
                   Line2D([], [], color="#555555", lw=0.7, ls=(0, (3, 1.6)), label="EKE, KE: bins following the monthly ADT"),
                   Line2D([], [], color="#555555", lw=0, marker="o", mfc="white", ms=2.5, mew=0.5, label="open: p ≥ 0.05")],
          loc="upper left", ncol=3, fontsize=4.9)
# ---------------------------------------------------------------- c half-power series
g2 = gR[1].subgridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.22)
ax = fig.add_subplot(g2[0])
keep = TB["keep"].astype(bool); a = TB["lat_hp"][:, keep]
y = KM * np.nanmean(a - np.nanmean(a, 0, keepdims=True), 1)
t = 2001 + (np.arange(300) + 0.5) / 12
ax.plot(t, y, color="#D0D0D0", lw=0.35, zorder=2, label="monthly")
yr = pd.Series(y).rolling(12, center=True).mean()
ax.plot(t, yr, color=C["pos"], lw=1.0, zorder=3, label="12-month running mean")
r = RT.loc["baseline|ACC|lat_hp"]; tc_ = t.mean(); xx = np.array([2001, 2026])
for b_, ls in ((r.slope, "-"),):
    ax.plot(xx, KM * b_ * (xx - tc_), color="k", lw=0.8, ls=(0, (3, 1.6)), zorder=4)
ax.fill_between(xx, KM * r.ci_lo * (xx - tc_), KM * r.ci_hi * (xx - tc_), color="k", alpha=0.10, lw=0, zorder=1)
S.zero(ax); S.ts_style(ax, (2000.8, 2026.2)); ax.set_ylim(-30, 30)
ax.set_ylabel("half-power latitude anomaly (km, + north)")
ax.text(0.02, 0.04, f"trend {10 * KM * r.slope:+.1f} km per decade ({10 * KM * r.ci_lo:+.1f} to {10 * KM * r.ci_hi:+.1f}); "
        f"25 km would be detected", transform=ax.transAxes, fontsize=5.0, ha="left", va="bottom")
ax.legend(loc="upper left", ncol=2, fontsize=5.0)
S.plabel(ax, "c", "Position of the energetic band, 900 tracked longitudes")
# ---------------------------------------------------------------- d meanders vs elsewhere
ax = fig.add_subplot(g2[1]); S.zero(ax)
R_ = P[P["bin"] == "region"]; names = list(MEANDERS) + ["Away from the five meanders"]; w_ = 0.24
def treg(nm):
    r = TA8.loc[nm] if nm in TA8.index else T18.set_index("key").loc[nm]
    return pd.Series(dict(pct_dec=r.pct, pct_dec_lo=r.pct_lo, pct_dec_hi=r.pct_hi, p=r.p))
for k, (v, col) in enumerate((("EKE", C["eke"]), ("KE", C["ke"]), ("T", CSPD))):
    for gi, nm in enumerate(names):
        rr = R_[(R_.sector == nm) & (R_["var"] == v)].iloc[0] if v != "T" else treg(nm); x = gi + (k - 1) * w_
        ax.bar(x, rr.pct_dec, w_ * 0.9, color=col if sigp(rr.p) else "white", edgecolor=col, lw=0.8, zorder=3)
        ax.plot([x, x], [rr.pct_dec_lo, rr.pct_dec_hi], color="#333333", lw=0.5, zorder=4)
ax.axvline(len(names) - 1.5, color="#9A9A9A", lw=0.5, ls=(0, (2, 2)))
ax.set_xticks(range(len(names))); ax.set_xticklabels(["KP", "SEIR", "MR", "PAR", "DP–SS", "away from\nmeanders"], fontsize=5.5); ax.tick_params(axis="x", length=0)
ax.set_ylabel("trend (% per decade)"); ax.set_ylim(-10, 30)
ax.legend(handles=[Patch(facecolor=C["eke"], label="EKE"), Patch(facecolor=C["ke"], label="KE"), Patch(facecolor=CSPD, label="surface transport proxy"),
                   Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper left", ncol=4, fontsize=5.0)
S.plabel(ax, "d", "Over and downstream of standing meanders, ACC band")
S.save(fig, FIGD, "Figure3"); plt.close(fig)
for nm in names:
    print(nm, R_[R_.sector == nm][["var", "mean", "pct_dec", "pct_dec_lo", "pct_dec_hi", "p", "area_frac"]].round(3).to_string(index=False))
