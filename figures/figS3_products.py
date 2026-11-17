"""Supplement Fig. S3 (final layout, 10 Nov 2026): altimetry products and constellation sensitivity, Product 1 (all-satellite,
1/8 deg) and Product 4 (two-satellite climate product, 1/4 deg), 2001-2025. The figure carries the whole product comparison
(product comparisons live only in the Supplement).

Panels (one message each):
  a  ACC-band EKE and KE, 12-month running means in % of the 2001-2025 mean, both products; 2001-2017 and 2018-2025 means of the
     monthly EKE (short segments) and the increase of the late over the early mean (P1 +13 %, P4 +6 %). Kept as panel a because
     the main text (Sect. 3.2) and SI Sect. S2.3 cite "Fig. S3a" for the 13 % increase.
     Source: outputs/a5/a5_series_{allsat,twosat}.csv (ACC_KE, ACC_EKE; the Product 4 file is reproduced to 1e-13 by p4_01).
  b, c  EKE trend maps on 1-deg cells (10 x slope, cm2 s-2 per decade), Product 1 (outputs/a1/a1_maps_bc.nc, as main Fig. 3a, no
     significance dots) and Product 4 (p4/p4_eke_trend_map.nc), same projection, extent and colour scale; identical overlays on
     both maps: mean SAF and PF (Product 1 mean ADT, levels of dk12), 3000 m isobath, sea ice (as main Fig. 3a). Pattern correlation
     r over the cells finite in both maps (Pearson, unweighted; recomputed here and asserted against p4_summary.csv).
  d  EKE trend between the mean SAF and PF in 2-deg bins with 95 % intervals, P1 (outputs/drake/dk12/dk12_band_trends_allsat.csv)
     and P4 (p4/p4_bins_trends.csv, same band: P1 mean ADT levels), offset -0.45 / +0.45 deg; mean seafloor depth between the
     fronts (as main Fig. 3d); feature keys style.KEYS; flow left to right.
  e  KE and EKE trends of the ACC band and its sectors, % of the mean per decade (outputs/retest/retest_summary.csv, keys
     <allsat|twosat>|<region>_<KE|EKE>|2001-2025: the paper's test).
  f  Position trends, km per decade (north +): half-power latitude (relative threshold, main-text tracks) and the fixed ADT contour
     through its mean position, circumpolar, by sector and 86-40 W (p4/p4_09_relthr_trends.csv: P1 = x7 tracks, P4 = p4_08
     relative-threshold tracks, raw monthly fields as x6).
  g  Lagged correlation of KE and EKE (ACC band + 100 km) with the wind stress, wind leading 0-36 months, P1 and P4
     (p4/p4_lagcorr.csv); maxima marked, filled where the family-wise p over the 37 lags is < 0.05 (p4/p4_lagcorr_summary.csv).
Conventions: Product 1 circles / solid, Product 4 squares / dashed; filled p < 0.05, open p >= 0.05; colours by quantity.
Usage: python figures/figS3_products.py [outputs dir] [figure dir]   (defaults: $OSR11_OUTPUTS, rendered/)
"""
import _sifig as F
import numpy as np, pandas as pd, xarray as xr
from matplotlib.lines import Line2D

S, Q, C, plt, OUT, FIGD = F.S, F.Q, F.C, F.plt, F.OUT, F.FIGD
P = F.P
P4D = OUT / "p4"   # Product 4 diagnostics (scripts/06_product4/)
cartopy = F.cartopy_local()
import cartopy.crs as ccrs, cartopy.feature as cfeature  # noqa: E402
PC = ccrs.PlateCarree()
SUM = pd.read_csv(P4D / "p4_summary.csv").set_index("quantity")
PRINTED = []


def say(*a):
    s = " ".join(str(x) for x in a); PRINTED.append(s); print(s)


H = 229.0
fig = S.figure(S.W2, H)
MS = 3.0
LS = {k: Q.PROD[k]["ls"] for k in ("P1", "P4")}
MK = {k: Q.PROD[k]["marker"] for k in ("P1", "P4")}

# ---------------------------------------------------------------- top legend (products, significance, quantities)
h1 = [Line2D([], [], ls=LS[k], lw=0.8, marker=MK[k], color="#222222", ms=MS, label=Q.PROD[k]["name"]) for k in ("P1", "P4")]
h1 += [Line2D([], [], ls="none", marker="o", color="#555555", ms=MS, label="filled: p < 0.05"),
       Line2D([], [], ls="none", marker="o", color="#555555", mfc="white", mew=0.8, ms=MS, label="open: p ≥ 0.05")]
h2 = [Line2D([], [], color=C["ke"], lw=1.2, label="KE"), Line2D([], [], color=C["eke"], lw=1.2, label="EKE"),
      Line2D([], [], color=C["pos"], lw=1.2, label="half-power latitude"), Line2D([], [], color=C["grey"], lw=1.2, label="fixed ADT contour")]
leg = fig.legend(handles=h1 + h2, loc="upper left", bbox_to_anchor=(2 / S.W2, (H - 0.8) / H), ncol=4, handlelength=2.4,
                 columnspacing=1.6, fontsize=F.FS)

# ---------------------------------------------------------------- a: band series
A = {k: pd.read_csv(OUT / "a5" / f"a5_series_{Q.PROD[k]['key']}.csv") for k in ("P1", "P4")}
yr = A["P1"].year.values; tA = A["P1"].x.values
assert (A["P4"].x.values == tA).all()
YA = 190.0
ax = S.axmm(fig, 16, YA, 122, 20)
inc = {}
for v, col in (("ACC_KE", C["ke"]), ("ACC_EKE", C["eke"])):
    for k in ("P1", "P4"):
        y = A[k][v].values; m = y.mean()
        ax.plot(tA, S.runfull(100 * (y / m - 1)), color=col, ls=LS[k], lw=S.LW["smooth"], zorder=3)
        if v == "ACC_EKE":
            e, l = y[yr <= 2017].mean(), y[yr >= 2018].mean(); inc[k] = 100 * (l / e - 1)
            for (x0, x1, mm) in ((2001, 2018, e), (2018, 2026, l)):
                ax.plot([x0, x1], [100 * (mm / m - 1)] * 2, color="#111111", ls=LS[k], lw=0.8, zorder=4, solid_capstyle="butt")
            say(f"a {k} EKE period means {e:.2f} / {l:.2f} cm2 s-2 -> {100 * (e / m - 1):+.2f} / {100 * (l / m - 1):+.2f} % of mean; increase {inc[k]:+.2f} %")
assert abs(inc["P1"] - SUM.loc["ACC band EKE: 2018-2025 mean over 2001-2017 mean", "P1"]) < 1e-9
assert abs(inc["P4"] - SUM.loc["ACC band EKE: 2018-2025 mean over 2001-2017 mean", "P4"]) < 1e-9
ax.axvline(2018, color="#BDBDBD", lw=0.5, zorder=0)
S.zero(ax); S.ts_axis(ax, 2000.8, 2026.2); ax.set_ylabel("% of 2001–2025 mean"); ax.set_ylim(-14, 16); ax.set_yticks([-10, 0, 10])
S.yticks_minus(ax)
hseg = [Line2D([], [], color="#111111", lw=0.8, ls=LS[k], label=S.um(f"{inc[k]:+.0f} %, {Q.PROD[k]['short']}")) for k in ("P1", "P4")]
fig.legend(handles=hseg, loc="upper left", bbox_to_anchor=(140 / S.W2, (YA + 21) / H), ncol=1, handlelength=2.4, fontsize=F.FS,
           title="EKE: 2018–2025 mean\nover 2001–2017 mean\n(black segments)", title_fontsize=F.FS, alignment="left", borderaxespad=0)
S.letter(fig, 2, YA + 24.5, "a", "ACC-band EKE and KE, 12-month running means")

# ---------------------------------------------------------------- b, c: EKE trend maps
bat = xr.open_dataset(P.BATHY).elevation; bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
clim = xr.open_dataset(P.CLIM)
ice_ds = xr.open_dataset(P.ICE_MASK); ice_g = ice_ds.ice_mask.values > 0.5
Z = np.load(OUT / "drake" / "dk12" / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z["lev_saf"]), float(Z["lev_pf"])
M1 = xr.open_dataset(OUT / "a1" / "a1_maps_bc.nc"); M4 = xr.open_dataset(P4D / "p4_eke_trend_map.nc")
assert np.allclose(M1.lat, M4.lat) and np.allclose(M1.lon, M4.lon)
e1, e4 = 10 * M1["EKE_slope"].values, 10 * M4["EKE_slope"].values; ok = np.isfinite(e1) & np.isfinite(e4)
r_map = np.corrcoef(e1[ok], e4[ok])[0, 1]
assert abs(r_map - SUM.loc["1-deg EKE trend maps: P1-P4 pattern correlation, all common cells (n = 8184), unweighted", "P1"]) < 1e-9 and ok.sum() == 8184
say(f"b/c pattern correlation r = {r_map:.4f} over {ok.sum()} common cells")
tc = plt.get_cmap(S.TREND_CMAP).copy(); tc.set_bad(C["masked"])
D = 54.0; YM = 118.0
for k, (x0, let, M) in {"P1": (11, "b", M1), "P4": (83, "c", M4)}.items():
    ax = S.axmm(fig, x0, YM, D, D, projection=ccrs.SouthPolarStereo())
    S.polar(ax, ccrs, cfeature)
    S.add_isobath(ax, bc.longitude, bc.latitude, -bc.values, transform=PC)
    S.add_fronts(ax, clim.longitude, clim.latitude, clim.adt.values, LEV_PF, LEV_SAF, transform=PC)
    S.add_ice(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_g)
    im = ax.pcolormesh(S.edges(M.lon.values), S.edges(M.lat.values), np.ma.masked_invalid(10 * M["EKE_slope"].values), transform=PC,
                       cmap=tc, vmin=-S.EKE_TREND_LIM, vmax=S.EKE_TREND_LIM, shading="flat", rasterized=True, zorder=0)
    S.letter(fig, 2 if let == "b" else 70, YM + D + 8, let, f"EKE trend, {Q.PROD[k]['name']}" + ("" if let == "b" else f"; pattern correlation with b: r = {r_map:.2f}"))
XR = 147.0
cax = S.axmm(fig, XR, YM + 24, 2.2, 28)
cb = fig.colorbar(im, cax=cax, orientation="vertical", extend="both", extendfrac=0.04, ticks=S.EKE_TREND_TICKS)
cb.ax.set_yticklabels([S.um(f"{t:g}") for t in S.EKE_TREND_TICKS]); cb.ax.tick_params(labelsize=F.FS, width=0.4, length=2, pad=1.5)
cb.outline.set_linewidth(0.4); cb.set_label("EKE trend (cm$^2$ s$^{-2}$ per decade)", fontsize=F.FS, labelpad=1.5)
fig.legend(handles=S.map_legend_handles(fdr=False), loc="upper left", bbox_to_anchor=(140.5 / S.W2, (YM + 17) / H), ncol=1,
           handlelength=1.4, fontsize=F.FS, borderaxespad=0)

# ---------------------------------------------------------------- d: EKE trend by longitude, both products
T1 = pd.read_csv(OUT / "drake" / "dk12" / "dk12_band_trends_allsat.csv").sort_values("lon").rename(columns={"slope_dec": "slope"})
T4 = pd.read_csv(P4D / "p4_bins_trends.csv").sort_values("lon")
assert (T1.lon.values == T4.lon.values).all()
YD = 67.0
ax = S.axmm(fig, 16, YD, 140, 38); axd = ax.twinx()
_ab = clim.adt.interp(latitude=bat.latitude, longitude=bat.longitude).values
_band = (_ab >= LEV_PF) & (_ab <= LEV_SAF) & (bat.values < 0)
with np.errstate(all="ignore"):
    BDEP = np.nanmean(np.where(_band, -bat.values, np.nan), 0) / 1000
axd.fill_between(bat.longitude.values, BDEP, 7, color=C["land"], lw=0, zorder=0)
axd.set_ylim(6.2, -5.6); axd.set_yticks([2, 4, 6]); axd.set_ylabel("Depth (km)", color="#333333")
axd.yaxis.set_label_coords(1.04, 0.17); axd.tick_params(axis="y", colors="#333333")
axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5); axd.spines["right"].set_bounds(6.2, 1.2)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
for k, T, dx in (("P1", T1, -0.45), ("P4", T4, 0.45)):
    T = T.dropna(subset=["slope"]); x = T.lon.values + dx
    ax.vlines(x, T.lo, T.hi, color=C["eke"] if k == "P1" else "#8C8C8C", lw=0.4, alpha=0.45 if k == "P1" else 0.6, zorder=2, linestyles="solid")
    S.profile(ax, x, T.slope, T.p, C["eke"], mk=MK[k], ms=1.7, lw=0.6, ls=LS[k])
    say(f"d {k}: {len(T)} bins, positive {(T.slope > 0).sum()} ({100 * (T.slope > 0).mean():.1f} %), p<0.05 {sum(S.sigp(p) for p in T.p)}, "
        f"median {T.slope.median():.2f}, max {T.slope.max():.2f} at {T.lon.values[np.argmax(T.slope.values)]}")
r_bins = np.corrcoef(T1.slope, T4.slope)[0, 1]
assert abs(r_bins - SUM.loc["2-deg bins: correlation of P1 and P4 EKE trends", "P1"]) < 1e-9
say(f"d r(P1, P4 bins) = {r_bins:.4f}")
S.zero(ax); S.lon_axis(ax); ax.set_ylim(-300, 215); ax.set_yticks([-100, 0, 100, 200]); ax.spines["left"].set_bounds(-110, 215)
S.yticks_minus(ax); ax.set_ylabel("EKE trend\n(cm$^2$ s$^{-2}$ per decade)")
ax.yaxis.set_label_coords(-0.065, 0.62)
S.keys_on_axis(ax, y=1.0); S.flow_arrow(ax, 166, 179, 0.91)
S.letter(fig, 2, YD + 41.5, "d", f"EKE trend between the mean SAF and PF (2° bins; r = {r_bins:.2f} between products); mean seafloor depth (bottom)")


# ---------------------------------------------------------------- e, f: forest rows; g: lagged correlation
def forest(ax, rows, xl, xt, xlabel, dodge=0.19):
    ys = -np.arange(len(rows), dtype=float)
    for (lab, ent), y in zip(rows, ys):
        n = len(ent)
        for i, (col, k, r) in enumerate(ent):
            yy = y + ((n - 1) / 2 - i) * dodge
            lo, hi = max(r["lo"], xl[0]), min(r["hi"], xl[1]); assert (lo, hi) == (r["lo"], r["hi"]), (lab, r)
            ax.plot([lo, hi], [yy, yy], color=col, lw=0.8, ls=LS[k], solid_capstyle="butt", zorder=2)
            sg = S.sigp(r["p"])
            ax.plot(r["val"], yy, MK[k], color=col, mfc=col if sg else "white", mec=col, mew=0.8, ms=MS, zorder=3)
    for y in ys[:-1]:
        ax.axhline(y - 0.5, color=C["light"], lw=0.5, zorder=0)
    ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0)
    ax.set_yticks(ys); ax.set_yticklabels([l for l, _ in rows]); ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False); ax.set_ylim(ys[-1] - 0.5, ys[0] + 0.5)
    ax.set_xlim(*xl); ax.set_xticks(xt); ax.set_xlabel(xlabel); S.xticks_minus(ax)


YE = 11; HE = 40
rr = lambda key: Q.retest_row(OUT, key)
rows = []
for reg, lab in (("ACC", "ACC band"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian")):
    ent = []
    for v, col in (("KE", C["ke"]), ("EKE", C["eke"])):
        for k in ("P1", "P4"):
            r = rr(f"{Q.PROD[k]['key']}|{reg}_{v}|2001-2025"); ent.append((col, k, r))
            q = SUM.loc[f"{reg + ' band' if reg == 'ACC' else reg + ' sector'} {v} trend 2001-2025"]
            q = q[q.unit == "% of mean per decade"].iloc[0] if isinstance(q, pd.DataFrame) else q
            assert abs(r["val"] - q[k]) < 1e-9 and abs(r["lo"] - q[f"{k}_lo"]) < 1e-9 and abs(r["p"] - q[f"{k}_p"]) < 1e-9
            say(f"e {reg} {v} {k}: {r['val']:+.2f} [{r['lo']:+.2f}, {r['hi']:+.2f}] p = {r['p']:.3f}")
    rows.append((lab, ent))
ax = S.axmm(fig, 17, YE, 38, HE)
forest(ax, rows, (-4, 18), [0, 5, 10, 15], "Trend (% of mean per decade)", dodge=0.19)
S.letter(fig, 2, YE + HE + 3.5, "e", "KE and EKE trends")

TP = pd.read_csv(P4D / "p4_09_relthr_trends.csv")
TP = TP[((TP["product"] == "P1") & (TP["set"] == "relthr")) | ((TP["product"] == "P4") & (TP["set"] == "relthr_none"))]
X7 = pd.read_csv(OUT / "x7" / "x7_trends.csv")
for key in ("lat_hp", "lat_contour"):
    a = TP[(TP["product"] == "P1") & (TP.region == "circumpolar") & (TP.series == key)].iloc[0]
    b = X7[(X7.region == "circumpolar") & (X7.series == key)].iloc[0]
    assert abs(a.slope - b.slope) < 1e-3 and abs(a.lo - b.lo) < 1e-3 and abs(a.hi - b.hi) < 1e-3 and a.p == b.p
rows = []
for reg, lab in (("circumpolar", "Circumpolar"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian"), ("86-40W", "86–40° W")):
    ent = []
    for key, col in (("lat_hp", C["pos"]), ("lat_contour", C["grey"])):
        for k in ("P1", "P4"):
            r = TP[(TP["product"] == k) & (TP.region == reg) & (TP.series == key)].iloc[0]
            ent.append((col, k, dict(val=r.slope, lo=r.lo, hi=r.hi, p=r.p)))
            say(f"f {reg} {key} {k}: {r.slope:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}] p = {r.p:.3f} (n = {r.n_tracks})")
    rows.append((lab, ent))
ax = S.axmm(fig, 79, YE, 38, HE)
forest(ax, rows, (-50, 45), [-40, -20, 0, 20, 40], "Trend (km per decade, + north)", dodge=0.17)
S.letter(fig, 59, YE + HE + 3.5, "f", "Position trends")

LC = pd.read_csv(P4D / "p4_lagcorr.csv"); LS_ = pd.read_csv(P4D / "p4_lagcorr_summary.csv")
ax = S.axmm(fig, 132, YE, 36, HE)
for v, col in (("KE", C["ke"]), ("EKE", C["eke"])):
    for k in ("P1", "P4"):
        d = LC[(LC["product"] == k) & (LC.series == v)].sort_values("lag_months")
        ax.plot(d.lag_months, d.r, color=col, ls=LS[k], lw=0.9, zorder=3)
        s = LS_[(LS_["product"] == k) & (LS_.series == v)].iloc[0]
        assert abs(d.r.max() - s.r_max) < 1e-12 and int(d.lag_months.values[np.argmax(d.r.values)]) == s.lag_max
        sg = S.sigp(s.p_fw)
        ax.plot(s.lag_max, s.r_max, MK[k], color=col, mfc=col if sg else "white", mec=col, mew=0.8, ms=MS, zorder=4)
        say(f"g {v} {k}: r_max {s.r_max:.3f} at {s.lag_max} months, family-wise p {s.p_fw:.4f}, pointwise p {s.p_pt_max:.4f}")
S.zero(ax); ax.set_xlim(0, 36); ax.set_xticks([0, 12, 24, 36]); ax.set_ylim(-0.3, 0.7); ax.set_yticks([-0.2, 0, 0.2, 0.4, 0.6])
S.yticks_minus(ax); ax.set_xlabel("Lag (months, wind leading)"); ax.set_ylabel("r with wind stress")
S.letter(fig, 121, YE + HE + 3.5, "g", "Lagged correlation with wind")

F.check(fig)
S.save(fig, FIGD, "FigureS3"); plt.close(fig)
