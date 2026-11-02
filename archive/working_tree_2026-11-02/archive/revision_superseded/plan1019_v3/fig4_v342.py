"""Figure 4 (v3.42): eddies, mean flow and transport of the ACC, 2001-2025 (altimetry, Product 1) and 2004-2025 (Argo, RG).
a  Trend of the 0-1975 dbar depth-mean eastward geostrophic velocity relative to 1975 dbar (Argo; dk42), with the mean SAF and PF.
b  Annual anomalies (% of the mean): wind stress over the ACC band (bars), EKE and KE (altimetry).
c  Annual anomalies: surface transport proxy (altimetry, int u dy over the fixed band) and baroclinic transport across the band
   (Argo, 0-1975 dbar relative to 1975 dbar), with the reported trends.
d  Trends by region (% of the mean per decade, 95 % intervals; open: p >= 0.05): wind stress, EKE, KE, surface transport proxy,
   baroclinic transport.
e  Trends across streamlines (bins of the 2001-2025 mean ADT): EKE in fixed bins (shading 95 %) and in bins following the monthly
   streamlines, surface transport proxy per bin (dk18) and Argo depth-mean baroclinic velocity per 0.1 m bin (dk42; bars 95 %).
Usage: python fig4_v342.py <revision dir> <drake out dir> <work out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, DO, WO, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
OUT = REV / "outputs"; ROOT_OUT = REV.parent / "outputs"
CA = "#8E44AD"      # Argo (baroclinic) colour
FEAT = [(75, -47.5, "Kerguelen\nPlateau"), (122, -46.0, "Southeast\nIndian Ridge"), (166, -49.0, "Campbell\nPlateau"),
        (-140, -60.5, "Pacific–Antarctic\nRidge"), (-66, -60.0, "Drake\nPassage")]
def sigp(p): return np.isfinite(p) and round(float(p), 3) < 0.05
def fmt(v, lo, hi): return f"{v:+.1f} ({lo:+.1f} to {hi:+.1f})".replace("-", "−")
# ---------------------------------------------------------------- data
clim = xr.open_dataset(OUT / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(Z12["lev_saf"]), float(Z12["lev_pf"])
ice_ds = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc"); bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation; bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
M42 = np.load(WO / "dk42" / "dk42_maps.npz"); T42 = pd.read_csv(WO / "dk42" / "dk42_trends.csv"); S42 = pd.read_csv(WO / "dk42" / "dk42_series.csv"); B42 = pd.read_csv(WO / "dk42" / "dk42_stream_bins.csv")
A5 = pd.read_csv(OUT / "a5" / "a5_series_allsat.csv"); A8 = pd.read_csv(OUT / "a8" / "a8_series_allsat.csv"); A8T = pd.read_csv(OUT / "a8" / "a8_trends_allsat.csv").set_index("series")
RT5 = pd.read_csv(OUT / "retest" / "retest_A5.csv").set_index("key")
WS = pd.read_csv(DO / "dk17" / "s1_wind_series.csv"); WT = pd.read_csv(DO / "dk17" / "s1_wind_trends.csv").set_index("series")
P5 = pd.read_csv(DO / "dk17" / "s5_profiles.csv"); D18 = pd.read_csv(DO / "dk18" / "dk18_trends.csv")
def ann(df, col):
    s = df.groupby("year")[col].mean(); return 100 * (s / df[col].mean() - 1)
def rt5(reg, var):
    r = RT5.loc[f"allsat|{reg}_{var}|2001-2025"]; mu = A5[f"{reg}_{var}"].mean()
    return dict(pct=r.pct_per_decade, lo=1000 * r.ci_lo / mu, hi=1000 * r.ci_hi / mu, p=r.p_new)
def t42(reg, var):
    r = T42[(T42.region == reg) & (T42["var"] == var)].iloc[0]; return dict(pct=r.pct, lo=r.pct_lo, hi=r.pct_hi, p=r.p)
def trend_line(ax, pct, p, col, x0, x1, y0=0.0, lw=0.9):
    xx = np.array([x0, x1]); ax.plot(xx, y0 + pct / 10 * (xx - (x0 + x1) / 2), color=col, lw=lw, ls="-" if sigp(p) else (0, (3, 1.6)), zorder=3)
# ---------------------------------------------------------------- layout
fig = plt.figure(figsize=(S.W2, 0.84 * S.W2))
gA = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.78], hspace=0.38)
g1 = gA[0].subgridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.12); gS = g1[1].subgridspec(2, 1, hspace=0.42)
g2 = gA[1].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.3)
# a: Argo depth-mean baroclinic velocity trend
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F2F2F2")
L1, A1 = np.meshgrid(M42["lon"], M42["lat"]); F0 = M42["u_mean_trend"]
# 3 x 3 running mean (cells with data only; longitude periodic) to show the pattern rather than cell noise
Fp = np.concatenate([F0[:, -1:], F0, F0[:, :1]], 1); Fp = np.concatenate([np.full((1, Fp.shape[1]), np.nan), Fp, np.full((1, Fp.shape[1]), np.nan)], 0)
F = np.nanmean(np.stack([Fp[i:i + F0.shape[0], j:j + F0.shape[1]] for i in range(3) for j in range(3)]), 0); F[~np.isfinite(F0)] = np.nan
im = ax.pcolormesh(L1, A1, np.ma.masked_invalid(F), transform=PC, cmap=tc, vmin=-0.3, vmax=0.3, shading="nearest", rasterized=True, zorder=0)
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[LEV_PF, LEV_SAF], colors="#1A1A1A", linewidths=0.6, linestyles="solid", transform=PC, zorder=6)
add_ice_mask(ax, ccrs, ice_ds.longitude.values, ice_ds.latitude.values, ice_ds.ice_mask.values > 0.5)
jj, ii = np.nonzero((M42["u_mean_p"] < 0.05) & M42["band"])
ax.scatter(M42["lon"][ii], M42["lat"][jj], s=1.6, c="k", lw=0, transform=PC, zorder=7)
for lo_, la_, s_ in FEAT:
    ax.text(lo_, la_, s_, transform=PC, fontsize=4.6, ha="center", va="center", color="#222222", zorder=9, linespacing=0.9, bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))
S.plabel(ax, "a", "Baroclinic velocity trend, 2004–2025 (Argo)")
S.hcbar(fig, im, ax, "0–1975 dbar mean eastward velocity trend (cm s$^{-1}$ per decade)", shrink=0.72, extend="both")
# b: energy and wind anomalies
ax = fig.add_subplot(gS[0]); S.zero(ax)
wa = ann(WS, "tau_accmask_ACC"); ax.bar(wa.index + 0.5, wa.values, 0.8, color="#D9D9D9", zorder=1, label="wind stress")
for col, cc, lab in (("ACC_KE", C["ke"], "KE"), ("ACC_EKE", C["eke"], "EKE")):
    s = ann(A5, col); ax.plot(s.index + 0.5, s.values, "-o", color=cc, lw=0.7, ms=1.8, zorder=4, label=lab)
    r = rt5("ACC", col.split("_")[1]); trend_line(ax, r["pct"], r["p"], cc, 2001, 2026)
rw = WT.loc["tau_accmask_ACC"]; trend_line(ax, rw.pct_dec, rw.p, "#7F7F7F", 2001, 2026)
S.ts_style(ax, (2000.8, 2026.2)); ax.set_ylim(-22, 30); ax.set_ylabel("Anomaly (% of mean)")
ax.legend(handles=[Patch(facecolor="#D9D9D9", label=f"wind stress, ACC band: {fmt(rw.pct_dec, rw.pct_dec_lo, rw.pct_dec_hi)} % per decade"),
                   Line2D([], [], color=C["eke"], marker="o", ms=1.8, lw=0.7, label="EKE: " + fmt(**{k: rt5('ACC', 'EKE')[k2] for k, k2 in (('v', 'pct'), ('lo', 'lo'), ('hi', 'hi'))})),
                   Line2D([], [], color=C["ke"], marker="o", ms=1.8, lw=0.7, label="KE: " + fmt(**{k: rt5('ACC', 'KE')[k2] for k, k2 in (('v', 'pct'), ('lo', 'lo'), ('hi', 'hi'))}))],
          loc="upper left", fontsize=4.8, ncol=1, handlelength=1.6)
S.plabel(ax, "b", "Eddy and total kinetic energy of the ACC band")
# c: transports
ax = fig.add_subplot(gS[1]); S.zero(ax)
s = ann(A8, "ACC_T"); ax.plot(s.index + 0.5, s.values, "-o", color=C["tr"], lw=0.7, ms=1.8, zorder=4)
r8 = A8T.loc["ACC_T"]; trend_line(ax, r8.pct_per_decade, r8.p, C["tr"], 2001, 2026)
yb = S42["ACC band|T_bc"]; sb = 100 * (yb / yb.mean() - 1); ax.plot(S42.year + 0.5, sb, "-s", color=CA, lw=0.7, ms=1.9, zorder=5)
rb = t42("ACC band", "T_bc"); trend_line(ax, rb["pct"], rb["p"], CA, 2004, 2026)
S.ts_style(ax, (2000.8, 2026.2)); ax.set_ylim(-4.5, 5.5); ax.set_ylabel("Anomaly (% of mean)")
ax.legend(handles=[Line2D([], [], color=C["tr"], marker="o", ms=1.8, lw=0.7, label=f"surface transport proxy (altimetry): {fmt(r8.pct_per_decade, r8.pct_lo, r8.pct_hi)} % per decade"),
                   Line2D([], [], color=CA, marker="s", ms=1.9, lw=0.7, label=f"baroclinic transport, 0–1975 dbar (Argo): {fmt(rb['pct'], rb['lo'], rb['hi'])} % per decade")],
          loc="upper left", fontsize=4.8, handlelength=1.6)
S.plabel(ax, "c", "Transport across the ACC band")
# d: trends by region
ax = fig.add_subplot(g2[0]); S.zero(ax)
REGS = [("ACC band", "ACC", "ACC band"), ("Pacific", "Pacific", "Pacific"), ("Indian", "Indian", "Indian"), ("Atlantic", "Atlantic", "Atlantic")]
SER = [("wind", "#9A9A9A", "wind stress"), ("EKE", C["eke"], "EKE"), ("KE", C["ke"], "KE"), ("T", C["tr"], "surface transport proxy"), ("Tbc", CA, "baroclinic transport (Argo)")]
w_ = 0.16
for gi, (lab, k5, k42) in enumerate(REGS):
    rw = WT.loc[f"tau_accmask_{k5}"]; vals = {"wind": dict(pct=rw.pct_dec, lo=rw.pct_dec_lo, hi=rw.pct_dec_hi, p=rw.p), "EKE": rt5(k5, "EKE"), "KE": rt5(k5, "KE")}
    r8 = A8T.loc[f"{k5}_T"]; vals["T"] = dict(pct=r8.pct_per_decade, lo=r8.pct_lo, hi=r8.pct_hi, p=r8.p); vals["Tbc"] = t42(k42, "T_bc")
    for k, (key, col, _) in enumerate(SER):
        r = vals[key]; x = gi + (k - 2) * w_
        ax.bar(x, r["pct"], w_ * 0.9, color=col if sigp(r["p"]) else "white", edgecolor=col, lw=0.8, zorder=3); ax.plot([x, x], [r["lo"], r["hi"]], color="#333333", lw=0.5, zorder=4)
ax.set_xticks(range(4)); ax.set_xticklabels([r[0] for r in REGS]); ax.tick_params(axis="x", length=0); ax.set_ylim(-3, 17); ax.set_ylabel("Trend (% of mean per decade)")
ax.legend(handles=[Patch(facecolor=c, edgecolor=c, label=l) for _, c, l in SER] + [Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper right", fontsize=4.6, ncol=2)
S.plabel(ax, "d", "Trends by region")
# e: across streamlines
ax = fig.add_subplot(g2[1]); S.zero(ax)
for lev, nm in ((LEV_SAF, "SAF"), (LEV_PF, "PF")):
    ax.axvline(lev, color="#BBBBBB", lw=0.6, zorder=0); ax.text(lev, 1.0, nm, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5, color="#555555")
pf = P5[(P5.bin == "madt") & (P5.sector == "ACC") & (P5["var"] == "EKE")].sort_values("mid"); pi_ = P5[(P5.bin == "iadt") & (P5.sector == "ACC") & (P5["var"] == "EKE")].sort_values("mid")
ax.fill_between(pf.mid, pf.pct_dec_lo, pf.pct_dec_hi, color=C["eke"], alpha=0.10, lw=0, zorder=1)
ax.plot(pf.mid, pf.pct_dec, "-o", color=C["eke"], lw=0.8, ms=1.8, zorder=4, mfc=np.where(pf.p.round(3) < 0.05, C["eke"], "white").tolist()[0])
for x_, y_, p_ in zip(pf.mid, pf.pct_dec, pf.p): ax.plot(x_, y_, "o", ms=1.8, color=C["eke"], mfc=C["eke"] if p_ < 0.05 else "white", mew=0.6, zorder=5)
ax.plot(pi_.mid, pi_.pct_dec, "-", color=C["eke"], lw=0.7, ls=(0, (2.5, 1.2)), zorder=4)
pk = P5[(P5.bin == "madt") & (P5.sector == "ACC") & (P5["var"] == "KE")].sort_values("mid")
ax.plot(pk.mid, pk.pct_dec, "-", color=C["ke"], lw=0.7, zorder=3)
d18 = D18[D18.key.str.startswith("bin|")].copy(); d18["mid"] = d18.key.str.split("|").str[1].astype(float); d18 = d18.sort_values("mid")
ax.plot(d18.mid, d18.pct, "-", color=C["tr"], lw=0.8, zorder=4)
for x_, y_, p_ in zip(d18.mid, d18.pct, d18.p): ax.plot(x_, y_, "o", ms=1.6, color=C["tr"], mfc=C["tr"] if p_ < 0.05 else "white", mew=0.6, zorder=5)
bb = B42[B42["var"] == "u_mean"].sort_values("mid")
ax.vlines(bb.mid, bb.pct_lo, bb.pct_hi, color=CA, lw=0.6, zorder=3); ax.plot(bb.mid, bb.pct, "-", color=CA, lw=0.8, zorder=4)
for x_, y_, p_ in zip(bb.mid, bb.pct, bb.p): ax.plot(x_, y_, "s", ms=2.0, color=CA, mfc=CA if p_ < 0.05 else "white", mew=0.6, zorder=5)
ax.set_xlim(-0.92, 0.52); ax.set_ylim(-5, 19); ax.set_xlabel("Mean ADT of the bin (m)   south ←  → north"); ax.set_ylabel("Trend (% of mean per decade)")
ax.legend(handles=[Line2D([], [], color=C["eke"], marker="o", ms=1.8, lw=0.8, label="EKE, fixed bins (shading: 95 %)"), Line2D([], [], color=C["eke"], lw=0.7, ls=(0, (2.5, 1.2)), label="EKE, bins following the streamlines"),
                   Line2D([], [], color=C["ke"], lw=0.7, label="KE"), Line2D([], [], color=C["tr"], marker="o", ms=1.6, lw=0.8, label="surface transport proxy"),
                   Line2D([], [], color=CA, marker="s", ms=2.0, lw=0.8, label="baroclinic velocity, 0–1975 dbar (Argo; bars: 95 %)"), Line2D([], [], color="#555555", marker="o", mfc="white", ls="", ms=2, label="open: p ≥ 0.05")],
          loc="upper left", fontsize=4.5, ncol=2, handlelength=1.6)
S.plabel(ax, "e", "Trends across streamlines")
S.save(fig, FIGD, "Figure4")
