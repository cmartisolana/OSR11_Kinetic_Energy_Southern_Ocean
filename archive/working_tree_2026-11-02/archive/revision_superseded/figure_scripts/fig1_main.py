"""Main Fig. 1 (v1 style): mean KE/EKE and their 2001-2025 trends; ACC-band series and regional trends.
(a, b) mean KE and EKE 2001-2025 (all-sat, 1/8°); (c, d) deseasonalised Theil-Sen trends on 1° cells,
dots = significant after Benjamini-Hochberg FDR (alpha 0.10); (e) ACC-band KE and EKE (all-sat: raw
monthly, 12-month running mean, trend; two-sat: 12-month running mean); (f) ACC-band along-stream
speed, both products; (g) trends by region (ACC band), both products, 95 % CI, calibrated p.
Usage: python fig1_main.py <revision/data dir> <cache dir>
"""
import sys
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

data, cache = sys.argv[1], sys.argv[2]
ROOT_OUT = Path("/mnt/user-data/uploads/repository/outputs")
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
ice = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc").ice_mask.values > 0.5
lat = np.load(f"{cache}/lat.npy"); lon = np.load(f"{cache}/lon.npy")
ke = np.load(f"{cache}/ke.npy", mmap_mode="r"); eke = np.load(f"{cache}/eke.npy", mmap_mode="r")
kem = np.nanmean(np.asarray(ke, np.float32), 0) * 1e4; ekem = np.nanmean(np.asarray(eke, np.float32), 0) * 1e4
kem[ice] = np.nan; ekem[ice] = np.nan
M = xr.open_dataset(OUT / "a1_maps_bc.nc")   # revised test (a1_maps_bc.py)
S = {p: pd.read_csv(OUT / "a5" / f"a5_series_{p}.csv") for p in ("allsat", "twosat")}
T = {p: pd.read_csv(OUT / "a5" / f"a5_trends_{p}.csv").query("period == '2001-2025'").set_index("series") for p in ("allsat", "twosat")}
T = {p: retest_override(T[p], lambda s, r, p=p: f"{p}|{s}|2001-2025") for p in T}   # revised test

proj = ccrs.SouthPolarStereo()
fig = plt.figure(figsize=(10.8, 12.6))
outer = fig.add_gridspec(2, 1, height_ratios=[2, 0.78], hspace=0.16)
gm = outer[0].subgridspec(2, 2, hspace=0.18, wspace=0.05)
gb = outer[1].subgridspec(1, 3, width_ratios=[1, 0.62, 0.85], wspace=0.38)
ax1 = fig.add_subplot(gm[0, 0], projection=proj); ax2 = fig.add_subplot(gm[0, 1], projection=proj)
ax3 = fig.add_subplot(gm[1, 0], projection=proj); ax4 = fig.add_subplot(gm[1, 1], projection=proj)
LON, LAT = np.meshgrid(lon, lat)
vmax_mean = np.nanpercentile(np.r_[kem[np.isfinite(kem)], ekem[np.isfinite(ekem)]], 99)
for ax, fld, ttl, lab in ((ax1, kem, "Mean KE (2001–2025)", "a"), (ax2, ekem, "Mean EKE (2001–2025)", "b")):
    polar_ax_style(ax, ccrs, cfeature)
    im1 = ax.pcolormesh(LON, LAT, fld, transform=ccrs.PlateCarree(), cmap="YlGnBu", vmin=0, vmax=vmax_mean,
                        shading="auto", rasterized=True)
    cs = add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
    add_ice_mask(ax, ccrs, lon, lat, ice); add_sector_boundaries(ax, ccrs)
    ax.set_title(ttl); panel_label(ax, lab, x=0.0, y=1.0)
ax2.clabel(cs, fmt="%.1f", fontsize=6, inline=True)
cb = fig.colorbar(im1, ax=[ax1, ax2], shrink=0.42, pad=0.02); cb.set_label(r"Kinetic Energy (cm$^2$ s$^{-2}$)"); cb.ax.tick_params(labelsize=8)
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F3F2EF")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
for ax, v, ttl, lab in ((ax3, "KE", "KE trend", "c"), (ax4, "EKE", "EKE trend", "d")):
    polar_ax_style(ax, ccrs, cfeature)
    im3 = ax.pcolormesh(L1, A1, M[f"{v}_slope"].values, transform=ccrs.PlateCarree(), cmap=tc, vmin=-8, vmax=8,
                        shading="nearest", rasterized=True, zorder=0)
    add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
    add_ice_mask(ax, ccrs, lon, lat, ice); add_sector_boundaries(ax, ccrs)
    jj, ii = np.nonzero(M[f"{v}_fdr"].values == 1)
    ax.scatter(M.lon.values[ii], M.lat.values[jj], s=7, c="#1F2933", lw=0, transform=ccrs.PlateCarree(), zorder=6)
    ax.set_title(ttl); panel_label(ax, lab, x=0.0, y=1.0)
cb = fig.colorbar(im3, ax=[ax3, ax4], shrink=0.42, pad=0.02, extend="both")
cb.set_label(r"Trend (cm$^2$ s$^{-2}$ yr$^{-1}$)"); cb.ax.tick_params(labelsize=8)

# (e) ACC-band KE and EKE
ax5 = fig.add_subplot(gb[0])
x = (S["allsat"].year + (S["allsat"].month - 0.5) / 12).values
for v, col in (("KE", STYLE["ke"]), ("EKE", STYLE["eke"])):
    y = S["allsat"][f"ACC_{v}"].values; b = T["allsat"].loc[f"ACC_{v}", "slope"]
    ts_panel(ax5, x, y, col, label=v, slope=b, intercept=theil_line(x, deseason(S["allsat"], f"ACC_{v}"), b))
for i, v in enumerate(("KE", "EKE")):
    a, t = T["allsat"].loc[f"ACC_{v}"], T["twosat"].loc[f"ACC_{v}"]
    ax5.text(0.98, 0.15 - 0.08 * i,
             f"{v}: {a.slope:+.2f} cm$^2$ s$^{{-2}}$ yr$^{{-1}}$ {stars(a.p_cal)} (two-sat {t.slope:+.2f})", ha="right",
             transform=ax5.transAxes, va="top", fontsize=8, color=STYLE["ke" if v == "KE" else "eke"])
ax5.set_xlabel("Year"); ax5.set_ylabel(r"Kinetic Energy (cm$^2$ s$^{-2}$)")
base_timeseries_style(ax5); ax5.legend(loc="upper left", ncol=2, fontsize=8, frameon=False); panel_label(ax5, "e", y=1.12)
yl = ax5.get_ylim(); ax5.set_ylim(yl[0] - 0.15 * (yl[1] - yl[0]), yl[1] + 0.08 * (yl[1] - yl[0]))

# (f) along-stream speed
ax6 = fig.add_subplot(gb[1])
for p, col, ls in (("allsat", STYLE["trend_sig"], "-"), ("twosat", STYLE["neutral"], (0, (1.5, 1.2)))):
    y = S[p]["ACC_SPD"].values; b = T[p].loc["ACC_SPD", "slope"]
    if p == "allsat":
        ts_panel(ax6, x, y, col, label="all-sat", slope=b, intercept=theil_line(x, deseason(S[p], "ACC_SPD"), b))
    else:
        ax6.plot(x, running(y), color=col, lw=1.2, ls=ls, label="two-sat")
a, t = T["allsat"].loc["ACC_SPD"], T["twosat"].loc["ACC_SPD"]
ax6.text(0.97, 0.05, f"{a.pct_per_decade:+.1f} % dec$^{{-1}}$ {stars(a.p_cal)}\n(two-sat {t.pct_per_decade:+.1f} %)",
         transform=ax6.transAxes, fontsize=8, color=STYLE["trend_sig"], ha="right", bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
ax6.set_xlabel("Year"); ax6.set_ylabel(r"Along-stream speed (cm s$^{-1}$)")
base_timeseries_style(ax6); ax6.legend(loc="upper left", fontsize=8, frameon=False); panel_label(ax6, "f", y=1.12)

# (g) regional bars
ax7 = fig.add_subplot(gb[2])
xpos = np.arange(len(REGION_ORDER)); w = 0.2
for k, (v, p, col, hatch) in enumerate((("KE", "allsat", STYLE["ke"], None), ("KE", "twosat", STYLE["ke"], "////"),
                                         ("EKE", "allsat", STYLE["eke"], None), ("EKE", "twosat", STYLE["eke"], "////"))):
    rows = [T[p].loc[f"{r}_{v}"] for r in REGION_ORDER]
    s_ = np.array([r_.slope for r_ in rows]); lo = np.array([r_.ci_cal_lo for r_ in rows]); hi = np.array([r_.ci_cal_hi for r_ in rows])
    pv = np.array([r_.p_cal for r_ in rows])
    xb = xpos + (k - 1.5) * w
    for xx, ss, l_, h_, pp in zip(xb, s_, lo, hi, pv):
        sig = pp < 0.05
        if hatch is None:
            ax7.bar(xx, ss, width=w, color=col if sig else STYLE["ns"], edgecolor="white", linewidth=0.6, zorder=3)
        else:
            ax7.bar(xx, ss, width=w, color="white", edgecolor=col if sig else "#9AA0A6", hatch=hatch, linewidth=0.6, zorder=3)
        ax7.errorbar(xx, ss, yerr=[[ss - l_], [h_ - ss]], color="#1F2933", lw=0.9, capsize=2.0, capthick=0.9, zorder=4)
        ax7.text(xx, h_ + 0.05, stars(pp), ha="center", va="bottom", fontsize=7)
for v, p, col, hatch in (("KE", "all-sat", STYLE["ke"], None), ("KE", "two-sat", STYLE["ke"], "////"),
                         ("EKE", "all-sat", STYLE["eke"], None), ("EKE", "two-sat", STYLE["eke"], "////")):
    ax7.bar([np.nan], [np.nan], color=col if hatch is None else "white", edgecolor="white" if hatch is None else col,
            hatch=hatch, label=f"{v} {p}")
ax7.axhline(0, color="#495057", lw=0.8)
ax7.set_xticks(xpos); ax7.set_xticklabels(REGION_ORDER, rotation=25, ha="right")
ax7.set_ylabel(r"Trend (cm$^2$ s$^{-2}$ yr$^{-1}$)"); base_timeseries_style(ax7)
ax7.legend(loc="upper left", fontsize=7, frameon=False, ncol=2); panel_label(ax7, "g", y=1.12)
yl = ax7.get_ylim(); ax7.set_ylim(min(0, yl[0]), yl[1] * 1.3)
fig.savefig(FIG / "Fig1_KE_EKE_trends.png"); fig.savefig(FIG / "Fig1_KE_EKE_trends.pdf")
print("saved")
