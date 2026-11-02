"""Fig. 6 (v1 style): standing meanders and EKE, 2001-2025 (all-sat).
Left (a, c, e): mean EKE (log scale) with the mean jet-core path (black; core-level contour of ADT relative
  to its regional mean, averaged over the three epochs of Fig. 3), the meander section (dashed) and the
  downstream EKE box (solid); grey: 3000 m isobath.
Right (b, d, f): downstream EKE response to the path length. Difference in downstream EKE anomaly between
  months whose path length was in the upper and in the lower tercile, as a function of lag (EKE after the
  path length for positive lags), in % of the mean downstream EKE. Series prepared as in b1_covar.py
  (deseasonalised, Theil-Sen detrended, 2-month running mean). Grey band: 95 % range of the same composite
  with phase-randomised path-length series (1000 surrogates). Shading: pre-specified lags 0-3 months.
  p: pre-specified test (b1_covar.py, maximum over lags 0-3, surrogates).
Bottom (g): lag correlation between path length and downstream EKE for each meander and their mean, with
  the 95 % surrogate range of the mean; pooled p of the pre-specified test.
Usage: OSR11_CACHE=<cache> python fig6_meander_eke.py <revision/data dir>
"""
import os, sys
import numpy as np, pandas as pd, xarray as xr
import contourpy
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.colors import LogNorm
from figstyle import *  # noqa: F401,F403
import b1_covar as BC

data = sys.argv[1]
CACHE = Path(os.environ.get("OSR11_CACHE", "/home/claude/w/cache"))
MEAN = [("PAR", "Pacific–Antarctic Ridge", (-156, -74, -64, -44)),
        ("SEIR150", "Southeast Indian Ridge (150°E)", (122, 166, -64, -44)),
        ("KERG", "Kerguelen Plateau (Fawn Trough)", (64, 100, -62, -42))]
BOX = {"PAR": dict(mea=(-135.13, -109.38), mlat=(-62.0, -45.0), down=(-109.38, -80.0), band=(-60.0, -48.0)),
       "SEIR150": dict(mea=(140.0, 150.0), mlat=(-60.0, -46.0), down=(150.0, 158.0), band=(-60.0, -47.0)),
       "KERG": dict(mea=(76.0, 84.0), mlat=(-58.0, -48.0), down=(84.0, 94.0), band=(-57.0, -48.0))}   # = b1_meander_run.py / b1_pilot.py
COL = {"PAR": "#E07A5F", "SEIR150": "#2A9D8F", "KERG": "#457B9D"}
T = pd.read_csv(OUT / "b1" / "b1_covar_KERG_PAR_SEIR150.csv")
P = pd.read_csv(OUT / "b1" / "b1_covar_pooled_KERG_PAR_SEIR150.csv")
p_pool = P[P.link == "L2_length->down_eke"].p_surr_maxlag.iloc[0]
LAGS = np.arange(-6, 13)
NS = 1000
rng = np.random.default_rng(20261013)
pc_ = ccrs.PlateCarree()

lat = np.load(CACHE / "lat.npy"); lon = np.load(CACHE / "lon.npy")
eke = np.load(CACHE / "eke.npy", mmap_mode="r")
ice = xr.open_dataset("/mnt/user-data/uploads/repository/outputs/sea_ice_mask.nc").ice_mask.values > 0.5


def lagged(x, y, L):
    """pairs (x_t, y_{t+L})"""
    if L >= 0:
        return x[:len(x) - L], y[L:]
    return x[-L:], y[:len(y) + L]


def composite(x, y):
    lo, hi = np.nanquantile(x, [1 / 3, 2 / 3])
    out = []
    for L in LAGS:
        xx, yy = lagged(x, y, L)
        ok = np.isfinite(xx) & np.isfinite(yy)
        out.append(np.mean(yy[ok & (xx >= hi)]) - np.mean(yy[ok & (xx <= lo)]))
    return np.array(out)


def lagr(x, y):
    out = []
    for L in LAGS:
        xx, yy = lagged(x, y, L); ok = np.isfinite(xx) & np.isfinite(yy)
        out.append(np.corrcoef(xx[ok], yy[ok])[0, 1])
    return np.array(out)


fig = plt.figure(figsize=(11.0, 12.6))
gs = fig.add_gridspec(5, 2, width_ratios=[1.75, 1], height_ratios=[1, 1, 1, 0.16, 0.95], hspace=0.42, wspace=0.18)
norm = LogNorm(vmin=20, vmax=2000)
R_all, Rs_all, R_mea = {}, {}, {}
for r, (name, title, (lo0, lo1, la0, la1)) in enumerate(MEAN):
    ax = fig.add_subplot(gs[r, 0], projection=pc_)
    ax.set_extent([lo0, lo1, la0, la1], crs=pc_); ax.set_aspect("auto")
    j0, j1 = np.searchsorted(lat, [la0, la1]); i0, i1 = np.searchsorted(lon, [lo0, lo1])
    E = np.nanmean(np.asarray(eke[:, j0:j1, i0:i1], np.float32), 0) * 1e4
    E[ice[j0:j1, i0:i1]] = np.nan
    im = ax.pcolormesh(lon[i0:i1], lat[j0:j1], E, cmap="YlGnBu", norm=norm, shading="nearest", rasterized=True, transform=pc_, zorder=0)
    bat = xr.open_dataset(f"{data}/bathymetry_so.nc").elevation.sel(longitude=slice(lo0, lo1), latitude=slice(la0, la1))
    ax.contour(bat.longitude, bat.latitude, -bat, levels=[3000], colors="#6b7280", linewidths=0.5, transform=pc_, zorder=1)
    ax.add_feature(cfeature.LAND, facecolor="#000000", edgecolor="#555555", linewidth=0.2, zorder=4)
    O = xr.open_dataset(OUT / "b1" / f"occ_{name}.nc")
    gen = contourpy.contour_generator(O.longitude.values, O.latitude.values,
                                      np.ma.masked_invalid(O.epoch_rel_adt.values.mean(0)), name="serial")
    for k, Lp in enumerate([Lp for Lp in gen.lines(O.attrs["level"]) if np.ptp(Lp[:, 0]) > 2]):
        ax.plot(Lp[:, 0], Lp[:, 1], color="#111111", lw=1.6, transform=pc_, zorder=3, label="mean jet-core path" if k == 0 else None)
    b = BOX[name]
    m0, m1 = b["mea"]; ml0, ml1 = b["mlat"]
    ax.plot([m0, m1, m1, m0, m0], [ml0, ml0, ml1, ml1, ml0], color="#111111", lw=1.0, ls=(0, (4, 2.5)), transform=pc_, zorder=5,
            label="meander section")
    d0, d1 = b["down"]; bl0, bl1 = b["band"]
    ax.plot([d0, d1, d1, d0, d0], [bl0, bl0, bl1, bl1, bl0], color="#C0392B", lw=1.6, transform=pc_, zorder=5,
            label="downstream EKE box")
    gl = ax.gridlines(crs=pc_, draw_labels=True, linewidth=0.3, color="#7A8288", alpha=0.45, linestyle=":")
    gl.top_labels = False; gl.right_labels = False
    gl.xlabel_style = {"size": 8, "color": "#4A4A4A"}; gl.ylabel_style = {"size": 8, "color": "#4A4A4A"}
    ax.set_title(title, loc="left"); panel_label(ax, "ace"[r], x=-0.09, y=1.13)
    if r == 0:
        ax.legend(loc="lower left", fontsize=7.5, frameon=True, framealpha=0.9, edgecolor="#cccccc")

    # composite
    S = pd.read_csv(OUT / "b1" / f"b1_series_{name}.csv")
    x = BC.prep(S, "core_s150_length_km"); y = BC.prep(S, "down_eke") / S.down_eke.mean() * 100
    ym = BC.prep(S, "mea_eke") / S.mea_eke.mean() * 100
    c = composite(x, y); cm = composite(x, ym)
    R_mea[name] = lagr(x, ym)
    cs = np.array([composite(BC.phase_surr(x, rng), y) for _ in range(NS)])
    R_all[name] = lagr(x, y)
    Rs_all[name] = np.array([lagr(BC.phase_surr(x, rng), y) for _ in range(NS)])
    a2 = fig.add_subplot(gs[r, 1])
    a2.axvspan(-0.3, 3.3, color="#E9ECEF", lw=0, zorder=0)
    a2.fill_between(LAGS, np.percentile(cs, 2.5, 0), np.percentile(cs, 97.5, 0), color="#BFC5CC", alpha=0.6, lw=0, zorder=1,
                    label="95 % surrogate range")
    a2.axhline(0, color="#495057", lw=0.8)
    a2.plot(LAGS, cm, "--", color=COL[name], lw=1.3, zorder=2, label="EKE in meander section")
    a2.plot(LAGS, c, "-o", ms=3.5, color=COL[name], lw=1.8, zorder=3, label="EKE downstream")
    pm = T[(T.meander == name) & (T.link == "L2_length->down_eke")].p_surr_maxlag.iloc[0]
    a2.text(0.03, 0.95, f"lags 0–3: p={pm:.3f} {stars(pm)}", transform=a2.transAxes, va="top", fontsize=STYLE["annot_size"],
            color=STYLE["trend_sig"])
    a2.set_ylabel("Δ EKE, upper − lower tercile (%)", fontsize=9); a2.set_xticks([-6, -3, 0, 3, 6, 9, 12])
    base_timeseries_style(a2); panel_label(a2, "bdf"[r], x=-0.24, y=1.13)
    if r == 0:
        a2.set_title("Path length → EKE", loc="left", fontsize=9.5)
        a2.legend(loc="upper right", fontsize=6.5, frameon=False)
    if r == 2:
        a2.set_xlabel("Lag (months; + = EKE after path length)", fontsize=9)
host = fig.add_subplot(gs[3, 0]); host.set_axis_off()
cax = host.inset_axes([0.1, 0.62, 0.8, 0.25])
cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label(r"Mean EKE 2001–2025 (cm$^2$ s$^{-2}$)", fontsize=9)
cb.ax.tick_params(labelsize=8)

# pooled lag correlation
ax = fig.add_subplot(gs[4, :])
ax.axvspan(-0.3, 3.3, color="#E9ECEF", lw=0, zorder=0)
Rm = np.mean([R_all[n] for n in R_all], 0)
Rsm = np.mean([Rs_all[n] for n in Rs_all], 0)
ax.fill_between(LAGS, np.percentile(Rsm, 2.5, 0), np.percentile(Rsm, 97.5, 0), color="#BFC5CC", alpha=0.6, lw=0, zorder=1,
                label="95 % surrogate range (mean)")
for n in R_all:
    ax.plot(LAGS, R_all[n], "-", color=COL[n], lw=1.1, alpha=0.9, zorder=2, label=n.replace("SEIR150", "SEIR 150°E").replace("KERG", "Kerguelen"))
ax.plot(LAGS, Rm, "-o", ms=4, color="#1F2933", lw=2.2, zorder=3, label="mean of the three (downstream EKE)")
ax.plot(LAGS, np.mean([R_mea[n] for n in R_mea], 0), "--", color="#1F2933", lw=1.6, zorder=3, label="mean of the three (meander-section EKE)")
ax.axhline(0, color="#495057", lw=0.8)
ax.text(0.01, 0.95, f"Pre-specified pooled test (lags 0–3, surrogates): p={p_pool:.3f} {stars(p_pool)}", transform=ax.transAxes,
        va="top", fontsize=STYLE["annot_size"], color=STYLE["trend_sig"])
ax.set_xlabel("Lag (months; + = EKE after path length)"); ax.set_ylabel("r (path length, EKE)")
ax.set_xticks(LAGS); base_timeseries_style(ax); ax.legend(loc="lower right", fontsize=7, frameon=False, ncol=3)
yl = ax.get_ylim(); ax.set_ylim(yl[0] - 0.12, yl[1] + 0.05)
panel_label(ax, "g", x=-0.05, y=1.12)
fig.savefig(FIG / "Fig6_meander_EKE.png"); fig.savefig(FIG / "Fig6_meander_EKE.pdf")
print("saved; mean r at lags 0-3:", np.round(Rm[6:10], 3))
