"""Figure 4 of v2.1 (python fig_eddies.py main) and SI eddy figure (v2 style): tracked mesoscale eddies (AVISO META4.0 DT) in the ACC band, 2001-2023.
(a) eddy births per year per 1° cell (first observation of each track), ACC contours, wake boxes (solid);
(b) ACC-band monthly series, 12-month running means, as % of the 2001-2023 mean: eddies present, mean
    amplitude, sum of squared amplitudes (eddy-energy proxy) and gridded EKE;
(c) trends 2001-01..2023-08 in % per decade with 95 % CI (revised test).
Usage: python fig_eddies.py
"""
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

FIGV3 = OUT / "figures_v21"
MAIN = __import__("sys").argv[1:] == ["main"]
D = xr.open_dataset(OUT / "eddies" / "e1_density.nc")
S = pd.read_csv(OUT / "eddies" / "e2_series.csv"); T = pd.read_csv(OUT / "eddies" / "e2_e3_trends.csv").set_index("series")
adt = xr.open_dataset("/mnt/user-data/uploads/repository/outputs/mean_adt.nc").mean_adt
BOXES = {"Kerguelen": (75, 95), "SEIR–Macquarie": (145, 175), "PAR": (-150, -125), "Drake": (-55, -35)}

if MAIN:
    fig = plt.figure(figsize=(10.8, 4.0)); gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1], wspace=0.3)
else:
    fig = plt.figure(figsize=(10.8, 9.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], width_ratios=[1.6, 1], hspace=0.3, wspace=0.3)
ax = fig.add_subplot(gs[0, :], projection=ccrs.SouthPolarStereo()) if not MAIN else None
if ax is not None:
  polar_ax_style(ax, ccrs, cfeature)
if ax is not None:
    L, A = np.meshgrid(D.lon.values, D.lat.values)
    im = ax.pcolormesh(L, A, D.births_per_year.values, transform=ccrs.PlateCarree(), cmap="YlGnBu", vmin=0, vmax=6, shading="nearest", rasterized=True)
    add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values); add_sector_boundaries(ax, ccrs)
    for name, (lo0, lo1) in BOXES.items():
        lo = np.r_[np.linspace(lo0, lo1, 40), np.full(20, lo1), np.linspace(lo1, lo0, 40), np.full(20, lo0)]
        la = np.r_[np.full(40, -64), np.linspace(-64, -42, 20), np.full(40, -42), np.linspace(-42, -64, 20)]
        ax.plot(lo, la, color="#111111", lw=1.3, transform=ccrs.PlateCarree(), zorder=7)
        ax.text(0.5 * (lo0 + lo1), -38.8, name, transform=ccrs.PlateCarree(), ha="center", va="center", fontsize=7.5, zorder=8,
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.8))
    cb = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.03); cb.set_label("Eddy births (tracks yr$^{-1}$ per 1° cell)", fontsize=9); cb.ax.tick_params(labelsize=8)
    ax.set_title("a) Eddy generation sites (META4.0 DT, 2001–2023)", fontsize=12, fontweight="bold")
    
ax = fig.add_subplot(gs[0, 0] if MAIN else gs[1, 0])
x = S.x.values
for col, lab, c in (("ACC_present_per_day", "eddies present", "#457B9D"), ("ACC_amp", "mean amplitude", "#E07A5F"),
                    ("ACC_amp2_per_day", "Σ amplitude² (eddy-energy proxy)", "#2A9D8F"), ("ACC_EKE_allsat", "gridded EKE", "#1F2933")):
    y = 100 * (S[col].values / np.nanmean(S[col].values) - 1)
    ax.plot(x, running(y), color=c, lw=1.8 if col != "ACC_EKE_allsat" else 2.2, ls="-" if col != "ACC_EKE_allsat" else (0, (3, 1.5)), label=lab)
ax.axhline(0, color="#6b7280", lw=0.8, ls="--"); ax.set_xlabel("Year"); ax.set_ylabel("Anomaly (% of 2001–2023 mean)")
base_timeseries_style(ax); ax.legend(loc="upper left", fontsize=7.5, frameon=False, ncol=2); ax.set_xlim(2000.5, 2024)
yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.3 * (yl[1] - yl[0]))
ax.set_title(("a" if MAIN else "b") + ") ACC band, 12-month running means", loc="left", fontsize=12, fontweight="bold")

ax = fig.add_subplot(gs[0, 1] if MAIN else gs[1, 1])
def row(k):
    r = T.loc[k]; f = 1000 / r["mean"]
    return dict(pct=r.pct_per_decade, pct_lo=r.ci_lo * f, pct_hi=r.ci_hi * f, p=r.p)
rows = [("Gridded EKE", {"t": row("ACC_EKE_allsat")}), ("Gridded KE", {"t": row("ACC_KE_allsat")}),
        ("Σ amplitude²", {"t": row("ACC_amp2_per_day")}), ("Mean amplitude", {"t": row("ACC_amp")}),
        ("Rotational speed", {"t": row("ACC_spd")}), ("Eddies present", {"t": row("ACC_present_per_day")}),
        ("Births per month", {"t": row("ACC_births")}), ("Mean radius", {"t": row("ACC_rad")})]
forest_panel(ax, rows, series_colors={"t": "#2A9D8F"}, xlim=(-6, 18), label_size=8.5)
ax.get_legend().remove() if ax.get_legend() else None
ax.set_title(("b" if MAIN else "c") + ") Trends 2001–2023", loc="left", fontsize=12, fontweight="bold")
nm = "Figure4" if MAIN else "FigureS_eddies"   # eddies = Figure 4 in v2.1; meanders = Figure 5
fig.savefig(FIGV3 / f"{nm}.png"); fig.savefig(FIGV3 / f"{nm}.pdf")
print("saved")
