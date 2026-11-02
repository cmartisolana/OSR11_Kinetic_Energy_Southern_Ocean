"""Fig. 5 (v1 style): is the ACC energy increase concentrated downstream of the standing meanders? (A7)
(a) relative EKE trend 2001-2025 (all-sat, 1° cells, % per decade of the cell mean; dots: FDR-significant
    absolute trend), with the downstream (WAKE) boxes (solid) and upstream boxes (dashed)
(b) difference in relative trend, WAKE minus REST of the ACC band (% per decade, 95 % CI, calibrated p),
    both products, 2001-2025 and 2001-2017
(c) EKE anomaly relative to the 2001-2025 mean (%), WAKE vs REST, 12-month running means, both products
(d) per feature: EKE relative-trend difference, downstream minus upstream box and downstream box minus REST
Usage: python fig5_main.py   (reads outputs/a1_maps.nc, outputs/a7, repository outputs/mean_adt.nc)
"""
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from figstyle import *  # noqa: F401,F403

ROOT_OUT = Path("/mnt/user-data/uploads/repository/outputs")
adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
M = xr.open_dataset(OUT / "a1_maps_bc.nc")   # revised test
T = pd.concat([pd.read_csv(OUT / "a7" / f"a7_tests_{p}.csv") for p in ("allsat", "twosat")], ignore_index=True)
T = retest_override(T, lambda i, r: f"{r['product']}|{r.period}|{r['var']}|{r.region_a}-{r.region_b}",
                    {"p_cal": "p_new", "diff_pct_dec": "slope", "diff_lo": "ci_lo", "diff_hi": "ci_hi"})
T["fdr_sig"] = False   # revised test: no secondary A7 test passes BH-FDR (retest_compare.py)
S = {p: pd.read_csv(OUT / "a7" / f"a7_series_{p}.csv") for p in ("allsat", "twosat")}
BOXES = {"Kerguelen": ((50, 60), (75, 95)), "SEIR–Macquarie": ((110, 130), (145, 175)),
         "PAR": ((-175, -160), (-150, -125)), "Drake": ((-100, -80), (-55, -35))}
KEYS = {"Kerguelen": "Kerguelen", "SEIR–Macquarie": "SEIR_Macquarie", "PAR": "PAR", "Drake": "Drake"}
C_PROD = {"allsat": "#1F2933", "twosat": "#7A8288"}
C_WAKE, C_REST = STYLE["eke"], "#7A8288"

fig = plt.figure(figsize=(11.0, 9.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.3], height_ratios=[1.1, 1], hspace=0.32, wspace=0.28)

# (a) map
ax = fig.add_subplot(gs[0, 0], projection=ccrs.SouthPolarStereo())
polar_ax_style(ax, ccrs, cfeature)
rel = 1000 * M.EKE_slope.values / M.EKE_mean.values
tc = plt.get_cmap("RdBu_r").copy(); tc.set_bad("#F3F2EF")
L1, A1 = np.meshgrid(M.lon.values, M.lat.values)
im = ax.pcolormesh(L1, A1, rel, transform=ccrs.PlateCarree(), cmap=tc, vmin=-25, vmax=25, shading="nearest",
                   rasterized=True, zorder=0)
add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
jj, ii = np.nonzero(M.EKE_fdr.values == 1)
ax.scatter(M.lon.values[ii], M.lat.values[jj], s=7, c="#1F2933", lw=0, transform=ccrs.PlateCarree(), zorder=6)


def box(ax, lo0, lo1, ls, col, lw):
    lo = np.r_[np.linspace(lo0, lo1, 40), np.full(20, lo1), np.linspace(lo1, lo0, 40), np.full(20, lo0)]
    la = np.r_[np.full(40, -64), np.linspace(-64, -42, 20), np.full(40, -42), np.linspace(-42, -64, 20)]
    ax.plot(lo, la, ls=ls, color=col, lw=lw, transform=ccrs.PlateCarree(), zorder=7)


for name, ((u0, u1), (d0, d1)) in BOXES.items():
    box(ax, d0, d1, "-", "#111111", 1.4); box(ax, u0, u1, (0, (3, 2)), "#555555", 0.9)
    lc = 0.5 * (d0 + d1)
    ax.text(lc, -38.8, name, transform=ccrs.PlateCarree(), ha="center", va="center", fontsize=7.5,
            color="#111111", zorder=8, bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.8))
ax.set_title("Relative EKE trend (all-sat)"); panel_label(ax, "a", x=0.0, y=1.02)
cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.75, pad=0.05, extend="both")
cb.set_label(r"% dec$^{-1}$", fontsize=9); cb.ax.tick_params(labelsize=8)
ax.plot([], [], "-", color="#111111", lw=1.4, label="downstream (WAKE)")
ax.plot([], [], ls=(0, (3, 2)), color="#555555", lw=0.9, label="upstream")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.27), fontsize=7.5, frameon=False, ncol=2)

# (b) WAKE - REST differences
ax = fig.add_subplot(gs[0, 1])
VARS = [("EKE", "EKE"), ("LF", "EKE, slow\n(>1 month)"), ("SUB", "EKE, fast\n(<1 month)"), ("KE", "KE"),
        ("SPD", "Along-stream\nspeed")]
off = {("allsat", "2001-2025"): -0.24, ("twosat", "2001-2025"): -0.08, ("allsat", "2001-2017"): 0.08,
       ("twosat", "2001-2017"): 0.24}
ax.axhline(0, color="#495057", lw=0.8)
for i, (v, lab) in enumerate(VARS):
    for (p, per), dx in off.items():
        r = T[(T["product"] == p) & (T.period == per) & (T["var"] == v) & (T.region_a == "WAKE") & (T.region_b == "REST")].iloc[0]
        mk = "o" if per == "2001-2025" else "^"
        ax.errorbar(i + dx, r.diff_pct_dec, yerr=[[r.diff_pct_dec - r.diff_lo], [r.diff_hi - r.diff_pct_dec]], fmt=mk,
                    ms=5, color=C_PROD[p], mfc=C_PROD[p] if per == "2001-2025" else "white", mew=1.1, elinewidth=1.0,
                    capsize=2, zorder=3)
        sig = "†" if r.fdr_sig else stars(r.p_cal)
        if r.kind == "primary":
            sig = stars(r.p_cal)
        if sig:
            ax.text(i + dx, r.diff_hi + 0.3, sig, ha="center", va="bottom", fontsize=8)
ax.set_xticks(range(len(VARS))); ax.set_xticklabels([l for _, l in VARS], fontsize=8.5)
ax.set_ylabel("WAKE − REST (% dec$^{-1}$)"); base_timeseries_style(ax)
ax.axvspan(-0.45, 0.45, color="#E9ECEF", lw=0, zorder=0)
ax.text(0, ax.get_ylim()[1] * 0.98, "primary", ha="center", va="top", fontsize=7.5, color="#495057")
for (p, per) in off:
    mk = "o" if per == "2001-2025" else "^"
    ax.plot([], [], mk, color=C_PROD[p], mfc=C_PROD[p] if per == "2001-2025" else "white", mew=1.1,
            label=f"{'all-sat' if p == 'allsat' else 'two-sat'}, {per.replace('-', '–')}")
ax.legend(loc="lower left", fontsize=7.5, frameon=False, ncol=2)
yl = ax.get_ylim(); ax.set_ylim(yl[0] - 2.5, yl[1])
panel_label(ax, "b", x=-0.1, y=1.06)

# (c) normalised EKE series
ax = fig.add_subplot(gs[1, 0])
x = S["allsat"].x.values
for p, ls in (("allsat", "-"), ("twosat", (0, (2, 1.5)))):
    for reg, col in (("WAKE", C_WAKE), ("REST", C_REST)):
        y = S[p][f"{reg}_EKE"].values; y = 100 * (y / y.mean() - 1)
        ax.plot(x, running(y), color=col, lw=STYLE["smooth_lw"] if p == "allsat" else 1.3, ls=ls,
                label=f"{reg}, {'all-sat' if p == 'allsat' else 'two-sat'}")
ax.axhline(0, color="#6b7280", lw=0.8, ls="--")
for p, y0 in (("allsat", 0.05), ("twosat", 0.13)):
    r = T[(T["product"] == p) & (T.period == "2001-2025") & (T["var"] == "EKE") & (T.region_a == "WAKE")].iloc[0]
    ax.text(0.98, y0, f"{'all-sat' if p == 'allsat' else 'two-sat'}: WAKE {r.pct_dec_a:+.1f}, REST {r.pct_dec_b:+.1f} % dec$^{{-1}}$ "
            f"(difference p={r.p_cal:.2f})", transform=ax.transAxes, ha="right", fontsize=7.8, color="#1f2937")
ax.set_xlabel("Year"); ax.set_ylabel("EKE anomaly (% of mean)"); base_timeseries_style(ax)
ax.legend(loc="upper left", fontsize=7.5, frameon=False, ncol=2); ax.set_xlim(2000.5, 2026)
yl = ax.get_ylim(); ax.set_ylim(yl[0] - 6, yl[1] + 5)
panel_label(ax, "c", x=-0.12, y=1.06)

# (d) per feature
ax = fig.add_subplot(gs[1, 1])
ax.axhline(0, color="#495057", lw=0.8)
offd = {("allsat", "up"): -0.24, ("twosat", "up"): -0.08, ("allsat", "REST"): 0.08, ("twosat", "REST"): 0.24}
for i, (name, key) in enumerate(KEYS.items()):
    for (p, cmp_), dx in offd.items():
        b = f"up_{key}" if cmp_ == "up" else "REST"
        r = T[(T["product"] == p) & (T.period == "2001-2025") & (T["var"] == "EKE") & (T.region_a == f"down_{key}") & (T.region_b == b)].iloc[0]
        mk = "s" if cmp_ == "up" else "D"
        ax.errorbar(i + dx, r.diff_pct_dec, yerr=[[r.diff_pct_dec - r.diff_lo], [r.diff_hi - r.diff_pct_dec]], fmt=mk,
                    ms=4.5, color=C_PROD[p], mfc=C_PROD[p] if cmp_ == "up" else "white", mew=1.1, elinewidth=1.0, capsize=2, zorder=3)
        if r.p_cal < 0.05:
            ax.text(i + dx, r.diff_hi + 0.5, stars(r.p_cal), ha="center", va="bottom", fontsize=8)
for (p, cmp_) in offd:
    mk = "s" if cmp_ == "up" else "D"
    ax.plot([], [], mk, color=C_PROD[p], mfc=C_PROD[p] if cmp_ == "up" else "white", mew=1.1,
            label=f"{'all-sat' if p == 'allsat' else 'two-sat'}: downstream − {'upstream' if cmp_ == 'up' else 'REST'}")
ax.set_xticks(range(4)); ax.set_xticklabels(list(KEYS), fontsize=9)
ax.set_ylabel("EKE trend difference (% dec$^{-1}$)"); base_timeseries_style(ax)
ax.legend(loc="lower left", fontsize=7.5, frameon=False, ncol=2)
yl = ax.get_ylim(); ax.set_ylim(yl[0] - 7, yl[1] + 1)
panel_label(ax, "d", x=-0.1, y=1.06)
fig.savefig(FIG / "Fig5_meander_localisation.png"); fig.savefig(FIG / "Fig5_meander_localisation.pdf")
print("saved")
