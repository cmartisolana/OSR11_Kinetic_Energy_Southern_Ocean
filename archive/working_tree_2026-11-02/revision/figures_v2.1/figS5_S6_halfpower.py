"""Supporting Information figures (v2.6), same style as the main figures (style_nc).
S1  two-satellite product: a mean KE map; b ACC-band KE and EKE in both products; c ACC-band trends of KE, EKE, its
    low-frequency and sub-monthly parts, along-stream speed and surface transport proxy; d sector trends, both products.
S2  half-power analysis: a sensitivity of the circumpolar position trend to the detection parameters; b trends in 5-degree
    longitude bins; c detection and d recovery of injected shifts.
S3  half-power latitude and peak KE trends in both products.
S4  tracked eddies: a births per year per 1-degree cell; b nonlinearity ratio U/c.
S5  month-to-month covariation of path length and EKE (b2): same-window and downstream links, SAF and PF levels.
Usage: OSR11_CACHE=<cache> python si_figures.py [S1 S2 S3 S4 S5]
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
import cartopy.crs as ccrs, cartopy.feature as cfeature
from matplotlib.colors import LogNorm
from style_nc import *  # noqa

FIGD = OUT / "figures_final"
ROOT_OUT = HERE.parent / "outputs"
WHICH = sys.argv[1:] or ["S1", "S2", "S3", "S4", "S5"]
PC = ccrs.PlateCarree()
A_, B_ = "All-satellite (Product 1)", "Two-satellite (Product 4)"
COLS = {A_: C["ink"], B_: C["two"]}; MK = {A_: "o", B_: "s"}
KM = 111.195
pa = lambda y: 100 * (y / np.nanmean(y) - 1)


def lon_labels(ax):
    t = np.arange(-180, 181, 60); ax.set_xticks(t)
    ax.set_xticklabels([("180°" if abs(v) == 180 else "0°" if v == 0 else f"{abs(v)}°{'W' if v < 0 else 'E'}") for v in t])


def a8_row(prod, s):
    r = pd.read_csv(OUT / "a8" / f"a8_trends_{prod}.csv").set_index("series").loc[s]
    return dict(pct=r.pct_per_decade, pct_lo=r.pct_lo, pct_hi=r.pct_hi, p=r.p)


# ---------------------------------------------------------------- S1
if "S1" in WHICH:
    cache = Path(os.environ.get("OSR11_CACHE", HERE / "cache"))
    lat = np.load(cache / "lat.npy"); lon = np.load(cache / "lon.npy")
    ice = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc").ice_mask.values > 0.5
    kem = np.nanmean(np.asarray(np.load(cache / "ke.npy", mmap_mode="r"), np.float32), 0) * 1e4; kem[ice] = np.nan
    S = {p: pd.read_csv(OUT / "a5" / f"a5_series_{p}.csv") for p in ("allsat", "twosat")}
    x = 2001 + (np.arange(300) + 0.5) / 12
    fig = plt.figure(figsize=(W2, 5.9))
    g = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.25], hspace=0.32)
    gt = g[0].subgridspec(1, 2, width_ratios=[0.75, 1.5], wspace=0.25); gb = g[1].subgridspec(1, 2, wspace=0.85)
    ax = fig.add_subplot(gt[0], projection=ccrs.SouthPolarStereo()); polar(ax, ccrs, cfeature)
    LON, LAT = np.meshgrid(lon, lat)
    im = ax.pcolormesh(LON, LAT, kem, transform=PC, cmap="YlGnBu", norm=LogNorm(20, 2000), shading="auto", rasterized=True)
    add_ice_mask(ax, ccrs, lon, lat, ice)
    hcbar(fig, im, ax, "KE (cm$^2$ s$^{-2}$)")
    plabel(ax, "a", "Mean KE, Product 1")
    ax = fig.add_subplot(gt[1])
    for v, col in (("KE", C["ke"]), ("EKE", C["eke"])):
        for p, ls in (("allsat", "-"), ("twosat", (0, (3, 1.5)))):
            ax.plot(x, running(pa(S[p][f"ACC_{v}"].values)), color=col, ls=ls, lw=1.0,
                    label=f"{v}, {'Product 1' if p == 'allsat' else 'Product 4'}")
    zero(ax); ts_style(ax); ax.set_ylabel("% of mean")
    ax.legend(loc="upper left", ncol=2)
    yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.25 * (yl[1] - yl[0]))
    plabel(ax, "b", "ACC band, 12-month running means", dx=-24)
    terms = [("KE", "KE"), ("EKE", "EKE"), ("LF", "EKE, low-frequency part"), ("SUB", "EKE, sub-monthly part"),
             ("SPD", "Along-stream speed")]
    rows = [(lab, {A_: retest_row(f"allsat|ACC_{t}|2001-2025"), B_: retest_row(f"twosat|ACC_{t}|2001-2025")}) for t, lab in terms]
    rows.append(("Surface transport proxy", {A_: a8_row("allsat", "ACC_T"), B_: a8_row("twosat", "ACC_T")}))
    ax = fig.add_subplot(gb[0])
    forest(ax, rows, COLS, MK, xlim=(-2, 25), xlabel="Trend 2001–2025 (% per decade)")
    plabel(ax, "c", "ACC band", dx=-80)
    rows = []
    for reg in ("Pacific", "Atlantic", "Indian"):
        for v in ("KE", "EKE", "SPD"):
            rows.append((f"{reg}: {'speed' if v == 'SPD' else v}",
                         {A_: retest_row(f"allsat|{reg}_{v}|2001-2025"), B_: retest_row(f"twosat|{reg}_{v}|2001-2025")}))
        rows.append((f"{reg}: transport", {A_: a8_row("allsat", f"{reg}_T"), B_: a8_row("twosat", f"{reg}_T")}))
    ax = fig.add_subplot(gb[1])
    forest(ax, rows, COLS, MK, xlim=(-6, 20), xlabel="Trend 2001–2025 (% per decade)", ms=3.2)
    for yy in (-3.5, -7.5): ax.axhline(yy, color=C["light"], lw=0.5)
    plabel(ax, "d", "Sectors", dx=-60)
    fig.legend(handles=sig_handles(COLS, MK), loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    save(fig, FIGD, "FigureS1")

# ---------------------------------------------------------------- S2, S3
if "S2" in WHICH or "S3" in WHICH:
    R = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")

    def kmrow(key):
        r = R.loc[key]; return dict(pct=r.slope * KM * 10, pct_lo=r.ci_lo * KM * 10, pct_hi=r.ci_hi * KM * 10, p=r.p_new)
if "S2" in WHICH:
    fig = plt.figure(figsize=(W2, 4.4))
    g = fig.add_gridspec(2, 2, width_ratios=[0.85, 1.25], height_ratios=[1.0, 0.85], hspace=0.6, wspace=0.55)
    sets = [("baseline", "Baseline (200 cm$^2$ s$^{-2}$, 25 %, 100 km)"), ("thr150", "Threshold 150 cm$^2$ s$^{-2}$"),
            ("thr250", "Threshold 250 cm$^2$ s$^{-2}$"), ("thr300", "Threshold 300 cm$^2$ s$^{-2}$"),
            ("thr400", "Threshold 400 cm$^2$ s$^{-2}$"), ("relthr", "Relative threshold"),
            ("frac20", "Peak fraction 20 %"), ("frac33", "Peak fraction 33 %"), ("frac50", "Peak fraction 50 %"),
            ("minw50", "Minimum width 50 km"), ("minw200", "Minimum width 200 km"), ("minw300", "Minimum width 300 km")]
    P = "circumpolar position"
    ax = fig.add_subplot(g[:, 0])
    forest(ax, [(lab, {P: kmrow(f"{s}|ACC|lat_hp")}) for s, lab in sets], {P: C["pos"]}, xlim=(-10, 12),
           xlabel="Trend (km per decade, positive northward)", gap=0.6)
    ax.axhline(-0.5, color=C["light"], lw=0.5)
    plabel(ax, "a", "Detection parameters", dx=-118)
    ax = fig.add_subplot(g[0, 1])
    L = sorted([(int(k.split("|")[1]) + 2.5, kmrow(k)) for k in R.index if k.startswith("lonbin|")])
    lo_ = np.array([l for l, _ in L]); d = pd.DataFrame([r for _, r in L]); sig = (d.p < 0.05).values
    zero(ax)
    ax.vlines(lo_, d.pct_lo, d.pct_hi, color=C["pos"], lw=0.6)
    ax.plot(lo_[sig], d.pct[sig], "o", color=C["pos"], ms=2.8)
    ax.plot(lo_[~sig], d.pct[~sig], "o", mfc="white", color=C["pos"], ms=2.8, mew=0.6)
    ax.set_xlim(-180, 180); lon_labels(ax); ax.set_ylabel("km per decade"); ax.set_xlabel("Longitude (bin centre)")
    plabel(ax, "b", "Trends in 5° longitude bins", dx=-26)
    gg = g[1, 1].subgridspec(1, 2, wspace=0.55)
    A3 = pd.read_csv(OUT / "a3" / "a3_summary.csv").set_index("scenario")
    sh = [("S0_null", 0), ("S1_shift5km", 5), ("S1_shift10km", 10), ("S1_shift25km", 25), ("S1_shift50km", 50)]
    xs = np.array([s for _, s in sh]); pw = np.array([A3.loc[k, "lat_hp_slope_pow"] for k, _ in sh]) * 100
    rec = np.array([A3.loc[k, "lat_hp_slope_med"] for k, _ in sh]) * KM * 25
    ax = fig.add_subplot(gg[0])
    ax.plot(xs, pw, "o-", color=C["pos"], ms=2.8, lw=0.8); ax.axhline(5, color=C["grey"], lw=0.5, ls="--")
    ax.set_xlabel("Injected shift (km)"); ax.set_ylabel("Detected (%)"); ax.set_ylim(0, 105)
    plabel(ax, "c", "Detection", dx=-26)
    ax = fig.add_subplot(gg[1])
    ax.plot(xs, rec, "o-", color=C["pos"], ms=2.8, lw=0.8); ax.plot([0, 50], [0, 50], color=C["grey"], lw=0.5, ls="--")
    ax.set_xlabel("Injected shift (km)"); ax.set_ylabel("Recovered (km)")
    plabel(ax, "d", "Recovery", dx=-26)
    save(fig, FIGD, "FigureS2")
    print("lonbins", len(d), "sig", int(sig.sum()))

if "S3" in WHICH:
    H = pd.read_csv(OUT / "a2" / "a2h_twosat_vs_allsat.csv")
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.1), gridspec_kw=dict(wspace=0.6))
    for ax, (fld, lim, xl, t, l) in zip(axs, (("pos", (-12, 16), "Trend (km per decade, positive northward)", "Half-power latitude", "a"),
                                             ("peak", (-4, 12), "Trend (% per decade)", "KE at the envelope peak", "b"))):
        rows = []
        for reg, lab in (("ACC", "Circumpolar"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian")):
            dd = {}
            for prod, nm in (("all-sat", A_), ("two-sat", B_)):
                r = H[(H["product"] == prod) & (H.region == reg)].iloc[0]
                dd[nm] = (dict(pct=r.pos_km_dec, pct_lo=r.pos_lo, pct_hi=r.pos_hi, p=r.pos_p) if fld == "pos" else
                          dict(pct=r.peak_pct_dec, pct_lo=r.peak_lo, pct_hi=r.peak_hi, p=r.peak_p))
            rows.append((lab, dd))
        forest(ax, rows, COLS, MK, xlim=lim, xlabel=xl)
        plabel(ax, l, t, dx=-50)
    fig.legend(handles=sig_handles(COLS, MK), loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12))
    save(fig, FIGD, "FigureS3")

# ---------------------------------------------------------------- S4
if "S4" in WHICH:
    D = xr.open_dataset(OUT / "eddies" / "e1_density.nc")
    N = pd.read_csv(OUT / "eddies" / "e5_nonlinearity_series_lag7.csv")
    adt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt
    fig = plt.figure(figsize=(W2, 2.7))
    g = fig.add_gridspec(1, 2, width_ratios=[0.75, 1.4], wspace=0.25)
    ax = fig.add_subplot(g[0], projection=ccrs.SouthPolarStereo()); polar(ax, ccrs, cfeature)
    Lo, La = np.meshgrid(D.lon.values, D.lat.values)
    im = ax.pcolormesh(Lo, La, D.births_per_year.values, transform=PC, cmap="magma_r", vmin=0, vmax=6, shading="auto", rasterized=True)
    add_acc_contours(ax, ccrs, adt.longitude, adt.latitude, adt.values)
    hcbar(fig, im, ax, "Births per year per 1° cell", extend="max")
    plabel(ax, "a", "Eddy generation, 2001–2023")
    ax = fig.add_subplot(g[1])
    xx = N.x.values
    ax.plot(xx, N.ratio_median, color=C["te"], lw=0.35, alpha=0.4)
    ax.plot(xx, running(N.ratio_median.values), color=C["te"], lw=1.0)
    ts_style(ax, (2000.8, 2023.8)); ax.set_ylabel("Median U/c")
    plabel(ax, "b", "Nonlinearity ratio of tracked eddies, ACC band", dx=-24)
    save(fig, FIGD, "FigureS4")

# ---------------------------------------------------------------- S5
if "S5" in WHICH:
    fig, axs = plt.subplots(2, 1, figsize=(W2, 3.4), sharex=True, gridspec_kw=dict(hspace=0.4))
    for ax, link, t, l in zip(axs, ("E", "D"), ("EKE within 150 km of the path, same window", "EKE 15–35° downstream"), "ab"):
        for sub, lab, col in (("b2", "Subantarctic Front level (+0.18 m)", C["ink"]), ("b2_pf", "Polar Front level (−0.20 m)", C["two"])):
            T = pd.read_csv(OUT / sub / "b2_window_tests_v2.csv")
            R0 = T[[f"r_{link}_lag{k}" for k in range(4)]].values
            ib = np.abs(R0).argmax(1); r = R0[np.arange(len(T)), ib]; ok = T[f"p_{link}"].values < 0.05
            ax.plot(T.centre, r, "-", color=col, lw=0.7, label=lab)
            ax.plot(T.centre[ok], r[ok], "o", color=col, ms=2.8); ax.plot(T.centre[~ok], r[~ok], "o", mfc="white", color=col, ms=2.8, mew=0.6)
        zero(ax); ax.set_ylabel("r"); ax.set_ylim(-0.45, 0.9)
        plabel(ax, l, f"Path length vs {t}", dx=-24)
    axs[0].legend(loc="upper right", ncol=2)
    axs[1].set_xlim(-185, 180); lon_labels(axs[1]); axs[1].set_xlabel("Window centre (30° windows)")
    save(fig, FIGD, "FigureS5")
