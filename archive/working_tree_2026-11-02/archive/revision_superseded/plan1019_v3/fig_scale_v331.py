"""Scale-aware analysis (dk20-dk22), one figure: map of the eddy-dominance scale l_e (dk21) and the eddy fraction per scale band over the ACC band, Product 1.
a l_e (km): scale below which the monthly eddy part dominates the surface KE of that band; mean SAF and PF; 3000 m isobath.
b eddy fraction EKE_band / KE_band against scale, ACC band mean, with the 25-75 % range over ACC-band pixels.
Usage: python figS_eddy_scale_v331.py <revision dir> <drake out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, xarray as xr
REV, DO, FIGD = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import cartopy.crs as ccrs, cartopy.feature as cfeature
from scipy.ndimage import gaussian_filter
from matplotlib.colors import BoundaryNorm
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
Z = np.load(DO / "dk21" / "dk21_eddy_scale.npz"); lat, lon = Z["lat"], Z["lon"]; le = Z["le"].astype(float); FB = Z["FB"]; cent = Z["cent"]
OUT = REV / "outputs"; ROOT_OUT = REV.parent / "outputs"
clim = xr.open_dataset(OUT / "clim_allsat.nc"); Z12 = np.load(DO / "dk12" / "dk12_band_series_allsat.npz")
madt = xr.open_dataset(ROOT_OUT / "mean_adt.nc").mean_adt.interp(latitude=lat, longitude=lon).values
ice = xr.open_dataset(ROOT_OUT / "sea_ice_mask.nc").ice_mask.values > 0.5
bat = xr.open_dataset(REV / "data" / "bathymetry_so.nc").elevation; bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
# l_e from band energies smoothed over ~2 deg (local filtering spectra are not positive definite at single pixels)
from scipy.ndimage import gaussian_filter as gf
E, M = Z["E"].astype(float), Z["M"].astype(float); K = E + M; val = Z["valid"]
def sm2(a): return gf(np.where(val, a, 0), (8, 16)) / np.maximum(gf(val.astype(float), (8, 16)), 1e-6)
nb = len(cent); FBs = []
for i in range(nb):
    eb = E[i] - (E[i + 1] if i + 1 < nb else 0); kb = K[i] - (K[i + 1] if i + 1 < nb else 0)
    FBs.append(sm2(eb) / np.maximum(sm2(kb), 1e-6))
FBs = np.array(FBs); lc = np.log(cent); sm = np.full(madt.shape, np.nan)
for i in range(nb - 1):
    a, b = FBs[i], FBs[i + 1]; cross = np.isnan(sm) & (a >= 0.5) & (b < 0.5)
    with np.errstate(all="ignore"): t = np.clip((a - 0.5) / (a - b), 0, 1)
    sm = np.where(cross, np.exp(lc[i] + t * (lc[i + 1] - lc[i])), sm)
sm = np.where(np.isnan(sm) & (FBs[0] < 0.5), 50.0, sm)          # mean flow dominates at all resolved scales
sm = np.where(np.isnan(sm), 1500.0, sm)                          # eddies dominate at all scales
band = np.isfinite(madt) & (madt >= -0.8) & (madt <= 0.4) & val
sm[~band] = np.nan; le = sm
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
D20 = pd.read_csv(DO / "dk20" / "dk20_trends_allsat.csv"); D20 = D20[D20.band_km != "all"]
D22 = pd.read_csv(DO / "dk22" / "dk22_trends.csv")
fig = plt.figure(figsize=(S.W2, 1.18 * S.W2))
gR = fig.add_gridspec(3, 1, height_ratios=[1.15, 0.72, 0.95], hspace=0.38)
g1 = gR[0].subgridspec(1, 2, width_ratios=[1.15, 1], wspace=0.25)
# a map
ax = fig.add_subplot(g1[0], projection=ccrs.SouthPolarStereo()); S.polar(ax, ccrs, cfeature)
lev = [0, 100, 150, 200, 250, 300, 400, 500, 700, 2000]; cm = plt.get_cmap("viridis", len(lev) - 1); nm = BoundaryNorm(lev, cm.N)
im = ax.pcolormesh(lon[::2], lat[::2], sm[::2, ::2], transform=PC, cmap=cm, norm=nm, shading="nearest", rasterized=True, zorder=0)
ax.contour(bc.longitude, bc.latitude, -bc, levels=[3000], colors="#8A8A8A", linewidths=0.25, transform=PC, zorder=2)
ax.contour(clim.longitude, clim.latitude, clim.adt.values, levels=[float(Z12["lev_pf"]), float(Z12["lev_saf"])], colors="white", linewidths=0.6, linestyles="solid", transform=PC, zorder=6)
add_ice_mask(ax, ccrs, lon, lat, ice)
S.plabel(ax, "a", "Scale below which eddies dominate the KE")
cb = S.hcbar(fig, im, ax, "eddy-dominance scale $\\ell_e$ (km)", shrink=0.75); cb.set_ticks(lev[:-1] + [2000]); cb.set_ticklabels([str(v) for v in lev[:-1]] + ["all"])
# b eddy fraction
ax = fig.add_subplot(g1[1])
acc = np.isfinite(madt) & (madt >= -0.6) & (madt <= 0.2) & Z["valid"]; wA = np.cos(np.radians(lat))[:, None] * acc
med = [float((wA * (E[i] - (E[i + 1] if i + 1 < nb else 0))).sum() / (wA * (K[i] - (K[i + 1] if i + 1 < nb else 0))).sum()) for i in range(nb)]
p25 = [np.nanpercentile(FBs[i][acc], 25) for i in range(nb)]; p75 = [np.nanpercentile(FBs[i][acc], 75) for i in range(nb)]
ax.fill_between(cent, p25, p75, color=C["eke"], alpha=0.2, lw=0); ax.plot(cent, med, "-o", color=C["eke"], ms=3, lw=1)
ax.axhline(0.5, color="#777777", lw=0.5, ls=(0, (3, 2))); ax.axvspan(30, 100, color="#F1F1F1", zorder=0)
ax.set_xscale("log"); ax.set_xlim(40, 1800); ax.set_ylim(0, 1); ax.set_xticks([50, 100, 200, 500, 1000]); ax.set_xticklabels(["50", "100", "200", "500", "1000"])
ax.set_xlabel("length scale (km, band centre)"); ax.set_ylabel("eddy fraction of band KE")
ax.text(70, 0.04, "below effective\nresolution", fontsize=4.8, color="#777777", ha="center")
ax.text(0.97, 0.53, "eddy dominated", transform=ax.transAxes, ha="right", fontsize=5, color="#555555")
ax.text(0.97, 0.45, "mean-flow dominated", transform=ax.transAxes, ha="right", va="top", fontsize=5, color="#555555")
S.plabel(ax, "b", "Eddy fraction by scale, ACC band")
# c, d trends by scale band
B = ["0-100", "100-200", "200-300", "300-500", "500-1000", ">1000"]; LAB = ["<100", "100–200", "200–300", "300–500", "500–1000", ">1000"]
g2 = gR[1].subgridspec(1, 2, wspace=0.25)
for j, (col, ylab, lims, lt, tt) in enumerate((("abs_dec", "trend (cm$^2$ s$^{-2}$ per decade)", (-0.5, 8.5), "c", "Trend by scale band, ACC band"),
                                              ("pct", "trend (% per decade)", (0, 11), "d", "Relative trend by scale band"))):
    a = fig.add_subplot(g2[j])
    for p, c, dx in (("KE", C["ke"], -0.18), ("EKEm", C["eke"], 0.18)):
        q = D20[D20["field"] == p].set_index("band_km").loc[B]; x = np.arange(len(B)) + dx
        lo, hi = (q.abs_lo, q.abs_hi) if col == "abs_dec" else (q.pct_lo, q.pct_hi); sg = q.p.round(3) < 0.05
        if col == "pct" and p == "EKEm": sg = sg & (np.arange(len(B)) < 4); keep = np.arange(len(B)) < 5
        else: keep = np.ones(len(B), bool)
        a.bar(x[sg & keep], q[col][sg & keep], 0.34, color=c, zorder=3); a.bar(x[~sg & keep], q[col][~sg & keep], 0.34, color="white", edgecolor=c, lw=0.8, zorder=3)
        a.vlines(x[keep], lo[keep], np.minimum(hi[keep], lims[1]), color="#333333", lw=0.5, zorder=4)
    a.axvspan(-0.5, 0.5, color="#F1F1F1", zorder=0); a.set_xticks(range(len(B))); a.set_xticklabels(LAB, fontsize=5.2)
    a.set_xlabel("length scale (km)"); a.set_ylabel(ylab); a.set_ylim(*lims); S.zero(a); S.plabel(a, lt, tt)
    if j == 0:
        a.legend(handles=[Patch(facecolor=C["ke"], label="KE"), Patch(facecolor=C["eke"], label="EKE (monthly-mean part)"),
                          Patch(facecolor="white", edgecolor="#555555", label="open: p ≥ 0.05")], loc="upper right", fontsize=4.9)
# e, f meander sectors
REG = ["Kerguelen Plateau", "Southeast Indian Ridge", "Macquarie Ridge", "Pacific–Antarctic Ridge", "Drake Passage–Scotia Sea"]
LABR = ["Kerguelen\nPlateau", "Southeast\nIndian Ridge", "Macquarie\nRidge", "Pacific–Antarctic\nRidge", "Drake Passage–\nScotia Sea", "away from\nmeanders"]
SEC = [("upstream", "#BFD7EA"), ("ridge", "#6A9FCB"), ("downstream", "#1F4E79")]
g3 = gR[2].subgridspec(1, 2, wspace=0.18)
for j, (fld, band, title, lt) in enumerate((("EKEm", "100-300", "EKE, 100–300 km (eddy dominated), at meanders", "e"),
                                           ("KE", ">1000", "KE above 1000 km (jets), at meanders", "f"))):
    a = fig.add_subplot(g3[j]); q = D22[(D22.field == fld) & (D22.band_km == band)]
    for gi, rg in enumerate(REG):
        for si, (sec, c) in enumerate(SEC):
            r = q[(q.region == rg) & (q.sector == sec)].iloc[0]; x = gi + (si - 1) * 0.26; sg = round(r.p, 3) < 0.05
            a.bar(x, r.pct, 0.24, color=c if sg else "white", edgecolor=c, lw=0.8, zorder=3)
            a.vlines(x, max(r.pct_lo, -20), min(r.pct_hi, 30), color="#333333", lw=0.45, zorder=4)
            if r.fdr10: a.plot(x, min(r.pct_hi, 30) + 1.3, "v", color="#111111", ms=2.2, zorder=5)
    r = q[q.region.str.startswith("Away")].iloc[0]; sg = round(r.p, 3) < 0.05
    a.bar(5, r.pct, 0.3, color="#9A9A9A" if sg else "white", edgecolor="#9A9A9A", lw=0.8, zorder=3); a.vlines(5, r.pct_lo, r.pct_hi, color="#333333", lw=0.45, zorder=4)
    a.axvline(4.5, color="#9A9A9A", lw=0.5, ls=(0, (2, 2))); S.zero(a); a.set_ylim(-20, 33)
    a.set_xticks(range(6)); a.set_xticklabels(LABR, fontsize=5.5); a.tick_params(axis="x", length=0)
    if j == 0: a.set_ylabel("trend (% per decade)")
    S.plabel(a, lt, title)
    if j == 0:
        a.legend(handles=[Patch(facecolor=c, edgecolor=c, label=s) for s, c in SEC] + [Line2D([], [], ls="", marker="v", color="#111111", ms=3, label="FDR")],
                 loc="upper left", ncol=4, fontsize=4.7, columnspacing=0.7)
S.save(fig, FIGD, "Figure_scale"); plt.close(fig)
