"""Supplement Fig. S6 (final layout, 10 Nov 2026): occurrence of the monthly half-power points and height of the energetic
pathway, Product 1, 2001-2025.

Panel c: same data and numbers as figures_old/figS6_occurrence_adt.py. Panels a,b: relative-threshold occurrence (10 Nov 2026,
see below). Layout: (a) and (b) share one row with
the colour bar (under a) and the legend (under b); the legend includes the sea-ice key of the main figures (ice mask drawn
in b only, as in main-text Figs. 2a and 3a); the two 10 % contours use the main-figure
occurrence-contour convention, distinguishable in colour, lightness and line style (2001-2017 light blue, solid, wider;
2018-2025 near-black, dashed, on top), over a grey scale whose darkest shade stays mid-grey so that the dark dashed contour
remains visible; the overlap fractions are written in each panel; the c annotations are inside the canvas (checked).
  a  Drake Passage region, 86-40 W: occurrence of the monthly half-power points of the RELATIVE threshold (20 % of each
     longitude's maximum mean KE; outputs/x7/x7_points_all.npy) in % of months per cell (grey, raw density 2001-2025), with
     the 10 % contours of the Gaussian-smoothed (sigma = 1 cell) densities of 2001-2017 and 2018-2025; the area overlap
     written in the panel is the cos-weighted intersection over union of exactly these smoothed 10 % regions, computed here
     and checked against figures/s6_occurrence_relthr.json (densities: figures/s6_occurrence_relthr.npz, made by
     figures/s6_occurrence_relthr.py with the method of review_calc/c06/c07).
  b  As a, circumpolar (same file, circ_* grids).
  c  ADT at the half-power point, mean anomaly over the 1826 retained tracks (outputs/x7/x7_tracks.npz; trend adt_hp of
     outputs/x7/x7_trends.csv), with the ACC-band mean ADT (outputs/x1/x1_band_adt.npy); implied displacement of a fixed
     contour from outputs/x7/x7_implied_shift.csv (circumpolar).
Usage: python figures/figS6_occurrence_adt.py [outputs dir] [figure dir]   (defaults: $OSR11_OUTPUTS, rendered/)
"""
import _sifig as F
import numpy as np, pandas as pd, xarray as xr
from scipy.ndimage import gaussian_filter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as _plt

S, C, plt, OUT, FIGD, P = F.S, F.C, F.plt, F.OUT, F.FIGD, F.P
F.cartopy_local()
import cartopy.crs as ccrs, cartopy.feature as cfeature   # noqa: E402
PC = ccrs.PlateCarree()
bat = xr.open_dataset(P.BATHY).elevation; bc = bat.coarsen(latitude=4, longitude=4, boundary="trim").mean()
ice_ds = xr.open_dataset(P.ICE_MASK); ice_g = ice_ds.ice_mask.values > 0.5; lat8, lon8 = ice_ds.latitude.values, ice_ds.longitude.values
from pathlib import Path
HERE = Path(__file__).resolve().parent
OC = np.load(OUT / "s6" / "s6_occurrence_relthr.npz")                    # relative-threshold occurrence (s6_occurrence_relthr.py)
DD, DC = ({k[len(d) + 1:]: OC[k] for k in OC.files if k.startswith(d + "_")} for d in ("drake", "circ"))
X7 = OUT / "x7"; Z = np.load(X7 / "x7_tracks.npz"); TR = pd.read_csv(X7 / "x7_trends.csv"); TR = TR[TR.region == "circumpolar"].set_index("series")
import json
OV = json.load(open(OUT / "s6" / "s6_occurrence_relthr.json"))["relative_threshold"]     # overlaps of the smoothed 10 % regions, 2001-2017 vs 2018-2025
IMP = pd.read_csv(X7 / "x7_implied_shift.csv").set_index("region").loc["circumpolar"]; band = 100 * np.load(OUT / "x1" / "x1_band_adt.npy")
box = dict(longitude=slice(-87, -33), latitude=slice(-66, -46)); bs = bat.sel(**box)
EXT = (-86, -40, -65, -50)
CE, CL = C["early"], C["late"]                      # main-figure occurrence-contour colours
LWE, LWL = 1.3, 0.85
LSD = (0, (3.0, 1.6))
OCC = LinearSegmentedColormap.from_list("greys_mid", _plt.get_cmap("Greys")(np.linspace(0.0, 0.78, 256)))


def occ(ax, DN, transform=None):
    kw = dict(transform=transform) if transform is not None else {}
    dd = np.where(DN["all"] > 0, DN["all"], np.nan)
    im = ax.pcolormesh(DN["lon_edges"], DN["lat_edges"], dd, cmap=OCC, vmin=0, vmax=S.OCC_LIM, shading="flat", rasterized=True, zorder=1, **kw)
    lc = 0.5 * (DN["lon_edges"][1:] + DN["lon_edges"][:-1]); ac = 0.5 * (DN["lat_edges"][1:] + DN["lat_edges"][:-1])
    SM = {k: DN[k + "_smooth"] for k in ("early", "late")}
    assert all(np.allclose(SM[k], gaussian_filter(DN[k], 1.0)) for k in SM)          # contours = the grids of the IoU
    ax.contour(lc, ac, SM["early"], levels=[S.LEVB], colors=CE, linewidths=LWE, linestyles=["-"], zorder=4, **kw)
    ax.contour(lc, ac, SM["late"], levels=[S.LEVB], colors=CL, linewidths=LWL, linestyles=[LSD], zorder=5, **kw)
    me, ml = SM["early"] >= S.LEVB, SM["late"] >= S.LEVB; w = np.cos(np.radians(ac))[:, None] * np.ones_like(me, float)
    return im, (w * (me & ml)).sum() / (w * (me | ml)).sum()


H = 142.5; MH = 52.0; MW = MH * 46 / 15 * 0.5          # Drake map: 46 deg lon x 15 deg lat at 1:2 aspect
fig = S.figure(S.W2, H)
TOPY = H - MH - 9.0
# a ---------------------------------------------------------------- Drake Passage occurrence
ax = S.axmm(fig, 14, TOPY, MW, MH)
im, iou_d = occ(ax, DD)
assert abs(iou_d - OV["drake"]["split2018"]) < 5e-4
S.add_isobath(ax, bs.longitude, bs.latitude, -bs.values)
ax.contourf(bs.longitude, bs.latitude, bs.values, levels=[0, 1e5], colors=[C["land"]], zorder=6)
ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
ax.set_xticks([-80, -70, -60, -50, -40]); ax.set_xticklabels([S.lonfmt(v) for v in (-80, -70, -60, -50, -40)])
ax.set_yticks([-65, -60, -55, -50]); ax.set_yticklabels([S.latfmt(v) for v in (-65, -60, -55, -50)])
ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
ax.text(0.985, 0.03, f"area overlap of the\n{S.LEVB} % regions: {100 * iou_d:.0f} %", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=F.FS, linespacing=1.15, bbox=dict(facecolor="white", edgecolor="none", pad=0.8, alpha=0.9), zorder=9)
S.letter(fig, 2, H - 4.5, "a", "Drake Passage region")
# b ---------------------------------------------------------------- circumpolar occurrence
D = 60.0; XB = 102.5
axb = S.axmm(fig, XB, TOPY - 3.0, D, D, projection=ccrs.SouthPolarStereo())
S.polar(axb, ccrs, cfeature); S.add_isobath(axb, bc.longitude, bc.latitude, -bc.values, transform=PC); S.add_ice(axb, ccrs, lon8, lat8, ice_g)
_, iou_c = occ(axb, DC, transform=PC)
assert abs(iou_c - OV["circ"]["split2018"]) < 5e-4
S.letter(fig, 96, H - 4.5, "b", "Circumpolar")
F.note(fig, XB - 1.5, TOPY - 3.5, f"area overlap of the\n{S.LEVB} % regions: {100 * iou_c:.0f} %", ha="left", va="top")
# shared colour bar (under a) and legend (under b) ------------------------------------------
YR = TOPY - 15.5
S.cbar(fig, im, 14, YR, 62, "Half-power points (% of months per cell)", extend="max", ticks=S.OCC_TICKS)
fig.legend(handles=[Line2D([], [], color=CE, lw=LWE, label=f"{S.LEVB} % contour, 2001–2017"),
                    Line2D([], [], color=CL, lw=LWL, ls=LSD, label=f"{S.LEVB} % contour, 2018–2025"),
                    Line2D([], [], color=C["isobath"], lw=0.6, label="3000 m isobath"),
                    Patch(facecolor=C["ice_flat"], edgecolor="none", label="sea ice (b)")],
           loc="upper left", bbox_to_anchor=(86 / S.W2, (YR + 3.2) / H), ncol=2, handlelength=2.4, columnspacing=1.4, labelspacing=0.45,
           borderaxespad=0)
# c ---------------------------------------------------------------- ADT at the pathway vs band-mean ADT
t = 2001 + (np.arange(300) + 0.5) / 12; tc = t.mean()


def series(a):
    return np.nanmean(a - np.nanmean(a, 0), 1)


CA = C["vel"]
YC, HC = 8.5, 40.0
ax = S.axmm(fig, 16, YC, 112, HC)
y = series(Z["adt_hp"]) * 100; r = TR.loc["adt_hp"]
ax.plot(t, y, color=CA, lw=S.LW["raw"], alpha=0.5, label="_")
ax.plot(t, S.runfull(y), color=CA, lw=S.LW["smooth"] + 0.1, label="ADT at the half-power point (12-month running mean)")
S.trend_line(ax, 2001, 2026, tc, 0.0, r.slope / 10, r.lo / 10, r.hi / 10, r.p, CA, wedge=True)
yb = band - band.mean(); ax.plot(t, S.runfull(yb), color="#222222", lw=0.8, ls=LSD, label="ACC-band mean ADT (12-month running mean)")
ax.text(2026.7, 4.0, S.um(f"{r.slope:+.1f} cm per decade\n({r.lo:+.1f} to {r.hi:+.1f}), {S.fp(r.p)}"), color=CA, fontsize=F.FS, ha="left",
        va="center", linespacing=1.15, clip_on=False)
ax.text(2026.7, -3.0, S.um(f"implied displacement\nof a fixed contour across\nthe median gradient\n({IMP.median_grad_cm_per_10km:.1f} cm per 10 km):\n"
                           f"{IMP.implied_contour_shift_km_per_decade:+.1f} km per decade"), color="#222222", fontsize=F.FS, ha="left", va="center",
        linespacing=1.15, clip_on=False)
ax.legend(loc="upper left", frameon=False, handlelength=2.2)
ax.set_ylabel("ADT anomaly (cm)"); ax.set_ylim(-9, 9); S.zero(ax); S.yticks_minus(ax); S.ts_axis(ax, 2000.8, 2026.2)
S.letter(fig, 2, YC + HC + 4.5, "c", f"Height of the energetic pathway ({Z['lat_hp'].shape[1]} tracks) and sea level of the ACC band")
F.check(fig); S.save(fig, FIGD, "FigureS6"); plt.close(fig)
print(f"S6a,b relative threshold, area overlap of the smoothed 10 % regions 2001-2017 vs 2018-2025 (as written): Drake {iou_d:.3f}, circumpolar {iou_c:.3f}")
print(f"S6c ADT at half-power point: {r.slope:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}] cm per decade, p = {r.p:.3f}; tracks {Z['lat_hp'].shape[1]}; "
      f"median gradient {IMP.median_grad_cm_per_10km:.2f} cm per 10 km; implied {IMP.implied_contour_shift_km_per_decade:+.2f} km per decade")
