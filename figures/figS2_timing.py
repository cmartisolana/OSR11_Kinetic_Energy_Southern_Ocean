"""Supplement Fig. S2 (final layout, 10 Nov 2026): timing of the EKE increase between the mean SAF and PF, 86-40 W, Product 1.

Same data and numbers as figures_old/figS2_timing.py; layout and presentation changed, one reference line added:
170 x 104 mm; a and b share the longitude axis (a above b, same width), vertical colour bar at the right of a; the domain
labels and values of b sit above the bars (no overlap); new in b: the ACC-band-wide increase of Product 1 as a light
reference line labelled "ACC band" (cited in the main text and SI Sect. S2.3).
  a  Annual EKE between the mean SAF and PF in the 2 deg bins of Figs. 2c and 3d (outputs/drake/dk12/dk12_band_annual_allsat.csv),
     anomaly in % of each bin's 2001-2025 mean; horizontal line: 2017/2018 split; dashed verticals: the passage (68-58 W).
  b  Increase of the 2018-2025 mean over the 2001-2017 mean per bin (ratio of period means, bars) and of the upstream /
     passage / downstream domain means (black lines and labels). ACC band (grey dashed line): the same ratio of period means
     for the monthly ACC-band EKE of Product 1 (outputs/a5/a5_series_allsat.csv, column ACC_EKE; mean of the 204 months of
     2001-2017 and of the 96 months of 2018-2025); Product 4 (a5_series_twosat.csv) is printed for reference.
Usage: python figures/figS2_timing.py [outputs dir] [figure dir]   (defaults: $OSR11_OUTPUTS, rendered/)
"""
import _sifig as F
import numpy as np, pandas as pd
from matplotlib.colors import TwoSlopeNorm

S, C, plt, OUT, FIGD = F.S, F.C, F.plt, F.OUT, F.FIGD
SPLIT = 2018
A = pd.read_csv(OUT / "drake" / "dk12" / "dk12_band_annual_allsat.csv", index_col=0); lon = A.columns.astype(float).values
sel = (lon > -86) & (lon < -40); A = A.loc[:, sel]; lon = lon[sel]; yrs = A.index.values.astype(int)
e, l = yrs < SPLIT, yrs >= SPLIT
DOM = [("upstream", "86–68°W", (-86, -68)), ("passage", "68–58°W", (-68, -58)), ("downstream", "58–40°W", (-58, -40))]
XT = [-80, -70, -60, -50, -40]


def ratio(ser, early, late):
    return 100 * (np.mean(ser[late]) / np.mean(ser[early]) - 1)


# ACC band, Product 1 (and Product 4 for reference): ratio of the period means of the monthly series
SER = {p: pd.read_csv(OUT / "a5" / f"a5_series_{p}.csv") for p in ("allsat", "twosat")}
ACC = {p: ratio(d.ACC_EKE.values, d.year.values < SPLIT, d.year.values >= SPLIT) for p, d in SER.items()}
NM = {p: (int((d.year < SPLIT).sum()), int((d.year >= SPLIT).sum())) for p, d in SER.items()}

H = 104.0; X0, AW = 18.0, 128.0
fig = S.figure(S.W2, H)
# a ---------------------------------------------------------------- Hovmoller
YA, HA = 49.0, 49.0
ax = S.axmm(fig, X0, YA, AW, HA)
im = ax.pcolormesh(S.edges(lon), np.r_[yrs - 0.5, yrs[-1] + 0.5], A.values, cmap=S.TREND_CMAP, norm=TwoSlopeNorm(0, -40, 40),
                   shading="flat", rasterized=True)
for x_ in (-68, -58): ax.axvline(x_, color="#444444", lw=0.6, ls=(0, (2, 2)))
ax.axhline(SPLIT - 0.5, color="#222222", lw=0.7)
ax.set_xlim(-86, -40); ax.set_xticks(XT); ax.set_xticklabels([])
ax.set_ylim(2025.5, 2000.5); ax.set_yticks([2001, 2005, 2010, 2015, 2020, 2025]); ax.set_ylabel("Year")
ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
S.letter(fig, 2, H - 5, "a", "Annual EKE anomaly between the SAF and PF, 2° bins")
cax = S.axmm(fig, X0 + AW + 3, YA, 2.4, HA)
cb = fig.colorbar(im, cax=cax, orientation="vertical", extend="both", extendfrac=0.04)
cb.set_ticks([-40, -20, 0, 20, 40]); cb.set_ticklabels([S.um(f"{t:g}") for t in (-40, -20, 0, 20, 40)])
cb.ax.tick_params(labelsize=F.FS, width=0.4, length=2, pad=1.5); cb.outline.set_linewidth(0.4)
cb.set_label("EKE anomaly (% of the bin's\n2001–2025 mean)", fontsize=F.FS, labelpad=2)
# b ---------------------------------------------------------------- change per bin
YB, HB = 11.0, 28.0
ax = S.axmm(fig, X0, YB, AW, HB); S.zero(ax)
frac = A.values / 100; ch = 100 * ((1 + frac[l].mean(0)) / (1 + frac[e].mean(0)) - 1)
ax.bar(lon, ch, 1.8, color=C["eke"], alpha=0.75, lw=0, zorder=2)
dom_v = {}
for nm, rng, (lo_, hi_) in DOM:
    m = (lon > lo_) & (lon < hi_); s_ = frac[:, m].mean(1); v = 100 * ((1 + s_[l].mean()) / (1 + s_[e].mean()) - 1); dom_v[nm] = v
    ax.plot([lo_, hi_], [v, v], color="#222222", lw=1.4, solid_capstyle="butt", zorder=4)
    xm = (lo_ + hi_) / 2
    ax.text(xm, 47.5, f"{nm} {rng}", ha="center", va="top", fontsize=F.FS, color=F.GREY_TXT)
    ax.text(xm, 40.5, S.um(f"{v:+.0f} %"), ha="center", va="top", fontsize=F.FS, color="#222222", fontweight="bold")
for x_ in (-68, -58): ax.axvline(x_, color="#444444", lw=0.6, ls=(0, (2, 2)), zorder=1)
# ACC-band reference (Product 1)
ax.axhline(ACC["allsat"], color="#8A8A8A", lw=0.8, ls=(0, (4, 2)), zorder=3)
ax.text(-39.4, ACC["allsat"], S.um(f"ACC band\n{ACC['allsat']:+.0f} %"), ha="left", va="center", fontsize=F.FS, color=F.GREY_TXT,
        linespacing=1.1, clip_on=False)
ax.set_xlim(-86, -40); ax.set_xticks(XT); ax.set_xticklabels([S.lonfmt(v) for v in XT])
ax.set_ylim(-3, 48); ax.set_yticks([0, 10, 20, 30]); S.yticks_minus(ax); ax.spines["left"].set_bounds(-3, 33)
ax.set_ylabel(f"{SPLIT}–2025 mean over\n2001–{SPLIT - 1} mean (%)"); ax.set_xlabel("Longitude (flow from left to right)")
S.letter(fig, 2, YB + HB + 4.5, "b", f"Increase of the {SPLIT}–2025 mean over the 2001–{SPLIT - 1} mean")
F.check(fig); S.save(fig, FIGD, "FigureS2"); plt.close(fig)
print("S2b domain increases (%):", {k: round(v, 1) for k, v in dom_v.items()})
print("S2b bins: min %.1f max %.1f, all positive: %s" % (ch.min(), ch.max(), bool((ch > 0).all())))
print("S2b ACC band (ratio of period means of monthly ACC_EKE, 2018-2025 over 2001-2017): Product 1 %+.2f %% (months %d/%d); "
      "Product 4 %+.2f %% (reference, not drawn)" % (ACC["allsat"], *NM["allsat"], ACC["twosat"]))
