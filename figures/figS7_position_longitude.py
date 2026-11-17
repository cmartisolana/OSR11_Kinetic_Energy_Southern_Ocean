"""Supplement Fig. S7 (final layout, 10 Nov 2026): position of the energetic pathway by 2 deg longitude bin, Product 1, 2001-2025.

Same data and numbers as figures_old/figS7_position_longitude.py; layout and presentation only: legend above the panel (no
longer on the bathymetry); feature labels (PAR, DP, NSR, SWIR, KP, SEIR, MR; osr11.style.KEYS, as main-text Fig. 3d) and the
flow arrow just above the trend axis so that they never sit on data; purple and grey markers offset by -0.5 / +0.5 deg in
longitude; trend axis -80..+80 km per decade (was -66..+62, which hid one half-power value of +66 and cut 41 half-power
intervals): every bin value is now on scale and the intervals that still leave the axis end in a small triangle at the axis
limit.
Tracks of the longitude-dependent threshold (1826 tracks; outputs/x6/x6_relthr_bins.csv). Per bin: Theil-Sen trend of the
bin-mean half-power latitude (purple) and of the latitude of the fixed ADT contour through its mean position (grey),
calibrated test (Sect. S1.7), 95 % intervals as bars; filled: p < 0.05. Mean seafloor depth between the fronts below (dark
fill; outputs/base/bathymetry_so.nc and clim_allsat.nc, band between the SAF and PF levels of outputs/drake/dk12/
dk12_band_series_allsat.npz), as main-text Fig. 3d. Light grey: longitude window of main-text Fig. 2 (86-40 W).
Usage: python figures/figS7_position_longitude.py [outputs dir] [figure dir]   (defaults: $OSR11_OUTPUTS, rendered/)
"""
import _sifig as F
import numpy as np, pandas as pd, xarray as xr
from matplotlib.lines import Line2D

S, C, plt, OUT, FIGD, P = F.S, F.C, F.plt, F.OUT, F.FIGD, F.P
X5 = pd.read_csv(OUT / "x6" / "x6_relthr_bins.csv").sort_values("lon").rename(columns={"hp_slope": "ke_slope", "hp_lo": "ke_lo", "hp_hi": "ke_hi", "hp_p": "ke_p"})
ZB = np.load(OUT / "drake" / "dk12" / "dk12_band_series_allsat.npz"); LEV_SAF, LEV_PF = float(ZB["lev_saf"]), float(ZB["lev_pf"])
bat = xr.open_dataset(P.BATHY).elevation; clim = xr.open_dataset(P.CLIM)
_ab = clim.adt.interp(latitude=bat.latitude, longitude=bat.longitude).values
_band = (_ab >= LEV_PF) & (_ab <= LEV_SAF) & (bat.values < 0)
with np.errstate(all="ignore"):
    BDEP = np.nanmean(np.where(_band, -bat.values, np.nan), 0) / 1000
YL = (-80.0, 80.0)
H = 74.0; AX, AY, AW, AH = 16.0, 11.0, 142.0, 52.0
fig = S.figure(S.W2, H)
ax = S.axmm(fig, AX, AY, AW, AH); axd = ax.twinx()
axd.axvspan(-86, -40, color="#EEEEEE", lw=0, zorder=-1)
axd.fill_between(bat.longitude.values, BDEP, 7, color=C["land"], lw=0, zorder=0)
axd.set_ylim(6.2, -7.0); axd.set_yticks([2, 3, 4, 5, 6]); axd.set_ylabel("Depth (km)", color="#333333")
axd.yaxis.set_label_coords(1.045, 0.19); axd.tick_params(axis="y", colors="#333333")
axd.spines["right"].set_visible(True); axd.spines["right"].set_linewidth(0.5); axd.spines["right"].set_bounds(6.2, 1.2)
ax.set_zorder(axd.get_zorder() + 1); ax.patch.set_visible(False)
YLO = -190.0
ncl = {}
for key, col, dx in (("contour", "#7F7F7F", 0.5), ("ke", C["pos"], -0.5)):
    q = X5.dropna(subset=[f"{key}_slope"])
    lo, hi = q[f"{key}_lo"].values, q[f"{key}_hi"].values; x = q.lon.values + dx
    ax.vlines(x, np.maximum(lo, YL[0]), np.minimum(hi, YL[1]), color=col, lw=0.55, alpha=0.75, zorder=2)
    up, dn = hi > YL[1], lo < YL[0]
    ax.plot(x[up], np.full(up.sum(), YL[1]), "^", color=col, ms=2.0, mew=0, alpha=0.9, zorder=2.5, clip_on=False)
    ax.plot(x[dn], np.full(dn.sum(), YL[0]), "v", color=col, ms=2.0, mew=0, alpha=0.9, zorder=2.5, clip_on=False)
    v = q[f"{key}_slope"].values
    assert ((v >= YL[0]) & (v <= YL[1])).all()                     # every bin value on scale
    S.profile(ax, x, v, q[f"{key}_p"].values, col, ms=1.9, lw=0, z=3)
    ncl[key] = (int(up.sum() + dn.sum()), len(q))
S.zero(ax, lim=(-180, 180)); S.lon_axis(ax)
ax.set_ylim(YLO, YL[1] + 2); ax.set_yticks([-80, -40, 0, 40, 80]); ax.spines["left"].set_bounds(*YL); S.yticks_minus(ax)
ax.set_ylabel("Latitude trend\n(km per decade, + north)"); ax.yaxis.set_label_coords(-0.07, 0.70)
# feature labels and flow arrow just above the axis (labels and positions of main-text Fig. 3d)
for x_, s_ in S.KEYS:
    ax.text(x_, 1.012, s_, ha="center", va="bottom", fontsize=F.FS, color="#222222", transform=ax.get_xaxis_transform(), clip_on=False)
ax.annotate("", xy=(179, 1.055), xytext=(163, 1.055), xycoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", lw=0.6, color="#555555"),
            annotation_clip=False)
ax.text(171, 1.075, "flow", ha="center", va="bottom", fontsize=F.FS, color="#555555", transform=ax.get_xaxis_transform(), clip_on=False)
fig.legend(handles=[Line2D([], [], color=C["pos"], marker="o", ms=2.6, lw=0.6, label="half-power latitude (bin mean of the tracks)"),
                    Line2D([], [], color="#7F7F7F", marker="o", ms=2.6, lw=0.6, label="fixed ADT contour"),
                    Line2D([], [], color="#555555", marker="o", mfc="white", ms=2.6, lw=0, label="open: p ≥ 0.05"),
                    Line2D([], [], color="#555555", marker="^", ms=2.6, mew=0, lw=0, label="interval beyond the axis")],
           loc="lower left", bbox_to_anchor=(AX / S.W2, (AY + AH + 6.5) / H), ncol=4, handlelength=1.6, columnspacing=1.6, borderaxespad=0)
F.check(fig); S.save(fig, FIGD, "FigureS7"); plt.close(fig)
_q = X5.dropna(subset=["ke_slope"])
print(f"S7: half-power trend: median {_q.ke_slope.median():+.1f} km/dec, significant {int((_q.ke_p < 0.05).sum())} of {len(_q)} ({int(((_q.ke_p < 0.05) & (_q.ke_slope > 0)).sum())} north); "
      f"contour: median {X5.contour_slope.median():+.1f}, significant {int((X5.contour_p < 0.05).sum())} of {int(X5.contour_slope.notna().sum())}, south {int(((X5.contour_p < 0.05) & (X5.contour_slope < 0)).sum())}")
print(f"S7: intervals reaching beyond {YL[0]:.0f}..{YL[1]:.0f} km per decade (ends marked): half-power {ncl['ke'][0]} ends of {ncl['ke'][1]} bins, "
      f"contour {ncl['contour'][0]} ends of {ncl['contour'][1]} bins; bin values range half-power {_q.ke_slope.min():+.1f}..{_q.ke_slope.max():+.1f}, "
      f"contour {X5.contour_slope.min():+.1f}..{X5.contour_slope.max():+.1f}")
