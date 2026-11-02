"""Figure 4 (v3.28): response of the surface ACC to the wind stress, an observational version of the eddy-saturation
diagnostics of model studies (transport and eddy energy against wind stress: Meredith et al. 2012; Constantinou and
Hogg 2019; Stewart et al. 2023; lagged EKE response: Meredith and Hogg 2006; Hogg et al. 2015). Product 1 unless stated;
wind stress = eastward stress of Product 2 averaged over the ACC band (dk17/s1, 'tau_accmask_ACC').
a  Annual anomalies (% of the 2001-2025 mean) of the wind stress and surface transport proxy;
b  of EKE (Product 1 and Product 4) and KE.
c  Lagged correlation with the wind stress of deseasonalised, detrended 12-month running means (positive lag: ocean
   lags wind); grey: pointwise 95 % range of 1000 phase-randomised wind series; dashed: 95 % level for the maximum |r|
   over all lags (dk17/s3).
d  Sensitivity d ln X / d ln tau of annual means: interannual (both detrended; lag 0 and 1 yr), and the ratio of the
   2001-2025 trends (dk17/s4), with 95 % intervals; references: eddy saturation (transport insensitive, EKE ~ tau,
   Marshall et al. 2017) and EKE ~ tau^0.7 in Stewart et al. (2023).
Usage: python fig4_v328.py <revision dir> <drake out dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
REV, DO, FIGD = (Path(a) for a in sys.argv[1:4])
sys.path.insert(0, str(REV))
import style_nc as S
from matplotlib.lines import Line2D
plt = S.plt; C = S.C
D17 = DO / "dk17"; CSPD = "#56B4E9"; CT = "#0072B2"; CTAU = "#555555"
A = pd.read_csv(D17 / "s3_all_series.csv")
LC = pd.read_csv(D17 / "s3_lagcorr_curves.csv"); LS = pd.read_csv(D17 / "s3_lagcorr_summary.csv")
E = pd.read_csv(D17 / "s4_elasticities.csv")
an = A.groupby("year").mean(numeric_only=True); yrs = an.index.values + 0.5
pct = lambda c: 100 * (an[c] / an[c].mean() - 1)
fig = plt.figure(figsize=(S.W2, 0.92 * S.W2))
gA = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.42], hspace=0.38)
g = gA[0].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.25)
gl = g[0].subgridspec(2, 1, hspace=0.12)
# ---------------------------------------------------------------- a, b
ax = fig.add_subplot(gl[0]); S.zero(ax)
ax.bar(yrs, pct("tau_accmask_ACC"), 0.75, color="#D9D9D9", edgecolor="none", zorder=1, label="wind stress")
ax.plot(yrs, pct("ACC_T"), "-o", color=CT, lw=0.8, ms=2.0, zorder=3, label="surface transport proxy")
S.ts_style(ax, (2000.6, 2025.4)); ax.set_ylim(-12, 12); ax.set_ylabel("anomaly (%)"); plt.setp(ax.get_xticklabels(), visible=False)
ax.legend(loc="lower right", ncol=2, fontsize=5.0)
S.plabel(ax, "a", "Wind stress and surface transport proxy, ACC band")
ax = fig.add_subplot(gl[1]); S.zero(ax)
ax.bar(yrs, pct("tau_accmask_ACC"), 0.75, color="#D9D9D9", edgecolor="none", zorder=1, label="wind stress")
ax.plot(yrs, pct("ACC_EKE"), "-o", color=C["eke"], lw=0.8, ms=2.0, zorder=3, label="EKE")
ax.plot(yrs, pct("ACC2s_EKE"), "--", color=C["two"], lw=0.8, zorder=3, label="EKE, two-satellite")
ax.plot(yrs, pct("ACC_KE"), "-^", color=C["ke"], lw=0.7, ms=1.8, zorder=3, label="KE")
S.ts_style(ax, (2000.6, 2025.4)); ax.set_ylim(-20, 20); ax.set_ylabel("anomaly (%)")
ax.legend(loc="upper left", ncol=2, fontsize=5.0)
S.plabel(ax, "b", "Wind stress and surface eddy energy, ACC band")
# ---------------------------------------------------------------- c lag correlation
gr = g[1].subgridspec(2, 1, hspace=0.42, height_ratios=[1.0, 1.0])
ax = fig.add_subplot(gr[0]); S.zero(ax); ax.axvline(0, color="#9A9A9A", lw=0.5)
sel = (LC.scale == "12mo") & LC.detrended & (LC.driver == "tau_accmask")
band = LC[sel & (LC.response == "ACC_EKE")]
ax.fill_between(band.lag_months, -band.r95_pointwise, band.r95_pointwise, color="#EDEDED", lw=0, zorder=0)
gq = LS[(LS.scale == "12mo") & LS.detrended & (LS.driver == "tau_accmask") & (LS.response == "ACC_EKE")].r_crit95_global.iloc[0]
for s_ in (-1, 1): ax.axhline(s_ * gq, color="#777777", lw=0.5, ls=(0, (3, 2)))
for resp, col, ls, lab in (("ACC_T", CT, "-", "surface transport proxy"),
                           ("ACC_EKE", C["eke"], "-", "EKE"), ("ACC2s_EKE", C["two"], (0, (3, 1.6)), "EKE, two-satellite")):
    q = LC[sel & (LC.response == resp)]
    ax.plot(q.lag_months, q.r, color=col, ls=ls, lw=0.9, label=lab)
    s = LS[(LS.scale == "12mo") & LS.detrended & (LS.driver == "tau_accmask") & (LS.response == resp)].iloc[0]
    if True:
        ax.plot(s.peakpos_lag_months, s.peakpos_r, "o", ms=3, color=col, mfc=col if s.p_peakpos_global < 0.05 else "white", mew=0.6, zorder=5)
ax.set_xlim(-24, 48); ax.set_xticks(range(-24, 49, 12)); ax.set_ylim(-0.7, 0.8)
ax.set_xlabel("lag (months; positive: ocean lags wind)"); ax.set_ylabel("correlation with wind stress")
ax.legend(loc="upper right", ncol=2, fontsize=4.9)
S.plabel(ax, "c", "Lagged response, 12-month means, detrended")
# ---------------------------------------------------------------- d sensitivity
ax = fig.add_subplot(gr[1]); S.zero(ax)
VARS = [("T", "transport\nproxy", CT), ("KE", "KE", C["ke"]), ("EKE", "EKE", C["eke"])]
YL = (-0.6, 1.8)
for i, (v, lab, col) in enumerate(VARS):
    q = E[(E.sector == "ACC") & (E.tau == "accmask") & (E["var"] == v)]
    for k, (L, ver, mk, dx) in enumerate(((0, "detrended", "o", -0.22), (1, "detrended", "s", 0.0))):
        r = q[(q.lag_yr == L) & (q.version == ver)].iloc[0]
        ax.plot([i + dx] * 2, [max(r.ci_lo, YL[0]), min(r.ci_hi, YL[1])], color=col, lw=0.7)
        ax.plot(i + dx, r.elasticity, mk, color=col, ms=3.2, mfc=col if r.p < 0.05 else "white", mew=0.6)
    r = q.iloc[0]
    lo_, hi_ = r.trend_ratio_lo, r.trend_ratio_hi
    ax.plot([i + 0.22] * 2, [max(lo_, YL[0]), min(hi_, YL[1])], color=col, lw=0.7)
    ax.plot(i + 0.22, r.trend_ratio, "D", color=col, ms=3.0, mfc=col if lo_ > 0 else "white", mew=0.6)
    if hi_ > YL[1]: ax.annotate("", xy=(i + 0.22, YL[1]), xytext=(i + 0.22, YL[1] - 0.12), arrowprops=dict(arrowstyle="-|>", lw=0.5, color=col, mutation_scale=4))
    if r.trend_ratio > YL[1]: ax.text(i + 0.30, YL[1] - 0.05, f"{r.trend_ratio:.1f}", fontsize=4.8, color=col, va="top")
ax.plot([1.55, 2.45], [1, 1], color="#777777", lw=0.5, ls=(0, (3, 2)), zorder=0)
ax.text(1.56, 1.02, "EKE ~ τ", fontsize=4.8, color="#555555", ha="left", va="bottom")
ax.plot([1.55, 2.45], [0.7, 0.7], color="#A2508A", lw=0.5, ls=(0, (1, 1.5)), zorder=0)
ax.text(1.56, 0.72, "EKE ~ τ$^{0.7}$", fontsize=4.8, color="#A2508A", ha="left", va="bottom")
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l, _ in VARS], fontsize=5.3); ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.5, 2.5); ax.set_ylim(*YL); ax.set_ylabel("d ln X / d ln τ")
ax.legend(handles=[Line2D([], [], ls="", marker="o", color="#444444", ms=3, label="interannual, lag 0"),
                   Line2D([], [], ls="", marker="s", color="#444444", ms=3, label="interannual, lag 1 yr"),
                   Line2D([], [], ls="", marker="D", color="#444444", ms=3, label="ratio of trends"),
                   Line2D([], [], ls="", marker="o", color="#444444", mfc="white", ms=3, label="open: interval includes 0")],
          loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=4, fontsize=4.8)
S.plabel(ax, "d", "Sensitivity to the wind stress, annual means")
# ---------------------------------------------------------------- e, f intensification and displacement
D19 = DO / "dk19"; SS = pd.read_csv(D19 / "dk19_series.csv"); TR19 = pd.read_csv(D19 / "dk19_trends.csv").set_index("series")
B19 = pd.read_csv(D19 / "dk19_bins.csv"); BC = pd.read_csv(D19 / "dk19_bin_corr.csv").set_index("pair")
gb = gA[1].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.42)
ax = fig.add_subplot(gb[0]); S.zero(ax)
an19 = SS.groupby("year").mean(numeric_only=True); y19 = an19.index.values + 0.5
axw = ax.twinx(); axw.spines["right"].set_visible(True); axw.spines["right"].set_linewidth(0.5)
axw.plot(y19, an19.westerly_lat_km, "-s", color=CTAU, lw=0.6, ms=1.6, label="latitude of the westerly maximum (right)")
ax.plot(y19, an19.hp_km, "-o", color=C["pos"], lw=0.9, ms=2.0, label="half-power latitude, 900 tracked longitudes (left)")
xx = np.array([2001, 2026])
for c_, col, a_ in (("westerly_lat_km", CTAU, axw), ("hp_km", C["pos"], ax)):
    r_ = TR19.loc[c_]; a_.plot(xx, r_.km_dec / 10 * (xx - 2013.5), color=col, lw=0.7, ls=(0, (3, 1.6)) if r_.p >= 0.05 else "-")
S.ts_style(ax, (2000.6, 2025.6)); ax.set_ylabel("half-power anomaly (km, + north)", color=C["pos"]); ax.set_ylim(-25, 25)
axw.set_ylim(-500, 500); axw.set_ylabel("westerly-maximum anomaly (km)", color=CTAU, fontsize=5.5); axw.tick_params(axis="y", colors=CTAU, labelsize=5)
ax.set_zorder(axw.get_zorder() + 1); ax.patch.set_visible(False)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axw.get_legend_handles_labels()
r1, r2 = TR19.loc["westerly_lat_km"], TR19.loc["hp_km"]
ax.text(0.02, 0.03, f"westerlies {r1.km_dec:+.0f} km per decade ({r1.lo:+.0f} to {r1.hi:+.0f}); half-power {r2.km_dec:+.1f} ({r2.lo:+.1f} to {r2.hi:+.1f})",
        transform=ax.transAxes, fontsize=4.8, va="bottom")
ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=4.8)
S.plabel(ax, "e", "Position of the westerlies and of the energetic band")
ax = fig.add_subplot(gb[1]); S.zero(ax); S.zero(ax, "x")
b = B19.dropna(subset=["speed_pct_dec", "hp_km_dec"])
ax.errorbar(b.hp_km_dec, b.speed_pct_dec, xerr=[b.hp_km_dec - b.hp_lo, b.hp_hi - b.hp_km_dec], fmt="none", ecolor="#CFCFCF", lw=0.4, zorder=1)
ax.scatter(b.hp_km_dec, b.speed_pct_dec, s=10, color=C["pos"], zorder=3, lw=0)
r = BC.loc["speed vs half-power latitude"]
ax.text(0.02, 0.97, f"r = {r.r:.2f} (p < 0.001, circular shifts; {int(r.n)} bins)\njet-path latitude: r = {BC.loc['speed vs path latitude'].r:.2f}",
        transform=ax.transAxes, fontsize=4.9, va="top")
ax.set_xlabel("half-power latitude trend (km per decade, + north)"); ax.set_ylabel("jet-speed trend (% per decade)")
ax.set_xlim(-35, 35)
S.plabel(ax, "f", "Intensification and apparent displacement, 5° bins")
S.save(fig, FIGD, "Figure4"); plt.close(fig)
pd.set_option("display.width", 200)
print(E[(E.sector == "ACC") & (E.tau == "accmask")][["var", "lag_yr", "version", "elasticity", "ci_lo", "ci_hi", "p", "trend_ratio", "trend_ratio_lo", "trend_ratio_hi"]].round(3).to_string(index=False))
print(LS[(LS.scale == "12mo") & LS.detrended & (LS.driver == "tau_accmask") & LS.response.str.startswith("ACC")].round(3).to_string(index=False))
