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
fig = plt.figure(figsize=(S.W2, 0.62 * S.W2))
g = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.25)
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
S.save(fig, FIGD, "Figure4"); plt.close(fig)
pd.set_option("display.width", 200)
print(E[(E.sector == "ACC") & (E.tau == "accmask")][["var", "lag_yr", "version", "elasticity", "ci_lo", "ci_hi", "p", "trend_ratio", "trend_ratio_lo", "trend_ratio_hi"]].round(3).to_string(index=False))
print(LS[(LS.scale == "12mo") & LS.detrended & (LS.driver == "tau_accmask") & LS.response.str.startswith("ACC")].round(3).to_string(index=False))
