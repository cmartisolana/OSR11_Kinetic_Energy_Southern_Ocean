"""Main Fig. 2 (v1 style): position of the ACC kinetic-energy envelope, 2001-2025.
(a) ACC-mean half-power latitude over the matched fixed sample (raw monthly, 12-month running mean,
    deseasonalised Theil-Sen trend); (b) trends by sector and ACC (95 % CI, sieve p);
(c) trends in 5° longitude bins (95 % CI) with Chambers (2018) regional estimates;
(d) detection frequency by longitude.
Usage: python fig2_main.py   (env OSR11_CACHE must point to the KE cache)
"""
import numpy as np, pandas as pd
from figstyle import *  # noqa: F401,F403
import a2_run as A

KM = 111.195
d = np.load(OUT / "a2" / "tracks_baseline.npz", allow_pickle=True)
lon = np.load(A.CACHE / "lon.npy")
keep = d["keep"] & d["acc"]
lat_hp = d["lat_hp"][:, keep].astype(float)
sec = A.sector(lon[d["j"][keep]])
x = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25)
S = pd.read_csv(OUT / "a2" / "a2_summary.csv"); S = S[S.set == "baseline"].set_index("region")
S = retest_override(S, lambda g, r: f"baseline|{g}|lat_hp", {"m_lat_hp_p": "p_new", "m_lat_hp_lo": "ci_lo", "m_lat_hp_hi": "ci_hi"})
LT = pd.read_csv(OUT / "a2" / "a2f_lon_trends.csv")
LT = retest_override(LT, lambda i, r: f"lonbin|{int(r.lon0)}", {"p": "p_new", "lo": "ci_lo", "hi": "ci_hi"})
DF = pd.read_csv(OUT / "a2" / "detection_frequency_lon.csv")

fig = plt.figure(figsize=(10.4, 9.0))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.85, 0.42], width_ratios=[1.75, 1], hspace=0.42, wspace=0.22)
ax = fig.add_subplot(gs[0, 0])
y = np.nanmean(lat_hp, 1)
r = S.loc["ACC"]
yd = pd.Series(y); yd = (yd - yd.groupby(month).transform("mean") + yd.mean()).values
ts_panel(ax, x, y, STYLE["lat"], slope=r.m_lat_hp_slope, intercept=theil_line(x, yd, r.m_lat_hp_slope))
ax.text(0.98, 0.04, f"{r.m_lat_hp_slope * KM * 10:+.1f} km dec$^{{-1}}$ [{r.m_lat_hp_lo * KM * 10:+.1f}, {r.m_lat_hp_hi * KM * 10:+.1f}]\n"
        f"p={r.m_lat_hp_p:.2f}   ({int(r.ntracks)} longitudes)", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=STYLE["annot_size"] - 1, color=STYLE["neutral"])
ax.set_xlabel("Year"); ax.set_ylabel("Latitude (°)"); base_timeseries_style(ax); panel_label(ax, "a", x=-0.1, y=1.08)

ax = fig.add_subplot(gs[0, 1])
regs = ["Indian", "Pacific", "Atlantic", "ACC"]
s_ = np.array([S.loc[g].m_lat_hp_slope for g in regs]) * KM * 10
lo = np.array([S.loc[g].m_lat_hp_lo for g in regs]) * KM * 10; hi = np.array([S.loc[g].m_lat_hp_hi for g in regs]) * KM * 10
pv = np.array([S.loc[g].m_lat_hp_p for g in regs])
ax.bar(np.arange(4), s_, width=0.62, color=[STYLE["lat"] if p < 0.05 else STYLE["ns"] for p in pv], edgecolor="white", zorder=3)
ax.errorbar(np.arange(4), s_, yerr=[s_ - lo, hi - s_], fmt="none", ecolor="#1F2933", lw=1.0, capsize=2.8, capthick=1.0, zorder=4)
for i in range(4):
    ax.text(i, hi[i] + 0.3, stars(pv[i]), ha="center", va="bottom", fontsize=10)
ax.axhline(0, color="#495057", lw=0.8); ax.set_xticks(np.arange(4)); ax.set_xticklabels(regs)
ax.set_ylabel("Trend (km dec$^{-1}$, + north)"); base_timeseries_style(ax); panel_label(ax, "b", x=-0.2, y=1.08)
ax.spines["bottom"].set_visible(False)

ax = fig.add_subplot(gs[1, :])
xc = LT.lon0.values + 2.5
s2 = LT.slope_deg_yr.values * KM * 10; lo2 = LT.lo.values * KM * 10; hi2 = LT.hi.values * KM * 10
ax.axhline(0, color="#495057", lw=0.8)
ax.errorbar(xc, s2, yerr=[s2 - lo2, hi2 - s2], fmt="o", ms=4, color=STYLE["lat"], ecolor=STYLE["lat"], elinewidth=1.0,
            capsize=2.0, label="this study, 5° bins (95 % CI)", zorder=3)
ax.hlines(-29, 90, 170, color="#E07A5F", lw=3, label="Chambers (2018), along-track, 1993–2016", zorder=2)
ax.hlines(46, -160, -90, color="#E07A5F", lw=3, zorder=2)
ax.hlines(-1.3, 90, 170, color="#2F3E46", lw=3, label="this study, same regions", zorder=2)
ax.hlines(3.7, -160, -90, color="#2F3E46", lw=3, zorder=2)
ax.set_xlim(-180, 180); ax.set_xticks(range(-180, 181, 60))
ax.set_xticklabels([f"{abs(v)}°{'W' if v < 0 else ('E' if v > 0 else '')}" for v in range(-180, 181, 60)])
ax.set_ylabel("Trend (km dec$^{-1}$, + north)"); base_timeseries_style(ax); panel_label(ax, "c", x=-0.06, y=1.1)
ax.legend(loc="lower left", fontsize=8, frameon=False, ncol=3)

ax = fig.add_subplot(gs[2, :])
f = DF.freq.rolling(8, center=True, min_periods=1).mean()
ax.fill_between(DF.lon, 0, f, color=STYLE["lat"], alpha=0.35, lw=0)
ax.plot(DF.lon, f, color=STYLE["lat"], lw=0.8)
ax.set_xlim(-180, 180); ax.set_ylim(0, 1.05); ax.set_xticks(range(-180, 181, 60))
ax.set_xticklabels([f"{abs(v)}°{'W' if v < 0 else ('E' if v > 0 else '')}" for v in range(-180, 181, 60)])
ax.set_ylabel("Detection\nfrequency"); base_timeseries_style(ax); panel_label(ax, "d", x=-0.06, y=1.2)
fig.savefig(FIG / "Fig2_frontal_position.png"); fig.savefig(FIG / "Fig2_frontal_position.pdf")
print("saved")
