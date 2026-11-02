"""Figures for the 19 Oct revision-plan analyses (style_nc.py: 180 mm, 5-7 pt, no statistics inside panels).
  FigureP1_drake_coverage   Drake Passage: coverage of the altimetric proxy and coverage-matched comparison (d1)
  FigureP2_comparison       EKE trends against earlier estimates; product and period sensitivity (m1, m2)
  FigureP3_halfpower_gates  Half-power meander analysis: coverage gate and synthetic tests (h1)
Usage: python fig_plan1019.py <revision dir> <plan1019 outputs dir> <GV2023 repo dir> <figure out dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio

REV, OUTS, GV, FIGD = (Path(a) for a in sys.argv[1:5])
sys.path.insert(0, str(REV))
import style_nc as S
plt = S.plt; C = S.C

# ------------------------------------------------------------------ P1 Drake
sec = pd.read_csv(OUTS/"d1"/"d1_section_profile.csv")
lim = sec.lat[sec.perm].min()
m = sio.loadmat(GV/"Datasets"/"TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
tr = np.array(m.transp, float); blat = np.array(m.lat, float); dist = np.array(m.dist, float)
tr[dist > 840, :] = np.nan
time = pd.to_datetime(np.array(m.time, float) - 719529, unit="D")
nvalid = np.isfinite(tr).sum(0)
north = blat > lim
obs = pd.DataFrame(dict(year=time.year, month=time.month, U_tot=np.nansum(tr, 0)/1e6,
                        U_north=np.nansum(np.where(north, tr, 0), 0)/1e6, n=nvalid))
obs = obs[obs.n >= 30]; obs["U_south"] = obs.U_tot - obs.U_north
alt = pd.read_csv(REV/"outputs"/"a6"/"a6c_transport_fixed_series.csv")[["year", "month", "Tperm_v2"]]
mm = obs.merge(alt, on=["year", "month"])
ann = obs.groupby("year")[["U_tot", "U_north", "U_south"]].mean().loc[2005:2019]
annA = alt.groupby("year").Tperm_v2.mean().loc[2005:2019]
z = lambda s: (s - s.mean())/s.std()
bin_lat = np.nanmean(np.where(np.isfinite(tr), blat, np.nan), 1); bin_T = np.nanmean(tr, 1)/1e6
ok = np.isfinite(bin_lat) & np.isfinite(bin_T)

fig = plt.figure(figsize=(S.W2, 0.52*S.W2))
gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.2, 0.95], height_ratios=[1, 1], wspace=0.85, hspace=0.6)
ax = fig.add_subplot(gs[:, 0])
ax.axvspan(sec.lat.min()-0.1, lim, color=C["ice"], lw=0, zorder=0)
ax.bar(bin_lat[ok], bin_T[ok], width=0.2, color=C["obs"], lw=0, label="in situ, 0–760 m (Sv per 25 km)")
ax.set_xlim(-55.0, -62.8); ax.set_xlabel("Latitude (°)"); ax.set_ylabel("In situ transport (Sv per 25-km bin)")
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True); ax2.spines["top"].set_visible(False)
ax2.plot(sec.lat, sec.vbar_cm_s, color=C["tr"], lw=1.0)
ax2.plot(sec.lat[sec.perm], sec.vbar_cm_s[sec.perm], "o", color=C["tr"], ms=1.6, mew=0)
ax2.set_ylabel("Surface cross-section velocity (cm s$^{-1}$)", color=C["tr"]); ax2.tick_params(axis="y", colors=C["tr"])
ax2.set_ylim(-5, 50); ax.set_ylim(-5*5.8/50, 5.8); S.zero(ax)
ax.axvline(lim, color=C["grey"], lw=0.6, ls="--")
ax.text(lim - 0.12, ax.get_ylim()[1]*0.97, "ice-affected in\nat least one month", fontsize=5.5, va="top", ha="left", color="#4A6A86")
S.plabel(ax, "a")
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=C["obs"], label="in situ, 0–760 m"), Line2D([], [], color=C["tr"], marker="o", ms=1.6, label="altimetry, mean (dots: fixed set)")],
          loc="lower left", bbox_to_anchor=(0, 1.04), ncol=1)

for k, (part, lab, letter) in enumerate([("U_north", "in situ north of 58.5°S vs altimetric proxy", "b"), ("U_tot", "in situ full section vs altimetric proxy", "c")]):
    a = fig.add_subplot(gs[k, 1])
    a.plot(ann.index, z(ann[part]), "-o", color=C["obs"], ms=2.2, lw=0.9)
    a.plot(annA.index, z(annA), "-o", color=C["tr"], ms=2.2, lw=0.9)
    S.zero(a); a.set_xlim(2004.5, 2019.5); a.set_xticks([2005, 2010, 2015]); a.set_ylabel("Annual mean\n(s.d. units)")
    a.set_ylim(-2.6, 2.6)
    if k == 0:
        a.legend(handles=[Line2D([], [], color=C["obs"], marker="o", ms=2.2, label="in situ, 0–760 m"),
                          Line2D([], [], color=C["tr"], marker="o", ms=2.2, label="altimetric proxy (fixed ice-free points)")],
                 loc="lower left", bbox_to_anchor=(0, 1.13), ncol=2, fontsize=5.5)
    S.plabel(a, letter, lab)
a = fig.add_subplot(gs[:, 2])
zt = z(mm.Tperm_v2)
a.plot(zt, z(mm.U_south), "o", color=C["grey"], ms=1.8, mew=0, alpha=0.6, label="south of 58.5°S")
a.plot(zt, z(mm.U_north), "o", color=C["obs"], ms=1.8, mew=0, alpha=0.8, label="north of 58.5°S")
a.set_xlabel("Altimetric proxy (s.d. units)"); a.set_ylabel("In situ transport per transect\n(s.d. units)")
a.set_aspect("equal"); a.set_xlim(-3.5, 3.5); a.set_ylim(-3.5, 3.5); S.zero(a); S.zero(a, "x")
a.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.5, handletextpad=0.1); S.plabel(a, "d", dy=14)
S.save(fig, FIGD, "FigureP1_drake_coverage"); plt.close(fig)

# ------------------------------------------------------------------ P2 comparison
t1 = pd.read_csv(OUTS/"m1"/"m1_trends.csv")
def row(prod, var, per):
    r = t1[(t1["product"] == prod) & (t1.variable == var) & (t1.period == per)].iloc[0]
    if prod.startswith("diff"):
        return dict(pct=r.slope_per_decade, pct_lo=r.ci_lo, pct_hi=r.ci_hi, p=r.p)
    return dict(pct=r.pct_per_decade, pct_lo=r.pct_ci_lo, pct_hi=r.pct_ci_hi, p=r.p)
REG = [("Hogg_Indian", "Indian sector box"), ("Hogg_Pacific", "Pacific sector box"), ("Hogg_Atlantic", "Atlantic sector box"),
       ("box35_65", "35–65°S"), ("ACC", "ACC band")]
fig = plt.figure(figsize=(S.W2, 0.62*S.W2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 0.9, 1.15], wspace=0.6)
ax = fig.add_subplot(gs[0])
ref = {"Hogg_Indian": ("Hogg et al. (2015), from 1993", 6.0, 1.6), "Hogg_Pacific": ("Hogg et al. (2015), from 1993", 6.4, 1.8),
       "box35_65": ("Martínez-Moreno et al. (2021), 1993–2020", 5.1, np.nan)}
y = 0.0; yt, yl, seps = [], [], []
def pt(ax, r, yy, c, mk="o"):
    ax.plot([r["pct_lo"], r["pct_hi"]], [yy, yy], color=c, lw=0.8, solid_capstyle="butt")
    sig = np.isfinite(r["p"]) and r["p"] < 0.05
    ax.plot(r["pct"], yy, mk, color=c, mfc=c if sig else "white", mec=c, mew=0.8, ms=3.4, zorder=3)
for i, (v, lab) in enumerate(REG):
    ax.text(-0.02, y + 0.55, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.3, fontweight="bold")
    for per in ("2001–2020", "2001–2025"):
        p_ = per.replace("–", "-")
        pt(ax, row("allsat", f"{v}_EKE", p_), y + 0.13, C["eke"]); pt(ax, row("twosat", f"{v}_EKE", p_), y - 0.13, C["two"])
        yt.append(y); yl.append(per); y -= 1
    if v in ref:
        nm, c0, e0 = ref[v]
        pt(ax, dict(pct=c0, pct_lo=c0 - (e0 if np.isfinite(e0) else 0), pct_hi=c0 + (e0 if np.isfinite(e0) else 0), p=0.01), y, C["ink"], "D")
        yt.append(y); yl.append(nm); y -= 1
    seps.append(y + 0.5); y -= 0.6
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=5.6); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
for sp in seps[:-1]:
    ax.axhline(sp - 0.3, color=C["light"], lw=0.5)
ax.set_ylim(y + 0.6, 1.2); ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0); ax.set_xlim(-12, 22); ax.set_xlabel("EKE trend (% per decade)")
h = S.sig_handles({"all-satellite": C["eke"], "two-satellite": C["two"]})
h.insert(2, Line2D([], [], ls="-", lw=0.8, marker="D", color=C["ink"], ms=3.2, label="earlier estimate (as published)"))
ax.legend(handles=h, loc="lower left", bbox_to_anchor=(-0.05, 1.01), ncol=2)
S.plabel(ax, "a", x=-0.62, dy=26)
ax = fig.add_subplot(gs[1])
rows = [(lab, {"2001–2020": row("diff_all-two", f"{v}_EKE", "2001-2020"), "2001–2025": row("diff_all-two", f"{v}_EKE", "2001-2025")}) for v, lab in REG]
S.forest(ax, rows, {"2001–2020": C["grey"], "2001–2025": C["ink"]}, xlabel="All − two-satellite\n(percentage points per decade)", xlim=(-15, 20))
ax.legend(handles=S.sig_handles({"2001–2020": C["grey"], "2001–2025": C["ink"]}), loc="lower left", bbox_to_anchor=(-0.05, 1.01), ncol=1)
S.plabel(ax, "b", x=-0.6, dy=26)
ax = fig.add_subplot(gs[2])
for prod, col in (("allsat", C["eke"]), ("twosat", C["two"])):
    s = pd.read_csv(OUTS/"m2"/f"m2_series_{prod}.csv")
    y = 100*(s.EKE_PJm/s.EKE_PJm.mean() - 1)
    ax.plot(s.x, y, color=col, lw=0.4, alpha=0.45)
    ax.plot(s.x, 100*(s.EKE_12mo/s.EKE_PJm.mean() - 1), color=col, lw=1.2, label="all-satellite" if prod == "allsat" else "two-satellite")
ax.axvspan(2020.0, 2026.0, color=C["light"], lw=0, zorder=0)
S.zero(ax); S.ts_style(ax); ax.set_ylabel("EKE integrated over 35–65°S\n(% of 2001–2025 mean)")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=2); S.plabel(ax, "c", dy=14)
S.save(fig, FIGD, "FigureP2_comparison"); plt.close(fig)

# ------------------------------------------------------------------ P3 half-power gates
g = pd.read_csv(OUTS/"h1"/"h1_gate.csv"); syn = pd.read_csv(OUTS/"h1"/"h1_synthetic.csv")
names = {"SEIR150": "SE Indian Ridge (140–150°E)", "KERG": "Kerguelen (76–84°E)", "PAR": "Pacific–Antarctic Ridge", "DRAKE": "Drake Passage (70–55°W)"}
fig = plt.figure(figsize=(S.W2, 0.40*S.W2))
gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 1, 1, 1], wspace=0.75)
ax = fig.add_subplot(gs[0])
x = np.arange(len(g))
ax.barh(x, 100*g.frac_usable, color=[C["pos"] if p else C["light"] for p in g.G1_pass], edgecolor=C["pos"], lw=0.5, height=0.6)
ax.axvline(90, color=C["grey"], lw=0.6, ls="--")
ax.set_yticks(x); ax.set_yticklabels([names[s_].replace(" (", "\n(") for s_ in g.site], fontsize=5.5); ax.invert_yaxis()
ax.set_xlabel("Longitudes with a continuous\nhalf-power branch (%)"); ax.tick_params(axis="y", length=0)
ax.set_xlim(0, 100); S.plabel(ax, "a", x=-0.75, dy=4)
SC = [("N", "no change"), ("T", "shift 25 km"), ("F", "flexing 25 km"), ("U", "uniform ×1.2"), ("A", "one flank ×1.2")]
for k, (col, lab, letter) in enumerate([("meanlat_km_s150", "Change in mean position\nover 25 yr (km)", "b"),
                                         ("amp_km_s150", "Change in meander\namplitude over 25 yr (km)", "c"),
                                         ("sinuosity_s150", "Change in sinuosity\nover 25 yr", "d")]):
    a = fig.add_subplot(gs[k+1]); s = syn[syn.site == "SEIR150"]
    data = [s[s.scenario == c][col].dropna().values for c, _ in SC]
    bp = a.boxplot(data, positions=range(len(SC)), widths=0.55, showfliers=False, patch_artist=True,
                   medianprops=dict(color=C["ink"], lw=0.8), boxprops=dict(facecolor=C["light"], edgecolor=C["grey"], lw=0.5),
                   whiskerprops=dict(color=C["grey"], lw=0.5), capprops=dict(color=C["grey"], lw=0.5))
    for i, d in enumerate(data):
        a.plot(np.full(d.size, i) + np.random.default_rng(i).uniform(-0.15, 0.15, d.size), d, "o", ms=1.2, mew=0, color=C["pos"], alpha=0.6)
    if col == "meanlat_km_s150":
        a.plot([0.7, 1.3], [25, 25], color=C["pos"], lw=0.8, ls="--")
    if col == "amp_km_s150":
        a.plot([1.7, 2.3], [25, 25], color=C["pos"], lw=0.8, ls="--")
    S.zero(a); a.set_xticks(range(len(SC))); a.set_xticklabels([l for _, l in SC], fontsize=5.5, rotation=40, ha="right", rotation_mode="anchor"); a.set_ylabel(lab)
    S.plabel(a, letter, dy=4)
S.save(fig, FIGD, "FigureP3_halfpower_gates"); plt.close(fig)
