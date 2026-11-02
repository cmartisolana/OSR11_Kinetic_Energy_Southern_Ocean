"""Figure D1: Drake Passage first -- altimetry checked against the LM Gould record (2005-2019) and extended to 2001-2025.
Usage: python fig_drake.py <revision dir> <drake outputs dir> <plan1019 outputs dir> <GV repo dir> <figure dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
REV, DO, PO, GV, FIGD = (Path(a) for a in sys.argv[1:6])
sys.path.insert(0, str(REV)); import style_nc as S
plt = S.plt; C = S.C
from matplotlib.lines import Line2D
m = sio.loadmat(GV/"Datasets"/"TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
tr = np.array(m.transp, float); blat = np.array(m.lat, float); dist = np.array(m.dist, float); tr[dist > 840, :] = np.nan
time = pd.to_datetime(np.array(m.time, float)-719529, unit="D"); nval = np.isfinite(tr).sum(0)
def ins_pos(a, b):
    out = np.full(tr.shape[1], np.nan)
    for k in range(tr.shape[1]):
        if nval[k] < 30: continue
        bb, qq = blat[:, k], tr[:, k]; s = np.isfinite(bb) & np.isfinite(qq) & (bb >= min(a, b)) & (bb <= max(a, b))
        if s.sum() < 3: continue
        o = np.argsort(-bb[s]); x, q = bb[s][o], np.clip(qq[s][o], 0, None); c = np.cumsum(q)
        if c[-1] > 0: out[k] = np.interp(0.5*c[-1], c, x)
    return out
ins = pd.DataFrame(dict(year=time.year, pos=ins_pos(-55.1, -59.33), n=nval))
ins = ins[ins.n >= 30]
split = pd.read_csv(PO/"d1"/"d1b_insitu_split.csv")
ekei = pd.read_csv(DO/"dk3"/"dk3_eke_transects_allsat.csv"); ekei = ekei[ekei.band == "north"]
ser = pd.read_csv(DO/"dk4"/"dk4_section_series.csv"); pos = pd.read_csv(DO/"dk2"/"dk2_positions.csv")
tr1 = pd.read_csv(DO/"dk1"/"dk1_trends.csv")
ov = (2005, 2019)
def rel(s, years):   # % anomaly relative to the 2005-2019 mean of that series
    mu = s[(years >= ov[0]) & (years <= ov[1])].mean(); return 100*(s/mu - 1)

fig = plt.figure(figsize=(S.W2, 0.78*S.W2))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], width_ratios=[1, 1, 1], wspace=0.55, hspace=0.62)
gs2 = gs[1, 1:].subgridspec(1, 2, width_ratios=[0.42, 1])
# a: mean velocity profile
ax = fig.add_subplot(gs[0, 0])
for p, col, lab in (("allsat", C["eke"], "altimetry, all-satellite"), ("twosat", C["two"], "altimetry, two-satellite")):
    pr = pd.read_csv(DO/"dk3"/f"dk3_velocity_profile_{p}.csv"); pr = pr[pr.uin_count >= 30]
    if p == "allsat":
        ax.plot(pr.latb, pr.uin_mean, color=C["obs"], lw=1.2, label="in situ ADCP, 70–118 m")
    ax.plot(pr.latb, pr.ualt_mean, color=col, lw=1.0, label=lab)
ax.set_xlim(-55, -63); S.zero(ax); ax.set_xlabel("Latitude (°)"); ax.set_ylabel("Mean cross-section velocity (cm s$^{-1}$)")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=1, fontsize=5.5); S.plabel(ax, "a", x=-0.3, dy=30)
# b, c, d: annual series, altimetry 2001-2025 vs in situ 2005-2019
def ann_alt(col):
    a = ser.groupby("year")[col].mean(); return a
panels = [("b", "EKE, north of 58.5°S\n(% of 2005–2019 mean)", ekei.groupby("year").eke_in.mean(), "allsat_EKE_north", "twosat_EKE_north", True),
          ("c", "Transport, north of 58.5°S\n(% of 2005–2019 mean)", split.groupby("year").U_north.mean(), "allsat_Tperm", "twosat_Tperm", True),
          ("d", "Position of the energetic band\n(latitude, °)", ins.groupby("year").pos.mean(), "allsat_front1_lat", "twosat_front1_lat", False)]
for k, (letter, ylab, insitu, ca, ct, relative) in enumerate(panels):
    ax = fig.add_subplot(gs[0, 1+k] if k < 2 else gs[1, 0])
    for c, col in ((ca, C["eke"] if letter == "b" else (C["tr"] if letter == "c" else C["pos"])), (ct, C["two"])):
        src = pos if letter == "d" else ser
        a = src.groupby("year")[c].mean()
        y = rel(a, a.index) if relative else a
        ax.plot(a.index, y, "-", color=col, lw=1.0, marker="o", ms=1.8)
    insitu = insitu.loc[ov[0]:ov[1]]
    y = rel(insitu, insitu.index) if relative else insitu
    ax.plot(insitu.index, y, "-o", color=C["obs"], lw=1.2, ms=2.4, zorder=4)
    ax.axvspan(ov[0]-0.5, ov[1]+0.5, color=C["light"], lw=0, zorder=0)
    if relative: S.zero(ax)
    else: ax.invert_yaxis() if False else None
    S.ts_style(ax, (2000.5, 2025.5)); ax.set_ylabel(ylab); S.plabel(ax, letter, dy=4)
handles = [Line2D([], [], color=C["obs"], marker="o", ms=2.4, label="in situ (LM Gould)"),
           Line2D([], [], color=C["ink"], marker="o", ms=1.8, label="altimetry, all-satellite (green, blue, purple)"),
           Line2D([], [], color=C["two"], marker="o", ms=1.8, label="altimetry, two-satellite")]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.66, 0.995), ncol=3, fontsize=5.6)
# e: trends forest (custom rows with region headers)
ax = fig.add_subplot(gs2[0, 1])
R_ = {"box": "Drake box (80–52°W, 63–54°S)", "up": "Upstream (100–80°W, ACC band)", "down": "Downstream (55–35°W, ACC band)"}
def g(name):
    r = tr1[(tr1.name == name) & (tr1.period == "2001-2025")].iloc[0]
    return dict(pct=r.pct, pct_lo=r.pct_lo, pct_hi=r.pct_hi, p=r.p)
def pt(r, yy, c):
    ax.plot([r["pct_lo"], r["pct_hi"]], [yy, yy], color=c, lw=0.8, solid_capstyle="butt")
    sig = r["p"] < 0.05
    ax.plot(r["pct"], yy, "o", color=c, mfc=c if sig else "white", mec=c, mew=0.8, ms=3.4, zorder=3)
y = 0.0; yt, yl = [], []
for reg, rl in R_.items():
    ax.text(-0.01, y + 0.75, rl, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.0, fontweight="bold")
    for var, vl in (("KE", "KE"), ("EKE", "EKE"), ("SPD", "along-stream speed")):
        pt(g(f"allsat_{reg}_{var}"), y + 0.14, C["eke"]); pt(g(f"twosat_{reg}_{var}"), y - 0.14, C["two"]); yt.append(y); yl.append(vl); y -= 1
    pt(g(f"wind_{reg}_TAUX"), y, C["wind"]); yt.append(y); yl.append("eastward wind stress"); y -= 1.6
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=5.8); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
ax.set_ylim(y + 0.9, 1.3); ax.axvline(0, color="#9A9A9A", lw=0.5, zorder=0); ax.set_xlim(-8, 28)
ax.set_xlabel("Trend 2001–2025 (% of mean per decade)")
ax.legend(handles=S.sig_handles({"all-satellite": C["eke"], "two-satellite": C["two"], "wind stress": C["wind"]}),
          loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, fontsize=5.5)
S.plabel(ax, "e", x=-0.36, dy=14)
S.save(fig, FIGD, "FigureD1_drake_validation_extension"); plt.close(fig)
