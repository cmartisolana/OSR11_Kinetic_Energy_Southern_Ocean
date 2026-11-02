import sys, numpy as np, pandas as pd
from pathlib import Path
REV=Path('/mnt/user-data/uploads/repository/revision'); sys.path.insert(0,str(REV)); import style_nc as S
plt=S.plt; C=S.C; DO=Path('/home/claude/work/drake/out')
K5=pd.read_csv(DO/'dk5'/'dk5_trends.csv')
def prof(coord,src,var,per): return K5[(K5.coord==coord)&(K5.src==src)&(K5['var']==var)&(K5.period==per)].sort_values('bin')
fig=plt.figure(figsize=(S.W1*1.4,S.W1*0.9))
ax = fig.add_subplot(111); ax2 = ax.twinx()
me = prof("passage", "allsat", "eke", "2001-2025")
ax2.fill_between(me.bin, 0, me["mean"], color="#EDEDED", lw=0, zorder=0); ax2.set_ylim(0, 700)
ax2.set_ylabel("Mean surface EKE (cm$^2$ s$^{-2}$)", color="#8A8A8A", labelpad=1); ax2.spines["right"].set_visible(True); ax2.tick_params(colors="#8A8A8A")
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.axvspan(408, 850, facecolor="none", hatch="....", edgecolor="#C9DDF0", lw=0, zorder=0)
q = prof("passage", "allsat", "eke", "2001-2025")
ax.fill_between(q.bin, 1000 * q.ci_lo / q["mean"], 1000 * q.ci_hi / q["mean"], color=C["eke"], alpha=0.15, lw=0, zorder=1)
INS = pd.read_csv(DO / "dk7" / "dk7_insitu_eke_passage.csv")
curves = [(INS.bin, INS.pct, INS.p, C["obs"], "in situ ADCP, 0–780 m, Oct 2005–Apr 2019", "-", 3.0)]
for per, ls, lab in (("2005-2019", (0, (3, 1.5)), "altimetry, Oct 2005–Apr 2019"), ("2001-2025", "-", "altimetry, 2001–2025 (shading: 95 % interval)")):
    q = prof("passage", "allsat", "eke", per); curves.append((q.bin, 1000 * q.slope / q["mean"], q.p, C["eke"], lab, ls, 2.4))
for x, y, pp, col, lab, ls, ms in curves:
    x, y, pp = np.asarray(x), np.asarray(y), np.asarray(pp); sg = np.round(pp, 3) < 0.05
    ax.plot(x, y, color=col, lw=0.9, ls=ls, label=lab, zorder=3)
    ax.plot(x[sg], y[sg], "o", color=col, ms=ms, zorder=4); ax.plot(x[~sg], y[~sg], "o", mfc="white", color=col, ms=ms, mew=0.6, zorder=4)
S.zero(ax); ax.set_xlim(0, 850); ax.set_ylim(-60, 60); ax.set_xlabel("Distance across the passage from the north (km)")
ax.set_ylabel("EKE trend (% of mean per decade)")
for x_, lab in ((40, "SAF"), (400, "PF"), (790, "SACCF")):
    ax.text(x_, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.5, color="#555555", fontweight="bold")
ax.legend(loc="lower left", fontsize=4.8, ncol=1, framealpha=0.8)
fig.savefig('figS/new_S6.png',dpi=300,bbox_inches='tight',pad_inches=0.02); print('ok')
