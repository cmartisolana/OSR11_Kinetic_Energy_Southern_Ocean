import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[1])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import sys, numpy as np, pandas as pd
from pathlib import Path
REV=Path(_REV); sys.path.insert(0,str(REV)); import style_nc as S
plt=S.plt; C=S.C; DO=Path(_DK_OUT)
K5=pd.read_csv(DO/'dk5'/'dk5_trends.csv')
def prof(coord,src,var,per): return K5[(K5.coord==coord)&(K5.src==src)&(K5['var']==var)&(K5.period==per)].sort_values('bin')
REC=((("insitu","2005-2019"),C["obs"],"in situ ADCP (70–118 m), Oct 2005–Apr 2019",3.0),(("allsat","2005-2019"),C["eke"],"altimetry, same period",2.4),(("allsat","2001-2025"),C["ink"],"altimetry, 2001–2025",2.4))
fig=plt.figure(figsize=(S.W1*1.4,S.W1*0.9))
ax = fig.add_subplot(111); ax2 = ax.twinx()
mu = prof("passage", "allsat", "u", "2001-2025")
ax2.fill_between(mu.bin, 0, mu["mean"], color="#EDEDED", lw=0, zorder=0); ax2.set_ylim(0, 45)
ax2.set_ylabel("Mean velocity (cm s$^{-1}$)", color="#8A8A8A", labelpad=1); ax2.spines["right"].set_visible(True); ax2.tick_params(colors="#8A8A8A")
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.axvspan(408, 850, facecolor="none", hatch="....", edgecolor="#C9DDF0", lw=0, zorder=0)
for (src, per), col, lab, ms in REC:
    q = prof("passage", src, "u", per); sg = q.p.round(3) < 0.05
    ax.plot(q.bin, 10 * q.slope, "-", color=col, lw=0.9, label=lab)
    ax.plot(q.bin[sg], 10 * q.slope[sg], "o", color=col, ms=ms); ax.plot(q.bin[~sg], 10 * q.slope[~sg], "o", mfc="white", color=col, ms=ms, mew=0.6)
S.zero(ax); ax.set_xlim(0, 850); ax.set_xlabel("Distance across the passage from the north (km)")
ax.set_ylabel("Velocity trend (cm s$^{-1}$ per decade)")
for x_, lab in ((40, "SAF"), (400, "PF"), (790, "SACCF")):
    ax.text(x_, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.5, color="#555555", fontweight="bold")
ax.annotate("faster between\nSAF and PF", xy=(215, 4.6), xytext=(300, 8.5), fontsize=5.2, color="#333333",
            arrowprops=dict(arrowstyle="-", lw=0.4, color="#555555"))
ax.legend(loc="upper left", bbox_to_anchor=(0, -0.18), ncol=3, fontsize=5.0, columnspacing=0.8)
fig.savefig('figS/new_S12.png',dpi=300,bbox_inches='tight',pad_inches=0.02); print('ok')
