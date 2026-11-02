import sys, numpy as np, pandas as pd
sys.path.insert(0,'rev'); import rstats as R
for p in sys.argv[1:]:
    t=pd.read_csv(p+'.csv'); g=t.dropna(subset=['h_noseas_cm'])
    o=R.trend_bcz(g.x.values,g.h_noseas_cm.values,month=g.month.values)
    hf=g.h_f_cm.values; T=g.T_hogg_model_calibrated_Sv.values
    hf_tr=10*np.polyfit(g.x,hf,1)[0]; hi=T+1.11*hf
    o2=R.trend_bcz(g.x.values,hi,month=g.month.values)
    A=g.groupby('year').T_hogg_model_calibrated_Sv.mean()
    J=g[g.month<=3].groupby('year').h_noseas_cm.mean()
    oj=R.trend_bcz(J.index.values+0.5,J.values)
    print(f"{p}: h trend {10*o['slope']:+.2f} cm/dec [{10*o['ci_lo']:+.2f},{10*o['ci_hi']:+.2f}] p={o['p']:.3f}; "
          f"h_f linear trend {hf_tr:+.2f} cm/dec -> -1.11*h_f {-1.11*hf_tr:+.2f} Sv/dec; "
          f"high-freq part -3.69(h-hf) trend {10*o2['slope']:+.2f} Sv/dec p={o2['p']:.2f}; "
          f"T monthly sd {np.nanstd(T):.1f} Sv, annual-mean sd {A.std():.1f} Sv; "
          f"JFM-only h trend {10*oj['slope']:+.2f} cm/dec [{10*oj['ci_lo']:+.2f},{10*oj['ci_hi']:+.2f}] p={oj['p']:.3f} -> x(-1.11) {-1.11*10*oj['slope']:+.1f} Sv/dec")
a=pd.read_csv(sys.argv[1]+'.csv'); b=pd.read_csv(sys.argv[2]+'.csv')
print('r allsat vs twosat T:', a.T_hogg_model_calibrated_Sv.corr(b.T_hogg_model_calibrated_Sv).round(3))
print(a.groupby('year').T_hogg_model_calibrated_Sv.mean().round(1).to_dict())
