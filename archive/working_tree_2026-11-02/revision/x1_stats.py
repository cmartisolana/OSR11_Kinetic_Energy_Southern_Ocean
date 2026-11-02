import numpy as np, pandas as pd, sys
sys.path.insert(0,'.'); import rstats as R
from scipy import stats
Z=np.load('outputs/x1/x1_tracks.npz'); KM=111.195
lat_hp,adt_hp,lat_c,lon=Z['lat_hp'],Z['adt_hp'],Z['lat_c'],Z['lon']
def series(a):
    an=a-np.nanmean(a,0); return np.nanmean(an,1)
t=2001+(np.arange(300)+0.5)/12; month=np.tile(np.arange(1,13),25); yr=np.repeat(np.arange(2001,2026),12)
rows=[]
for nm,a,sc,u in (('lat_hp',lat_hp,KM,'km'),('lat_contour',lat_c,KM,'km'),('adt_hp',adt_hp,100,'cm')):
    y=series(a)*sc; r=R.trend_bcz(t,y,month=month,monthly=True,n_outer=300,n_bias=100,seed=0)
    rows.append(dict(series=nm,unit=u+' per decade',slope=10*r['slope'],lo=10*r['ci_lo'],hi=10*r['ci_hi'],p=r['p']))
    print(f"{nm:12s} {10*r['slope']:+.2f} {u}/decade ({10*r['ci_lo']:+.2f} to {10*r['ci_hi']:+.2f}) p={r['p']:.3f}", flush=True)
pd.DataFrame(rows).to_csv('outputs/x1/x1_circumpolar_trends.csv',index=False)
def tr(a):
    out=[]
    for k in range(a.shape[1]):
        y=a[:,k]; ok=np.isfinite(y); out.append(stats.linregress(t[ok],y[ok]).slope*10 if ok.sum()>100 else np.nan)
    return np.array(out)
th,tc,ta=tr(lat_hp)*KM,tr(lat_c)*KM,tr(adt_hp)*100
print('per-track median trend: hp %.1f km/dec, contour %.1f km/dec, adt_hp %.2f cm/dec'%(np.nanmedian(th),np.nanmedian(tc),np.nanmedian(ta)))
print('fraction of tracks southward: contour %.2f, hp %.2f'%((tc<0).mean(),(th<0).mean()))
ok=np.isfinite(th)&np.isfinite(tc); print('corr per-track trends hp vs contour r=%.2f'%np.corrcoef(th[ok],tc[ok])[0,1])
np.savez('outputs/x1/x1_pertrack.npz',th=th,tc=tc,ta=ta,lon=lon)
sec=np.where((lon>=-70)&(lon<20),'Atlantic',np.where((lon>=20)&(lon<147),'Indian','Pacific'))
for s in ('Atlantic','Indian','Pacific'):
    m=sec==s; line=f'  {s:9s} n={m.sum():3d}'
    for nm,a,sc in (('hp',lat_hp,KM),('contour',lat_c,KM),('adt_hp',adt_hp,100)):
        y=series(a[:,m])*sc; ann=pd.Series(y).groupby(yr).mean(); rr=stats.linregress(ann.index,ann.values)
        line+=f' | {nm} {10*rr.slope:+.2f} p={rr.pvalue:.2f}'
    print(line)
