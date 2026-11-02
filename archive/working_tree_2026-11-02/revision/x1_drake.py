"""Drake Passage (70-55W) version of x1: half-power latitude vs fixed-ADT-contour latitude, 2001-2025, on the retained
"""NOTE (30 Oct 2026, run in a cloud container from the Drake columns 90-35W of data/monthly_allsat extracted with extract_drake_cols.py;
run from the revision directory). Uses the 62 retained circumpolar tracks (>= 90 % rule) in 70-55W, NOT the Drake-rule envelopes of Fig. 2b/Table S9."""
circumpolar tracks (tracks_baseline.npz, keep=True) that fall in 70-55W. Same method as x1_contour_vs_energy.py / x1_stats.py."""
import numpy as np, pandas as pd, xarray as xr, sys, glob
sys.path.insert(0,'.'); import rstats as R
from scipy import stats
KM=111.195
T=np.load('outputs/a2/tracks_baseline.npz')
keep=T['keep']; j=T['j'][keep]; lat_hp=T['lat_hp'][:,keep].astype(float); lhc=T['lat_hp_clim'][keep]
madt=xr.open_dataset('../outputs/mean_adt.nc')['mean_adt']
lat_full=madt.latitude.values.astype(float); lon_full=madt.longitude.values.astype(float)
lon_tr=lon_full[j]
files=sorted(glob.glob('outputs/x1/drake_cols/*.npz')); assert len(files)==300
c0=np.load(files[0]); lat=c0['lat'].astype(float); lon_c=c0['lon'].astype(float)
assert np.allclose(lat,lat_full)
def run(lo,hi,label):
    sel=np.where((lon_tr>=lo)&(lon_tr<=hi))[0]
    jj=j[sel]; jc=np.array([np.argmin(np.abs(lon_c-lon_full[k])) for k in jj]); assert np.allclose(lon_c[jc],lon_full[jj])
    hp=lat_hp[:,sel]; lhc_s=lhc[sel]
    M=madt.values
    lev=np.array([np.interp(lhc_s[k],lat,M[:,jj[k]]) for k in range(sel.size)])
    W=3.0; adt_hp=np.full(hp.shape,np.nan); lat_c=np.full(hp.shape,np.nan)
    for m,f in enumerate(files):
        A=np.load(f)['adt'].astype(float)
        for k in range(sel.size):
            col=A[:,jc[k]]
            if np.isfinite(hp[m,k]):
                v=np.interp(hp[m,k],lat,col); adt_hp[m,k]=v if np.isfinite(v) else np.nan
            w=(lat>=lhc_s[k]-W)&(lat<=lhc_s[k]+W); la=lat[w]; c=col[w]-lev[k]; ok=np.isfinite(c)
            if ok.sum()<3: continue
            la,c=la[ok],c[ok]; s=np.where(np.diff(np.sign(c))!=0)[0]
            if s.size==0: continue
            x=la[s]+(0-c[s])*(la[s+1]-la[s])/(c[s+1]-c[s]); lat_c[m,k]=x[np.argmin(np.abs(x-lhc_s[k]))]
    t=2001+(np.arange(300)+0.5)/12; month=np.tile(np.arange(1,13),25); yr=np.repeat(np.arange(2001,2026),12)
    def series(a): an=a-np.nanmean(a,0); return np.nanmean(an,1)
    print(f'\n== {label}: {sel.size} tracks, lon {lon_tr[sel].min():.1f} to {lon_tr[sel].max():.1f}; contour found in {np.isfinite(lat_c).mean()*100:.0f} % of track-months, hp in {np.isfinite(hp).mean()*100:.0f} %')
    rows=[]
    for nm,a,sc,u in (('lat_hp',hp,KM,'km'),('lat_contour',lat_c,KM,'km'),('adt_hp',adt_hp,100,'cm')):
        y=series(a)*sc; r=R.trend_bcz(t,y,month=month,monthly=True,n_outer=300,n_bias=100,seed=0)
        rows.append(dict(region=label,series=nm,unit=u+' per decade',slope=10*r['slope'],lo=10*r['ci_lo'],hi=10*r['ci_hi'],p=r['p']))
        print(f"  {nm:12s} {10*r['slope']:+.1f} {u}/decade ({10*r['ci_lo']:+.1f} to {10*r['ci_hi']:+.1f}) p={r['p']:.3f}")
    # early/late period means
    e=(yr<=2015); l=(yr>=2016)
    for nm,a in (('hp',hp),('contour',lat_c)):
        an=a-np.nanmean(a,0); d=(np.nanmean(an[l],0)-np.nanmean(an[e],0))*KM
        print(f"  {nm:8s} late-minus-early mean position: {np.nanmean(d):+.1f} km (median {np.nanmedian(d):+.1f}; southward at {(d<0).mean()*100:.0f} % of tracks)")
    g=np.array([np.nanmedian(np.abs(np.gradient(M[:,k],lat*KM))) for k in jj])*100  # cm per km
    print(f"  median |dADT/dlat| {np.nanmedian(g)*10:.1f} cm per 10 km -> implied contour shift from adt_hp trend: {-(rows[2]['slope']/(np.nanmedian(g)*10))*10:+.1f} km/decade")
    return rows, dict(hp=hp,lat_c=lat_c,adt_hp=adt_hp,lon=lon_tr[sel],lhc=lhc_s,lev=lev)
out=[]; 
r1,D1=run(-70,-55,'Drake 70-55W'); out+=r1
r2,D2=run(-86,-40,'Drake region 86-40W'); out+=r2
pd.DataFrame(out).to_csv('x1_drake_trends.csv',index=False)
np.savez_compressed('x1_drake_tracks.npz',**{k+'_7055':v for k,v in D1.items()},**{k+'_8640':v for k,v in D2.items()})
