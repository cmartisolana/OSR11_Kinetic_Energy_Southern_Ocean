"""D1: Drake Passage partial-section label and coverage-matched comparison with the in situ series.
(i) Share of the full-line mean surface proxy and of the mean ADT drop carried by the T_perm points
    (53 of 121, north of ~58.47S) on the v2 line.
(ii) In situ TOTAL 0-760 m ADCP transport (Gutierrez-Villanueva et al. 2023, per 25-km bin) restricted to
    bins north of the T_perm southern limit, i.e. the same latitudinal coverage as T_perm; annual-mean
    trends 2005-2019 (trend_bcz, annual AR(1)) and agreement with T_perm (same months, annual means).
Usage: python d1_drake_coverage.py <revision dir> <GV repo dir> <out dir>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, scipy.io as sio
REV, GV, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
sys.path.insert(0, str(REV)); import rstats as R
OUT.mkdir(parents=True, exist_ok=True)
R_E = 6.371e6
lon0, lat0, lon1, lat1, n = -65.0, -55.0, -60.0, -63.0, 121
tt = np.linspace(0, 1, n); lo = lon0 + tt*(lon1-lon0); la = lat0 + tt*(lat1-lat0)
latm = np.radians(0.5*(la[:-1]+la[1:]))
seg = np.hypot(R_E*np.cos(latm)*np.radians(np.diff(lo)), R_E*np.radians(np.diff(la)))
dsx = np.r_[seg[0]/2, 0.5*(seg[:-1]+seg[1:]), seg[-1]/2]
latr = np.radians(0.5*(lat0+lat1)); DX = R_E*np.cos(latr)*np.radians(lon1-lon0); DY = R_E*np.radians(lat1-lat0)
L = np.hypot(DX, DY); nx, ny = DY/L, -DX/L
if nx < 0: nx, ny = -nx, -ny
files = [REV/"data"/"monthly_allsat"/f"{y}"/f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
VN, ICE, ADT = [], [], []
pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
for f in files:
    with xr.open_dataset(f) as ds:
        d = ds.sel(longitude=slice(-70, -55), latitude=slice(-65, -53)).load()
    VN.append(d.u_mean.interp(**pts).values*nx + d.v_mean.interp(**pts).values*ny)
    ICE.append(d.ice_frac.interp(**pts).values); ADT.append(d.adt_mean.interp(**pts).values)
VN, ICE, ADT = map(np.array, (VN, ICE, ADT))
perm = (np.nan_to_num(ICE).max(0) <= 0.10) & np.all(np.isfinite(VN), 0)
south_lim = la[perm].min()
vbar = np.nanmean(VN, 0); abar = np.nanmean(ADT, 0)
share_T = np.sum(vbar[perm]*dsx[perm]) / np.sum(vbar*dsx)
iN, iS = 0, n-1; ilim = np.where(perm)[0].max()
share_dADT = (abar[iN]-abar[ilim]) / (abar[iN]-abar[iS])
# months each southern point is ice-affected
ice_months = (np.nan_to_num(ICE) > 0.10).sum(0)
sec = pd.DataFrame(dict(lat=la, lon=lo, perm=perm, vbar_cm_s=100*vbar, adt_mean_m=abar, months_ice_gt10=ice_months))
sec.to_csv(OUT/"d1_section_profile.csv", index=False)
summ = dict(n_points=n, n_perm=int(perm.sum()), south_limit_lat=float(south_lim),
            share_mean_surface_proxy=float(share_T), share_mean_dADT=float(share_dADT),
            dADT_full_m=float(abar[iN]-abar[iS]), dADT_perm_m=float(abar[iN]-abar[ilim]))
# in situ, per-bin total transport
m = sio.loadmat(GV/"Datasets"/"TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
tr = np.array(m.transp, float); blat = np.array(m.lat, float); dist = np.array(m.dist, float)
tr[dist > 840, :] = np.nan
def dn(x): return pd.Timestamp("0000-12-31") if False else (pd.to_datetime(x-719529, unit="D"))
time = pd.to_datetime(np.array(m.time, float) - 719529, unit="D")
north = blat > south_lim
U_tot = np.nansum(tr, 0)/1e6
U_n = np.nansum(np.where(north, tr, 0), 0)/1e6
U_s = U_tot - U_n
obs = pd.DataFrame(dict(date=time, U_tot=U_tot, U_north=U_n, U_south=U_s,
                        frac_bins_north=np.nanmean(north & np.isfinite(tr), 0)))
obs = obs[np.isfinite(obs.U_tot) & (obs.U_tot > 0)]
obs["year"] = obs.date.dt.year; obs["month"] = obs.date.dt.month
summ.update(insitu_n_transects=len(obs), insitu_mean_total_Sv=float(obs.U_tot.mean()),
            insitu_mean_north_Sv=float(obs.U_north.mean()), insitu_share_north=float(obs.U_north.mean()/obs.U_tot.mean()))
# altimetric T_perm / T_all series
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
t["Tperm"] = np.nansum(np.where(perm[None, :], VN, 0)*dsx, 1)
t["Tall"] = np.nansum(np.where(np.all(np.isfinite(VN), 0)[None, :], VN, 0)*dsx, 1)
mm = obs.merge(t, on=["year", "month"], how="left")
ann = obs.groupby("year")[["U_tot", "U_north", "U_south"]].mean()
ann = ann.loc[2005:2019]
annA = t.groupby("year")[["Tperm", "Tall"]].mean().loc[2005:2019]
rows = []
for c in ["U_tot", "U_north", "U_south"]:
    o = R.trend_bcz(ann.index.values.astype(float), ann[c].values, monthly=False)
    rows.append(dict(series=c, period="2005-2019 annual", slope_per_decade=10*o["slope"], ci_lo=10*o["ci_lo"], ci_hi=10*o["ci_hi"],
                     p=o["p"], mean=float(ann[c].mean()), pct_per_decade=1000*o["slope"]/float(ann[c].mean()),
                     r_transects_Tperm=float(np.corrcoef(mm[c], mm.Tperm)[0, 1]),
                     r_annual_Tperm=float(np.corrcoef(ann[c], annA.Tperm)[0, 1])))
for c in ["Tperm", "Tall"]:
    o = R.trend_bcz(annA.index.values.astype(float), annA[c].values, monthly=False)
    rows.append(dict(series=c, period="2005-2019 annual", slope_per_decade=10*o["slope"], ci_lo=10*o["ci_lo"], ci_hi=10*o["ci_hi"],
                     p=o["p"], mean=float(annA[c].mean()), pct_per_decade=1000*o["slope"]/float(annA[c].mean())))
pd.DataFrame(rows).to_csv(OUT/"d1_trends.csv", index=False)
pd.Series(summ).to_csv(OUT/"d1_summary.csv")
print(pd.Series(summ).to_string()); print(pd.DataFrame(rows).round(3).to_string())
