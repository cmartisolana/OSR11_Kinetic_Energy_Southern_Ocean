"""D1b: surface proxy along the ACTUAL mean LM Gould ADCP track (bin-mean positions of GV2023, 0-840 km),
fixed ice-free point set; coverage-matched in situ total transport; trends and agreement.
Usage: python d1b_gvtrack.py <revision dir> <GV repo dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, scipy.io as sio
REV, GV, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
sys.path.insert(0, str(REV)); import rstats as R
R_E = 6.371e6
m = sio.loadmat(GV/"Datasets"/"TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
tr = np.array(m.transp, float); blat = np.array(m.lat, float); blon = np.array(m.lon, float)-360; dist = np.array(m.dist, float)
keep = dist <= 840
mlat = np.nanmean(blat[keep], 1); mlon = np.nanmean(blon[keep], 1)
# 121 points along the bin-mean track (piecewise linear in distance)
s = np.r_[0, np.cumsum(np.hypot(np.diff(mlon)*np.cos(np.radians(mlat[:-1])), np.diff(mlat)))]
ss = np.linspace(0, s[-1], 121); lo = np.interp(ss, s, mlon); la = np.interp(ss, s, mlat)
latm = np.radians(0.5*(la[:-1]+la[1:]))
dxs = R_E*np.cos(latm)*np.radians(np.diff(lo)); dys = R_E*np.radians(np.diff(la)); seg = np.hypot(dxs, dys)
dsx = np.r_[seg[0]/2, 0.5*(seg[:-1]+seg[1:]), seg[-1]/2]
# local normal (eastward positive) from local track direction
tx = np.gradient(lo)*np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty/nn, tx/nn
flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
files = [REV/"data"/"monthly_allsat"/f"{y}"/f"{mo:02d}.nc" for y in range(2001, 2026) for mo in range(1, 13)]
pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
VN, ICE = [], []
for f in files:
    with xr.open_dataset(f) as ds:
        d = ds.sel(longitude=slice(-70, -55), latitude=slice(-65, -53)).load()
    VN.append(d.u_mean.interp(**pts).values*nx + d.v_mean.interp(**pts).values*ny); ICE.append(d.ice_frac.interp(**pts).values)
VN = np.array(VN); ICE = np.nan_to_num(np.array(ICE))
perm = (ICE.max(0) <= 0.10) & np.all(np.isfinite(VN), 0)
lim = la[perm].min()
t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
t["x"] = t.year + (t.month-0.5)/12
t["Tperm_gv"] = np.nansum(np.where(perm[None, :], VN, 0)*dsx, 1)
t["Tall_gv"] = np.nansum(np.where(np.all(np.isfinite(VN), 0)[None, :], VN, 0)*dsx, 1)
trk = tr.copy(); trk[~keep, :] = np.nan
time = pd.to_datetime(np.array(m.time, float)-719529, unit="D")
north = blat > lim
obs = pd.DataFrame(dict(year=time.year, month=time.month, U_tot=np.nansum(trk, 0)/1e6,
                        U_north=np.nansum(np.where(north, trk, 0), 0)/1e6, nvalid=np.isfinite(trk).sum(0)))
obs = obs[obs.nvalid >= 30]; obs["U_south"] = obs.U_tot-obs.U_north
mm = obs.merge(t, on=["year", "month"])
ann = obs.groupby("year")[["U_tot", "U_north", "U_south"]].mean().loc[2005:2019]
annA = t.groupby("year")[["Tperm_gv", "Tall_gv"]].mean().loc[2005:2019]
res = dict(n_perm=int(perm.sum()), south_limit=float(lim), n_transects=len(obs),
           share_north_insitu=float(obs.U_north.mean()/obs.U_tot.mean()),
           share_perm_surface=float(np.nanmean(t.Tperm_gv)/np.nanmean(t.Tall_gv)))
for a in ["U_tot", "U_north", "U_south"]:
    for b in ["Tperm_gv", "Tall_gv"]:
        res[f"r_tr_{a}_{b}"] = float(np.corrcoef(mm[a], mm[b])[0, 1]); res[f"r_ann_{a}_{b}"] = float(np.corrcoef(ann[a], annA[b])[0, 1])
res["r_tr_Unorth_Usouth"] = float(np.corrcoef(obs.U_north, obs.U_south)[0, 1])
# significance of the annual U_north-Tperm correlation from phase-randomised surrogates of Tperm annual series
rng = np.random.default_rng(1); y = annA.Tperm_gv.values - annA.Tperm_gv.mean(); x = ann.U_north.values
F = np.fft.rfft(y); rs = []
for _ in range(5000):
    ph = np.exp(1j*rng.uniform(0, 2*np.pi, F.size)); ph[0] = 1
    if y.size % 2 == 0: ph[-1] = 1
    rs.append(np.corrcoef(x, np.fft.irfft(F*ph, n=y.size))[0, 1])
res["p_ann_Unorth_Tperm_surrogate"] = float(np.mean(np.abs(rs) >= abs(res["r_ann_U_north_Tperm_gv"])))
rows = []
for c, df in [("U_tot", ann), ("U_north", ann), ("U_south", ann), ("Tperm_gv", annA), ("Tall_gv", annA)]:
    o = R.trend_bcz(df.index.values.astype(float), df[c].values, monthly=False)
    rows.append(dict(series=c, period="2005-2019 annual", pct_per_decade=1000*o["slope"]/df[c].mean(),
                     ci_lo_pct=1000*o["ci_lo"]/df[c].mean(), ci_hi_pct=1000*o["ci_hi"]/df[c].mean(), p=o["p"]))
for c in ["Tperm_gv", "Tall_gv"]:
    o = R.trend_bcz(t.x.values, t[c].values, month=t.month.values)
    rows.append(dict(series=c, period="2001-2025 monthly", pct_per_decade=1000*o["slope"]/t[c].mean(),
                     ci_lo_pct=1000*o["ci_lo"]/t[c].mean(), ci_hi_pct=1000*o["ci_hi"]/t[c].mean(), p=o["p"]))
pd.Series(res).to_csv(OUT/"d1b_summary.csv"); pd.DataFrame(rows).to_csv(OUT/"d1b_trends.csv", index=False)
t.to_csv(OUT/"d1b_series.csv", index=False); obs.to_csv(OUT/"d1b_insitu_split.csv", index=False)
pd.DataFrame(dict(lon=lo, lat=la, perm=perm)).to_csv(OUT/"d1b_track.csv", index=False)
print(pd.Series(res).round(3).to_string()); print(pd.DataFrame(rows).round(3).to_string())
