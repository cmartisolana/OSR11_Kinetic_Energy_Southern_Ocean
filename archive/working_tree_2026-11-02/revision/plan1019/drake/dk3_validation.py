"""DK3: does altimetry reproduce the in situ Drake Passage record? (basis for extending the record with altimetry)
In situ (LM Gould ADCP, Gutierrez-Villanueva et al. 2023; 2005-2019, 248 transects):
  u   across-transect velocity per 25-km bin at 70-118 m (mean of the 70, 94, 118 m levels; below the Ekman layer)
  eke depth-integrated EKE 0-780 m per bin (their eke_uvprime_all_os38.mat)
Altimetry (each product): monthly-mean surface geostrophic velocity normal to the local track direction (east-positive)
and monthly EKE = 0.5|U_m - ubar|^2 + 0.5 var_intra, interpolated to each bin position, same calendar month.
Points with monthly ice fraction > 0.10 are excluded (as in the proxies).
Outputs: (1) time-mean velocity profiles by 0.25 deg latitude; (2) per-bin correlation of velocity (all matched bins,
and by latitude band); (3) section-mean EKE north of 58.5S and over the ice-free part: annual-mean correlation and trends
(trend_bcz, annual) 2005-2019, in % of mean per decade.
Note: a monthly mean field is compared with a 2-day transect, so mesoscale noise lowers transect-level correlations.
Usage: python dk3_validation.py <revision dir> <GV repo dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, scipy.io as sio
REV, GV, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R
v = sio.loadmat(GV/"Datasets"/"lmgvelos38nb_50_300.mat", squeeze_me=True, struct_as_record=False)["lmgdata"]
e = sio.loadmat(GV/"Datasets"/"eke_uvprime_all_os38.mat", squeeze_me=True, struct_as_record=False)["lmgdata"]
time = pd.to_datetime(np.array(v.time, float)-719529, unit="D")
uin = np.nanmean(np.array(v.u)[2:5], 0)                       # (48 bins, 248)
vlat, vlon = np.array(v.lat, float), np.array(v.lon, float)-360
elat, elon, eke_in = np.array(e.lat, float), np.array(e.lon, float), np.array(e.eke, float)
# local transect direction from the bin-mean track (south->north order irrelevant for an east-positive normal)
def normal(lo, la):
    tx = np.gradient(lo)*np.cos(np.radians(la)); ty = np.gradient(la); n = np.hypot(tx, ty)
    nx, ny = -ty/n, tx/n; f = nx < 0; nx[f] *= -1; ny[f] *= -1; return nx, ny
LIM = -58.5
res = {}
for p in ("allsat", "twosat"):
    # 2001-2025 mean velocity of the product (for EKE)
    fs = {(y, m): REV/"data"/f"monthly_{p}"/f"{y}"/f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)}
    acc = None
    for f in fs.values():
        with xr.open_dataset(f) as d:
            dd = d.sel(longitude=slice(-68, -58), latitude=slice(-64, -54))
            if acc is None: acc = [np.zeros(dd.u_mean.shape), np.zeros(dd.u_mean.shape), 0]
            acc[0] += np.nan_to_num(dd.u_mean.values); acc[1] += np.nan_to_num(dd.v_mean.values); acc[2] += 1
            glat, glon = dd.latitude.values, dd.longitude.values
    ub = xr.DataArray(acc[0]/acc[2], coords=dict(latitude=glat, longitude=glon), dims=("latitude", "longitude"))
    vb = xr.DataArray(acc[1]/acc[2], coords=dict(latitude=glat, longitude=glon), dims=("latitude", "longitude"))
    U_alt = np.full(uin.shape, np.nan); E_alt = np.full(eke_in.shape, np.nan)
    for k, t in enumerate(time):
        if (t.year, t.month) not in fs: continue
        with xr.open_dataset(fs[(t.year, t.month)]) as d:
            dd = d.sel(longitude=slice(-68, -58), latitude=slice(-64, -54)).load()
        for (la, lo, out, kind) in ((vlat[:, k], vlon[:, k], U_alt, "u"), (elat[:, k], elon[:, k], E_alt, "e")):
            ok = np.isfinite(la) & np.isfinite(lo)
            if ok.sum() < 3: continue
            pts = dict(longitude=xr.DataArray(lo[ok], dims="p"), latitude=xr.DataArray(la[ok], dims="p"))
            uu = dd.u_mean.interp(**pts).values; vv = dd.v_mean.interp(**pts).values
            ice = np.nan_to_num(dd.ice_frac.interp(**pts).values); bad = ice > 0.10
            if kind == "u":
                nx, ny = normal(lo[ok], la[ok]); val = uu*nx + vv*ny
            else:
                val = (0.5*((uu-ub.interp(**pts).values)**2 + (vv-vb.interp(**pts).values)**2) + 0.5*dd.var_intra.interp(**pts).values)*1e4
            val[bad] = np.nan; tmp = np.full(la.shape, np.nan); tmp[ok] = val; out[:, k] = tmp
    # (1)+(2) velocity
    ok = np.isfinite(uin) & np.isfinite(U_alt)
    latb = np.round(vlat/0.25)*0.25
    prof = pd.DataFrame(dict(latb=latb[ok], uin=uin[ok]*100, ualt=U_alt[ok]*100)).groupby("latb").agg(["mean", "count"])
    prof.columns = ["_".join(c) for c in prof.columns]; prof.to_csv(OUT/f"dk3_velocity_profile_{p}.csv")
    band = {"all": np.ones_like(ok), "north": vlat > LIM, "south": vlat <= LIM}
    for b, msk in band.items():
        mm = ok & msk
        res[(p, f"u_r_bins_{b}")] = np.corrcoef(uin[mm], U_alt[mm])[0, 1]
        res[(p, f"u_n_bins_{b}")] = int(mm.sum())
        # anomalies relative to each latitude bin's mean (removes the shared mean profile)
        dfa = pd.DataFrame(dict(l=latb[mm], a=uin[mm], b=U_alt[mm]))
        dfa["a"] -= dfa.groupby("l").a.transform("mean"); dfa["b"] -= dfa.groupby("l").b.transform("mean")
        res[(p, f"u_r_anom_{b}")] = np.corrcoef(dfa.a, dfa.b)[0, 1]
    res[(p, "u_rms_diff_cm_s")] = 100*np.sqrt(np.nanmean((uin[ok]-U_alt[ok])**2))
    res[(p, "u_mean_profile_r")] = np.corrcoef(prof.uin_mean, prof.ualt_mean)[0, 1]
    # (3) EKE section means, annual
    rows = []
    for k, t in enumerate(time):
        for b, msk in (("north", elat[:, k] > LIM), ("icefree", np.isfinite(E_alt[:, k]))):
            sel = msk & np.isfinite(eke_in[:, k]) & np.isfinite(E_alt[:, k])
            if sel.sum() < 10: continue
            rows.append(dict(year=t.year, month=t.month, band=b, eke_in=np.mean(eke_in[sel, k]), eke_alt=np.mean(E_alt[sel, k])))
    D = pd.DataFrame(rows); D.to_csv(OUT/f"dk3_eke_transects_{p}.csv", index=False)
    for b in ("north", "icefree"):
        A = D[D.band == b].groupby("year")[["eke_in", "eke_alt"]].mean().loc[2005:2019]
        res[(p, f"eke_r_transects_{b}")] = np.corrcoef(D[D.band == b].eke_in, D[D.band == b].eke_alt)[0, 1]
        res[(p, f"eke_r_annual_{b}")] = np.corrcoef(A.eke_in, A.eke_alt)[0, 1]
        for c in ("eke_in", "eke_alt"):
            o = R.trend_bcz(A.index.values.astype(float), A[c].values, monthly=False)
            mu = A[c].mean()
            res[(p, f"{c}_trend_pct_dec_{b}")] = 1000*o["slope"]/mu
            res[(p, f"{c}_trend_ci_{b}")] = (round(1000*o["ci_lo"]/mu, 1), round(1000*o["ci_hi"]/mu, 1))
            res[(p, f"{c}_trend_p_{b}")] = o["p"]
        A.to_csv(OUT/f"dk3_eke_annual_{b}_{p}.csv")
S = pd.Series(res); S.index = [f"{a}|{b}" for a, b in S.index]; S.to_csv(OUT/"dk3_summary.csv")
print(S.to_string())
