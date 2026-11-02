"""DK2: half-power front positions ALONG the LM Gould section, compared with the in situ jet positions.
Section: bin-mean track of the GV2023 ADCP transects (0-840 km, 55.1S-62.3S), 121 points (~6.9 km).
Altimetry (each product, each month): KE = 0.5|U_m|^2 interpolated to the section; points with monthly ice fraction
> 0.10 set to NaN; fronts.detect_profile with the D1 baseline (200 cm2 s-2, 100 km, 25 %) on along-section distance.
Climatological envelopes from the 2001-2025 mean KE profile; a monthly envelope is assigned to a climatological one if its
half-power point lies inside the climatological segment (strongest if several) -> one position series per front.
In situ (each transect): within each climatological segment, the latitude that splits the eastward transport of the
0-760 m total transport per 25-km bin into two equal halves (the transport analogue of the half-power point).
Comparison on matching months (2005-2019): correlation, mean offset and RMS difference of positions; trends (trend_bcz)
of the altimetric positions 2001-2025 and 2005-2019 and of the in situ annual positions 2005-2019.
Usage: python dk2_section_fronts.py <revision dir> <GV repo dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, scipy.io as sio
REV, GV, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV)); import rstats as R, fronts as F
FI = {k: i for i, k in enumerate(F.FIELDS)}
import os; THR = float(os.environ.get("DK2_THR", "0.02"))   # baseline 200 cm2 s-2; S3 sensitivity range 150-400
m = sio.loadmat(GV/"Datasets"/"TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
tr = np.array(m.transp, float); blat = np.array(m.lat, float); blon = np.array(m.lon, float)-360; dist = np.array(m.dist, float)
keep = dist <= 840; tr[~keep, :] = np.nan
mlat = np.nanmean(blat[keep], 1); mlon = np.nanmean(blon[keep], 1)
s = np.r_[0, np.cumsum(111.2*np.hypot(np.diff(mlon)*np.cos(np.radians(mlat[:-1])), np.diff(mlat)))]
ss = np.linspace(0, s[-1], 121); lo = np.interp(ss, s, mlon); la = np.interp(ss, s, mlat)
pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month-0.5)/12
time = pd.to_datetime(np.array(m.time, float)-719529, unit="D")
nval = np.isfinite(tr).sum(0)
def insitu_pos(lat_lo, lat_hi):
    out = np.full(tr.shape[1], np.nan)
    for k in range(tr.shape[1]):
        if nval[k] < 30: continue
        b = blat[:, k]; q = tr[:, k]
        sel = np.isfinite(b) & np.isfinite(q) & (b >= min(lat_lo, lat_hi)) & (b <= max(lat_lo, lat_hi))
        if sel.sum() < 3: continue
        bb, qq = b[sel], np.clip(q[sel], 0, None); o = np.argsort(-bb); bb, qq = bb[o], qq[o]   # north -> south
        c = np.cumsum(qq)
        if c[-1] <= 0: continue
        out[k] = np.interp(0.5*c[-1], c, bb)
    return out
rows, comp = [], []
for p in ("allsat", "twosat"):
    KE, ICE = [], []
    for y in range(2001, 2026):
        for mo in range(1, 13):
            with xr.open_dataset(REV/"data"/f"monthly_{p}"/f"{y}"/f"{mo:02d}.nc") as d:
                dd = d.sel(longitude=slice(-68, -58), latitude=slice(-64, -54)).load()
            u = dd.u_mean.interp(**pts).values; v = dd.v_mean.interp(**pts).values
            KE.append(0.5*(u**2+v**2)); ICE.append(np.nan_to_num(dd.ice_frac.interp(**pts).values))
    KE = np.array(KE); ICE = np.array(ICE); KEm = np.where(ICE > 0.10, np.nan, KE)
    LATMIN = float(os.environ.get("DK2_LATMIN", "-90"))   # SI S3 sensitivity: truncate the profile at the ice-free limit (-58.55)
    KEm[:, la < LATMIN] = np.nan
    clim = np.nanmean(KEm, 0)
    cc, ct = F.detect_field(clim[:, None], la, ss, THR, 100.0, 0.25, 8)
    env = [(la[int(ct[0, c, FI["i_start"]])], la[int(ct[0, c, FI["i_end"]])], ct[0, c, FI["lat_hp"]]) for c in range(cc[0])]
    print(p, "climatological envelopes (lat start, lat end, half-power lat):", [tuple(round(v, 2) for v in e) for e in env], flush=True)
    pos = np.full((len(KE), len(env)), np.nan)
    for k in range(len(KE)):
        c1, t1 = F.detect_field(KEm[k][:, None], la, ss, THR, 100.0, 0.25, 8)
        for e, (a, b, _) in enumerate(env):
            lo_, hi_ = min(a, b), max(a, b)
            cand = [(t1[0, j, FI["peak_val"]], t1[0, j, FI["lat_hp"]]) for j in range(c1[0]) if lo_ <= t1[0, j, FI["lat_hp"]] <= hi_]
            if cand: pos[k, e] = max(cand)[1]
    for e, (a, b, lhp) in enumerate(env):
        name = f"front{e+1}"
        T[f"{p}_{name}_lat"] = pos[:, e]
        cov = np.isfinite(pos[:, e]).mean()
        for per, (y0, y1) in {"2001-2025": (2001, 2025), "2005-2019": (2005, 2019)}.items():
            sel = (T.year >= y0) & (T.year <= y1) & np.isfinite(pos[:, e])
            o = R.trend_bcz(T.x[sel].values, pos[sel.values, e]*111.2, month=T.month[sel].values) if sel.sum() > 120 else dict(slope=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan)
            rows.append(dict(product=p, front=name, clim_lat_hp=lhp, seg=(round(a, 2), round(b, 2)), coverage=cov, period=per,
                             trend_km_dec=10*o["slope"], ci_lo=10*o["ci_lo"], ci_hi=10*o["ci_hi"], p=o["p"]))
        ins = insitu_pos(a, b)
        I = pd.DataFrame(dict(year=time.year, month=time.month, ins=ins)).dropna()
        M = I.merge(T[["year", "month", f"{p}_{name}_lat"]], on=["year", "month"]).dropna()
        d_ = M[f"{p}_{name}_lat"] - M.ins
        ann = M.groupby("year")[["ins", f"{p}_{name}_lat"]].mean()
        oi = R.trend_bcz(ann.index.values.astype(float), ann.ins.values*111.2, monthly=False) if len(ann) >= 10 else dict(slope=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan)
        comp.append(dict(product=p, front=name, clim_lat_hp=lhp, n_matched=len(M), r_transects=np.corrcoef(M.ins, M[f"{p}_{name}_lat"])[0, 1] if len(M) > 5 else np.nan,
                         r_annual=np.corrcoef(ann.ins, ann[f"{p}_{name}_lat"])[0, 1] if len(ann) > 5 else np.nan,
                         mean_offset_km=111.2*d_.mean(), rms_diff_km=111.2*np.sqrt((d_**2).mean()), sd_alt_km=111.2*M[f"{p}_{name}_lat"].std(), sd_ins_km=111.2*M.ins.std(),
                         insitu_trend_km_dec_2005_2019=10*oi["slope"], insitu_ci=(10*oi["ci_lo"], 10*oi["ci_hi"]), insitu_p=oi["p"]))
    np.savez_compressed(OUT/f"dk2_profiles_{p}.npz", lat=la, lon=lo, dist=ss, ke_clim=clim, ice_months=(ICE > 0.10).sum(0), pos=pos)
T.to_csv(OUT/"dk2_positions.csv", index=False)
pd.DataFrame(rows).to_csv(OUT/"dk2_trends.csv", index=False); pd.DataFrame(comp).to_csv(OUT/"dk2_vs_insitu.csv", index=False)
pd.set_option("display.width", 250)
print(pd.DataFrame(rows).round(3).to_string()); print(pd.DataFrame(comp).round(3).to_string())
