"""DK4: monthly altimetric series along the mean LM Gould track for 2001-2025, both products, for the 'validate on the
overlap, then extend' figure: section-mean EKE north of 58.5S (ice-free points) and the surface transport proxy T_perm
(points ice-free in every month of that product). Usage: python dk4_section_series.py <revision dir> <GV repo dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr, scipy.io as sio
REV, GV, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
m = sio.loadmat(GV/"Datasets"/"TotalTransport_os38nbmaxz_760_2.mat", squeeze_me=True, struct_as_record=False)["totaltransp"]
blat = np.array(m.lat, float); blon = np.array(m.lon, float)-360; dist = np.array(m.dist, float); keep = dist <= 840
mlat = np.nanmean(blat[keep], 1); mlon = np.nanmean(blon[keep], 1)
s = np.r_[0, np.cumsum(np.hypot(np.diff(mlon)*np.cos(np.radians(mlat[:-1])), np.diff(mlat)))]
ss = np.linspace(0, s[-1], 121); lo = np.interp(ss, s, mlon); la = np.interp(ss, s, mlat)
R_E = 6.371e6
seg = np.hypot(R_E*np.cos(np.radians(0.5*(la[:-1]+la[1:])))*np.radians(np.diff(lo)), R_E*np.radians(np.diff(la)))
dsx = np.r_[seg[0]/2, 0.5*(seg[:-1]+seg[1:]), seg[-1]/2]
tx = np.gradient(lo)*np.cos(np.radians(la)); ty = np.gradient(la); n = np.hypot(tx, ty); nx, ny = -ty/n, tx/n
f = nx < 0; nx[f] *= -1; ny[f] *= -1
pts = dict(longitude=xr.DataArray(lo, dims="p"), latitude=xr.DataArray(la, dims="p"))
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month-0.5)/12
for p in ("allsat", "twosat"):
    U, V, S2, I = [], [], [], []
    for y in range(2001, 2026):
        for mo in range(1, 13):
            with xr.open_dataset(REV/"data"/f"monthly_{p}"/f"{y}"/f"{mo:02d}.nc") as d:
                dd = d.sel(longitude=slice(-68, -58), latitude=slice(-64, -54)).load()
            U.append(dd.u_mean.interp(**pts).values); V.append(dd.v_mean.interp(**pts).values)
            S2.append(dd.var_intra.interp(**pts).values); I.append(np.nan_to_num(dd.ice_frac.interp(**pts).values))
    U, V, S2, I = map(np.array, (U, V, S2, I))
    eke = (0.5*((U-np.nanmean(U, 0))**2 + (V-np.nanmean(V, 0))**2) + 0.5*S2)*1e4
    north = la > -58.5
    T[f"{p}_EKE_north"] = np.nanmean(np.where((I <= 0.10) & north[None], eke, np.nan), 1)
    perm = (I.max(0) <= 0.10) & np.all(np.isfinite(U), 0)
    T[f"{p}_Tperm"] = np.nansum(np.where(perm[None], U*nx+V*ny, 0)*dsx, 1)
    T[f"{p}_n_perm"] = perm.sum()
T.to_csv(OUT/"dk4_section_series.csv", index=False); print(T.describe().T.round(2).to_string())
