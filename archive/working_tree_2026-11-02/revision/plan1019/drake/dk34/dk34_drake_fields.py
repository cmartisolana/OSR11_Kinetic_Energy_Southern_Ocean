"""dk34: Drake Passage fields for Figure 2 (rebuilt): EKE trend map 2001-2025 (Product 1, annual means, least squares,
per decade; 2x2 pixel average; pixels with sea-ice fraction > 0.10 in any month masked) and annual EKE between the SAF and
PF per 2 deg longitude (band as dk8: mean ADT between the SAF and PF streamlines, pixels ice-free in every month), % of the
2001-2025 mean of each bin. Usage: python dk34_drake_fields.py <revision dir> <dk5 stream edges csv> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
REV, EDG, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
mf = sorted((REV / "data" / "monthly_allsat").glob("*/*.nc")); mf = [f for f in mf if 2001 <= int(f.parent.name) <= 2025]
box = dict(longitude=slice(-87, -33), latitude=slice(-66, -46))
U, V, VI, IC, Y = [], [], [], [], []
for f in mf:
    with xr.open_dataset(f) as d:
        s = d.sel(**box); U.append(s.u_mean.values.astype(np.float32)); V.append(s.v_mean.values.astype(np.float32))
        VI.append(s.var_intra.values.astype(np.float32)); IC.append(np.nan_to_num(s.ice_frac.values).astype(np.float32)); Y.append(int(f.parent.name))
        lat, lon = s.latitude.values.astype(float), s.longitude.values.astype(float)
U, V, VI, IC, Y = np.array(U), np.array(V), np.array(VI), np.array(IC), np.array(Y)
E = (0.5 * ((U - np.nanmean(U, 0)) ** 2 + (V - np.nanmean(V, 0)) ** 2) + 0.5 * VI) * 1e4; del U, V, VI
icemax = IC.max(0)
YY = np.arange(2001, 2026); EA = np.array([np.nanmean(E[Y == y], 0) for y in YY])
xa = YY - YY.mean(); sl = np.tensordot(xa, EA - EA.mean(0), axes=(0, 0)) / (xa ** 2).sum() * 10; sl[icemax > 0.10] = np.nan
cs = 2; ny, nx = (lat.size // cs) * cs, (lon.size // cs) * cs
sm = np.nanmean(sl[:ny, :nx].reshape(ny // cs, cs, nx // cs, cs), axis=(1, 3))
_e = pd.read_csv(EDG).set_index("psi").adt_edge
LEV_SAF, LEV_PF = (float(np.interp(p, _e.index[::-1], _e.values[::-1])) for p in (-0.5, -1.125))
a = xr.open_dataset(REV / "outputs" / "clim_allsat.nc").adt.interp(latitude=lat, longitude=lon).values
band = np.isfinite(a) & (a >= LEV_PF) & (a <= LEV_SAF) & (icemax <= 0.10); w = np.cos(np.radians(lat))[:, None] * band
lb = np.arange(-86, -34, 2.0); H = np.full((YY.size, lb.size), np.nan)
for i, l0 in enumerate(lb):
    m = (lon >= l0) & (lon < l0 + 2); ww = w[:, m]
    if ww.sum() == 0: continue
    H[:, i] = np.nansum(EA[:, :, m] * ww[None], axis=(1, 2)) / ww.sum()
Hp = 100 * (H / np.nanmean(H, 0) - 1)
np.savez(OUT / "dk34_drake_fields.npz", slope_map=sm, slat=lat[:ny].reshape(-1, cs).mean(1), slon=lon[:nx].reshape(-1, cs).mean(1),
         hov=Hp, hov_years=YY, hov_lon=lb + 1, lev_saf=LEV_SAF, lev_pf=LEV_PF)
print("ok", sm.shape, np.nanpercentile(sm, [5, 50, 95]).round(1), np.nanmin(Hp).round(), np.nanmax(Hp).round())
