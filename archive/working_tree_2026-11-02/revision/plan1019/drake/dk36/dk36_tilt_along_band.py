"""dk36: trend of the isopycnal tilt between the SAF and PF per longitude bin in the Drake Passage domain (RG Argo 1 deg,
annual 2004-2025, dk32). Tilt = |grad z| of the sigma0 surface (m per 100 km); band cells: mean ADT (altimetry, on the RG
grid) between the PF and SAF streamline levels (dk34/dk8 definition); bin mean of the tilt, cos-lat weighted; trend: least
squares on annual values, % of the bin mean per decade, 95 % with n_eff from lag-1 autocorrelation.
Usage: python dk36_tilt_along_band.py <work out dir> <out dir> [bin width deg=4]"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
WO, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True); BW = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
I = np.load(WO / "dk32" / "dk32_iso.npz"); F34 = np.load(WO / "dk34" / "dk34_drake_fields.npz")
yrs, lat, lon, A = I["years"], I["lat"], I["lon"], I["adt"]; t = yrs - yrs.mean()
L = np.where(lon > 180, lon - 360, lon)
band = np.isfinite(A) & (A >= float(F34["lev_pf"])) & (A <= float(F34["lev_saf"]))
rows = []
for s in (27.3, 27.5, 27.7):
    Z = I[f"p{s}"]
    o = np.argsort(L); LO = L[o]; Zs = Z[:, :, o]; Bs = band[:, o]
    dy = np.gradient(Zs, lat, axis=1) / 111.195 * 100; dx = np.gradient(Zs, LO, axis=2) / (111.195 * np.cos(np.radians(lat))[None, :, None]) * 100
    G = np.hypot(dx, dy); w = np.cos(np.radians(lat))[:, None] * np.ones(LO.size)[None, :]
    for l0 in np.arange(-88, -36, BW):
        m = Bs & (LO[None, :] >= l0) & (LO[None, :] < l0 + BW) & np.isfinite(G).all(0)
        if m.sum() < 2: continue
        y = (G[:, m] * w[m]).sum(1) / w[m].sum(); yp = 100 * (y / y.mean() - 1)
        b, a = np.polyfit(t, yp, 1); r = yp - a - b * t; r1 = np.corrcoef(r[1:], r[:-1])[0, 1]; ne = max(3, len(t) * (1 - r1) / (1 + r1))
        se = np.sqrt((r ** 2).sum() / (ne - 2) / (t ** 2).sum()); q = stats.t.ppf(0.975, ne - 2)
        rows.append(dict(sigma0=s, lon=l0 + BW / 2, ncells=int(m.sum()), tilt_mean=y.mean(), pct=10 * b, lo=10 * (b - q * se), hi=10 * (b + q * se),
                         p=2 * stats.t.sf(abs(b / se), ne - 2)))
D = pd.DataFrame(rows); D.to_csv(OUT / "dk36_tilt_along_band.csv", index=False); print(D[D.sigma0 == 27.5].round(2).to_string(index=False))
