"""DK5c: zone-level trends along the LM Gould transect (passage zones of GV23 Figs. 3-4; streamwise frontal zones of GV23
Figs. 5, 7), altimetry (both products) vs in situ; uses dk5_section_fields.npz and the GV23 files.
Usage: python dk5c_zones.py <revision dir> <GV repo dir> <dk5 out dir>"""
import sys, numpy as np, pandas as pd, scipy.io as sio, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
REV, GV, OUT = (Path(a) for a in sys.argv[1:4]); sys.path.insert(0, str(REV)); import rstats as R
exec(open(Path(__file__).with_name("dk5_geometry.py")).read())        # lo, la, dkm, dsx
Z = np.load(OUT/"dk5_section_fields.npz"); E = pd.read_csv(OUT/"dk5_stream_edges.csv"); psi2adt = dict(zip(E.psi.round(3), E.adt_edge))
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - .5)/12
ov = ((T.x >= 2005 + 9/12) & (T.x <= 2019 + 4/12)).values
def tr(x, y, monthly, month=None):
    ok = np.isfinite(y)
    o = R.trend_bcz(x[ok], y[ok], month=None if month is None else month[ok], monthly=monthly); return o
def ann(years, v, minn):
    df = pd.DataFrame(dict(year=years, v=v)).dropna(); g = df.groupby("year").v; a = g.mean()[g.count() >= minn]; return a.index.values.astype(float), a.values
rows = []
PZ = {"SAF-PF zone (150-300 km)": (150, 300), "PF southern flank (350-550 km)": (350, 550), "SACCF (750-850 km)": (750, 850)}
SZ = {"SAF (-0.40 to -0.60 m)": (-0.40, -0.60), "SAF-PF zone (-0.60 to -1.00 m)": (-0.60, -1.00), "PF (-1.00 to -1.25 m)": (-1.00, -1.25)}
for p in ("allsat", "twosat"):
    ice = Z[f"{p}_ice"] > 0.10; un = np.where(ice, np.nan, Z[f"{p}_un"]) * 100; ek = np.where(ice, np.nan, Z[f"{p}_eke"]) * 1e4; adt = Z[f"{p}_adt"]
    for name, (a, b) in PZ.items():
        sel = (dkm >= a) & (dkm < b); y = np.nanmean(un[:, sel], 1)
        xa, ya = ann(T.year[ov].values, y[ov], 3); o1 = tr(xa, ya, False)
        if np.isfinite(y).all(): o2 = tr(T.x.values, y, True, T.month.values)
        else: xa2, ya2 = ann(T.year.values, y, 6); o2 = tr(xa2, ya2, False)
        for per, o in (("2005-2019", o1), ("2001-2025", o2)):
            rows.append(dict(coord="passage", zone=name, src=p, var="u (cm/s)", period=per, mean=np.nanmean(y), slope_dec=10*o["slope"], lo=10*o["ci_lo"], hi=10*o["ci_hi"], p=o["p"]))
    for name, (a, b) in SZ.items():
        lo_e, hi_e = sorted((psi2adt[round(b, 3)], psi2adt[round(a, 3)]))
        inb = (adt >= lo_e) & (adt < hi_e); icy = (inb & ice).any(1)
        q = np.where(icy, np.nan, np.nansum(np.where(inb, Z[f"{p}_un"]*dsx[None, :], 0), 1))
        e = np.array([np.nanmean(ek[k, inb[k]]) if inb[k].any() else np.nan for k in range(len(adt))])
        for var, y in (("transport (m2/s)", q), ("EKE (cm2/s2)", e)):
            xa, ya = ann(T.year[ov].values, y[ov], 3); o1 = tr(xa, ya, False)
            xa2, ya2 = ann(T.year.values, y, 6); o2 = tr(xa2, ya2, False)
            for per, o in (("2005-2019", o1), ("2001-2025", o2)):
                rows.append(dict(coord="stream", zone=name, src=p, var=var, period=per, mean=np.nanmean(y), slope_dec=10*o["slope"], lo=10*o["ci_lo"], hi=10*o["ci_hi"], p=o["p"]))
# in situ
LV = sio.loadmat(GV/"Datasets"/"lmgvelos38nb_50_300.mat", squeeze_me=True, struct_as_record=False)["lmgdata"]
TH = sio.loadmat(GV/"Datasets"/"trandsmiss_theta.mat", squeeze_me=True)
uoffs = np.mean(np.abs(np.r_[np.atleast_1d(TH["uoffs_north"]), np.atleast_1d(TH["uoffs_south"])]))
U = np.array(LV.u, float); z = np.array(LV.z, float); d = np.array(LV.dist, float); nb = np.array(LV.northbound) == 1
U[:, :, nb] -= uoffs; U[:, :, ~nb] += uoffs
for k in np.where(nb)[0]:
    dd = np.abs(np.nanmax(d[:, k]) - d[:, k]) + 12.5; o = np.argsort(dd); d[:, k] = dd[o]; U[:, :, k] = U[:, o, k]
zs = (z >= 70) & (z <= 118); us = np.nanmean(U[zs], 0) * 100          # (bin, transect)
t = pd.to_datetime(np.array(LV.time, float) - 719529, unit="D"); yr = t.year.values; inov = (yr >= 2005) & (yr <= 2019)
for name, (a, b) in PZ.items():
    y = np.array([np.nanmean(us[(d[:, k] >= a) & (d[:, k] < b), k]) for k in range(us.shape[1])])
    xa, ya = ann(yr[inov], y[inov], 3); o = tr(xa, ya, False)
    rows.append(dict(coord="passage", zone=name, src="insitu", var="u (cm/s)", period="2005-2019", mean=np.nanmean(y), slope_dec=10*o["slope"], lo=10*o["ci_lo"], hi=10*o["ci_hi"], p=o["p"]))
L = lambda f: sio.loadmat(GV/"Datasets"/f, squeeze_me=True, struct_as_record=False)["stream"]
ST, SE = L("Streamwise_binning_os38nb.mat"), L("Streamwise_binning_eke.mat"); psic = np.array(ST.psic, float)
for name, (a, b) in SZ.items():
    j = (psic < a) & (psic > b)
    for var, arr, f, tm in (("transport (Sv)", np.array(ST.transp, float), 1e-6, ST.time), ("EKE (J/m, 0-780 m)", np.array(SE.eke, float), 1.0, SE.time)):
        X = arr[:, j] * f; y = np.nansum(X, 1) if var.startswith("transport") else np.nanmean(X, 1)
        y[np.isnan(X).all(1)] = np.nan
        yy = pd.to_datetime(np.array(tm, float) - 719529, unit="D").year.values
        xa, ya = ann(yy, y, 3); sel = (xa >= 2005) & (xa <= 2019); o = tr(xa[sel], ya[sel], False)
        rows.append(dict(coord="stream", zone=name, src="insitu", var=var, period="2005-2019", mean=np.nanmean(y), slope_dec=10*o["slope"], lo=10*o["ci_lo"], hi=10*o["ci_hi"], p=o["p"]))
D = pd.DataFrame(rows); D["pct"] = 100*D.slope_dec/D["mean"]; D["pct_lo"] = 100*D.lo/D["mean"]; D["pct_hi"] = 100*D.hi/D["mean"]
D.to_csv(OUT/"dk5_zones.csv", index=False)
pd.set_option("display.width", 220); print(D.round(3).to_string())
