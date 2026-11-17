"""Step 1: ACC-band and sector KE/EKE monthly series of Product 4, recomputed with the method of 02a_energy.py, validated against
outputs/a5/a5_series_twosat.csv; period means 2001-2017 vs 2018-2025 for both products.

Purpose : reproduce the Product 4 band series that the paper uses (a5_series_twosat.csv) from the raw monthly fields with the exact
          02a_energy.py code path, so that every other Product 4 diagnostic here is built on validated inputs; then the late-versus-
          early increase of the band means (SI S2.3: +13 % Product 1, +6 % Product 4).
Inputs  : Product 4 monthly fields (p4_common.FILES); outputs/a5/a5_series_allsat.csv, a5_series_twosat.csv
Outputs : outputs/p4/p4_band_series.csv     (year, month, x, ACC/Atlantic/Indian/Pacific _KE and _EKE, cm2 s-2; Product 4, recomputed)
          outputs/p4/p4_band_validation.csv (column, max |diff|, max relative diff vs a5_series_twosat.csv)
          outputs/p4/p4_band_periods.csv    (product, series, mean 2001-2017, mean 2018-2025, % increase, mean 2001-2025)
          outputs/p4/p4_band_series_both.csv (year, month, x, P1_ACC_KE, P4_ACC_KE, P1_ACC_EKE, P4_ACC_EKE; P1 copied from a5_series_allsat.csv)
          outputs/p4/p4_band_mask.npz       (acc mask, sector labels, lat, lon, P4 2001-2025 mean ADT on its own grid)
Method  : as 02a_energy.py (lines 'pass 1' and the monthly loop): ubar = unweighted mean of the 300 monthly-mean velocities;
          valid = finite velocity in all 300 months and maximum monthly ice fraction <= 0.10; ACC band = valid pixels whose 2001-2025
          mean ADT OF THE SAME PRODUCT (Product 4, on its own 1/4 deg grid; unweighted mean of adt_mean) lies in [-0.6, +0.2] m;
          KE = 0.5|U_m|^2, EKE = 0.5|U_m - ubar|^2 + 0.5 var_intra; cos(lat)-weighted area means; x 1e4 -> cm2 s-2.
          Sectors: Indian 20-147 E, Pacific 147 E-70 W, Atlantic 70 W-20 E (lon mod 360). Period increase = 100 (mean of the monthly
          values 2018-2025 / mean 2001-2017 - 1).
Run     : python scripts/06_product4/p4_01_band_series.py
"""
import p4_common as C
import numpy as np, pandas as pd

lat, lon = C.grid()
U, V, S2, A, I = (C.stack(v) for v in ("u_mean", "v_mean", "var_intra", "adt_mean", "ice_frac"))
nt = U.shape[0]
su = np.zeros(U.shape[1:]); sv = np.zeros_like(su); sa = np.zeros_like(su); imax = np.zeros_like(su); nvalid = np.zeros_like(su)
for k in range(nt):                                   # pass 1, as 02a
    u, v, a, ice = (np.asarray(Z[k], np.float64) for Z in (U, V, A, I))
    ok = np.isfinite(u) & np.isfinite(v)
    su += np.where(ok, u, 0); sv += np.where(ok, v, 0); nvalid += ok
    sa += np.nan_to_num(a); imax = np.fmax(imax, np.nan_to_num(ice))
valid = (nvalid == nt) & ~(imax > 0.10)
ub = su / np.maximum(nvalid, 1); vb = sv / np.maximum(nvalid, 1); adt_mean = sa / nt
acc = valid & (adt_mean >= -0.6) & (adt_mean <= 0.2)
lon360 = np.mod(lon, 360)
sec = np.where((lon360 >= 20) & (lon360 < 147), "Indian", np.where((lon360 >= 147) & (lon360 < 290), "Pacific", "Atlantic"))
Sx = np.broadcast_to(sec[None, :], acc.shape)
masks = {"ACC": acc, "Atlantic": acc & (Sx == "Atlantic"), "Indian": acc & (Sx == "Indian"), "Pacific": acc & (Sx == "Pacific")}
w = np.broadcast_to(np.cos(np.radians(lat))[:, None], acc.shape)
rows = []
for k in range(nt):
    Uk, Vk, S2k = (np.asarray(Z[k], np.float64) for Z in (U, V, S2))
    ke = 0.5 * (Uk ** 2 + Vk ** 2); eke = 0.5 * ((Uk - ub) ** 2 + (Vk - vb) ** 2) + 0.5 * S2k
    rec = {}
    for name, m in masks.items():
        ww = w[m]
        for vn, F in (("KE", ke), ("EKE", eke)):
            aa = F[m]; ok = np.isfinite(aa)
            rec[f"{name}_{vn}"] = 1e4 * float(np.sum(aa[ok] * ww[ok]) / np.sum(ww[ok]))
    rows.append(rec)
ts = pd.concat([C.T, pd.DataFrame(rows)], axis=1)
ts.to_csv(C.HERE / "p4_band_series.csv", index=False)
np.savez_compressed(C.HERE / "p4_band_mask.npz", acc=acc, sector=sec, lat=lat, lon=lon, adt_mean=adt_mean.astype(np.float32))
print("P4 ACC band pixels", int(acc.sum()))

ref = pd.read_csv(C.OUT / "a5" / "a5_series_twosat.csv")
val = []
for c in [c for c in ts.columns if c.endswith(("_KE", "_EKE"))]:
    d = np.abs(ts[c].values - ref[c].values)
    val.append(dict(column=c, max_abs_diff=float(d.max()), max_rel_diff=float(np.max(d / np.abs(ref[c].values)))))
val = pd.DataFrame(val); val.to_csv(C.HERE / "p4_band_validation.csv", index=False); print(val.to_string(index=False))

per = []
P1 = pd.read_csv(C.OUT / "a5" / "a5_series_allsat.csv")
for prod, D in (("P1 (a5_series_allsat)", P1), ("P4 (a5_series_twosat)", ref), ("P4 (recomputed here)", ts)):
    e, l = D.year.values <= 2017, D.year.values >= 2018
    for c in [f"{r}_{v}" for r in ("ACC", "Pacific", "Atlantic", "Indian") for v in ("KE", "EKE")]:
        y = D[c].values
        per.append(dict(product=prod, series=c, mean_2001_2017=y[e].mean(), mean_2018_2025=y[l].mean(),
                        pct_increase=100 * (y[l].mean() / y[e].mean() - 1), mean_2001_2025=y.mean()))
per = pd.DataFrame(per); per.to_csv(C.HERE / "p4_band_periods.csv", index=False)
both = C.T.copy()
for c in ("ACC_KE", "ACC_EKE"):
    both[f"P1_{c}"] = P1[c].values; both[f"P4_{c}"] = ts[c].values
both.to_csv(C.HERE / "p4_band_series_both.csv", index=False)
print(per[per.series.str.startswith("ACC")].round(2).to_string(index=False))
