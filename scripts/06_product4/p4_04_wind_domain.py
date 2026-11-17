"""Step 4: Product 4 KE and EKE averaged over the domain of the wind relations (ACC band + 100 km, main Fig. 4c-e), as
outputs/f4e/f4e_series.csv is for Product 1.

Purpose : monthly domain-mean wind stress, KE, EKE and isopycnal slope with Product 4 altimetry, input of p4_05/p4_06.
Inputs  : Product 4 monthly fields; outputs/x13/x13_cells.npz (Product 1 1-deg cell fields WC, KC, EC,
          SL, distance dmin of each cell centre from the 1/8-deg P1 band, band, ice, slope_ok; validated by x14 against
          outputs/f4e/f4e_series.csv to < 1e-7); outputs/f4e/f4e_series.csv (check)
Outputs : outputs/p4/p4_wind_series.csv          year, month, x, wind, KE, EKE, slope  (P4 KE/EKE; domain 'common', the primary)
          outputs/p4/p4_wind_series_p4band.csv   same on the sensitivity domain built from Product 4's own band
          outputs/p4/p4_wind_series_p1check.csv  Product 1 series rebuilt from x13 on the P1 domain (must equal f4e_series.csv)
          outputs/p4/p4_wind_domains.npz         cell masks (dom_p1, dom_common, dom_p4band), dmin_p4, 1-deg P4 cell fields KC4, EC4
          outputs/p4/p4_wind_domain_log.txt      cell counts and checks
Method  : How the P1 domain was built (04c_slope_fields.py -> 04g_common_cells.py; reproduced by x13/x14, and prep_band.py takes its
          band-only variant from the same x13 cells): the domain is a set of 1-DEGREE CELLS, not native pixels: cells whose centre is
          within 100 km (great-circle distance, cKDTree on unit vectors) of a 1/8-deg pixel of the ACC band (P1 mean ADT -0.6..+0.2 m,
          outputs/base/mean_adt.nc), not ice-covered (f4_tilt_maps 'ice'), with the sigma0 = 27.3 slope defined in all 264 months and
          wind, KE and EKE finite in every month: 2830 cells. KE and EKE of a cell are the plain (unweighted) mean of the native-grid
          pixels whose centres fall in the cell (np.digitize on the cell edges), with the permanent ice mask applied to the pixels;
          KE = 0.5|U_m|^2, EKE = 0.5|U_m - ubar|^2 + 0.5 var_intra, ubar = unweighted mean of the 300 monthly velocities; domain means
          are cos(lat)-weighted over the cells. For Product 4 the native grid is 1/4 deg (16 pixels per cell); the permanent ice mask
          is the Product 4 one (any JJA month with ice_frac > 0, as p4_03).
          Primary domain 'common' = the P1 domain cells in which the P4 KE and EKE are also finite in every month: the wind and the
          slope then come from the same cells as for P1 (identical if no cell is lost).
          Sensitivity domain 'p4band' = the same construction with the 100 km distance measured from Product 4's own 1/4-deg ACC band
          (its 2001-2025 mean ADT, p4_band_mask.npz), all other conditions unchanged (slope defined, f4 ice = 0, all finite).
Run     : python scripts/06_product4/p4_04_wind_domain.py
"""
import p4_common as C
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

Z = np.load(C.OUT / "x13" / "x13_cells.npz")
lat, lon = Z["lat"], Z["lon"]; YY, MO = Z["year"], Z["month"]
LON, LAT = np.meshgrid(lon, lat); w = np.cos(np.radians(LAT))
la_e = np.r_[lat - 0.5, lat[-1] + 0.5]; lo_e = np.r_[lon - 0.5, lon[-1] + 0.5]; NC = lat.size * lon.size
log = []; P = lambda s: (print(s, flush=True), log.append(s))

# ---- Product 4 fields binned to the 1-deg cells (as x13 / 04g)
la4, lo4 = C.grid()
U, V, S2, I = (C.stack(v) for v in ("u_mean", "v_mean", "var_intra", "ice_frac"))
ice4 = np.zeros(U.shape[1:], bool)
for k in np.nonzero(np.isin(C.MON, (6, 7, 8)))[0]:
    ice4 |= np.nan_to_num(np.asarray(I[k], float)) > 0
ub = np.asarray(U, np.float64).sum(0) / 300; vb = np.asarray(V, np.float64).sum(0) / 300   # as 04g/x13: su / 300


def binner(la, lo):
    i = np.digitize(la, la_e) - 1; j = np.digitize(lo, lo_e) - 1
    II, J = np.meshgrid(i, j, indexing="ij"); ok = (II >= 0) & (II < lat.size) & (J >= 0) & (J < lon.size)
    return np.where(ok, II * lon.size + J, NC).ravel()


def cellmean(F, idx):
    f = F.ravel(); g = np.isfinite(f); s = np.bincount(idx[g], f[g], NC + 1)[:NC]; n = np.bincount(idx[g], None, NC + 1)[:NC]
    with np.errstate(invalid="ignore"):
        return (s / np.where(n > 0, n, np.nan)).reshape(lat.size, lon.size)


i4 = binner(la4, lo4)
EC4 = np.empty((300, lat.size, lon.size), np.float32); KC4 = np.empty_like(EC4)
for k in range(300):
    u, v, vi = (np.asarray(F[k], np.float64) for F in (U, V, S2))
    e = (0.5 * ((u - ub) ** 2 + (v - vb) ** 2) + 0.5 * vi) * 1e4; kk = 0.5 * (u ** 2 + v ** 2) * 1e4
    e[ice4] = np.nan; kk[ice4] = np.nan; EC4[k] = cellmean(e, i4); KC4[k] = cellmean(kk, i4)

# ---- domains
fin1 = np.isfinite(Z["EC"]).all(0) & np.isfinite(Z["KC"]).all(0) & np.isfinite(Z["WC"]).all(0) & np.isfinite(Z["SL"]).all(0) & (Z["ice"] == 0)
dom_p1 = (Z["dmin"] <= 100 + 1e-6) & fin1 & Z["slope_ok"]                       # x14 'buf100' = 04g domain
fin4 = np.isfinite(EC4).all(0) & np.isfinite(KC4).all(0)
dom_common = dom_p1 & fin4


def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo); return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


B = np.load(C.HERE / "p4_band_mask.npz"); A4 = B["adt_mean"].astype(float)
b4 = np.isfinite(A4) & (A4 >= -0.6) & (A4 <= 0.2)          # as 04c on the P1 mean ADT (no ice condition in the band itself)
L4, B4 = np.meshgrid(lo4, la4)
dch, _ = cKDTree(xyz(B4[b4], L4[b4])).query(xyz(LAT.ravel(), LON.ravel()))
dmin4 = (2 * 6371.0 * np.arcsin(np.clip(dch / 2, 0, 1))).reshape(LAT.shape)
fin_other = np.isfinite(Z["WC"]).all(0) & np.isfinite(Z["SL"]).all(0) & (Z["ice"] == 0) & Z["slope_ok"]
dom_p4band = (dmin4 <= 100 + 1e-6) & fin_other & fin4
P(f"P1 domain (x14 buf100 = 04g): {int(dom_p1.sum())} cells; P4 KE/EKE finite in all months in {int(dom_common.sum())} of them "
  f"(lost: {int((dom_p1 & ~fin4).sum())}); P4-own-band domain: {int(dom_p4band.sum())} cells, overlap with P1 domain {int((dom_p4band & dom_p1).sum())}")


def dm(F, d): return np.array([np.sum(F[q][d] * w[d]) / np.sum(w[d]) for q in range(F.shape[0])])


def series(d, KC, EC):
    return pd.DataFrame({"year": YY, "month": MO, "x": YY + (MO - 0.5) / 12, "wind": dm(Z["WC"], d), "KE": dm(KC, d), "EKE": dm(EC, d),
                         "slope": np.r_[np.full(36, np.nan), dm(Z["SL"], d)]})


ref = pd.read_csv(C.OUT / "f4e" / "f4e_series.csv")
S1 = series(dom_p1, Z["KC"], Z["EC"]); S1.to_csv(C.HERE / "p4_wind_series_p1check.csv", index=False)
for c in ("wind", "KE", "EKE", "slope"):
    a, b = ref[c].values, S1[c].values; g = np.isfinite(a)
    P(f"check P1 series rebuilt from x13 vs outputs/f4e/f4e_series.csv ({c}): max rel diff {np.nanmax(np.abs(b[g] / a[g] - 1)):.2e}")
S4 = series(dom_common, KC4, EC4)
for c in ("wind", "slope"):
    a, b = ref[c].values, S4[c].values; g = np.isfinite(a)
    P(f"P4 'common' domain {c} vs P1 domain {c}: max rel diff {np.nanmax(np.abs(b[g] / a[g] - 1)):.2e}")
if (dom_common == dom_p1).all():       # same cells: take wind and slope bit-identical from f4e_series.csv (x13 is float32)
    S4["wind"] = ref.wind.values; S4["slope"] = ref.slope.values
    P("common domain = P1 domain: wind and slope columns copied from outputs/f4e/f4e_series.csv (identical inputs for both products)")
S4.to_csv(C.HERE / "p4_wind_series.csv", index=False)
series(dom_p4band, KC4, EC4).to_csv(C.HERE / "p4_wind_series_p4band.csv", index=False)
for nm, S in (("P1", S1), ("P4 common", S4)):
    P(f"{nm}: mean KE {S.KE[S.year >= 2004].mean():.1f}, EKE {S.EKE[S.year >= 2004].mean():.1f} cm2 s-2 (2004-2025)")
np.savez_compressed(C.HERE / "p4_wind_domains.npz", dom_p1=dom_p1, dom_common=dom_common, dom_p4band=dom_p4band, dmin_p4=dmin4,
                    KC4=KC4, EC4=EC4, lat=lat, lon=lon)
(C.HERE / "p4_wind_domain_log.txt").write_text("\n".join(log) + "\n")
