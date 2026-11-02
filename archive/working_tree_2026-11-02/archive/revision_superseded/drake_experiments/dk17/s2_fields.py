"""Stream the 300 monthly all-sat fields and build binned monthly series of EKE, KE and along-stream speed.
EKE = 0.5|U - Ubar|^2 + 0.5 var_intra (Ubar = 2001-2025 mean), KE = 0.5|U|^2, SPD = U.Ubar/|Ubar| (cm2/s2, cm/s).
Valid pixels: finite in all 300 months and not in the sea-ice mask (ice_mask > 0.5).
Binnings (cos-lat area weights), each circumpolar + Atlantic/Indian/Pacific:
  'madt' : time-mean ADT bins (fixed mean streamlines), -0.9..+0.5 m, 0.05 m
  'iadt' : monthly ADT bins after removing that month's domain-mean ADT anomaly (moving streamlines:
           follows meridional frontal displacements; the domain-mean rise is removed so contours are not
           pushed south by the basin-scale sea-level trend)
  'lat'  : 1 deg latitude bins 65-35S (all valid points, Martinez-Moreno style)
Output: out/s2_binned.npz
"""
import numpy as np, xarray as xr
from common import REV, OUT, acc_mask_eighth, sector_of

files = [REV / "data" / "monthly_allsat" / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
accm, lat, lon, madt = acc_mask_eighth()
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
shape = madt.shape
su = np.zeros(shape); sv = np.zeros(shape); sa = np.zeros(shape); nv = np.zeros(shape)
for f in files:
    with xr.open_dataset(f) as d:
        u = d.u_mean.values.astype(float); v = d.v_mean.values.astype(float); a = d.adt_mean.values.astype(float)
    ok = np.isfinite(u) & np.isfinite(v) & np.isfinite(a)
    su += np.where(ok, u, 0); sv += np.where(ok, v, 0); sa += np.where(ok, a, 0); nv += ok
valid = (nv == len(files)) & ~ice & np.isfinite(madt)
ub = su / np.maximum(nv, 1); vb = sv / np.maximum(nv, 1); abar = sa / np.maximum(nv, 1)
sp = np.hypot(ub, vb); eu = np.where(sp > 0, ub / np.where(sp > 0, sp, 1), 0); ev = np.where(sp > 0, vb / np.where(sp > 0, sp, 1), 0)
w = np.cos(np.deg2rad(lat))[:, None] * np.ones(lon.size)[None, :]
secid = np.broadcast_to(np.select([sector_of(lon) == "Atlantic", sector_of(lon) == "Indian"], [1, 2], 3)[None, :], shape)
madt_edges = np.round(np.arange(-0.9, 0.5001, 0.05), 3); nbA = madt_edges.size - 1
lat_edges = np.arange(-65, -34.99, 1.0); nbL = lat_edges.size - 1
V = valid.ravel(); wv = w.ravel()[V]; sv_ = secid.ravel()[V]
lab_madt = np.digitize(madt.ravel()[V], madt_edges) - 1
lab_lat = np.digitize(np.broadcast_to(lat[:, None], shape).ravel()[V], lat_edges) - 1
mean_dom = np.sum(abar.ravel()[V] * wv) / wv.sum()
# five standing meanders (dk13/b5d regions): ridge + downstream longitudes inside the ACC band; 'away' = ACC band outside
# the upstream, ridge and downstream ranges of all five (Thompson and Naveira Garabato 2014 Fig. 1b; meander vs elsewhere)
MEANDERS = {"Kerguelen Plateau": ((55, 65), (65, 80), (80, 95)), "Southeast Indian Ridge": ((120, 130), (130, 140), (140, 150)),
            "Macquarie Ridge": ((150, 155), (155, 160), (160, 175)), "Pacific–Antarctic Ridge": ((-155, -145), (-145, -115), (-115, -100)),
            "Drake Passage–Scotia Sea": ((-80, -70), (-70, -60), (-50, -35))}
inacc = (madt >= -0.6) & (madt <= 0.2) & valid
LON2 = np.broadcast_to(lon[None, :], shape); reg = np.full(shape, -1); anym = np.zeros(lon.size, bool)
for i, (nm, (u_, r_, d_)) in enumerate(MEANDERS.items()):
    m = ((lon > r_[0]) & (lon < r_[1])) | ((lon > d_[0]) & (lon < d_[1]))
    if nm.startswith("Drake"): m &= ~((lon > -60) & (lon < -50))
    reg[inacc & m[None, :]] = i
    anym |= (lon > u_[0]) & (lon < d_[1])
reg[inacc & ~anym[None, :] & (reg < 0)] = len(MEANDERS)
REGNAMES = list(MEANDERS) + ["Away from the five meanders"]
regv = reg.ravel()[V]; NR = len(REGNAMES)
regA = np.bincount(regv[regv >= 0], weights=wv[regv >= 0], minlength=NR)
regser = {v: np.zeros((len(files), NR)) for v in ("EKE", "KE", "SPD")}

nt = len(files); NS = 4
out = {k: np.zeros((nt, NS, nb)) for k, nb in (("madt", nbA), ("iadt", nbA), ("lat", nbL))}
res = {f"{v}_{b}": np.zeros((nt, NS, out[b].shape[2])) for v in ("EKE", "KE", "SPD") for b in out}
area = {b: np.zeros((nt, NS, out[b].shape[2])) for b in out}
domain_adt_anom = np.zeros(nt)


def accum(labels, nb, vals, key, t):
    ok = (labels >= 0) & (labels < nb)
    for s in range(NS):
        sel = ok if s == 0 else ok & (sv_ == s)
        A = np.bincount(labels[sel], weights=wv[sel], minlength=nb)
        if key.startswith("EKE"):
            area[key.split("_", 1)[1]][t, s] = A
        for vn, x in vals.items():
            res[f"{vn}_{key.split('_', 1)[1]}"][t, s] = np.bincount(labels[sel], weights=(wv * x)[sel], minlength=nb) / np.where(A > 0, A, np.nan)


for t, f in enumerate(files):
    with xr.open_dataset(f) as d:
        u = d.u_mean.values.astype(float); v = d.v_mean.values.astype(float)
        s2 = d.var_intra.values.astype(float); a = d.adt_mean.values.astype(float)
    eke = (0.5 * ((u - ub) ** 2 + (v - vb) ** 2) + 0.5 * s2) * 1e4
    ke = 0.5 * (u ** 2 + v ** 2) * 1e4
    spd = (u * eu + v * ev) * 1e2
    vals = {"EKE": eke.ravel()[V], "KE": ke.ravel()[V], "SPD": spd.ravel()[V]}
    av = a.ravel()[V]
    domain_adt_anom[t] = np.sum(av * wv) / wv.sum() - mean_dom
    accum(lab_madt, nbA, vals, "EKE_madt", t)
    for vn, x in vals.items():
        regser[vn][t] = np.bincount(regv[regv >= 0], weights=(wv * x)[regv >= 0], minlength=NR) / regA
    accum(lab_lat, nbL, vals, "EKE_lat", t)
    lab_i = np.digitize(av - domain_adt_anom[t], madt_edges) - 1
    accum(lab_i, nbA, vals, "EKE_iadt", t)
    if t % 60 == 0:
        print(t, flush=True)

np.savez(OUT / "s2_binned.npz", madt_edges=madt_edges, lat_edges=lat_edges, domain_adt_anom=domain_adt_anom,
         **res, **{f"area_{b}": area[b] for b in area},
         regnames=np.array(REGNAMES), reg_area=regA, **{f"reg_{v}": regser[v] for v in regser})
# sanity: ACC band (-0.6..0.2) EKE from madt bins vs a5 series
import pandas as pd
inb = (madt_edges[:-1] >= -0.6 - 1e-9) & (madt_edges[1:] <= 0.2 + 1e-9)
A = area["madt"][:, 0, inb]; E = res["EKE_madt"][:, 0, inb]
acc_eke = np.nansum(A * E, 1) / A.sum(1)
a5 = pd.read_csv(REV / "outputs" / "a5" / "a5_series_allsat.csv")
print("ACC EKE check: mean", acc_eke.mean(), "a5", a5.ACC_EKE.mean(), "corr", np.corrcoef(acc_eke, a5.ACC_EKE)[0, 1])
print("domain-mean ADT anomaly first/last yr (m):", domain_adt_anom[:12].mean(), domain_adt_anom[-12:].mean())
