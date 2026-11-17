"""Figure 1: the half-power method, one example month (Product 1, January 2001).
a  KE of the monthly-mean geostrophic velocity, January 2001, with the half-power points detected in the ACC band (purple);
   star = example longitude of b.
b  meridional KE profile at that longitude, January 2001 (solid) and its 2001-2025 mean (dashed): threshold 200 cm2 s-2,
   envelope (shaded) bounded at 25 % of the peak, envelope limits phi_S and phi_N, half-power latitude phi_hp (Eq. 2).
Baseline parameters (decision D1): threshold 200 cm2 s-2, minimum width 100 km, envelope boundary 25 % of peak.
Inputs: revision/data/monthly_allsat/2001/01.nc (map, profile); all 300 monthly files 2001-2025 for the dashed mean profile
(cached in outputs/fig1/clim_profile_<lon>.npy, then reused); outputs/mean_adt.nc, outputs/sea_ice_mask.nc;
revision/fronts.py (detection).

Publication repository: ported from figures_final/scripts/fig1_methods.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : $OSR11_DATA/monthly_allsat/2001/01.nc (one month, raw data), outputs/fig1/clim_profile_147.9.npy (2001-2025 mean KE
          profile at 147.9 E; recomputed from all 300 monthly files if absent), outputs/base/mean_adt.nc, sea_ice_mask.nc
Outputs : rendered/{png,pdf}/Figure1
Run     : python figures/fig1_halfpower_method.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, xarray as xr
from osr11 import style as S, paths as PA  # noqa: E402  (P is the half-power parameter dict below)
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else PA.FIG; LON_EX = float(sys.argv[2]) if len(sys.argv) > 2 else 147.9
from osr11 import fronts as F  # noqa: E402
import cartopy.crs as ccrs, cartopy.feature as cfeature  # noqa: E402
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
P = dict(thr=0.20, min_w=100.0, frac=0.25, relative=True); MAX_ENV = 8; FI = {k: i for i, k in enumerate(F.FIELDS)}   # 5 Nov 2026: longitude-dependent threshold
CLIM_MAX = np.load(PA.OUT / "x6" / "x6_relthr_tracks.npz")["clim_max"]            # 2001-2025 maximum KE in the ACC band per longitude (m2 s-2; x6)

_f = PA.monthly_dir("allsat") / "2001" / "01.nc"
if not _f.exists():
    _f = PA.OUT / "fig1" / "monthly_allsat_2001_01_uv.nc"   # distributed copy of the January 2001 monthly-mean u, v
with xr.open_dataset(_f) as d:
    lat, lon = d.latitude.values.astype(float), d.longitude.values.astype(float)
    fld = (0.5 * (d.u_mean.values ** 2 + d.v_mean.values ** 2)).astype(float)          # m2 s-2
JEX = int(np.argmin(np.abs(lon - LON_EX)))
cache = PA.OUT / "fig1" / f"clim_profile_{LON_EX:.1f}.npy"   # distributed
if cache.exists():
    clim = np.load(cache)
else:                                                                                  # 2001-2025 mean KE at JEX (300 files)
    mf = PA.monthly_files("allsat")
    assert len(mf) == 300, f"{len(mf)} monthly files"
    clim = np.zeros(lat.size)
    for f in mf:
        with xr.open_dataset(f) as d:
            clim += np.nan_to_num(0.5 * (d.u_mean.values[:, JEX] ** 2 + d.v_mean.values[:, JEX] ** 2))
    clim = clim / len(mf) * 1e4; cache.parent.mkdir(parents=True, exist_ok=True); np.save(cache, clim)
adt = xr.open_dataset(PA.MEAN_ADT)["mean_adt"].values
acc = np.isfinite(adt) & (adt >= -0.6) & (adt <= 0.2)
ice = xr.open_dataset(PA.ICE_MASK).ice_mask.values > 0.5
dist = F.meridional_dists(lat)
cnt, tab = F.detect_field(fld / CLIM_MAX[None, :], lat, dist, P["thr"], P["min_w"], P["frac"], MAX_ENV)
THR_EX = P["thr"] * CLIM_MAX[JEX] * 1e4   # threshold at the example longitude, cm2 s-2
jj, kk = np.nonzero(np.arange(MAX_ENV)[None, :] < cnt[:, None])
lhp = tab[jj, kk, FI["lat_hp"]]
inacc = acc[np.clip(np.searchsorted(lat, lhp), 0, lat.size - 1), jj]
ex = tab[JEX, :cnt[JEX]]
ex = ex[acc[np.clip(np.searchsorted(lat, ex[:, FI["lat_hp"]]), 0, lat.size - 1), JEX]][0]

H1 = 122.0; fig = S.figure(S.W2, H1)
# a ---------------------------------------------------------------- map
D = 57.0
ax = S.axmm(fig, 9, H1 - 66.0, D, D, projection=ccrs.SouthPolarStereo())
S.polar(ax, ccrs, cfeature, keys=False)
LON, LAT = np.meshgrid(lon, lat)
imk = ax.pcolormesh(LON, LAT, np.where(ice, np.nan, fld * 1e4), transform=PC, cmap="Greys", vmin=0, vmax=1500, shading="auto",
                    rasterized=True, zorder=0)
S.add_ice(ax, ccrs, lon, lat, ice)
ax.scatter(lon[jj[inacc]], lhp[inacc], s=0.6, c=C["pos"], lw=0, transform=PC, zorder=5)
ax.plot(lon[JEX], ex[FI["lat_hp"]], "*", ms=7, color=C["obs"], mec="k", mew=0.5, transform=PC, zorder=7)
S.letter(fig, 2, H1 - 3.5, "a", "Half-power points, January 2001")
S.cbar(fig, imk, 13, H1 - 72.5, 48, "KE of the monthly-mean flow (cm$^2$ s$^{-2}$)", extend="max", ticks=[0, 500, 1000, 1500])
# b ---------------------------------------------------------------- profile
ax = S.axmm(fig, 90, H1 - 70, 76, 62)
prof = fld[:, JEX] * 1e4
s_, e_ = int(ex[FI["i_start"]]), int(ex[FI["i_end"]])
ls_, ln_, pv, hp = ex[FI["lat_s"]], ex[FI["lat_n"]], ex[FI["peak_val"]] * CLIM_MAX[JEX] * 1e4, ex[FI["lat_hp"]]   # peak_val is relative to CLIM_MAX
ax.plot(lat, clim, color=C["grey"], lw=0.8, ls=(0, (3, 1.5)), label="2001–2025 mean")
ax.plot(lat, prof, color=C["ink"], lw=1.0, label="January 2001")
sel = (lat >= ls_) & (lat <= ln_)
ax.fill_between(lat[sel], 0, prof[sel], color=C["pos"], alpha=0.18, lw=0, label="envelope")
ax.axhline(THR_EX, color=C["ke"], lw=0.7, ls=":", label=f"threshold: 20 % of the 2001–2025\nmaximum KE here ({THR_EX:.0f} cm$^2$ s$^{{-2}}$)")
ax.plot([lat[s_] - 0.5, lat[e_] + 0.5], [P["frac"] * pv] * 2, color=C["ke"], lw=0.7, ls="-.", label="25 % of peak")
ax.axvline(hp, color=C["pos"], lw=1.0, label="half-power latitude")
for v_, lab in ((ls_, r"$\varphi_S$"), (ln_, r"$\varphi_N$"), (hp, r"$\varphi_{hp}$")):
    if v_ != hp: ax.axvline(v_, color=C["pos"], lw=0.5, ls="--")
    ax.text(v_, pv * 1.02, lab, ha="center", va="bottom", fontsize=S.FS_LAB, color=C["pos"],
            bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"), zorder=6)
ax.set_xlim(hp - 5.5, hp + 11); ax.set_ylim(0, pv * 1.12)
xt = np.arange(np.ceil(hp - 5.5), hp + 11, 2); ax.set_xticks(xt); ax.set_xticklabels([S.latfmt(v) for v in xt])
ax.set_xlabel("Latitude"); ax.set_ylabel("KE (cm$^2$ s$^{-2}$)")
ax.legend(loc="upper right", fontsize=S.FS_ANN)
lab = f"{lon[JEX]:.1f}°E" if lon[JEX] >= 0 else f"{-lon[JEX]:.1f}°W"
S.letter(fig, 80, H1 - 3.5, "b", f"KE profile at {lab}, January 2001")
# c ---------------------------------------------------------------- threshold by longitude
ax = S.axmm(fig, 16, 10, 144, 23)
ax.axvspan(-86, -40, color="#EEEEEE", lw=0, zorder=-1)
ax.plot(lon, P["thr"] * CLIM_MAX * 1e4, color=C["ke"], lw=0.8, label="20 % of the 2001–2025 maximum KE in the ACC band")
ax.axhline(200, color=C["grey"], lw=0.7, ls=(0, (3, 1.5)), label="fixed 200 cm$^2$ s$^{-2}$ (sensitivity set, Supporting Information S4)")
ax.plot(lon[JEX], THR_EX, "*", ms=7, color=C["obs"], mec="k", mew=0.5, zorder=7)
S.lon_axis(ax); ax.set_yscale("log"); ax.set_ylim(20, 1500); ax.set_yticks([30, 100, 300, 1000]); ax.set_yticklabels(["30", "100", "300", "1000"])
ax.set_ylabel("Threshold\n(cm$^2$ s$^{-2}$)"); S.keys_on_axis(ax)
ax.legend(loc="lower left", ncol=2, fontsize=S.FS_ANN, frameon=False, handlelength=2.0, bbox_to_anchor=(0.0, -0.02))
S.letter(fig, 2, 38.5, "c", "KE threshold by longitude")
S.save(fig, FIGD, "Figure1")
print(dict(lon=round(float(lon[JEX]), 3), lat_hp=round(float(hp), 3), lat_s=float(ls_), lat_n=float(ln_), peak=round(float(pv), 1), n_acc_points=int(inacc.sum())))
