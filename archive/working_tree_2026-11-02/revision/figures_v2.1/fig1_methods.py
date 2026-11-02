"""Methods figure (Figure 1 of manuscript v3.26): the half-power method, one example (after panels a-b of fig2.py, v2.6).
a  KE of the monthly-mean geostrophic velocity, January 2001, with all half-power points detected in the ACC band;
   star = example longitude of b.
b  meridional KE profile at that longitude, January 2001: threshold, envelope, 25 % of peak, phi_S, phi_N and the
   half-power latitude phi_hp (Eq. 2); dashed: 2001-2025 mean KE at the same longitude.
Baseline parameters (decision D1): threshold 200 cm2 s-2, minimum width 100 km, envelope boundary 25 % of peak.
Usage: python fig_methods_v326.py <revision dir> <figure dir> [example lon, default 147.9]
"""
import sys
from pathlib import Path
import numpy as np, xarray as xr
REV, FIGD = Path(sys.argv[1]), Path(sys.argv[2]); LON_EX = float(sys.argv[3]) if len(sys.argv) > 3 else 147.9
sys.path.insert(0, str(REV))
import style_nc as S
from style_nc import add_ice_mask
import fronts as F
import cartopy.crs as ccrs, cartopy.feature as cfeature
plt = S.plt; C = S.C; PC = ccrs.PlateCarree()
P = dict(thr=0.02, min_w=100.0, frac=0.25); MAX_ENV = 8; FI = {k: i for i, k in enumerate(F.FIELDS)}

mf = sorted((REV / "data" / "monthly_allsat").glob("*/*.nc")); mf = [f for f in mf if 2001 <= int(f.parent.name) <= 2025]
with xr.open_dataset(REV / "data" / "monthly_allsat" / "2001" / "01.nc") as d:
    lat, lon = d.latitude.values.astype(float), d.longitude.values.astype(float)
    fld = (0.5 * (d.u_mean.values ** 2 + d.v_mean.values ** 2)).astype(float)          # m2 s-2
JEX = int(np.argmin(np.abs(lon - LON_EX)))
clim = np.zeros(lat.size); n = 0
for f in mf:                                                                             # 2001-2025 mean KE at JEX
    with xr.open_dataset(f) as d:
        clim += np.nan_to_num(0.5 * (d.u_mean.values[:, JEX] ** 2 + d.v_mean.values[:, JEX] ** 2)); n += 1
clim = clim / n * 1e4
adt = xr.open_dataset(REV.parent / "outputs" / "mean_adt.nc")["mean_adt"].values
acc = np.isfinite(adt) & (adt >= -0.6) & (adt <= 0.2)
ice = xr.open_dataset(REV.parent / "outputs" / "sea_ice_mask.nc").ice_mask.values > 0.5
dist = F.meridional_dists(lat)
cnt, tab = F.detect_field(fld, lat, dist, P["thr"], P["min_w"], P["frac"], MAX_ENV)
jj, kk = np.nonzero(np.arange(MAX_ENV)[None, :] < cnt[:, None])
lhp = tab[jj, kk, FI["lat_hp"]]
inacc = acc[np.clip(np.searchsorted(lat, lhp), 0, lat.size - 1), jj]
ex = tab[JEX, :cnt[JEX]]
ex = ex[acc[np.clip(np.searchsorted(lat, ex[:, FI["lat_hp"]]), 0, lat.size - 1), JEX]][0]

fig = plt.figure(figsize=(S.W2, 2.75))
gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.2], wspace=0.25)
# a
ax = fig.add_subplot(gs[0], projection=ccrs.SouthPolarStereo())
S.polar(ax, ccrs, cfeature, extent_lat=-70)
LON, LAT = np.meshgrid(lon, lat)
imk = ax.pcolormesh(LON, LAT, np.where(ice, np.nan, fld * 1e4), transform=PC, cmap="Greys", vmin=0, vmax=1500, shading="auto", rasterized=True, zorder=0)
add_ice_mask(ax, ccrs, lon, lat, ice)
ax.scatter(lon[jj[inacc]], lhp[inacc], s=0.5, c=C["pos"], lw=0, transform=PC, zorder=5)
ax.plot(lon[JEX], ex[FI["lat_hp"]], "*", ms=7, color=C["obs"], mec="k", mew=0.5, transform=PC, zorder=7)
S.plabel(ax, "a", "Half-power points, January 2001")
cbk = fig.colorbar(imk, ax=ax, orientation="horizontal", shrink=0.6, pad=0.04, extend="max", aspect=30)
cbk.set_label("KE of the monthly-mean flow (cm$^2$ s$^{-2}$)", fontsize=5.5); cbk.ax.tick_params(labelsize=5, length=2, width=0.4); cbk.outline.set_linewidth(0.4)
# b
ax = fig.add_subplot(gs[1])
prof = fld[:, JEX] * 1e4
s_, e_ = int(ex[FI["i_start"]]), int(ex[FI["i_end"]])
ls_, ln_, pv, hp = ex[FI["lat_s"]], ex[FI["lat_n"]], ex[FI["peak_val"]] * 1e4, ex[FI["lat_hp"]]
ax.plot(lat, clim, color=C["grey"], lw=0.8, ls=(0, (3, 1.5)), label="2001–2025 mean")
ax.plot(lat, prof, color=C["ink"], lw=1.0, label="January 2001")
sel = (lat >= ls_) & (lat <= ln_)
ax.fill_between(lat[sel], 0, prof[sel], color=C["pos"], alpha=0.18, lw=0, label="Envelope")
ax.axhline(P["thr"] * 1e4, color=C["ke"], lw=0.7, ls=":", label="Threshold (200 cm$^2$ s$^{-2}$)")
ax.plot([lat[s_] - 0.5, lat[e_] + 0.5], [P["frac"] * pv] * 2, color=C["ke"], lw=0.7, ls="-.", label="25 % of peak")
ax.axvline(hp, color=C["pos"], lw=1.0, label="Half-power latitude")
for v_, lab in ((ls_, r"$\varphi_S$"), (ln_, r"$\varphi_N$"), (hp, r"$\varphi_{hp}$")):
    if v_ != hp: ax.axvline(v_, color=C["pos"], lw=0.5, ls="--")
    ax.text(v_, pv * 1.02, lab, ha="center", va="bottom", fontsize=6, color=C["pos"])
ax.set_xlim(hp - 5.5, hp + 11); ax.set_ylim(0, pv * 1.12)
ax.set_xlabel("Latitude (°)"); ax.set_ylabel("KE (cm$^2$ s$^{-2}$)")
ax.legend(loc="upper right", fontsize=5.6)
lab = f"{lon[JEX]:.1f}°E" if lon[JEX] >= 0 else f"{-lon[JEX]:.1f}°W"
S.plabel(ax, "b", f"KE profile at {lab}, January 2001", dx=-30)
S.save(fig, FIGD, "Figure1")
print(dict(lon=lon[JEX], lat_hp=hp, lat_s=ls_, lat_n=ln_, peak=pv, n_acc_points=int(inacc.sum())))
