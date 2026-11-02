"""Standing-meander paths over time from altimetry (descriptive maps).

Path definition: ADT minus its regional mean (each month; removes uniform sea-level rise,
Gille 2014) contoured at the climatological jet-core level of the meander box
(b1_meander.core_level, the same level as b1_meander_run.py), expressed relative to the
regional mean.
Panels per meander:
 (a) annual-mean paths (contour of the annual-mean relative ADT), coloured by year, on GEBCO depth
 (b) monthly path occurrence (fraction of months the monthly contour crosses each 0.25° cell),
     2001-2025, with epoch-mean paths 2001-08 / 2009-16 / 2017-25
 (c) Hovmöller: annual-mean path latitude minus its 2001-2025 mean, per 0.5° longitude
 (d) within-year SD of monthly path latitude per 0.5° longitude (as the Gulf Stream
     destabilisation diagnostic, Andres 2016; Sánchez-Román et al. 2024, OSR8), epoch means
Path latitude at a longitude: crossing of the core level closest to the 25-yr-mean path.
Usage: python b1_maps.py <revision/data dir> <PAR|SEIR150|KERG>
-> outputs/figures/figS_b1_paths_<name>.png, outputs/b1/paths_<name>.nc
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
import contourpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b1_meander as B

REG = {"PAR": dict(reg=(-156.0, -74.0, -64.0, -44.0), mea=(-135.13, -109.38), mlat=(-62.0, -45.0),
                   secs=(-150.0, -135.13, -109.38, -80.0), title="Pacific–Antarctic Ridge"),
       "SEIR150": dict(reg=(122.0, 166.0, -64.0, -44.0), mea=(140.0, 150.0), mlat=(-60.0, -46.0),
                       secs=(125.0, 138.0, 140.0, 150.0, 158.0), title="SE Indian Ridge (150°E)"),
       "KERG": dict(reg=(64.0, 100.0, -62.0, -42.0), mea=(76.0, 84.0), mlat=(-58.0, -48.0),
                    secs=(70.0, 77.0, 76.0, 84.0, 94.0), title="Kerguelen (Fawn Trough branch)")}
EPOCHS = [(2001, 2008), (2009, 2016), (2017, 2025)]
ECOL = ["#1f5fa8", "#7f7f7f", "#c0392b"]

data, name = Path(sys.argv[1]), sys.argv[2]
c = REG[name]
lo0, lo1, la0, la1 = c["reg"]
years = np.repeat(np.arange(2001, 2026), 12)
A, U, V = [], [], []
for y in range(2001, 2026):
    for m in range(1, 13):
        with xr.open_dataset(data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc") as ds:
            s = ds.sel(longitude=slice(lo0, lo1), latitude=slice(la0, la1))
            A.append(s.adt_mean.values); U.append(s.u_mean.values); V.append(s.v_mean.values)
lat = s.latitude.values.astype(float); lon = s.longitude.values.astype(float)
A = np.array(A, float); U = np.array(U, float); V = np.array(V, float)
Ab, ub, vb = np.nanmean(A, 0), np.nanmean(U, 0), np.nanmean(V, 0)
mbox = dict(lon0=c["mea"][0], lon1=c["mea"][1], lat0=c["mlat"][0], lat1=c["mlat"][1])
core_box, spd = B.core_level(Ab, ub, vb, lat, lon, mbox)
la_m = (lat >= mbox["lat0"]) & (lat <= mbox["lat1"]); lo_m = (lon >= mbox["lon0"]) & (lon <= mbox["lon1"])
regmean = np.array([np.nanmean(a) for a in A]); boxmean = np.array([np.nanmean(a[np.ix_(la_m, lo_m)]) for a in A])
level = core_box + float(np.mean(boxmean - regmean))          # same contour, relative to the regional mean
R = A - regmean[:, None, None]
print(name, "core level: box-relative", round(core_box, 3), "region-relative", round(level, 3), flush=True)


def lines_of(field):
    gen = contourpy.contour_generator(lon, lat, np.ma.masked_invalid(field), name="serial")
    return [L for L in gen.lines(level) if np.ptp(L[:, 0]) > 2.0]


# path latitude per 0.5° longitude bin
bins = np.arange(lo0, lo1 + 0.01, 0.5); bc = 0.5 * (bins[:-1] + bins[1:])
def crossings(field):
    """latitudes where `field` crosses `level` in each grid column -> list per column"""
    f = field - level
    out = []
    for j in range(f.shape[1]):
        col = f[:, j]; ok = np.isfinite(col[:-1]) & np.isfinite(col[1:])
        idx = np.nonzero(ok & (np.sign(col[:-1]) != np.sign(col[1:])))[0]
        out.append(lat[idx] - col[idx] * (lat[idx + 1] - lat[idx]) / (col[idx + 1] - col[idx]))
    return out
ref = crossings(np.nanmean(R, 0))
refbin = np.full(bc.size, np.nan)
def pathlat(field, refcol=None):
    cr = crossings(field); y = np.full(lon.size, np.nan)
    for j, cj in enumerate(cr):
        if cj.size == 0:
            continue
        r0 = refcol[j] if refcol is not None else np.nan
        y[j] = cj[np.argmin(np.abs(cj - r0))] if np.isfinite(r0) else (cj[0] if cj.size == 1 else np.nan)
    return pd.Series(y).groupby(np.digitize(lon, bins) - 1).mean().reindex(range(bc.size)).values
refcol = np.array([cj[0] if cj.size == 1 else (np.median(cj) if cj.size else np.nan) for cj in ref])
mpath = np.array([pathlat(R[k], refcol) for k in range(R.shape[0])])          # month x lon-bin
ann = pd.DataFrame(mpath).groupby(years).mean().values                         # annual mean of monthly path lat
annsd = pd.DataFrame(mpath).groupby(years).std().values                        # within-year SD
meanpath = np.nanmean(mpath, 0)
xr.Dataset({"path_lat": (("time", "lon"), mpath), "annual_sd": (("year", "lon"), annsd)},
           coords=dict(time=np.arange(300), year=np.arange(2001, 2026), lon=bc),
           attrs=dict(level_region_relative=level)).to_netcdf(HERE / "outputs" / "b1" / f"paths_{name}.nc")

# occurrence: monthly contour rasterised on 0.25° cells
g_lo = np.arange(lo0, lo1 + 0.01, 0.25); g_la = np.arange(la0, la1 + 0.01, 0.25)
occ = np.zeros((g_la.size - 1, g_lo.size - 1))
for k in range(R.shape[0]):
    hit = np.zeros_like(occ, bool)
    for L in lines_of(R[k]):
        # densify
        seg = np.hypot(np.diff(L[:, 0]), np.diff(L[:, 1])); n = np.maximum(1, (seg / 0.05).astype(int))
        pts = np.vstack([np.linspace(L[i], L[i + 1], n[i], endpoint=False) for i in range(len(L) - 1)])
        ii = np.digitize(pts[:, 1], g_la) - 1; jj = np.digitize(pts[:, 0], g_lo) - 1
        ok = (ii >= 0) & (ii < occ.shape[0]) & (jj >= 0) & (jj < occ.shape[1])
        hit[ii[ok], jj[ok]] = True
    occ += hit
occ /= R.shape[0]
# save occurrence and epoch-mean relative ADT for the main-text figure
EP = np.array([np.nanmean(R[(years >= e0) & (years <= e1)], 0) for e0, e1 in EPOCHS])
xr.Dataset({"occ": (("glat", "glon"), occ), "epoch_rel_adt": (("epoch", "latitude", "longitude"), EP.astype(np.float32))},
           coords=dict(glat=0.5 * (g_la[:-1] + g_la[1:]), glon=0.5 * (g_lo[:-1] + g_lo[1:]), latitude=lat, longitude=lon,
                       epoch=[f"{a}-{b}" for a, b in EPOCHS]),
           attrs=dict(level=level)).to_netcdf(HERE / "outputs" / "b1" / f"occ_{name}.nc")

bat = xr.open_dataset(data / "bathymetry_so.nc").elevation.sel(longitude=slice(lo0, lo1), latitude=slice(la0, la1))
fig = plt.figure(figsize=(11, 10))
ax1 = fig.add_axes([0.07, 0.56, 0.37, 0.38]); ax2 = fig.add_axes([0.57, 0.56, 0.37, 0.38])
ax3 = fig.add_axes([0.07, 0.07, 0.37, 0.38]); ax4 = fig.add_axes([0.57, 0.07, 0.37, 0.38])
for ax in (ax1, ax2):
    ax.pcolormesh(bat.longitude, bat.latitude, -bat, cmap="Greys", vmin=0, vmax=7000, shading="auto", rasterized=True)
    ax.contour(bat.longitude, bat.latitude, bat, levels=[-3000], colors="0.35", linewidths=0.5)
    for x in c["secs"]:
        ax.axvline(x, color="k", lw=0.5, ls=":")
    ax.set_ylim(la0, la1); ax.set_xlim(lo0, lo1); ax.set_ylabel("latitude")
cmap = plt.get_cmap("viridis")
for i, y in enumerate(range(2001, 2026)):
    for L in lines_of(np.nanmean(R[years == y], 0)):
        ax1.plot(L[:, 0], L[:, 1], color=cmap(i / 24), lw=0.9)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(2001, 2025)); plt.colorbar(sm, ax=ax1, pad=0.01, label="year")
ax1.set_title(f"(a) {c['title']}: annual-mean paths", fontsize=9)
oc = np.where(occ > 0, occ, np.nan)
pc = ax2.pcolormesh(g_lo, g_la, oc, cmap="magma_r", vmin=0, vmax=0.6, shading="flat", rasterized=True)
plt.colorbar(pc, ax=ax2, pad=0.01, label="fraction of months")
for (e0, e1), col in zip(EPOCHS, ECOL):
    for L in lines_of(np.nanmean(R[(years >= e0) & (years <= e1)], 0)):
        ax2.plot(L[:, 0], L[:, 1], color=col, lw=1.4, label=f"{e0}–{e1}")
h, l = ax2.get_legend_handles_labels(); d = dict(zip(l, h)); ax2.legend(d.values(), d.keys(), fontsize=7, loc="lower left")
ax2.set_title("(b) monthly path occurrence and epoch-mean paths", fontsize=9)
an = ann - meanpath[None, :]
pc3 = ax3.pcolormesh(bc, np.arange(2001, 2026), an, cmap="RdBu_r", vmin=-1.5, vmax=1.5, shading="nearest")
plt.colorbar(pc3, ax=ax3, pad=0.01, label="path latitude anomaly (°; + = north)")
ax3.set_xlim(lo0, lo1); ax3.set_ylabel("year"); ax3.set_xlabel("longitude")
ax3.set_title("(c) annual-mean path latitude anomaly", fontsize=9)
for (e0, e1), col in zip(EPOCHS, ECOL):
    sel = (np.arange(2001, 2026) >= e0) & (np.arange(2001, 2026) <= e1)
    ax4.plot(bc, np.nanmean(annsd[sel], 0), color=col, lw=1.3, label=f"{e0}–{e1}")
ax4.set_xlim(lo0, lo1); ax4.set_xlabel("longitude"); ax4.set_ylabel("within-year SD of path latitude (°)")
ax4.legend(fontsize=7); ax4.grid(alpha=0.3, lw=0.4)
for x in c["secs"]:
    ax4.axvline(x, color="k", lw=0.5, ls=":")
ax4.set_title("(d) month-to-month path variability (OSR8-style)", fontsize=9)
fig.savefig(HERE / "outputs" / "figures" / f"figS_b1_paths_{name}.png", dpi=170)
print("saved", name)
