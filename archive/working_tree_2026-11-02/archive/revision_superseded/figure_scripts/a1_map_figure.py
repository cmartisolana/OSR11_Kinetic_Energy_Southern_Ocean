"""Fig. 1 candidate: KE and EKE trend maps (2001-2025, 1 deg cells) with field significance.

Slope (cm2 s-2 yr-1) of deseasonalised Theil-Sen; dots = cells significant after
Benjamini-Hochberg FDR (alpha_FDR = 0.10) on Gaussian-tail sieve p-values (a1_maps.py);
grey contour = 3000 m isobath (GEBCO 2026); black contours = ACC band (mean ADT -0.6, +0.2 m).
Usage: python a1_map_figure.py <revision/data dir>
-> outputs/figures/fig1_trend_maps_fdr.png
"""
import sys
from pathlib import Path
import numpy as np, xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
data = Path(sys.argv[1])
M = xr.open_dataset(HERE / "outputs" / "a1_maps.nc")
bat = xr.open_dataset(data / "bathymetry_so.nc").elevation.coarsen(latitude=4, longitude=4, boundary="trim").mean()
# mean ADT of the all-sat product for the ACC-band contours
s = None
for y in range(2001, 2026):
    for m in range(1, 13):
        with xr.open_dataset(data / "monthly_allsat" / f"{y}" / f"{m:02d}.nc") as ds:
            a = ds.adt_mean.coarsen(latitude=4, longitude=4, boundary="trim").mean().load()
        s = a if s is None else s + a
adt = s / 300

fig, axs = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True)
for ax, v, lab in zip(axs, ("KE", "EKE"), ("(a) KE of monthly-mean flow", "(b) EKE")):
    sl = M[f"{v}_slope"]; fdr = M[f"{v}_fdr"]
    lim = float(np.nanpercentile(np.abs(sl), 98))
    pc = ax.pcolormesh(M.lon, M.lat, sl, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="nearest")
    ax.contour(bat.longitude, bat.latitude, bat, levels=[-3000], colors="0.45", linewidths=0.5)
    ax.contour(adt.longitude, adt.latitude, adt, levels=[-0.6, 0.2], colors="k", linewidths=0.7)
    jj, ii = np.nonzero(fdr.values == 1)
    ax.plot(M.lon.values[ii], M.lat.values[jj], ".", ms=1.6, color="k")
    nf = int(np.isfinite(M[f"{v}_p"]).sum()); nsig = int((fdr == 1).sum())
    ax.set_title(f"{lab}: trend 2001–2025; dots = FDR-significant ({nsig} of {nf} cells)", fontsize=9)
    ax.set_ylabel("latitude"); ax.set_ylim(-65, -35)
    plt.colorbar(pc, ax=ax, label="cm² s⁻² yr⁻¹", pad=0.01)
axs[-1].set_xlabel("longitude")
fig.tight_layout()
fig.savefig(HERE / "outputs" / "figures" / "fig1_trend_maps_fdr.png", dpi=200)
print("saved")
