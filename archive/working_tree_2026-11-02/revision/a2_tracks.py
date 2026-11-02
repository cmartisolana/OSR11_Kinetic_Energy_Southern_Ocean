"""Matched half-power tracks for the baseline parameters (900 longitudes) -> outputs/a2/tracks_baseline.npz (used by fig2.py).
Usage: OSR11_OUTPUTS=../outputs OSR11_CACHE=<cache> python a2_tracks.py"""
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import a2_run as A
ke, lat, lon, acc, t = A.load()
kemean = np.nanmean(np.asarray(ke, np.float64), axis=0); clim_max = np.nanmax(kemean, axis=0)
p = A.BASE
counts, tab = A.run_set(ke, lat, lon, p, clim_max)
cc, ct = A.clim_envelopes(kemean, lat, lon, p, clim_max)
arr, meta = A.matched_tracks(counts, tab, cc, ct, lat, lon, acc)
np.savez_compressed(A.OUT / "tracks_baseline.npz", **arr, j=meta.j.values, c=meta.c.values, keep=meta.keep.values,
                    coverage=meta.coverage.values, lat_hp_clim=meta.lat_hp_clim.values, acc=meta.acc.values)
print("tracks kept", int(meta.keep.sum()), "of acc", int(meta.acc.sum()))
