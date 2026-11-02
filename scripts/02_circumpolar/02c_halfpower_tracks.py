"""Matched half-power tracks for the baseline parameters (900 longitudes) -> outputs/a2/tracks_baseline.npz (used by fig2.py).

Purpose : matched half-power tracks for the baseline parameters only (900 ACC longitudes); same result as the baseline set of
          02b, written separately because later steps read it.
Inputs  : stored monthly KE fields via $OSR11_CACHE (01b, osr11.halfpower.load), outputs/base/mean_adt.nc
Outputs : outputs/a2/tracks_baseline.npz (19 MB; distributed: read by 02e, 02j, 03t and 05a A2)
Run     : python scripts/02_circumpolar/02c_halfpower_tracks.py
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np
from osr11 import halfpower as A
A.OUT.mkdir(parents=True, exist_ok=True)
ke, lat, lon, acc, t = A.load()
kemean = np.nanmean(np.asarray(ke, np.float64), axis=0); clim_max = np.nanmax(kemean, axis=0)
p = A.BASE
counts, tab = A.run_set(ke, lat, lon, p, clim_max)
cc, ct = A.clim_envelopes(kemean, lat, lon, p, clim_max)
arr, meta = A.matched_tracks(counts, tab, cc, ct, lat, lon, acc)
np.savez_compressed(A.OUT / "tracks_baseline.npz", **arr, j=meta.j.values, c=meta.c.values, keep=meta.keep.values,
                    coverage=meta.coverage.values, lat_hp_clim=meta.lat_hp_clim.values, acc=meta.acc.values)
print("tracks kept", int(meta.keep.sum()), "of acc", int(meta.acc.sum()))
