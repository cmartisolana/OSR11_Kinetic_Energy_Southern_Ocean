"""Shared helpers for the Product 4 (two-satellite, 1/4 deg) counterparts of the main-text diagnostics.

Purpose : one place for paths, the stacked Product 4 monthly fields and the time axis used by the p4_*.py scripts.
Inputs  : $OSR11_DATA/monthly_twosat/YYYY/MM.nc (u_mean, v_mean, adt_mean, var_intra, ice_frac, n_days; 300 files,
          2001-2025; written by scripts/00_extract/00a_extract_monthly.py --product twosat)
Outputs : $OSR11_CACHE/p4/p4_<var>.npy (float32 stacks 300 x 120 x 1440, rebuilt on first use; not distributed)
Method  : no computation beyond stacking; the pipeline package osr11 (src/) is imported unchanged.
Run     : imported by the p4_*.py scripts (run them as python scripts/06_product4/p4_NN_*.py).
Origin  : written for the 10 Nov 2026 Supplement; moved here on 17 Nov 2026 with the paths routed through osr11.paths.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)

import numpy as np            # noqa: E402
import pandas as pd           # noqa: E402
import xarray as xr           # noqa: E402

from osr11 import paths as PA  # noqa: E402,F401
from osr11 import rstats as R  # noqa: E402,F401

DATA4 = PA.monthly_dir("twosat")                     # Product 4 monthly fields (scripts/00_extract/00a_extract_monthly.py --product twosat)
OUT = PA.OUT                                         # Product 1 / existing outputs (read)
HERE = PA.out("p4")                                  # Product 4 diagnostics are written to outputs/p4
CACHE = PA.CACHE / "p4"
FILES = [DATA4 / f"{y}" / f"{m:02d}.nc" for y in range(2001, 2026) for m in range(1, 13)]
VARS = ("u_mean", "v_mean", "adt_mean", "var_intra", "ice_frac", "n_days")

T = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
T["x"] = T.year + (T.month - 0.5) / 12
X, MON, YR = T.x.values, T.month.values, T.year.values


def grid():
    with xr.open_dataset(FILES[0]) as d:
        return d.latitude.values.astype(float), d.longitude.values.astype(float)


def stack(var: str, dtype=np.float32) -> np.ndarray:
    """(300, 120, 1440) stack of one variable of the Product 4 monthly files (cached as .npy, memory-mapped)."""
    assert all(f.exists() for f in FILES), "missing Product 4 monthly files"
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"p4_{var}.npy"
    if not f.exists():
        lat, lon = grid()
        A = np.empty((len(FILES), lat.size, lon.size), dtype)
        for k, fn in enumerate(FILES):
            with xr.open_dataset(fn) as d:
                A[k] = d[var].values.astype(dtype)
        np.save(f, A)
    return np.load(f, mmap_mode="r")
