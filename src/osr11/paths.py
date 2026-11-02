"""Single place for every path used by the pipeline.

All scripts import their locations from here; nothing else in the repository builds paths from
``__file__`` or the current working directory. Every location can be overridden by an
environment variable (documented in README.md):

    OSR11_DATA      raw and large derived inputs (not distributed)   default <repo>/data
    OSR11_OUTPUTS   processed outputs (small; distributed)            default <repo>/outputs
    OSR11_CACHE     stacked monthly KE/EKE fields (1.6 GB, rebuilt)   default <repo>/cache
    OSR11_FIGURES   rendered figures                                  default <repo>/rendered
    GV23_DIR        unpacked Gutierrez-Villanueva et al. (2023) archive
                                                                      default $OSR11_DATA/TrendsDrakePassageTransp
    ARGO_RG_DIR     Roemmich-Gilson fields converted by 00d_argo_rg_extract.py
                                                                      default $OSR11_DATA/argo_rg

Layout under OUT (one folder per analysis step, names kept from the working repository):
    base/      mean_adt.nc, sea_ice_mask.nc, clim_allsat.nc, bathymetry_so.nc
    a1 a2 a3 a5 retest stats x1 m1 d1 f3w f4 f4m f4e f4s ...
    drake/dkN  Drake Passage steps (dk1 ... dk40)
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _env(name: str, default: Path) -> Path:
    v = os.environ.get(name)
    return Path(v).expanduser().resolve() if v else default


DATA = _env("OSR11_DATA", ROOT / "data")
OUT = _env("OSR11_OUTPUTS", ROOT / "outputs")
CACHE = _env("OSR11_CACHE", ROOT / "cache")
FIG = _env("OSR11_FIGURES", ROOT / "rendered")
GV23_DIR = _env("GV23_DIR", DATA / "TrendsDrakePassageTransp")
ARGO_RG = _env("ARGO_RG_DIR", DATA / "argo_rg")

# base fields (small, distributed in outputs/base)
BASE = OUT / "base"
MEAN_ADT = BASE / "mean_adt.nc"            # 2001-2025 mean ADT (01a_preprocessing.py)
ICE_MASK = BASE / "sea_ice_mask.nc"        # permanent JJA ice mask (01a_preprocessing.py)
CLIM = BASE / "clim_allsat.nc"             # 2001-2025 mean of the monthly fields (01d_climatology.py)
BATHY = BASE / "bathymetry_so.nc"          # GEBCO block means (00c_prepare_bathymetry.py)

# raw inputs (not distributed)
WIND = DATA / "wind_taux.nc"               # 00b_extract_wind.py
HP_POINTS = DATA / "monthly_half-power_points"              # 01b_halfpower_fields.py (all-sat)
HP_POINTS_TWOSAT = DATA / "monthly_half-power_points_two_sat"  # 01c_halfpower_fields_twosat.py
PREPROC = DATA / "preprocessing.nc"        # 01a_preprocessing.py (long-term mean u, v; used by 01b/01c only)


def monthly_dir(product: str = "allsat") -> Path:
    """Directory of the monthly fields written by 00a_extract_monthly.py."""
    return DATA / f"monthly_{product}"


def monthly_files(product: str = "allsat", y0: int = 2001, y1: int = 2025) -> list[Path]:
    """The 12*(y1-y0+1) monthly files in time order (no existence check)."""
    d = monthly_dir(product)
    return [d / f"{y}" / f"{m:02d}.nc" for y in range(y0, y1 + 1) for m in range(1, 13)]


def out(*parts: str, mkdir: bool = True) -> Path:
    """OUT/<parts...>; creates the directory (or the parent of a file path) when mkdir is True."""
    p = OUT.joinpath(*parts)
    if mkdir:
        (p.parent if p.suffix else p).mkdir(parents=True, exist_ok=True)
    return p


def dk(n: int | str, mkdir: bool = True) -> Path:
    """Output folder of Drake Passage step dkN: OUT/drake/dkN."""
    return out("drake", f"dk{n}", mkdir=mkdir)
