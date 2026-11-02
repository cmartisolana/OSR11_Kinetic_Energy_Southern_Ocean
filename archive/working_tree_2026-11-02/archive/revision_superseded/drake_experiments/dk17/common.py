"""Shared paths/helpers for the eddy-saturation quick-look prototypes."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

REV = Path("/mnt/user-data/uploads/OSR11/repository/revision")
REPO_OUT = REV.parent / "outputs"
OUT = Path("/home/claude/work/drake/out/dk17"); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV))
import rstats as R  # noqa: E402

SECTORS = ("ACC", "Atlantic", "Indian", "Pacific")


def sector_of(lon):
    l = np.mod(lon, 360)
    return np.where((l >= 20) & (l < 147), "Indian", np.where((l >= 147) & (l < 290), "Pacific", "Atlantic"))


def acc_mask_eighth():
    """ACC band on the 1/8 deg altimetry grid: mean ADT in [-0.6, 0.2] m and not ice (ice_mask <= 0.5)."""
    a = xr.open_dataset(REPO_OUT / "mean_adt.nc").mean_adt
    ice = xr.open_dataset(REPO_OUT / "sea_ice_mask.nc").ice_mask
    m = (a >= -0.6) & (a <= 0.2) & ~(ice > 0.5)
    return m.values, a.latitude.values.astype(float), a.longitude.values.astype(float), a.values


def time_frame():
    yr = np.repeat(np.arange(2001, 2026), 12); mo = np.tile(np.arange(1, 13), 25)
    return pd.DataFrame({"year": yr, "month": mo, "x": yr + (mo - 0.5) / 12})


def load_sam():
    s = pd.read_csv(REV / "data" / "sam_marshall_monthly.txt", sep=r"\s+", comment="#", header=None, na_values="NA")
    s.columns = ["year"] + list(range(1, 13))
    s = s.melt(id_vars="year", var_name="month", value_name="SAM")
    s["month"] = s.month.astype(int)
    return time_frame().merge(s, on=["year", "month"], how="left")


def pct_dec(slope, mean):
    return 100 * 10 * slope / mean
