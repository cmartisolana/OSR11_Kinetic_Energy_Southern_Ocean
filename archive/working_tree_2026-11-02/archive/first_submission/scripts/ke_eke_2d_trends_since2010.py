"""Compute pixel-wise KE/EKE 2D trends from monthly half-power outputs (since 2010).

This is a lightweight variant of `ke_eke_2d_trends.py` that restricts
the analysis to monthly samples from 2010.0 onwards and writes a
separate NetCDF file so figures can be generated from this subset.

Usage
-----
python scripts/ke_eke_2d_trends_since2010.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr
from scipy import stats


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
INPUT_HP_DIR = REPO_ROOT / "outputs" / "monthly_half-power_points"
OUTPUT_DIR = REPO_ROOT / "outputs" / "trends"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "ke_eke_2d_trends_since2010.nc"


def theil_sen_with_ci(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return {"slope": np.nan, "intercept": np.nan, "lo": np.nan, "hi": np.nan}

    res = stats.theilslopes(y[m], x[m], alpha=0.95)
    return {
        "slope": float(res.slope),
        "intercept": float(res.intercept),
        "lo": float(res.low_slope),
        "hi": float(res.high_slope),
    }


def _mk_s_var(y: np.ndarray) -> tuple[float, float]:
    n = len(y)
    s = 0.0
    for i in range(n - 1):
        s += np.sign(y[i + 1 :] - y[i]).sum()

    _, counts = np.unique(y, return_counts=True)
    tie_term = np.sum(counts * (counts - 1) * (2 * counts + 5))
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0
    return s, var_s


def mann_kendall_yue_wang(y: np.ndarray, alpha: float = 0.05, max_lag: int = 12) -> dict[str, float]:
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    n = len(y)
    if n < 5:
        return {"tau": np.nan, "p": np.nan, "n_eff": np.nan}

    s, var_s = _mk_s_var(y)
    tau = s / (0.5 * n * (n - 1))

    x = np.arange(n, dtype=float)
    ts = theil_sen_with_ci(x, y)
    if np.isfinite(ts["slope"]) and np.isfinite(ts["intercept"]):
        y_for_ac = y - (ts["slope"] * x + ts["intercept"])
    else:
        y_for_ac = y - np.nanmean(y)

    ranks = stats.rankdata(y_for_ac)
    r = ranks - np.mean(ranks)
    denom = np.sum(r**2)
    if denom <= 0:
        return {"tau": float(tau), "p": np.nan, "n_eff": np.nan}

    max_lag = min(max_lag, n - 2)
    sig_rhos: list[tuple[int, float]] = []
    for k in range(1, max_lag + 1):
        num = np.sum(r[:-k] * r[k:])
        rho_k = num / denom
        crit = stats.norm.ppf(1 - alpha / 2) / np.sqrt(n)
        if np.abs(rho_k) > crit:
            sig_rhos.append((k, rho_k))

    if len(sig_rhos) == 0:
        n_eff = float(n)
    else:
        corr_sum = np.sum([(1 - k / n) * rho for k, rho in sig_rhos])
        n_eff = float(np.clip(n / (1 + 2 * corr_sum), 2, n))

    var_s_adj = var_s * (n / n_eff)
    if s > 0:
        z = (s - 1) / np.sqrt(var_s_adj)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s_adj)
    else:
        z = 0.0

    p = 2 * stats.norm.sf(abs(z))
    return {"tau": float(tau), "p": float(p), "n_eff": n_eff}


def load_monthly_ke_eke(input_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    files = sorted(input_dir.glob("*/*.nc"))
    if not files:
        raise FileNotFoundError(f"No monthly NetCDF files found in {input_dir}")

    ke_stack: list[np.ndarray] = []
    eke_stack: list[np.ndarray] = []
    t_stack: list[float] = []
    lat_1d: np.ndarray | None = None
    lon_1d: np.ndarray | None = None

    for fp in files:
        ds = xr.open_dataset(fp)

        if ("ke" not in ds) or ("eke" not in ds):
            ds.close()
            continue

        yr = int(fp.parent.name)
        mo = int(fp.stem)
        t = yr + (mo - 0.5) / 12.0

        ke_stack.append(ds["ke"].values.astype(float))
        eke_stack.append(ds["eke"].values.astype(float))
        t_stack.append(t)

        if lat_1d is None:
            lat_1d = ds["latitude"].values.astype(float).copy()
            lon_1d = ds["longitude"].values.astype(float).copy()

        ds.close()

    if not ke_stack or lat_1d is None or lon_1d is None:
        raise RuntimeError("No valid KE/EKE fields were found in monthly files.")

    ke_3d = np.asarray(ke_stack, dtype=float)
    eke_3d = np.asarray(eke_stack, dtype=float)
    time_decimal = np.asarray(t_stack, dtype=float)

    order = np.argsort(time_decimal)
    time_decimal = time_decimal[order]
    ke_3d = ke_3d[order, :, :]
    eke_3d = eke_3d[order, :, :]

    return time_decimal, lat_1d, lon_1d, ke_3d, eke_3d


def compute_pixelwise_trends(
    time_decimal: np.ndarray,
    ke_3d: np.ndarray,
    eke_3d: np.ndarray,
    min_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    _, nlat, nlon = ke_3d.shape
    ke_slope = np.full((nlat, nlon), np.nan, dtype=float)
    eke_slope = np.full((nlat, nlon), np.nan, dtype=float)
    ke_p = np.full((nlat, nlon), np.nan, dtype=float)
    eke_p = np.full((nlat, nlon), np.nan, dtype=float)

    ke_counts = np.sum(np.isfinite(ke_3d), axis=0)
    eke_counts = np.sum(np.isfinite(eke_3d), axis=0)

    for j in range(nlat):
        for i in range(nlon):
            if ke_counts[j, i] >= min_points:
                y_ke = ke_3d[:, j, i]
                m_ke = np.isfinite(y_ke)
                ts_ke = theil_sen_with_ci(time_decimal[m_ke], y_ke[m_ke])
                mk_ke = mann_kendall_yue_wang(y_ke[m_ke])
                ke_slope[j, i] = ts_ke["slope"] * 1e4
                ke_p[j, i] = mk_ke["p"]

            if eke_counts[j, i] >= min_points:
                y_eke = eke_3d[:, j, i]
                m_eke = np.isfinite(y_eke)
                ts_eke = theil_sen_with_ci(time_decimal[m_eke], y_eke[m_eke])
                mk_eke = mann_kendall_yue_wang(y_eke[m_eke])
                eke_slope[j, i] = ts_eke["slope"] * 1e4
                eke_p[j, i] = mk_eke["p"]

    return ke_slope, eke_slope, ke_p, eke_p


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute 2D KE/EKE trends (since 2010).")
    parser.add_argument("--stride", type=int, default=1, help="Spatial decimation stride for trend calculation.")
    parser.add_argument("--min-points", type=int, default=24, help="Minimum finite monthly points per cell.")
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_FILE,
        help="Output NetCDF path (default: outputs/trends/ke_eke_2d_trends_since2010.nc).",
    )
    args = parser.parse_args()

    if args.stride < 1:
        raise ValueError("--stride must be >= 1")
    if args.min_points < 3:
        raise ValueError("--min-points must be >= 3")

    print("=" * 70)
    print("COMPUTE 2D KE/EKE TRENDS (SINCE 2010)")
    print("=" * 70)
    print(f"Input monthly fields: {INPUT_HP_DIR}")
    print(f"Output file: {args.output}")

    time_decimal, lat_1d, lon_1d, ke_3d, eke_3d = load_monthly_ke_eke(INPUT_HP_DIR)

    # Restrict to >= 2010.0
    mask_time = time_decimal >= 2010.0
    if mask_time.sum() == 0:
        raise RuntimeError("No monthly samples >= 2010 found in input directory.")

    time_decimal = time_decimal[mask_time]
    ke_3d = ke_3d[mask_time, :, :]
    eke_3d = eke_3d[mask_time, :, :]

    lat_sel = lat_1d[:: args.stride]
    lon_sel = lon_1d[:: args.stride]
    ke_sel = ke_3d[:, :: args.stride, :: args.stride]
    eke_sel = eke_3d[:, :: args.stride, :: args.stride]

    ke_mean = np.nanmean(ke_sel, axis=0) * 1e4
    eke_mean = np.nanmean(eke_sel, axis=0) * 1e4

    print(f"Monthly samples: {ke_sel.shape[0]}")
    print(f"Grid shape used for trends: {ke_sel.shape[1]} x {ke_sel.shape[2]}")
    print(f"Time span: {time_decimal.min():.2f} to {time_decimal.max():.2f}")

    ke_slope, eke_slope, ke_p, eke_p = compute_pixelwise_trends(
        time_decimal=time_decimal,
        ke_3d=ke_sel,
        eke_3d=eke_sel,
        min_points=args.min_points,
    )

    ke_slope_sig = np.where((ke_p < 0.05) & np.isfinite(ke_p), ke_slope, np.nan)
    eke_slope_sig = np.where((eke_p < 0.05) & np.isfinite(eke_p), eke_slope, np.nan)

    ds_out = xr.Dataset(
        data_vars={
            "ke_mean_cm2_s2": (("latitude", "longitude"), ke_mean),
            "eke_mean_cm2_s2": (("latitude", "longitude"), eke_mean),
            "ke_slope_cm2_s2_yr": (("latitude", "longitude"), ke_slope),
            "eke_slope_cm2_s2_yr": (("latitude", "longitude"), eke_slope),
            "ke_pvalue": (("latitude", "longitude"), ke_p),
            "eke_pvalue": (("latitude", "longitude"), eke_p),
            "ke_slope_sig_cm2_s2_yr": (("latitude", "longitude"), ke_slope_sig),
            "eke_slope_sig_cm2_s2_yr": (("latitude", "longitude"), eke_slope_sig),
        },
        coords={
            "latitude": lat_sel,
            "longitude": lon_sel,
        },
        attrs={
            "description": "Pixel-wise KE/EKE means and Theil-Sen trends (since 2010)",
            "units_mean": "cm2 s-2",
            "units_trend": "cm2 s-2 yr-1",
            "significance_threshold": "p < 0.05",
            "trend_method": "Theil-Sen",
            "significance_method": "Mann-Kendall Yue-Wang",
            "min_points": int(args.min_points),
            "stride": int(args.stride),
            "time_start_decimal_year": float(time_decimal.min()),
            "time_end_decimal_year": float(time_decimal.max()),
            "n_time": int(ke_sel.shape[0]),
        },
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    ds_out.to_netcdf(args.output)

    sig_ke = np.isfinite(ke_slope_sig).sum()
    sig_eke = np.isfinite(eke_slope_sig).sum()
    print(f"Saved: {args.output}")
    print(f"Significant KE cells: {sig_ke}")
    print(f"Significant EKE cells: {sig_eke}")


if __name__ == "__main__":
    main()
