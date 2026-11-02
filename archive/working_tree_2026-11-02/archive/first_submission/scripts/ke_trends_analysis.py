# KE / EKE / frontal-position trend analysis from monthly half-power outputs
# Cristina Marti-Solana, 2025-2026
#
# TREND ANALYSIS METHODS
# ──────────────────────────────────────────────────────────────────────
# • Theil–Sen:                        Robust non-parametric median slope estimator
# • Modified Mann–Kendall (Yue–Wang): Significance test (via pymannkendall when available)
#
# Scope:
#   - ACC-wide and basin-level (Atlantic, Indian, Pacific) diagnostics
#
# The script:
#   1) Loads all monthly half-power point NetCDF outputs
#   2) Assigns each front to a basin region by longitude
#   3) Computes monthly time series of KE/EKE and frontal metrics for:
#        - ACC (all longitudes)
#        - Atlantic basin
#        - Indian basin
#        - Pacific basin
#   4) Fits linear trends (Theil–Sen + Modified Mann–Kendall)
#   5) Saves CSV outputs for plotting and downstream analyses
#
# Usage:
#   cd scripts/
#   python ke_trends_analysis.py

import sys
import warnings
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path

# Allow running both as:
#   - python scripts/ke_trends_analysis.py   (repo root)
#   - cd scripts && python ke_trends_analysis.py
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from utils.trends import theil_sen_slope, mann_kendall

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ──────────────────────────────────────────────────────────────────────
# 1. CONFIGURATION
# ──────────────────────────────────────────────────────────────────────

# Path to the monthly half-power point outputs
_REPO_ROOT = _SCRIPT_DIR.parent
INPUT_DIR = _REPO_ROOT / "outputs" / "monthly_half-power_points"

# Output directory for trend results
OUTPUT_DIR = _REPO_ROOT / "outputs" / "trends"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Southern Ocean basin partitions (longitude in degrees east, wrapped to [0, 360)).
# Atlantic spans two ranges, so it is handled explicitly in assign_region().
BASIN_REGIONS = {
    "Indian": (20.0, 147.0),
    "Pacific": (147.0, 290.0),
}

REGION_ORDER = ["ACC", "Atlantic", "Indian", "Pacific"]

# Latitude range considered for the ACC
LAT_MIN, LAT_MAX = -65.0, -35.0


# ──────────────────────────────────────────────────────────────────────
# 2. HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────

def lon_in_range(lon, lo, hi):
    """Check if longitude falls inside [lo, hi) on a wrapped [0, 360) axis."""
    # Normalise to [0, 360)
    lon360 = lon % 360
    lo360 = lo % 360
    hi360 = hi % 360
    if lo360 < hi360:
        return lo360 <= lon360 < hi360
    else:  # wraps around 0/360
        return lon360 >= lo360 or lon360 < hi360


def load_all_monthly_outputs(input_dir: Path):
    """Load every monthly half-power point NetCDF into a single DataFrame.

    Returns
    -------
    df : pd.DataFrame
        Columns: year, month, time, front_lon, front_lat,
                 envelope_lat_south, envelope_lat_north,
                 peak_lat, peak_val, width_km
    eke_fields : dict
        Mapping (year, month) -> (eke_2d, lat_1d, lon_1d)
    ke_fields : dict
        Mapping (year, month) -> (ke_2d, lat_1d, lon_1d)
    coords : dict
        'latitude' and 'longitude' arrays from the grids
    """
    records = []
    eke_fields = {}
    ke_fields = {}
    coords = {}

    years = sorted([y for y in input_dir.iterdir()
                    if y.is_dir() and y.name != ".DS_Store"])

    for year_path in years:
        months = sorted([m for m in year_path.iterdir()
                         if m.is_file() and m.suffix == ".nc"])
        for month_path in months:
            ds = xr.open_dataset(month_path)
            yr = int(year_path.name)
            mo = int(month_path.stem)
            t = pd.Timestamp(yr, mo, 15)  # mid-month representative time

            n = len(ds.front)
            if n == 0:
                continue

            for i in range(n):
                records.append({
                    "year": yr,
                    "month": mo,
                    "time": t,
                    "front_lon": float(ds.front_lon.values[i]),
                    "front_lat": float(ds.front_lat.values[i]),
                    "envelope_lat_south": float(ds.envelope_lat_south.values[i]),
                    "envelope_lat_north": float(ds.envelope_lat_north.values[i]),
                    "peak_lat": float(ds.peak_lat.values[i]),
                    "peak_val": float(ds.peak_val.values[i]),
                    "width_km": float(ds.width_km.values[i]),
                })

            # Store 2-D energy fields alongside their own coordinates
            lat = ds.latitude.values.copy()
            lon = ds.longitude.values.copy()
            eke_fields[(yr, mo)] = (ds["eke"].values.copy(), lat, lon)
            ke_fields[(yr, mo)]  = (ds["ke"].values.copy(),  lat, lon)
            if not coords:
                coords["latitude"] = lat
                coords["longitude"] = lon
            ds.close()

    df = pd.DataFrame(records)
    if not df.empty:
        df["time"] = pd.to_datetime(df["time"])
    return df, eke_fields, ke_fields, coords


def assign_region(lon):
    """Return basin region for a given longitude.

    Regions
    -------
    Atlantic : [290, 360) U [0, 20)
    Indian   : [20, 147)
    Pacific  : [147, 290)
    """
    lon360 = lon % 360
    if lon360 >= 290.0 or lon360 < 20.0:
        return "Atlantic"
    if lon_in_range(lon360, *BASIN_REGIONS["Indian"]):
        return "Indian"
    if lon_in_range(lon360, *BASIN_REGIONS["Pacific"]):
        return "Pacific"
    return "Unknown"


def _mk_sigstars(pv: float) -> str:
    if not np.isfinite(pv):
        return ""
    if pv < 0.001:
        return "***"
    if pv < 0.01:
        return "**"
    if pv < 0.05:
        return "*"
    return ""


# ──────────────────────────────────────────────────────────────────────
# 3. COMPUTE MONTHLY REGION STATISTICS
# ──────────────────────────────────────────────────────────────────────

def compute_region_timeseries(df):
    """Aggregate front detections into monthly region-level statistics.

    Returns
    -------
    ts : pd.DataFrame
        Indexed by (region, year, month) with columns:
        n_fronts, mean_lat, median_lat, std_lat,
        mean_peak_val, median_peak_val,
        mean_width_km, mean_envelope_span,
        decimal_year
    """
    df = df.copy()
    df["region"] = df["front_lon"].apply(assign_region)
    df = df[df["region"].isin(["Atlantic", "Indian", "Pacific"])].copy()
    df["envelope_span"] = df["envelope_lat_north"] - df["envelope_lat_south"]

    groups = df.groupby(["region", "year", "month"])

    rows = []
    for (region, yr, mo), grp in groups:
        rows.append({
            "region": region,
            "year": yr,
            "month": mo,
            "decimal_year": yr + (mo - 0.5) / 12.0,
            "n_fronts": len(grp),
            "mean_lat": grp["front_lat"].mean(),
            "median_lat": grp["front_lat"].median(),
            "std_lat": grp["front_lat"].std(),
            "mean_peak_val": grp["peak_val"].mean(),
            "median_peak_val": grp["peak_val"].median(),
            "mean_width_km": grp["width_km"].mean(),
            "mean_envelope_span": grp["envelope_span"].mean(),
        })

    ts = pd.DataFrame(rows)
    if ts.empty:
        return ts
    ts["region"] = pd.Categorical(ts["region"], categories=REGION_ORDER[1:], ordered=True)
    return ts.sort_values(["region", "year", "month"]).reset_index(drop=True)


def compute_whole_SO_timeseries(df):
    """Aggregate all fronts per month regardless of basin."""
    df = df.copy()
    df["envelope_span"] = df["envelope_lat_north"] - df["envelope_lat_south"]
    groups = df.groupby(["year", "month"])

    rows = []
    for (yr, mo), grp in groups:
        rows.append({
            "region": "ACC",
            "year": yr,
            "month": mo,
            "decimal_year": yr + (mo - 0.5) / 12.0,
            "n_fronts": len(grp),
            "mean_lat": grp["front_lat"].mean(),
            "median_lat": grp["front_lat"].median(),
            "std_lat": grp["front_lat"].std(),
            "mean_peak_val": grp["peak_val"].mean(),
            "median_peak_val": grp["peak_val"].median(),
            "mean_width_km": grp["width_km"].mean(),
            "mean_envelope_span": grp["envelope_span"].mean(),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["year", "month"]).reset_index(drop=True)


def compute_energy_field_timeseries(fields, field_name):
    """Compute ACC and basin means of a 2-D energy field.

    Parameters
    ----------
    fields : dict
        Mapping (yr, mo) -> (field_2d, lat_1d, lon_1d).
    field_name : str
        Label used for the output columns, e.g. 'eke' or 'ke'.
    Returns
    -------
    pd.DataFrame with columns: region, year, month, decimal_year,
        mean_{field_name}, max_{field_name}.
    """
    col_mean = f"mean_{field_name}"
    col_max  = f"max_{field_name}"

    rows = []
    for (yr, mo), (data2d, lat, lon) in sorted(fields.items()):
        cos_lat = np.cos(np.radians(lat))
        valid = np.isfinite(data2d)
        weights = np.broadcast_to(cos_lat[:, None], data2d.shape)
        w = np.where(valid, weights, 0.0)
        mean_val = np.nansum(data2d * w) / np.nansum(w) if np.nansum(w) > 0 else np.nan
        max_val  = np.nanmax(data2d) if valid.any() else np.nan

        rows.append({
            "region": "ACC", "year": yr, "month": mo,
            "decimal_year": yr + (mo - 0.5) / 12.0,
            col_mean: mean_val, col_max: max_val,
        })

        region_by_lon = np.array([assign_region(l) for l in lon], dtype=object)
        for rname in ["Atlantic", "Indian", "Pacific"]:
            mask_lon = region_by_lon == rname
            region_data = data2d[:, mask_lon]
            region_weights = np.broadcast_to(cos_lat[:, None], region_data.shape)
            valid_r = np.isfinite(region_data)
            wr = np.where(valid_r, region_weights, 0.0)
            rmean = np.nansum(region_data * wr) / np.nansum(wr) if np.nansum(wr) > 0 else np.nan
            rmax  = np.nanmax(region_data) if valid_r.any() else np.nan

            rows.append({
                "region": rname, "year": yr, "month": mo,
                "decimal_year": yr + (mo - 0.5) / 12.0,
                col_mean: rmean, col_max: rmax,
            })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["region"] = pd.Categorical(out["region"], categories=REGION_ORDER, ordered=True)
    return out.sort_values(["region", "year", "month"]).reset_index(drop=True)


# ──────────────────────────────────────────────────────────────────────
# 4. TREND FITTING
# ──────────────────────────────────────────────────────────────────────

def fit_all_trends(ts, variable, groupby="region"):
    """Fit Theil–Sen + Mann–Kendall trends per group.

    Parameters
    ----------
    ts : pd.DataFrame
        Must contain 'decimal_year', `variable`, and `groupby` columns.
    variable : str
        Column name to analyse (e.g. 'mean_peak_val', 'median_lat').
    groupby : str
        Grouping column (default 'region').

    Returns
    -------
    trends : pd.DataFrame
        One row per group with Theil–Sen slope (+ 95% CI) and MK p-value.
    """
    results = []
    for name, grp in ts.groupby(groupby):
        x = grp["decimal_year"].values
        y = grp[variable].values

        ts_sl, ts_ic, ts_lo, ts_hi = theil_sen_slope(x, y)
        mk_tau, mk_p = mann_kendall(y)

        results.append({
            groupby: name,
            "variable": variable,
            "n_points": int(np.isfinite(y).sum()),
            # Theil-Sen
            "ts_slope": ts_sl,
            "ts_intercept": ts_ic,
            "ts_slope_lo95": ts_lo,
            "ts_slope_hi95": ts_hi,
            # Mann-Kendall (Modified with autocorrelation correction)
            "mk_tau": mk_tau,
            "mk_pvalue": mk_p,
        })

    return pd.DataFrame(results)


# ──────────────────────────────────────────────────────────────────────
# 5. SAVE RESULTS
# ──────────────────────────────────────────────────────────────────────

def save_results(ts, energy_ts, trend_tables, output_dir):
    """Save analysis results to CSV."""
    ts.to_csv(output_dir / "region_timeseries.csv", index=False)
    energy_ts.to_csv(output_dir / "energy_field_timeseries.csv", index=False)

    # Trend tables (includes mean_eke_gridded and mean_ke_gridded)
    for name, tbl in trend_tables.items():
        tbl.to_csv(output_dir / f"trends_{name}.csv", index=False)

    print(f"  All results saved to {output_dir}/")


# ──────────────────────────────────────────────────────────────────────
# 6. MAIN
# ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("KE/EKE/FRONTAL TRENDS — ACC and basin-level")
    print("=" * 70)

    # ── Load data ────────────────────────────────────────────────────
    print("\n[1/4] Loading monthly half-power point outputs...")
    df, eke_fields, ke_fields, coords = load_all_monthly_outputs(INPUT_DIR)

    if df.empty:
        print("ERROR: No data found. Run half_power_EKE.py first.")
        sys.exit(1)

    years = sorted(df["year"].unique())
    n_months = df.groupby(["year", "month"]).ngroups
    print(f"  Loaded {len(df):,} frontal detections over {n_months} months "
          f"({years[0]}–{years[-1]})")

    # ── Region time series (ACC + Atlantic/Indian/Pacific) ─────────
    print("\n[2/4] Computing region-level monthly time series (ACC + basins)...")
    ts_region = compute_region_timeseries(df)
    ts_whole = compute_whole_SO_timeseries(df)
    ts = pd.concat([ts_whole, ts_region], ignore_index=True)
    print(f"  Regions: {', '.join(ts['region'].astype(str).unique())}")

    # ── Gridded KE & EKE time series ─────────────────────────────────
    print("\n[3/4] Computing area-weighted mean KE and EKE from gridded fields...")
    eke_ts = compute_energy_field_timeseries(eke_fields, "eke")
    ke_ts  = compute_energy_field_timeseries(ke_fields,  "ke")
    energy_ts = eke_ts.merge(ke_ts, on=["region", "year", "month", "decimal_year"], how="outer")
    print(f"  {len(energy_ts)} region×month records")

    # ── Trend fitting ────────────────────────────────────────────────
    print("\n[4/4] Fitting linear trends...")
    variables_to_trend = [
        "median_lat",        # frontal position
        "mean_peak_val",     # peak EKE within envelopes
        "mean_width_km",     # envelope width
        "mean_envelope_span",# lat span of envelopes
        "n_fronts",          # number of detected fronts
    ]

    trend_tables = {}
    for var in variables_to_trend:
        tdf = fit_all_trends(ts, var)
        trend_tables[var] = tdf
        print(f"\n  ── {var} ──")
        for _, row in tdf.iterrows():
            sig = _mk_sigstars(float(row["mk_pvalue"]))
            print(
                f"    {str(row['region']):20s}  "
                f"Theil-Sen slope = {row['ts_slope']:+.4e}/yr  "
                f"MK p = {row['mk_pvalue']:.3f} {sig}  "
                f"MK τ = {row['mk_tau']:+.3f}"
            )

    # Also fit trends on gridded EKE and KE
    for energy_var, label in [("mean_eke", "EKE"), ("mean_ke", "KE")]:
        energy_trend = fit_all_trends(energy_ts, energy_var)
        trend_tables[f"{energy_var}_gridded"] = energy_trend
        print(f"\n  ── {energy_var} (gridded field, area-weighted) — {label} ──")
        for _, row in energy_trend.iterrows():
            sig = _mk_sigstars(float(row["mk_pvalue"]))
            print(
                f"    {str(row['region']):20s}  "
                f"Theil-Sen slope = {row['ts_slope']:+.4e}/yr  "
                f"MK p = {row['mk_pvalue']:.3f} {sig}"
            )

    # ── Save ─────────────────────────────────────────────────────────
    print("\nSaving results...")
    save_results(ts, energy_ts, trend_tables, OUTPUT_DIR)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
