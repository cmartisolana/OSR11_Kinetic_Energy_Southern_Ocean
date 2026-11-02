#!/usr/bin/env python3
"""Extract observed Drake Passage transport from TrendsDrakePassageTransp.

Source repository:
https://github.com/manuelogtzv/TrendsDrakePassageTransp

By default, reads `Datasets/TotalTransport_os38nbmaxz_760_2.mat` and exports:
- date (UTC-like calendar date from MATLAB datenum)
- total_transport_m3s
- total_transport_sv

Optional mode (`--source geostrophic`) preserves legacy extraction from
`Datasets/geostrophic_noaver.mat` with:
- net_transport_m3s
- net_transport_sv
"""

from __future__ import annotations

import argparse
import subprocess
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import scipy.io as sio

SOURCE_REPO_URL = "https://github.com/manuelogtzv/TrendsDrakePassageTransp"
DEFAULT_CLONE_DIR = Path("external/TrendsDrakePassageTransp")
TOTAL_MAT_REL_PATH = Path("Datasets/TotalTransport_os38nbmaxz_760_2.mat")
GEOS_MAT_REL_PATH = Path("Datasets/geostrophic_noaver.mat")
DEFAULT_OUTPUT_TOTAL = Path("outputs/trends/drake_passage_total_transport_observations.csv")
DEFAULT_OUTPUT_GEOS = Path("outputs/trends/drake_passage_geostrophic_transport_observations.csv")


def matlab_datenum_to_datetime(datenum: float) -> datetime:
    """Convert MATLAB datenum to Python datetime."""
    day = datetime.fromordinal(int(datenum))
    frac = timedelta(days=float(datenum) % 1)
    return day + frac - timedelta(days=366)


def ensure_repo(local_repo: Path) -> None:
    if local_repo.exists():
        return
    local_repo.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["git", "clone", "--depth", "1", SOURCE_REPO_URL, str(local_repo)],
        check=True,
    )


def build_timeseries_total(mat_path: Path) -> pd.DataFrame:
    """Extract observed TOTAL transport from totaltransp.nettransp."""
    mat = sio.loadmat(str(mat_path), squeeze_me=True, struct_as_record=False)
    totaltransp = mat["totaltransp"]

    time = pd.Series(totaltransp.time).astype(float)
    total_m3s = pd.Series(totaltransp.nettransp).astype(float)

    df = pd.DataFrame(
        {
            "date": time.map(matlab_datenum_to_datetime),
            "total_transport_m3s": total_m3s,
        }
    )
    df["total_transport_sv"] = df["total_transport_m3s"] / 1e6
    return df.sort_values("date").reset_index(drop=True)


def build_timeseries_geostrophic(mat_path: Path) -> pd.DataFrame:
    """Extract observed geostrophic transport from geosbaroc.nettransp."""
    mat = sio.loadmat(str(mat_path), squeeze_me=True, struct_as_record=False)
    geosbaroc = mat["geosbaroc"]

    time = pd.Series(geosbaroc.time).astype(float)
    net_m3s = pd.Series(geosbaroc.nettransp).astype(float)

    df = pd.DataFrame(
        {
            "date": time.map(matlab_datenum_to_datetime),
            "net_transport_m3s": net_m3s,
        }
    )
    df["net_transport_sv"] = df["net_transport_m3s"] / 1e6
    return df.sort_values("date").reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        choices=["total", "geostrophic"],
        default="total",
        help="Which transport product to extract.",
    )
    parser.add_argument(
        "--repo-path",
        type=Path,
        default=DEFAULT_CLONE_DIR,
        help="Local path to TrendsDrakePassageTransp repository.",
    )
    parser.add_argument(
        "--mat-path",
        type=Path,
        default=None,
        help="Path to MAT file (overrides --repo-path and --source defaults).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV path.",
    )
    parser.add_argument(
        "--no-clone",
        action="store_true",
        help="Do not clone if --repo-path is missing.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    default_mat = TOTAL_MAT_REL_PATH if args.source == "total" else GEOS_MAT_REL_PATH
    default_output = DEFAULT_OUTPUT_TOTAL if args.source == "total" else DEFAULT_OUTPUT_GEOS

    if args.mat_path is not None:
        mat_path = args.mat_path
    else:
        if not args.no_clone:
            ensure_repo(args.repo_path)
        mat_path = args.repo_path / default_mat

    if not mat_path.exists():
        raise FileNotFoundError(f"MAT file not found: {mat_path}")

    if args.source == "total":
        df = build_timeseries_total(mat_path)
        mean_col = "total_transport_sv"
        summary_label = "Mean total observed transport"
    else:
        df = build_timeseries_geostrophic(mat_path)
        mean_col = "net_transport_sv"
        summary_label = "Mean net geostrophic transport"

    output_path = args.output if args.output is not None else default_output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    print(f"Wrote {len(df)} rows to {output_path}")
    print(f"Date range: {df['date'].min().date()} to {df['date'].max().date()}")
    print(f"{summary_label}: {df[mean_col].mean():.2f} Sv")


if __name__ == "__main__":
    main()
