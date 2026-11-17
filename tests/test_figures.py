"""The figure scripts run end to end from the shipped outputs and write a PNG and a PDF each."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MAIN = ["fig1_halfpower_method", "fig2_drake", "fig3_circumpolar", "fig4_wind_slope"]
SI = ["figS1_drake_section", "figS2_timing", "figS3_products", "figS4_previous_studies", "figS5_pathway_robustness",
      "figS6_occurrence_adt", "figS7_position_longitude"]


@pytest.mark.parametrize("name", MAIN + SI)
def test_figure_runs(name, tmp_path):
    if not (ROOT / "outputs" / "MANIFEST.csv").exists():
        pytest.skip("shipped outputs not present")
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONPYCACHEPREFIX=str(tmp_path / "pyc"))
    args = [str(ROOT / "outputs"), str(tmp_path)] if name.startswith("figS") else [str(tmp_path)]
    r = subprocess.run([sys.executable, str(ROOT / "figures" / f"{name}.py"), *args], capture_output=True, text=True, env=env, timeout=900)
    assert r.returncode == 0, r.stderr[-2000:]
    pngs = list((tmp_path / "png").glob("Figure*.png")); pdfs = list((tmp_path / "pdf").glob("Figure*.pdf"))
    assert len(pngs) == 1 and pngs[0].stat().st_size > 50_000 and len(pdfs) == 1


def test_manifest_matches_outputs():
    """every file listed in outputs/MANIFEST.csv exists with the listed size"""
    import pandas as pd
    M = pd.read_csv(ROOT / "outputs" / "MANIFEST.csv")
    bad = [r.path for r in M.itertuples() if not (ROOT / "outputs" / r.path).exists() or (ROOT / "outputs" / r.path).stat().st_size != r.bytes]
    assert not bad, bad[:10]
