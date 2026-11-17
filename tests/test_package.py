"""Smoke tests: the package imports, its constants are the ones quoted in the paper, paths resolve."""
from pathlib import Path


def test_import_and_version():
    import osr11
    from osr11 import params, paths, rstats, fronts, halfpower, drake_halfpower, style  # noqa: F401
    assert osr11.__version__ == "2.0.0"


def test_params_match_manuscript():
    from osr11 import params as P
    assert (P.HP_REL_THR, P.HP_MIN_W, P.HP_FRAC, P.HP_MAX_ENV, P.HP_MIN_COVERAGE) == (0.20, 100.0, 0.25, 8, 0.90)   # 20 % of the local maximum, 100 km, 25 %, 90 % of months
    assert P.HP_THR == 0.02                                                                      # fixed 200 cm2 s-2 (sensitivity set)
    assert (P.ADT_SOUTH, P.ADT_NORTH) == (-0.6, 0.2)
    assert (P.Y0, P.Y1, P.NMONTHS) == (2001, 2025, 300)
    assert P.PERIOD_EARLY == (2001, 2017) and P.PERIOD_LATE == (2018, 2025)


def test_paths_point_into_repository():
    from osr11 import paths as PA
    root = Path(__file__).resolve().parents[1]
    assert PA.ROOT == root
    assert (PA.OUT / "MANIFEST.csv").exists() or PA.OUT != root / "outputs"   # shipped outputs unless redirected
    assert PA.dk(12, mkdir=False) == PA.OUT / "drake" / "dk12"
