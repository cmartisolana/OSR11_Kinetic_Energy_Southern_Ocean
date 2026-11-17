"""pytest configuration: make the osr11 package importable without installation and keep matplotlib headless."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("PYTHONPYCACHEPREFIX", str(Path(os.environ.get("TMPDIR", "/tmp")) / "osr11_pyc"))
