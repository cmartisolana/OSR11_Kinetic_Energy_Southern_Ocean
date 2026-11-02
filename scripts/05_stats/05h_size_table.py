"""Table of the false-positive (size) checks for SI S4, parsed from the size-check logs (written for the publication
repository; the logs are text and no table was produced before).

Purpose : take the last progress line of each log in outputs/stats/size_bc/ and tabulate replicates, rejection rate at 5 %
          (calibrated/bias-corrected test and plain sieve where printed) and CI coverage (annual runs). Logs that ended
          without a progress line (z_ar1_0.95: non-stationary AR fit; z_real_ACC_EKE: no result line) are listed as incomplete.
Inputs  : outputs/stats/size_bc/*.log (05f, 05g)
Outputs : outputs/stats/size_summary.csv
Run     : python scripts/05_stats/05h_size_table.py      (no raw data needed)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import re
import pandas as pd
from osr11 import paths as PA

rows = []
for f in sorted((PA.OUT / "stats" / "size_bc").glob("*.log")):
    lines = f.read_text().splitlines()
    ann = [l for l in lines if re.match(r"^\d+ [\d.]+ \d+ size ", l)]
    mon = [l for l in lines if re.match(r"^\d+ bc [\d.]+ sieve", l)]
    if ann:
        n, phi, k, _, size, _, cov = ann[-1].split()
        rows.append(dict(log=f.name, kind="annual trend_bcz", n_values=int(n), phi=float(phi), replicates=int(k),
                         size_test=float(size), size_sieve=None, coverage=float(cov), complete=True))
    elif mon:
        p = mon[-1].split()
        rows.append(dict(log=f.name, kind="monthly " + ("trend_bcz" if f.name.startswith(("z_", "ord12_z")) else "trend_bc"),
                         n_values=300, phi=None, replicates=int(p[0]), size_test=float(p[2]), size_sieve=float(p[4]),
                         coverage=None, complete=True))
    else:
        rows.append(dict(log=f.name, kind="", replicates=0, complete=False))
T = pd.DataFrame(rows); T.to_csv(PA.out("stats") / "size_summary.csv", index=False); print(T.to_string(index=False))
