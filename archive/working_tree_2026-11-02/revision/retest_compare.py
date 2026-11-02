"""Summarise retest_all.py: old vs new p for every reported trend; lists results whose 5 % significance changed.
Writes outputs/retest/retest_summary.csv and prints the changes. FDR families are re-evaluated with the new p:
A7 secondary tests (per product and period) and B1 path statistics (18 tests).
Usage: python retest_compare.py
"""
import json
import sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

RT = HERE / "outputs" / "retest"
D = pd.concat([pd.read_csv(f) for f in sorted(RT.glob("retest_*.csv")) if f.stem != "retest_summary"], ignore_index=True)
D["old_sig"] = D.old_p < 0.05
D["new_sig"] = D.p_new < 0.05
D["change"] = np.where(D.old_p.isna(), "new test", np.where(D.old_sig & ~D.new_sig, "LOST", np.where(~D.old_sig & D.new_sig, "GAINED", "")))
D["fdr_new"] = pd.Series([None] * len(D), dtype=object)
meta = D.meta.apply(json.loads)
a7 = D.family == "A7"
for (p, per) in {(m.get("product"), m.get("period")) for m in meta[a7]}:
    sel = a7 & meta.apply(lambda m: m.get("product") == p and m.get("period") == per and m.get("kind") == "secondary")
    D.loc[sel, "fdr_new"] = list(R.bh_fdr(D.loc[sel, "p_new"].values, 0.10))
pth = D.key.str.startswith("paths|")
D.loc[pth, "fdr_new"] = list(R.bh_fdr(D.loc[pth, "p_new"].values, 0.10))
D.to_csv(RT / "retest_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print(D.groupby("family").agg(n=("key", "size"), old_sig=("old_sig", "sum"), new_sig=("new_sig", "sum"),
                              lost=("change", lambda c: (c == "LOST").sum()), gained=("change", lambda c: (c == "GAINED").sum())))
cols = ["family", "key", "old_test", "old_p", "p_new", "z", "z_crit", "pct_per_decade", "phi_sum", "phi_sum_bc"]
print(D[D.change.isin(["LOST", "GAINED"])][cols].round(4).to_string())
