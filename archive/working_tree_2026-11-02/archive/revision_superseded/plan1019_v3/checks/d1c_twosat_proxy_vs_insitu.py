"""Product 4 (two-sat) surface transport proxy (GV track, T_perm) vs in situ 0-760 m transport north of 58.55 S.
Reproduces the Product 1 values of d1b_summary.csv (0.535 per transect, 0.670 annual) with the same procedure:
per transect = altimetric month of each transect; annual = in situ annual mean vs altimetric 12-month annual mean.
Usage: python d1c_twosat_proxy_vs_insitu.py out/d1/d1b_insitu_split.csv drake/out/dk4/dk4_section_series.csv"""
import sys, numpy as np, pandas as pd
s = pd.read_csv(sys.argv[1]); a = pd.read_csv(sys.argv[2])
s = s[(s.year >= 2005) & (s.year <= 2019)]
m = s.merge(a, on=["year", "month"]); ins = s.groupby("year").U_north.mean()
for p in ["allsat", "twosat"]:
    r_tr = np.corrcoef(m.U_north, m[f"{p}_Tperm"])[0, 1]
    r_an = np.corrcoef(ins, a.groupby("year")[f"{p}_Tperm"].mean().loc[ins.index])[0, 1]
    print(f"{p}: n={len(m)} r_transect={r_tr:.3f} r_annual={r_an:.3f}")
