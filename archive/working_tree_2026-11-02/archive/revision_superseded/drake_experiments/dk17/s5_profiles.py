"""4. Streamline-averaged (and latitude-averaged) profiles of mean EKE/KE and their 2001-2025 trends.
Input: out/s2_binned.npz. Trend per bin: the paper's test on the monthly series (rstats.trend_bcz, deseasonalised
Theil-Sen, calibrated AR noise, n_outer=300). Bins with < 0.3 % of the domain area are dropped. Circumpolar profiles in
fixed (mean-ADT) and moving (monthly ADT minus the domain-mean anomaly) streamline bins; meander regions (s2).
"""
import numpy as np, pandas as pd
from multiprocessing import Pool
from common import OUT, R

Z = np.load(OUT / "s2_binned.npz")
xa = 2001 + (np.arange(300) + 0.5) / 12; MON = np.tile(np.arange(1, 13), 25)
SEC = ["ACC", "Atlantic", "Indian", "Pacific"]
jobs = []
for b, secs, vars_ in (("madt", (0,), ("EKE", "KE", "SPD")), ("iadt", (0,), ("EKE", "KE", "SPD"))):
    edges = Z["madt_edges"] if b != "lat" else Z["lat_edges"]
    area = Z[f"area_{b}"]
    tot = area[:, 0].sum(1).mean()
    for s in secs:
        for k in range(edges.size - 1):
            am = area[:, s, k].mean()
            if am < 0.003 * tot:
                continue
            for v in vars_:
                y = Z[f"{v}_{b}"][:, s, k]
                jobs.append(dict(bin=b, sector=SEC[s], k=k, lo=edges[k], hi=edges[k + 1], var=v, area_frac=am / tot, y=y))
for i, nm in enumerate(Z["regnames"]):
    for v in ("EKE", "KE", "SPD"):
        jobs.append(dict(bin="region", sector=str(nm), k=i, lo=np.nan, hi=np.nan, var=v, area_frac=float(Z["reg_area"][i] / Z["reg_area"].sum()), y=Z[f"reg_{v}"][:, i]))


def run(j):
    y = j.pop("y")
    o = R.trend_bcz(xa, y, month=MON, n_outer=300)
    m = y.mean()
    return dict(**j, mean=m, slope=o["slope"], ci_lo=o["ci_lo"], ci_hi=o["ci_hi"], p=o["p"],
                pct_dec=1000 * o["slope"] / m, pct_dec_lo=1000 * o["ci_lo"] / m, pct_dec_hi=1000 * o["ci_hi"] / m)


if __name__ == "__main__":
    print(len(jobs), "trend calls", flush=True)
    with Pool(2) as p:
        rows = p.map(run, jobs, chunksize=4)
    P = pd.DataFrame(rows); P["mid"] = 0.5 * (P.lo + P.hi); P.to_csv(OUT / "s5_profiles.csv", index=False)
    pd.set_option("display.width", 220)
    print(P[["bin", "sector", "var", "mid", "area_frac", "mean", "pct_dec", "pct_dec_lo", "pct_dec_hi", "p"]].round(3).to_string())
