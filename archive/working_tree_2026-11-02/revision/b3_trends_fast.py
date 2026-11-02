"""B3 trends per 5-degree bin with the statistic of rstats.trend_bcz, calibrated per bin (300 null series simulated
from the bin's own bias-corrected AR model, both noise models, larger p reported, CI = union of slope +- q95(|z_null|) sd).
Same statistic as trend_bcz and a1_maps_bc; differs from trend_bcz only in the CI (single, not double, bootstrap).
Usage: python b3_trends_fast.py <b3 dir> [vars]"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import a1_maps_bc as M
d = Path(sys.argv[1]); VARS = sys.argv[2:] or ["lat", "sinu", "speed", "eke"]
S = dict(np.load(d / "b3_bin_series.npz"))
x = 2001 + (np.arange(300) + 0.5) / 12; month = np.tile(np.arange(1, 13), 25)
CEN = np.arange(-180, 180, 5.0) + 2.5
NN = 300


def one(args):
    v, k = args; y = S[v][:, k]
    if np.isfinite(y).mean() < 0.8:
        return dict(lon=CEN[k], var=v, cov=float(np.isfinite(y).mean()))
    res = {}
    for order in ("aic", "12"):
        b, z, ps, sig, sd, (phi, m) = M.cell((x, month, y.copy(), order))
        zz = []
        for s in range(NN // 40 + 1):
            zz += M.null_z((x, month, phi, sig, m, 1000 * k + s, order))
        zz = np.abs(np.array(zz[:NN])); zz = zz[np.isfinite(zz)]
        res[order] = dict(b=b, z=z, sd=sd, p=(1 + np.sum(zz >= abs(z))) / (1 + zz.size), q95=np.quantile(zz, 0.95))
    b = res["aic"]["b"]; p = max(res["aic"]["p"], res["12"]["p"])
    lo = min(b - res[o]["q95"] * res[o]["sd"] for o in res); hi = max(b + res[o]["q95"] * res[o]["sd"] for o in res)
    mval = float(np.nanmean(y)); f = 111.195 * 10 if v == "lat" else 1000 / mval
    return dict(lon=CEN[k], var=v, cov=float(np.isfinite(y).mean()), mean=mval, slope=b, tr=b * f, lo=lo * f, hi=hi * f, p=p)


if __name__ == "__main__":
    for v in VARS:
        with ProcessPoolExecutor(max_workers=2) as ex:
            rows = list(ex.map(one, [(v, k) for k in range(72)]))
        D = pd.DataFrame(rows)
        pv = D.p.values; ok = np.isfinite(pv); q = np.zeros(len(D), bool)
        if ok.any():
            ps = np.sort(pv[ok]); m = ok.sum(); thr = ps[ps <= 0.10 * np.arange(1, m + 1) / m]
            if thr.size: q[ok] = pv[ok] <= thr.max()
        D["fdr10"] = q; D.to_csv(d / f"b3_trends_{v}.csv", index=False)
        print(v, "tested", int(ok.sum()), "p<0.05", int((pv[ok] < 0.05).sum()), "FDR", int(q.sum()),
              "median tr", np.nanmedian(D.tr).round(2), flush=True)
