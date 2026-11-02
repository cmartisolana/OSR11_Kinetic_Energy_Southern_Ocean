"""B1 (revised design): TNG14 covariability test at standing meanders.

Specification fixed on 2026-10-08, after the trend pilot and BEFORE any of the
cross-correlations below were computed (see R2_revision_results.md, B1).
Rationale: in Thompson & Naveira Garabato (2014) the curvature-EKE link is a
month-to-month covariation (2-month smoothing, r = 0.34-0.46), not a trend; the
pilot's 12-month low-pass removed that time scale.

Series (monthly, per meander; from b1_pilot.py / b1_meander_run.py):
  curvature : core-level contour length (primary, TNG14's measure) and
              integral |kappa| ds (secondary), both at 150 km along-arc smoothing
  up_core_speed, mea_eke, down_eke
Preparation: remove the monthly climatology, remove a Theil-Sen trend, 2-month running mean.
Links (x leads y by L = 0..3 months):
  L1 up_core_speed -> curvature
  L2 curvature -> down_eke    (primary: spatially separate from the contour)
  L3 curvature -> mea_eke     (secondary: an eddy that bends the contour also raises local EKE)
Significance: (i) t-test with n_eff = n (1 - r1x r1y) / (1 + r1x r1y);
  (ii) phase-randomised surrogates of x (2000), per lag and for the maximum |r| over the
  4 lags (multiplicity across lags).
Pooled (several meanders): Stouffer z of Fisher-z values weighted by sqrt(n_eff - 3), per lag;
  surrogate p for the max over lags with independent surrogates per meander; sign consistency.
Usage: python b1_covar.py <out_dir> name1=series1.csv [name2=series2.csv ...]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

LAGS = range(4)
CURV = {"length": "core_s150_length_km", "kappa": "core_s150_int_abs_kappa"}
NSURR = 2000


def prep(P, col):
    y = R.deseasonalise(P[col].values.astype(float), P.month.values)
    ok = np.isfinite(y)
    b, a = R.theil_sen(P.x.values[ok], y[ok])
    r = y - (a + b * P.x.values)
    r = pd.Series(r).interpolate(limit=2).values
    return pd.Series(r).rolling(2, min_periods=2).mean().values[1:]


def lagcorr(x, y, L):
    xx, yy = (x[:len(x) - L], y[L:]) if L else (x, y)
    ok = np.isfinite(xx) & np.isfinite(yy)
    return float(np.corrcoef(xx[ok], yy[ok])[0, 1]), int(ok.sum())


def r1(z):
    z = z[np.isfinite(z)]
    return float(np.corrcoef(z[:-1], z[1:])[0, 1])


def phase_surr(x, rng):
    ok = np.isfinite(x); xs = x.copy(); v = x[ok] - x[ok].mean()
    F = np.fft.rfft(v); ph = rng.uniform(0, 2 * np.pi, F.size); ph[0] = 0
    if v.size % 2 == 0:
        ph[-1] = 0
    xs[ok] = np.fft.irfft(np.abs(F) * np.exp(1j * ph), n=v.size)
    return xs


def links(P):
    up = prep(P, "up_core_speed"); dn = prep(P, "down_eke"); me = prep(P, "mea_eke")
    out = []
    for cname, ccol in CURV.items():
        cv = prep(P, ccol)
        out += [("L1_up_speed->" + cname, up, cv), ("L2_" + cname + "->down_eke", cv, dn),
                ("L3_" + cname + "->mea_eke", cv, me)]
    return out


def run_one(name, P, rng):
    rows, surr = [], {}
    for lname, x, y in links(P):
        rx, ry = r1(x), r1(y)
        obs = []
        for L in LAGS:
            r, n = lagcorr(x, y, L)
            neff = n * (1 - rx * ry) / (1 + rx * ry)
            t = r * np.sqrt(max(neff - 2, 1) / max(1 - r * r, 1e-12))
            rows.append(dict(meander=name, link=lname, lag=L, r=r, n=n, n_eff=neff,
                             p_neff=float(2 * stats.t.sf(abs(t), max(neff - 2, 1)))))
            obs.append(r)
        S = np.empty((NSURR, len(LAGS)))
        for k in range(NSURR):
            xs = phase_surr(x, rng)
            S[k] = [lagcorr(xs, y, L)[0] for L in LAGS]
        surr[lname] = S
        for i, L in enumerate(LAGS):
            rows[-len(LAGS) + i]["p_surr"] = float((1 + np.sum(np.abs(S[:, i]) >= abs(obs[i]))) / (NSURR + 1))
        pmax = float((1 + np.sum(np.max(np.abs(S), 1) >= np.max(np.abs(obs)))) / (NSURR + 1))
        for i in range(len(LAGS)):
            rows[-len(LAGS) + i]["p_surr_maxlag"] = pmax
    return pd.DataFrame(rows), surr


def pooled(T, surr):
    out = []
    for lname in T.link.unique():
        sub = T[T.link == lname]
        mz = sorted(sub.meander.unique())
        zobs = []
        for L in LAGS:
            s = sub[sub.lag == L].set_index("meander").loc[mz]
            w = np.sqrt(np.maximum(s.n_eff.values - 3, 1))
            z = float(np.sum(w * np.arctanh(s.r.values)) / np.sqrt(len(mz)))
            zobs.append(z)
            out.append(dict(link=lname, lag=L, n_meanders=len(mz), z=z, p=float(2 * stats.norm.sf(abs(z))),
                            r_by_meander=";".join(f"{m}:{v:+.2f}" for m, v in zip(mz, s.r.values)),
                            same_sign=bool(np.all(np.sign(s.r.values) == np.sign(s.r.values[0])))))
        # surrogate max-over-lags for the pooled statistic
        W = {m: np.sqrt(np.maximum(sub[(sub.meander == m)].sort_values("lag").n_eff.values - 3, 1)) for m in mz}
        Zs = sum(W[m][None, :] * np.arctanh(np.clip(surr[m][lname], -0.999, 0.999)) for m in mz) / np.sqrt(len(mz))
        pmax = float((1 + np.sum(np.max(np.abs(Zs), 1) >= np.max(np.abs(zobs)))) / (NSURR + 1))
        for row in out[-len(LAGS):]:
            row["p_surr_maxlag"] = pmax
    return pd.DataFrame(out)


if __name__ == "__main__":
    outdir = Path(sys.argv[1]); outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261008)
    allT, allS = [], {}
    for arg in sys.argv[2:]:
        name, f = arg.split("=", 1)
        T, S = run_one(name, pd.read_csv(f), rng)
        allT.append(T); allS[name] = S
        print(T.round(3).to_string(), flush=True)
    T = pd.concat(allT)
    tag = "_".join(sorted(allS))
    T.to_csv(outdir / f"b1_covar_{tag}.csv", index=False)
    if len(allS) > 1:
        Pd = pooled(T, allS)
        Pd.to_csv(outdir / f"b1_covar_pooled_{tag}.csv", index=False)
        print(Pd.round(3).to_string())
