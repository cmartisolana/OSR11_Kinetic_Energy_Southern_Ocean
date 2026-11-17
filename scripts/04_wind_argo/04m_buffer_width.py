"""x14 (10 Nov 2026): sensitivity of the Sect. 3.3 / Fig. 4 diagnostics to (A) the width of the extension of the ACC band
(band only, +100 km = manuscript, +200 km, +300 km) and (B) the period of the EKE/KE response (from 2002 instead of 2004; the
wind extraction starts in 2001, so a one-year lag allows responses from 2002 only).
Methods as in the manuscript: trends rstats.trend_bcz (monthly, 2004-2025, % of mean per decade); lagged correlation as x12
(deseasonalised, detrended, 12-month running mean, wind leading 0-L months, phase-randomised wind surrogates, family-wise p);
exponent b = OLS slope of ln X (annual) on ln tau of the previous year, trends retained and both detrended, intervals as x10
(fitted model + phase-randomised wind, factor 1.07), p from wind surrogates (x9, fixed lag 1 yr).
Inputs: x13/x13_cells.npz. Outputs: x14/x14_series_<dom>.csv, x14/x14_results.csv, x14/x14_log.txt

Origin  : x14_buffer_and_period.py (10 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : SI S5.2, independence of the Sect. 3.3 results from the width of the extension of the ACC band.
Run     : python scripts/04_wind_argo/04m_buffer_width.py   (no raw data; about 15 min)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd
from scipy import stats
from osr11 import rstats as R, paths as PA
OUT = PA.out("x14")
Z = np.load(PA.OUT / "x13" / "x13_cells.npz"); lat = Z["lat"]; YY, MO = Z["year"], Z["month"]
w = np.cos(np.radians(np.meshgrid(Z["lon"], lat)[1])); x = YY + (MO - 0.5) / 12; c4 = YY >= 2004
fin = np.isfinite(Z["EC"]).all(0) & np.isfinite(Z["KC"]).all(0) & np.isfinite(Z["WC"]).all(0) & np.isfinite(Z["SL"]).all(0) & (Z["ice"] == 0)
DOMS = {"band": (Z["band"] > 0) & fin, "buf100": (Z["dmin"] <= 100 + 1e-6) & fin & Z["slope_ok"],
        "buf100_noslopeok": (Z["dmin"] <= 100 + 1e-6) & fin, "buf200": (Z["dmin"] <= 200 + 1e-6) & fin, "buf300": (Z["dmin"] <= 300 + 1e-6) & fin}
log = []; P = lambda s: (print(s, flush=True), log.append(s))
def dm(F, d): return np.array([np.sum(F[q][d] * w[d]) / np.sum(w[d]) for q in range(F.shape[0])])
SER = {}
for k, d in DOMS.items():
    s = pd.DataFrame({"year": YY, "month": MO, "x": x, "wind": dm(Z["WC"], d), "KE": dm(Z["KC"], d), "EKE": dm(Z["EC"], d),
                      "slope": np.r_[np.full(36, np.nan), dm(Z["SL"], d)]})
    s.to_csv(OUT / f"x14_series_{k}.csv", index=False); SER[k] = s; P(f"{k}: {int(d.sum())} cells")
ref = pd.read_csv(PA.OUT / "f4e" / "f4e_series.csv")
for c in ("wind", "KE", "EKE", "slope"):
    a, b = ref[c].values, SER["buf100"][c].values; g = np.isfinite(a)
    P(f"check vs f4e_series ({c}): max rel diff {np.nanmax(np.abs(b[g] / a[g] - 1)):.2e}")
rng = np.random.default_rng(20261110)
def surr(y, m):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, (m, f.size)); ph[:, 0] = 0
    if y.size % 2 == 0: ph[:, -1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size, axis=-1) + y.mean()
def prep(v, sel):
    out = np.full(v.size, np.nan); a = 100 * (v[sel] / np.nanmean(v[sel]) - 1); m = MO[sel]
    a = a - pd.Series(a).groupby(m).transform("mean").values; t = np.arange(a.size); a = a - np.polyval(np.polyfit(t, a, 1), t)
    out[sel] = pd.Series(a).rolling(12, center=True).mean().values; return out
def lagcorr(S, k, y0, Lmax, NS=1000):
    Wp = prep(S.wind.values, np.ones(len(S), bool)); Xp = prep(S[k].values, YY >= y0); L = np.arange(Lmax + 1)
    first_w = np.where(np.isfinite(Wp))[0][0]
    valid = np.where(np.isfinite(Xp))[0]; valid = valid[valid - Lmax >= first_w]
    r = np.array([np.corrcoef(Wp[valid - l], Xp[valid])[0, 1] for l in L])
    wv = S.wind.values; wa = 100 * (wv / wv.mean() - 1); wa = wa - pd.Series(wa).groupby(MO).transform("mean").values
    mx = np.empty(NS)
    for i, s_ in enumerate(surr(wa, NS)):
        t = np.arange(s_.size); s_ = s_ - np.polyval(np.polyfit(t, s_, 1), t); s_ = pd.Series(s_).rolling(12, center=True).mean().values
        mx[i] = np.max(np.abs([np.corrcoef(s_[valid - l], Xp[valid])[0, 1] for l in L]))
    km = int(np.argmax(r)); return dict(lag=int(L[km]), r=r[km], p_fw=(1 + np.sum(mx >= abs(r[km]))) / (NS + 1), n_months=valid.size,
                                        first=f"{YY[valid[0]]}-{MO[valid[0]]:02d}")
yrs = np.arange(2001, 2026)
def ols(a, b): a = a - a.mean(-1, keepdims=True); return (a * (b - b.mean(-1, keepdims=True))).sum(-1) / (a * a).sum(-1)
def detr(y):
    t = np.arange(y.shape[-1]); tc = t - t.mean(); return y - y.mean(-1, keepdims=True) - (((y - y.mean(-1, keepdims=True)) @ tc) / (tc @ tc))[..., None] * tc
def exponent(S, k, y0, m=2000):
    ann = lambda c: np.array([S[c][S.year == q].mean() for q in yrs]); lw = np.log(ann("wind") / ann("wind").mean())
    sy = yrs >= y0; sw = (yrs >= y0 - 1) & (yrs <= 2024); v = ann(k); yy = np.log(v / np.nanmean(v[sy]))[sy]; wx = lw[sw]; n = wx.size; t = np.arange(n)
    b = float(ols(wx, yy)); bd = float(ols(detr(wx), detr(yy)))
    SW = surr(lw, m)[:, sw]; p = (1 + np.sum(np.abs(ols(SW, np.broadcast_to(yy, SW.shape))) >= abs(b))) / (m + 1)
    pd_ = (1 + np.sum(np.abs(ols(detr(SW), np.broadcast_to(detr(yy), SW.shape))) >= abs(bd))) / (m + 1)
    A = np.c_[np.ones(n), wx, t]; coef, *_ = np.linalg.lstsq(A, yy, rcond=None); e = yy - A @ coef
    phi = np.corrcoef(e[1:], e[:-1])[0, 1]; sig = e.std(ddof=3) * np.sqrt(1 - min(phi, 0.95) ** 2); phic = float(np.clip(phi + (1 + 3 * phi) / n, -0.9, 0.95))
    E = np.empty((m, n)); E[:, 0] = rng.normal(0, sig / np.sqrt(1 - phic ** 2), m)
    for i in range(1, n): E[:, i] = phic * E[:, i - 1] + rng.normal(0, sig, m)
    TF = stats.t.ppf(0.975, n - 3) / stats.norm.ppf(0.975)
    bs = ols(SW, b * SW + coef[2] * t + E); lo, hi = b - TF * np.percentile(bs - b, [97.5, 2.5])
    bsd = ols(detr(SW), detr(bd * SW + coef[2] * t + E)); lod, hid = bd - TF * np.percentile(bsd - bd, [97.5, 2.5])
    return dict(n_pairs=n, b=b, b_lo=lo, b_hi=hi, p_b=p, b_detr=bd, bd_lo=lod, bd_hi=hid, p_bd=pd_)
rows = []
for dname in ("band", "buf100", "buf200", "buf300"):
    S = SER[dname]
    for k in ("wind", "KE", "EKE", "slope"):
        y = S[k].values[c4]; r = R.trend_bcz(x[c4], y, month=MO[c4], n_outer=300); mu = y.mean()
        row = dict(domain=dname, test="trend 2004-2025", series=k, pct=1000 * r["slope"] / mu, lo=1000 * r["ci_lo"] / mu, hi=1000 * r["ci_hi"] / mu, p=r["p"])
        rows.append(row); P(f"{dname} trend {k}: {row['pct']:+.1f} [{row['lo']:+.1f}, {row['hi']:+.1f}] % /dec, p {row['p']:.3f}")
    for k in ("KE", "EKE", "slope"):
        q = lagcorr(S, k, 2004, 36); rows.append(dict(domain=dname, test="lagcorr 2004 L0-36", series=k, **q))
        P(f"{dname} lag {k}: r {q['r']:.2f} at {q['lag']} mo, p_fw {q['p_fw']:.3f} ({q['n_months']} months from {q['first']})")
        e = exponent(S, k, 2004); rows.append(dict(domain=dname, test="exponent 2004-2025", series=k, **e))
        P(f"{dname} b {k}: {e['b']:.2f} [{e['b_lo']:.2f}, {e['b_hi']:.2f}] p {e['p_b']:.3f} | detr {e['b_detr']:.2f} [{e['bd_lo']:.2f}, {e['bd_hi']:.2f}] p {e['p_bd']:.3f}")
S = SER["buf100"]
for k in ("KE", "EKE"):
    e = exponent(S, k, 2002); rows.append(dict(domain="buf100", test="exponent 2002-2025", series=k, **e))
    P(f"PERIOD b {k} 2002-2025 ({e['n_pairs']} pairs): {e['b']:.2f} [{e['b_lo']:.2f}, {e['b_hi']:.2f}] p {e['p_b']:.3f} | detr {e['b_detr']:.2f} [{e['bd_lo']:.2f}, {e['bd_hi']:.2f}] p {e['p_bd']:.3f}")
    for y0, L in ((2001, 36), (2001, 24), (2004, 24)):
        q = lagcorr(S, k, y0, L); rows.append(dict(domain="buf100", test=f"lagcorr {y0} L0-{L}", series=k, **q))
        P(f"PERIOD lag {k} resp from {y0}, L0-{L}: r {q['r']:.2f} at {q['lag']} mo, p_fw {q['p_fw']:.3f} ({q['n_months']} months from {q['first']})")
    y = S[k].values; r = R.trend_bcz(x, y, month=MO, n_outer=300); mu = y.mean()
    rows.append(dict(domain="buf100", test="trend 2001-2025", series=k, pct=1000 * r["slope"] / mu, lo=1000 * r["ci_lo"] / mu, hi=1000 * r["ci_hi"] / mu, p=r["p"]))
    P(f"PERIOD trend {k} 2001-2025: {rows[-1]['pct']:+.1f} [{rows[-1]['lo']:+.1f}, {rows[-1]['hi']:+.1f}] p {rows[-1]['p']:.3f}")
pd.DataFrame(rows).to_csv(OUT / "x14_results.csv", index=False); (OUT / "x14_log.txt").write_text("\n".join(log)); P("done")
