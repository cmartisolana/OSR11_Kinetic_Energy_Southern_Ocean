"""x11 (9 Nov 2026): lag between the wind and EKE following the earlier observational studies.
Hogg et al. (2015): EKE averaged in the sector boxes of Meredith and Hogg (2006) (Indian 40-150E 57-44S; Pacific 150E-72W 62-48S;
Atlantic 35W-10E 56-46S) correlated with the 1-year running mean of the Marshall (2003) SAM index, SAM leading; no detrending stated;
reported peak correlation at a 3-year lag (Pacific 0.50, Indian 0.32; Atlantic n.s.). Meredith and Hogg (2006): annual means, lag from
the timing of the SAM and EKE peaks (2-3 yr).
Here: Product 1 monthly EKE in the same boxes (outputs/m1/m1_series_allsat.csv, 2001-2025) and over the ACC band; Marshall SAM index
(BAS station index, monthly 1997-2025, x11/marshall_sam_monthly_1997_2025.csv). Both as 12-month running means; Pearson r with SAM
leading by 0-48 months (monthly steps; EKE 2001-2025 fixed, SAM shifted). Significance: 2000 phase-randomised surrogates of the monthly
SAM (1997-2025), smoothed and lagged identically; pointwise p and family-wise p (max |r| over lags 0-48 months). As a sensitivity,
the same with both series linearly detrended.
Outputs: x11/x11_lag_corr.csv, x11/x11_summary.csv, x11/x11_log.txt

Origin  : x11_lag_previous_studies.py (9 Nov 2026); moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Purpose : Fig. S4d and SI S3.2 (lag of EKE behind the SAM index, method of Hogg et al., 2015).
Inputs  : outputs/m1/m1_series_allsat.csv (02l), outputs/x11/marshall_sam_monthly_1997_2025.csv (Marshall, 2003; British Antarctic Survey)
Run     : python scripts/05_stats/05l_sam_lag.py   (no raw data; minutes)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import numpy as np, pandas as pd
from osr11 import paths as PA
OUT = PA.out("x11")
E = pd.read_csv(PA.OUT / "m1" / "m1_series_allsat.csv")
S = pd.read_csv(OUT / "marshall_sam_monthly_1997_2025.csv")
S["t"] = S.year * 12 + S.month - 1; E["t"] = E.year * 12 + E.month - 1
sam = S.sort_values("t").sam.values; ts = S.sort_values("t").t.values
run = lambda x: pd.Series(x).rolling(12, center=True).mean().values
def detr(x):
    m = np.isfinite(x); t = np.arange(x.size); y = x.copy(); y[m] -= np.polyval(np.polyfit(t[m], x[m], 1), t[m]); return y
REG = {"Pacific": "Hogg_Pacific_EKE", "Indian": "Hogg_Indian_EKE", "Atlantic": "Hogg_Atlantic_EKE", "ACC band": "ACC_EKE"}
LAGS = np.arange(0, 49); te = E.t.values
def lagcorr(samser, eke):
    out = []
    for L in LAGS:
        idx = np.searchsorted(ts, te - L); x = samser[idx]; m = np.isfinite(x) & np.isfinite(eke)
        out.append(np.corrcoef(x[m], eke[m])[0, 1])
    return np.array(out)
rng = np.random.default_rng(20261109)
def surr(y):
    f = np.fft.rfft(y - y.mean()); ph = rng.uniform(0, 2 * np.pi, f.size); ph[0] = 0
    if y.size % 2 == 0: ph[-1] = 0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=y.size) + y.mean()
SUR = [surr(sam) for _ in range(2000)]
rows, summ, log = [], [], []
P = lambda s: (print(s, flush=True), log.append(s))
for variant in ("as Hogg et al. (2015): not detrended", "both detrended"):
    dt = variant.startswith("both")
    s_main = run(sam); s_main = detr(s_main) if dt else s_main
    SS = [detr(run(s)) if dt else run(s) for s in SUR]
    for reg, col in REG.items():
        e = run(E[col].values); e = detr(e) if dt else e
        r = lagcorr(s_main, e); nul = np.array([lagcorr(s, e) for s in SS])
        p_pt = (1 + np.sum(np.abs(nul) >= np.abs(r)[None, :], 0)) / (len(SS) + 1)
        p_fw = (1 + np.sum(np.abs(nul).max(1)[:, None] >= np.abs(r)[None, :], 0)) / (len(SS) + 1)
        for L, a, b, c in zip(LAGS, r, p_pt, p_fw): rows.append(dict(variant=variant, region=reg, lag_months=L, r=a, p_pointwise=b, p_familywise=c))
        k = int(np.argmax(r)); yr = {y: r[12 * y] for y in range(5)}
        summ.append(dict(variant=variant, region=reg, lag_max_months=LAGS[k], r_max=r[k], p_pt_max=p_pt[k], p_fw_max=p_fw[k], **{f"r_lag{y}yr": v for y, v in yr.items()},
                         **{f"p_pt_lag{y}yr": p_pt[12 * y] for y in range(5)}))
        P(f"[{variant}] {reg:9s}: r at 0/1/2/3/4 yr = " + ", ".join(f"{yr[y]:+.2f}" for y in range(5)) +
          f" | max r {r[k]:+.2f} at {LAGS[k]} months (pointwise p {p_pt[k]:.3f}, family-wise p {p_fw[k]:.3f})")
pd.DataFrame(rows).to_csv(OUT / "x11_lag_corr.csv", index=False); pd.DataFrame(summ).to_csv(OUT / "x11_summary.csv", index=False)
(OUT / "x11_log.txt").write_text("\n".join(log))
