"""A3: injection-recovery experiments for the half-power metrics.

Synthetic monthly KE fields = template(t) + noise(t), processed with the
identical pipeline as the observations (baseline D1 parameters, ACC mask,
matching to climatological envelopes, per-track anomalies averaged over the
ACC, Theil-Sen slope).

  template : the observed 2001-2025 mean monthly KE (real, asymmetric,
             multi-envelope profiles), on every 4th longitude (0.5 deg) to
             keep the cost manageable
  noise    : observed KE anomaly fields (monthly climatology removed, linear
             trend removed per pixel), resampled as whole fields in circular
             36-month blocks -> realistic spatial coherence and sub-3-yr
             autocorrelation (longer memory is NOT represented)
  scenarios: S0 none | S1 meridional shift | S2 broadening at fixed peak |
             S3 amplitude-only change | D resolution drift (meridional
             Gaussian smoothing decreasing from sigma_s to 0 over the record)

Null spread of slopes comes from S0; power = fraction of realisations outside
the central 95 % of S0; spurious trend = median slope under S3 / D.

Purpose : injection-recovery experiments for the half-power metrics (synthetic shift, broadening, amplitude, resolution drift).
Inputs  : stored monthly KE fields via $OSR11_CACHE (01b); outputs/base/mean_adt.nc
Outputs : outputs/a3/a3_realisations.csv; then outputs/a3/a3_summary.csv via 05c_detectability_summary.py (called at the end)
Run     : python scripts/02_circumpolar/02d_halfpower_detectability.py [NREAL=100]   (the shipped run: NREAL=41, i.e. 41 noisy realisations + the noise-free run)
Seeds   : numpy default_rng(123) for the block resampling.
"""
from __future__ import annotations
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from osr11 import fronts as F          # noqa: E402
from osr11 import halfpower as A       # noqa: E402
from osr11 import paths as PA           # noqa: E402

OUT = PA.out("a3")
OUT.mkdir(parents=True, exist_ok=True)
STRIDE = 4
NREAL = int(sys.argv[1]) if len(sys.argv) > 1 else 100
BLOCK = 36
DLAT_KM = 6371.0 * np.pi / 180.0 * 0.125      # grid step (km)


def shift_field(T, delta_km):
    """Shift every profile northward by delta_km (linear interpolation)."""
    g = delta_km / DLAT_KM
    n = int(np.floor(g)); f = g - n
    a = np.roll(T, n, axis=0); b = np.roll(T, n + 1, axis=0)
    out = (1 - f) * a + f * b
    if n + 1 > 0:
        out[: n + 1] = np.nan
    elif n < 0:
        out[n:] = np.nan
    return out


def _fft_gauss(A, sigma_grid):
    """Exact Gaussian smoothing along axis 0 in Fourier space (any sigma >= 0)."""
    n = A.shape[0]
    npad = 2 * n
    k = np.fft.rfftfreq(npad) * 2 * np.pi
    H = np.exp(-0.5 * (k * sigma_grid) ** 2)
    P = np.concatenate([A, A[::-1]], axis=0)       # even extension limits edge effects
    S = np.fft.irfft(np.fft.rfft(P, axis=0) * H[:, None], n=npad, axis=0)
    return S[:n]


def nan_smooth(T, sigma_grid):
    if sigma_grid <= 0:
        return T
    m = np.isfinite(T)
    num = _fft_gauss(np.where(m, T, 0.0), sigma_grid)
    den = _fft_gauss(m.astype(float), sigma_grid)
    out = num / np.where(den > 1e-6, den, np.nan)
    out[~m] = np.nan
    return out


def broaden(T, frac):
    """Broaden profiles by ~frac in width at fixed per-longitude peak:
    Gaussian convolution with sigma_b = sqrt((1+frac)^2-1) * sigma_jet,
    sigma_jet ~ 36 km (median second-moment width of observed envelopes)."""
    if frac <= 0:
        return T
    sb = np.sqrt((1 + frac) ** 2 - 1) * 36.0 / DLAT_KM
    S = nan_smooth(T, sb)
    return S * (np.nanmax(T, axis=0) / np.nanmax(S, axis=0))[None, :]


def main():
    ke, lat, lon, acc, t = A.load()
    lonS = lon[::STRIDE]
    accS = acc[:, ::STRIDE]
    keS = np.asarray(ke[:, :, ::STRIDE], np.float64)
    nt = keS.shape[0]
    month = t.month.values
    clim = np.stack([np.nanmean(keS[month == m], axis=0) for m in range(1, 13)])
    anom = keS - clim[month - 1]
    # remove linear trend per pixel
    xx = t.x.values - t.x.values.mean()
    good = np.isfinite(anom)
    a0 = np.where(good, anom, 0.0)
    cnt = good.sum(0)
    b = np.einsum("t,tij->ij", xx, a0) / np.maximum(np.einsum("t,tij->ij", xx ** 2, good.astype(float)), 1e-12)
    anom = anom - xx[:, None, None] * b[None]
    T0 = np.nanmean(keS, axis=0)                 # template
    tc = (t.x.values - t.x.values.mean())         # years from mid-record
    yrs = t.x.values.max() - t.x.values.min()

    scen = {
        "S0_null": {},
        "S1_shift5km": dict(shift=5.0), "S1_shift10km": dict(shift=10.0), "S1_shift25km": dict(shift=25.0), "S1_shift50km": dict(shift=50.0),
        "S2_broad0.25pct": dict(broad=0.0025), "S2_broad0.5pct": dict(broad=0.005), "S2_broad1pct": dict(broad=0.01), "S2_broad2pct": dict(broad=0.02), "S2_broad5pct": dict(broad=0.05),
        "S3_amp0.15": dict(amp=0.0015), "S3_amp0.30": dict(amp=0.0030), "S3_amp0.60": dict(amp=0.0060),
        "D_res5km": dict(drift=5.0), "D_res10km": dict(drift=10.0),
    }
    p = A.BASE
    d = F.meridional_dists(lat)
    rows = []
    rng = np.random.default_rng(123)
    for r in range(NREAL + 1):          # r = 0 is the noise-free run
        if r == 0:
            noise = np.zeros_like(anom)
        else:
            nb = int(np.ceil(nt / BLOCK))
            starts = rng.integers(0, nt, nb)
            idx = ((starts[:, None] + np.arange(BLOCK)[None]) % nt).ravel()[:nt]
            noise = np.nan_to_num(anom[idx], nan=0.0)
        for name, s in scen.items():
            fields = np.empty_like(keS)
            for k in range(nt):
                frac_t = (k + 0.5) / nt
                Tk = T0
                if "shift" in s:
                    Tk = shift_field(Tk, s["shift"] * (frac_t - 0.5))
                if "broad" in s:
                    Tk = broaden(Tk, s["broad"] * frac_t)
                if "amp" in s:
                    Tk = Tk * (1 + s["amp"] * tc[k])
                fk = Tk + noise[k]
                if "drift" in s:
                    fk = nan_smooth(fk, s["drift"] * (1 - frac_t) / DLAT_KM)
                fk = np.where(np.isfinite(T0), np.maximum(fk, 0.0), np.nan)
                fields[k] = fk
            counts = np.zeros((nt, lonS.size), np.int16)
            tab = np.full((nt, lonS.size, A.MAX_ENV, F.NF), np.nan, np.float32)
            for k in range(nt):
                c, tb = F.detect_field(fields[k], lat, d, p["thr"], p["min_w"], p["frac"], A.MAX_ENV)
                counts[k] = c; tab[k] = tb
            ccounts, ctab = F.detect_field(np.nanmean(fields, 0), lat, d, p["thr"], p["min_w"], p["frac"], A.MAX_ENV)
            arr, meta = A.matched_tracks(counts, tab, ccounts, ctab, lat, lonS, accS)
            keep = meta.keep.values
            rec = {"real": r, "scenario": name, "ntracks": int(keep.sum())}
            for f in ["lat_hp", "seg_km", "fwhm_km", "w25_km", "sig_km", "peak_val"]:
                a_ = arr[f][:, keep].astype(float)
                y = np.nanmean(a_ - np.nanmean(a_, 0, keepdims=True), axis=1)
                from osr11.rstats import theil_sen
                rec[f + "_slope"] = theil_sen(t.x.values, y)[0]
            # unmatched (original-style) metrics
            um = A.unmatched_series(counts, tab, lat, lonS, accS, t)
            u = um[um.region == "ACC"].sort_values("i")
            from osr11.rstats import theil_sen as TS
            rec["u_median_lat_slope"] = TS(u.x.values, u.median_lat.values)[0]
            rec["u_seg_km_slope"] = TS(u.x.values, u.seg_km.values)[0]
            rec["u_n_slope"] = TS(u.x.values, u.n.values.astype(float))[0]
            rows.append(rec)
        pd.DataFrame(rows).to_csv(OUT / "a3_realisations.csv", index=False)
        print("realisation", r, "done", flush=True)


if __name__ == "__main__":
    main()
    import runpy
    runpy.run_path(str(PA.ROOT / "scripts" / "05_stats" / "05c_detectability_summary.py"), run_name="__main__")
