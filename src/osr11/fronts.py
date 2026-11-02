"""Half-power envelope detection for the reviewer-2 revision (plan items A0, A2).

`detect_profile` reproduces the original algorithm in
scripts/utils/chambers_fronts.py exactly (segment above an absolute threshold,
minimum segment width, envelope bounds at `frac` x peak searched inside the
segment, trapezoid integral, interpolated half-power latitude) and adds:

  * seg_km     width of the above-threshold segment  (= original `width_km`)
  * env_km     width between the integration bounds (isouth..inorth)
  * edge       True when both integration bounds are the segment edges
  * fwhm_km    interpolated full width at half maximum, searched outward from
               the peak on the full profile and stopped at a local minimum;
               invariant to a pure amplitude change of a fixed-shape profile
  * w25_km     same at 25 % of peak (local-minimum bounded)
  * sig_km     KE-weighted second-moment width (std of latitude in km) between
               the w25 bounds
  * peak_lat, peak_val

Parameters are always passed explicitly (no hidden defaults).
"""

from __future__ import annotations

import numpy as np
from numba import njit

R_EARTH_KM = 6371.0
NF = 13  # fields per envelope (see FIELDS)
FIELDS = ("lat_hp", "peak_lat", "peak_val", "seg_km", "env_km", "edge",
          "fwhm_km", "w25_km", "sig_km", "i_start", "i_end", "lat_s", "lat_n")


@njit(cache=True)
def _cross(prof, lats, ipk, level, step):
    """Walk from ipk in direction step until prof <= level (interpolate) or a
    local minimum / NaN / edge is met. Returns (latitude, hit_level flag)."""
    n = prof.shape[0]
    k = ipk
    while True:
        k2 = k + step
        if k2 < 0 or k2 >= n or not np.isfinite(prof[k2]):
            return lats[k], False
        if prof[k2] <= level:
            # interpolate between k and k2
            d = prof[k] - prof[k2]
            t = (prof[k] - level) / d if d != 0.0 else 0.0
            return lats[k] + t * (lats[k2] - lats[k]), True
        if prof[k2] > prof[k]:          # local minimum at k
            return lats[k], False
        k = k2


@njit(cache=True)
def detect_profile(prof, lats, dists, thr, min_w, frac, out):
    """Detect envelopes on one south->north profile. Writes rows into `out`
    (max_env x NF) and returns the number of envelopes."""
    n = prof.shape[0]
    deg2km = R_EARTH_KM * np.pi / 180.0
    found = 0
    i = 0
    while i < n:
        if not (np.isfinite(prof[i]) and prof[i] >= thr):
            i += 1
            continue
        s = i
        j = i
        while j + 1 < n and np.isfinite(prof[j + 1]) and prof[j + 1] >= thr:
            j += 1
        e = j
        i = e + 1
        seg = dists[e] - dists[s]
        if seg < min_w:
            continue
        ipk = s
        pv = prof[s]
        for kk in range(s, e + 1):
            if prof[kk] > pv:
                pv = prof[kk]
                ipk = kk
        lev = pv * frac
        isouth = s
        for kk in range(s, ipk + 1):
            if prof[kk] < lev:
                isouth = kk
        inorth = e
        for kk in range(ipk, e + 1):
            if prof[kk] < lev:
                inorth = kk
                break
        if inorth <= isouth:
            continue
        total = 0.0
        for kk in range(isouth, inorth):
            total += 0.5 * (prof[kk] + prof[kk + 1]) * (dists[kk + 1] - dists[kk])
        if total <= 0.0:
            continue
        half = 0.5 * total
        cum = 0.0
        lat_hp = np.nan
        for kk in range(isouth, inorth):
            area = 0.5 * (prof[kk] + prof[kk + 1]) * (dists[kk + 1] - dists[kk])
            if cum + area >= half:
                f = (half - cum) / area if area > 0 else 0.0
                lat_hp = lats[kk] + f * (lats[kk + 1] - lats[kk])
                break
            cum += area
        if not np.isfinite(lat_hp):
            continue
        # amplitude-invariant widths
        ls5, _ = _cross(prof, lats, ipk, 0.5 * pv, -1)
        ln5, _ = _cross(prof, lats, ipk, 0.5 * pv, 1)
        ls2, _ = _cross(prof, lats, ipk, 0.25 * pv, -1)
        ln2, _ = _cross(prof, lats, ipk, 0.25 * pv, 1)
        # second moment between w25 bounds (grid points inside)
        sw = 0.0
        sy = 0.0
        syy = 0.0
        for kk in range(n):
            if lats[kk] >= ls2 and lats[kk] <= ln2 and np.isfinite(prof[kk]):
                w = prof[kk]
                sw += w
                sy += w * lats[kk]
                syy += w * lats[kk] * lats[kk]
        sig = np.nan
        if sw > 0:
            var = syy / sw - (sy / sw) ** 2
            sig = np.sqrt(max(var, 0.0)) * deg2km
        if found < out.shape[0]:
            o = out[found]
            o[0] = lat_hp
            o[1] = lats[ipk]
            o[2] = pv
            o[3] = seg
            o[4] = dists[inorth] - dists[isouth]
            o[5] = 1.0 if (isouth == s and inorth == e) else 0.0
            o[6] = (ln5 - ls5) * deg2km
            o[7] = (ln2 - ls2) * deg2km
            o[8] = sig
            o[9] = s
            o[10] = e
            o[11] = lats[isouth]
            o[12] = lats[inorth]
            found += 1
    return found


@njit(cache=True)
def detect_field(field, lats, dists, thr, min_w, frac, max_env):
    """Run detect_profile on every longitude of a (lat, lon) field.
    Returns (counts[lon], table[lon, max_env, NF])."""
    nlat, nlon = field.shape
    counts = np.zeros(nlon, np.int32)
    table = np.full((nlon, max_env, NF), np.nan)
    buf = np.empty((max_env, NF))
    for j in range(nlon):
        prof = field[:, j].astype(np.float64)
        c = detect_profile(prof, lats, dists, thr, min_w, frac, buf)
        counts[j] = c
        for k in range(c):
            table[j, k, :] = buf[k, :]
    return counts, table


def meridional_dists(lats: np.ndarray) -> np.ndarray:
    lats = np.asarray(lats, float)
    d = np.zeros_like(lats)
    d[1:] = np.cumsum(R_EARTH_KM * np.abs(np.diff(np.radians(lats))))
    return d
