"""Chambers front detection utilities.

Functions:
- meridional_distances_km_from_lats
- detect_envelopes_on_profile
- detect_fronts_from_dataset

Defaults follow Chambers (2018): threshold=0.02 m^2/s^2, min_width_km=100 km.
"""

from __future__ import annotations

import numpy as np
import xarray as xr
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Optional


# Optional numba acceleration
try:
    import numba
    from numba import njit
    _HAVE_NUMBA = True
except Exception:
    _HAVE_NUMBA = False


if _HAVE_NUMBA:
    @njit
    def _numba_detect_envelopes(profile, lats, lons, dists, threshold, min_width_km, peak_fraction):
        """Numba-accelerated inner loop to detect envelopes on a 1D profile.

        This helper is compiled with Numba (njit) and implements the same
        envelope / half-power calculation as the pure-Python version. It
        operates entirely on 1D numpy arrays and returns a compact set of
        parallel arrays describing found envelopes.

        Parameters
        ----------
        profile, lats, lons, dists : 1D float64 arrays
            Arrays of equal length describing the meridional profile (profile),
            the latitudes and longitudes of each grid point, and cumulative
            meridional distances `dists` (in km).
        threshold : float
            Minimum profile value to consider (Chambers-style threshold).
        min_width_km : float
            Minimum envelope width (km) to keep.
        peak_fraction : float
            Fraction of the peak used to determine envelope boundaries
            (e.g. 0.25 for 25%% of peak).

        Returns
        -------
        found : int
            Number of envelopes found.
        lat_mid, lon_mid, isouth_arr, inorth_arr, peak_lat_arr, peak_val_arr,
        width_arr : arrays
            Parallel arrays of length >= found containing envelope metadata:
            mid-point latitude/longitude (half-power point), southern and
            northern boundary indices, peak latitude/value and envelope width.

        Notes
        -----
        - The function is deliberately low-level and uses explicit loops to
            be compatible with Numba's nopython mode.
        - The first call to the Numba-compiled function includes JVM/LLVM
            compile overhead; subsequent calls are very fast.
        """
        # using numba for acceleration (prints removed in compiled mode)
        # profile, lats, lons, dists: 1D float64 arrays
        n = profile.shape[0]
        max_segs = n // 2 + 2
        lat_mid = np.empty(max_segs, dtype=np.float64)
        lon_mid = np.empty(max_segs, dtype=np.float64)
        isouth_arr = np.empty(max_segs, dtype=np.int64)
        inorth_arr = np.empty(max_segs, dtype=np.int64)
        peak_lat_arr = np.empty(max_segs, dtype=np.float64)
        peak_val_arr = np.empty(max_segs, dtype=np.float64)
        width_arr = np.empty(max_segs, dtype=np.float64)

        found = 0
        # mask invalid
        valid = np.isfinite(profile)
        # find indices above threshold
        above = np.empty(n, dtype=np.bool_)
        for i in range(n):
            above[i] = valid[i] and (profile[i] >= threshold)

        # find contiguous segments
        i = 0
        while i < n:
            # skip until above
            if not above[i]:
                i += 1
                continue
            start = i
            j = i
            while j + 1 < n and above[j+1]:
                j += 1
            end = j

            # check width
            width_km = dists[end] - dists[start]
            if width_km >= min_width_km:
                # segment arrays
                # find relative peak
                rel_peak = start
                pval = profile[start]
                for kk in range(start, end+1):
                    if profile[kk] > pval:
                        pval = profile[kk]
                        rel_peak = kk

                peak_global_idx = rel_peak
                peak_val = pval

                south_thresh = peak_val * peak_fraction

                # south boundary: last index before peak where profile < south_thresh
                isouth = start
                for kk in range(start, peak_global_idx+1):
                    if profile[kk] < south_thresh:
                        isouth = kk
                # if no such point found above, isouth remains start

                # north boundary: first index after peak where profile < south_thresh
                inorth = end
                for kk in range(peak_global_idx, end+1):
                    if profile[kk] < south_thresh:
                        inorth = kk
                        break

                if inorth > isouth:
                    # integrate y over x using trapezoid
                    # compute total = trapz(y, x)
                    total = 0.0
                    for kk in range(isouth, inorth):
                        dx = dists[kk+1] - dists[kk]
                        total += 0.5 * (profile[kk] + profile[kk+1]) * dx

                    if total > 0.0:
                        # cumulative
                        cum0 = 0.0
                        # We will search for half along segments
                        half = total * 0.5
                        cum = 0.0
                        placed = False
                        for kk in range(isouth, inorth):
                            dx = dists[kk+1] - dists[kk]
                            area = 0.5 * (profile[kk] + profile[kk+1]) * dx
                            if cum + area >= half and not placed:
                                # interpolate between kk and kk+1
                                if area <= 0.0:
                                    frac = 0.0
                                else:
                                    frac = (half - cum) / area
                                latm = lats[kk] + frac * (lats[kk+1] - lats[kk])
                                lonm = lons[kk] + frac * (lons[kk+1] - lons[kk])
                                lat_mid[found] = latm
                                lon_mid[found] = lonm
                                isouth_arr[found] = isouth
                                inorth_arr[found] = inorth
                                peak_lat_arr[found] = lats[peak_global_idx]
                                peak_val_arr[found] = peak_val
                                width_arr[found] = width_km
                                found += 1
                                placed = True
                            cum += area

            i = end + 1

        return found, lat_mid, lon_mid, isouth_arr, inorth_arr, peak_lat_arr, peak_val_arr, width_arr


def meridional_distances_km_from_lats(lats: np.ndarray) -> np.ndarray:
    """Compute cumulative meridional distances (km) for a 1D latitude array.

    Vectorized estimate assuming longitude roughly constant for a meridional
    column. Distance between adjacent latitudes is approximated by R*|dlat|.
    """
    lats = np.asarray(lats, dtype=float)
    R = 6371.0
    latr = np.radians(lats)
    dlat = np.diff(latr)
    step = R * np.abs(dlat)
    dists = np.empty(lats.shape, dtype=float)
    dists[0] = 0.0
    if step.size > 0:
        dists[1:] = np.cumsum(step)
    return dists


def detect_envelopes_on_profile(profile: np.ndarray,
                                lats: np.ndarray,
                                lons: np.ndarray,
                                dists: Optional[np.ndarray] = None,
                                threshold: float = 0.02,
                                min_width_km: float = 100.0,
                                peak_fraction: float = 0.25,
                                use_numba: bool = False) -> List[Dict]:
    """Detect contiguous high-CKE envelopes and return half-power points.

    Parameters
    ----------
    profile : ndarray
        1D EKE profile ordered south -> north.
    lats, lons : ndarray
        Coordinates matching profile.
    dists : ndarray or None
        Precomputed cumulative distances (km) matching profile.
    """
    profile = np.asarray(profile, dtype=float)
    lats = np.asarray(lats, dtype=float)
    lons = np.asarray(lons, dtype=float)

    if dists is None:
        dists = meridional_distances_km_from_lats(lats)

    # If requested and available, use the numba-accelerated inner routine.
    if use_numba and _HAVE_NUMBA:
        # Ensure float64 arrays for numba
        prof = np.asarray(profile, dtype=np.float64)
        lat_a = np.asarray(lats, dtype=np.float64)
        lon_a = np.asarray(lons, dtype=np.float64)
        d_a = np.asarray(dists, dtype=np.float64)
        found, lat_mid, lon_mid, isouth_arr, inorth_arr, peak_lat_arr, peak_val_arr, width_arr = _numba_detect_envelopes(
            prof, lat_a, lon_a, d_a, float(threshold), float(min_width_km), float(peak_fraction)
        )
        results: List[Dict] = []
        for ii in range(found):
            results.append({
                'lat_mid': float(lat_mid[ii]),
                'lon_mid': float(lon_mid[ii]),
                'isouth': int(isouth_arr[ii]),
                'inorth': int(inorth_arr[ii]),
                'peak_lat': float(peak_lat_arr[ii]),
                'peak_val': float(peak_val_arr[ii]),
                'width_km': float(width_arr[ii])
            })
        return results

    valid = np.isfinite(profile)
    if not np.any(valid):
        return []

    above = (profile >= threshold) & valid
    if not np.any(above):
        return []

    idx = np.where(above)[0]
    if idx.size == 0:
        return []

    gaps = np.where(np.diff(idx) > 1)[0]
    if gaps.size == 0:
        segments = [(int(idx[0]), int(idx[-1]))]
    else:
        starts = np.concatenate(([idx[0]], idx[gaps+1])).astype(int)
        ends = np.concatenate((idx[gaps], [idx[-1]])).astype(int)
        segments = list(zip(starts, ends))

    results: List[Dict] = []
    for iseg, jseg in segments:
        width_km = float(dists[jseg] - dists[iseg])
        if width_km < min_width_km:
            continue

        seg_profile = profile[iseg:jseg+1]
        seg_lats = lats[iseg:jseg+1]
        seg_lons = lons[iseg:jseg+1]

        rel_peak_idx = int(np.nanargmax(seg_profile))
        peak_val = float(seg_profile[rel_peak_idx])
        peak_global_idx = iseg + rel_peak_idx

        south_thresh = peak_val * peak_fraction

        sub = profile[iseg:peak_global_idx+1]
        below_south = np.where(sub < south_thresh)[0]
        isouth = iseg + int(below_south[-1]) if below_south.size > 0 else iseg

        sup = profile[peak_global_idx:jseg+1]
        below_north = np.where(sup < south_thresh)[0]
        inorth = peak_global_idx + int(below_north[0]) if below_north.size > 0 else jseg

        if inorth <= isouth:
            continue

        x = dists[isouth:inorth+1]
        y = profile[isouth:inorth+1]
        # lat/lon arrays aligned with x/y for correct interpolation
        lat_for_x = lats[isouth:inorth+1]
        lon_for_x = lons[isouth:inorth+1]
        total = float(np.trapezoid(y, x))
        if total <= 0.0:
            continue

        cum = np.concatenate(([0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))))
        half = total * 0.5
        if cum[-1] < half:
            continue

        k = int(np.searchsorted(cum, half))
        if k <= 0:
            lat_mid = float(lat_for_x[0])
            lon_mid = float(lon_for_x[0])
        elif k >= len(cum):
            lat_mid = float(lat_for_x[-1])
            lon_mid = float(lon_for_x[-1])
        else:
            denom = (cum[k] - cum[k-1])
            t = (half - cum[k-1]) / denom if denom != 0.0 else 0.0
            lat_mid = float(lat_for_x[k-1] + t * (lat_for_x[k] - lat_for_x[k-1]))
            lon_mid = float(lon_for_x[k-1] + t * (lon_for_x[k] - lon_for_x[k-1]))

        results.append({
            'lat_mid': lat_mid,
            'lon_mid': lon_mid,
            'isouth': int(isouth),
            'inorth': int(inorth),
            'peak_lat': float(lats[peak_global_idx]),
            'peak_val': float(peak_val),
            'width_km': float(width_km)
        })

    return results


def detect_fronts_from_dataset(ds: xr.Dataset,
                               eke_var: str = 'eke',
                               lon_list: Optional[List[float]] = None,
                               lon_window: int = 0,
                               threshold: float = 0.03,
                               min_width_km: float = 100.0,
                               peak_fraction: float = 0.25,
                               sample_every: int = 1,
                               use_numba: bool = False,
                               n_jobs: int = 1) -> Dict[float, List[Dict]]:
    """Scan a dataset and detect front envelopes for requested longitudes.

    Returns a dict mapping evaluation longitude -> list of envelope dicts.

    Parameters
    ----------
    n_jobs : int
        Number of threads for parallel longitude processing.
        -1 uses all available CPUs. 1 (default) disables threading.
        Threads are safe here because both numba and numpy release the GIL.
    """
    if eke_var not in ds:
        raise ValueError(f"Dataset does not contain variable '{eke_var}'")

    if 'time' in ds.dims:
        mean_2d = ds[eke_var].mean(dim='time', skipna=True)
    else:
        mean_2d = ds[eke_var]

    # Materialise once to avoid repeated lazy evaluation inside the loop
    eke_vals = np.asarray(mean_2d.values)

    lats = mean_2d['latitude'].values
    lons = mean_2d['longitude'].values

    if lats.size == 0 or lons.size == 0:
        return {}

    # Pre-orient south→north so every worker gets a consistent array
    if lats[0] > lats[-1]:
        eke_vals = eke_vals[::-1, :]
        lats = lats[::-1]

    dists = meridional_distances_km_from_lats(lats)

    if lon_list is None:
        eval_idxs = list(range(0, len(lons), sample_every))
        lon_list_eval = [float(lons[i]) for i in eval_idxs]
    else:
        lon_list_eval = lon_list

    def _process_one(lon_val: float):
        j = int(np.argmin(np.abs(lons - lon_val)))
        i0 = max(0, j - lon_window)
        i1 = min(len(lons) - 1, j + lon_window)
        lon_mid = float(np.mean(lons[i0:i1+1]))

        profile_2d = eke_vals[:, i0:i1+1]
        if profile_2d.size == 0:
            return lon_mid, []
        if profile_2d.ndim == 2:
            finite = np.isfinite(profile_2d)
            counts = finite.sum(axis=1)
            summed = np.nansum(profile_2d, axis=1)
            profile = np.full(profile_2d.shape[0], np.nan, dtype=float)
            np.divide(summed, counts, out=profile, where=counts > 0)
        else:
            profile = profile_2d

        if not np.any(np.isfinite(profile)):
            return lon_mid, []

        lon_arr = np.full_like(lats, lon_mid, dtype=float)

        envelopes = detect_envelopes_on_profile(
            profile, lats, lon_arr,
            dists=dists,
            threshold=threshold,
            min_width_km=min_width_km,
            peak_fraction=peak_fraction,
            use_numba=use_numba,
        )
        return lon_mid, envelopes

    import os as _os
    workers = _os.cpu_count() if n_jobs == -1 else n_jobs

    results: Dict[float, List[Dict]] = {}
    if workers == 1:
        for lon_val in lon_list_eval:
            lon_mid, envs = _process_one(lon_val)
            results[lon_mid] = envs
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_process_one, lon_val): lon_val for lon_val in lon_list_eval}
            for fut in as_completed(futures):
                lon_mid, envs = fut.result()
                results[lon_mid] = envs

    return results
