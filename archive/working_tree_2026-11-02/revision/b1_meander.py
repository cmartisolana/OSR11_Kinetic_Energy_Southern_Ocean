"""B1: standing-meander pathway metrics (pilot: Pacific-Antarctic Ridge).

Pathway = ADT contour relative to the box-mean ADT (removes spatially uniform
sea-level rise, cf. Gille 2014). Levels are fixed from the CLIMATOLOGY before
any trend is computed: the relative level whose contour carries the highest
mean geostrophic speed in the box (jet core), and +/- 0.10 m around it.

Per month and level:
  * the longest connected contour piece that spans the box west->east
  * projected to km (point-wise equirectangular, x = R cos(lat) dlon, y = R dlat)
  * resampled every DS km along the arc and Gaussian-smoothed along the arc
    with scale SMOOTH km (>= effective resolution; 50/100/150 tested)
  * metrics: length, sinuosity (length / end-to-end distance), meridional
    amplitude (max-min latitude, km), RMS curvature, integral |kappa| ds
Months with no spanning contour are flagged (NaN), never filled.
"""
from __future__ import annotations

import numpy as np
import contourpy
from scipy.ndimage import gaussian_filter1d

R_KM = 6371.0


def box_relative(adt, lat, lon, box):
    la = (lat >= box["lat0"]) & (lat <= box["lat1"])
    lo = (lon >= box["lon0"]) & (lon <= box["lon1"])
    sub = adt[np.ix_(la, lo)]
    return sub - np.nanmean(sub), lat[la], lon[lo]


def spanning_contour(field, lat, lon, level):
    """Longest contour piece at `level`; require it to span >= 90 % of the box width."""
    gen = contourpy.contour_generator(lon, lat, np.ma.masked_invalid(field), name="serial")
    lines = gen.lines(level)
    if not lines:
        return None
    best = max(lines, key=lambda L: np.ptp(L[:, 0]))
    if np.ptp(best[:, 0]) < 0.9 * (lon.max() - lon.min()):
        return None
    # orient west -> east
    if best[0, 0] > best[-1, 0]:
        best = best[::-1]
    return best


def path_metrics(line, ds=10.0, smooth=100.0):
    lon, lat = line[:, 0], line[:, 1]
    lat0 = np.radians(np.mean(lat))
    x = R_KM * np.cos(np.radians(lat)) * np.radians(lon - lon[0])
    y = R_KM * np.radians(lat - lat[0])
    seg = np.hypot(np.diff(x), np.diff(y))
    s = np.r_[0, np.cumsum(seg)]
    keep = np.r_[True, seg > 0]
    s, x, y = s[keep], x[keep], y[keep]
    if s[-1] < 5 * ds:
        return None
    si = np.arange(0, s[-1], ds)
    xi = np.interp(si, s, x); yi = np.interp(si, s, y)
    sg = smooth / ds
    xs = gaussian_filter1d(xi, sg, mode="nearest"); ys = gaussian_filter1d(yi, sg, mode="nearest")
    dx = np.gradient(xs, ds); dy = np.gradient(ys, ds)
    ddx = np.gradient(dx, ds); ddy = np.gradient(dy, ds)
    kappa = (dx * ddy - dy * ddx) / np.power(dx * dx + dy * dy, 1.5)
    L = float(np.sum(np.hypot(np.diff(xs), np.diff(ys))))
    chord = float(np.hypot(xs[-1] - xs[0], ys[-1] - ys[0]))
    return dict(length_km=L, sinuosity=L / chord, amp_km=float(np.ptp(ys)),
                rms_kappa=float(np.sqrt(np.mean(kappa ** 2))),
                int_abs_kappa=float(np.sum(np.abs(kappa)) * ds))


def core_level(adt_mean, u_mean, v_mean, lat, lon, box, levels=np.arange(-0.6, 0.61, 0.02)):
    """Relative ADT level with the highest mean speed along its contour (climatology)."""
    rel, la, lo = box_relative(adt_mean, lat, lon, box)
    spd, _, _ = box_relative(np.hypot(u_mean, v_mean), lat, lon, box)
    spd = spd + np.nanmean(np.hypot(u_mean, v_mean)[np.ix_((lat >= box["lat0"]) & (lat <= box["lat1"]),
                                                            (lon >= box["lon0"]) & (lon <= box["lon1"]))])
    best, bestv = None, -np.inf
    for lev in levels:
        line = spanning_contour(rel, la, lo, lev)
        if line is None:
            continue
        ii = np.clip(np.searchsorted(la, line[:, 1]), 0, la.size - 1)
        jj = np.clip(np.searchsorted(lo, line[:, 0]), 0, lo.size - 1)
        v = np.nanmean(spd[ii, jj])
        if v > bestv:
            best, bestv = lev, v
    return best, bestv
