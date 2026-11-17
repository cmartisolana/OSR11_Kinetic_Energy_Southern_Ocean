"""Occurrence of the monthly half-power points for SI Fig. S6a,b with the RELATIVE threshold (20 % of each longitude's
climatological maximum mean KE, the detection used throughout the paper; 10 Nov 2026).

Input : <outputs>/x7/x7_points_all.npy (month index 0-299 from Jan 2001, lon, lat of every monthly half-power point in the ACC
        band; producer scripts/02_circumpolar/02p_pathway_products.py).
Method: exactly review_calc/c06_overlap_null.py (Drake) and c07_overlap_null_circ.py (circumpolar): counts per whole year on
        the grids below, density = % of months per grid longitude of the cell column (1/8 deg grid: 4 longitudes per 0.5 deg
        Drake column, 8 per 1 deg circumpolar column), Gaussian smoothing sigma = 1 cell, 10 % level, cos(lat)-weighted
        intersection over union (IoU) of the two smoothed 10 % regions.  Periods 2001-2017 / 2018-2025 (and 2001-2015 /
        2016-2025 for reference).
Domain: the fixed-threshold Drake points (c06a) exclude the Malvinas Current (-62 <= lon < -56 and lat > -52) and keep
        ADT in [-0.6, 0.2] m; the x7 points already satisfy the ADT window (ACC-band mask) and the same Malvinas exclusion is
        applied here to both domains (467 Drake-box points, none of them inside a 10 % region: the IoU is unchanged).  The
        histogram grids bound the domains (Drake: 86-34 W, 64-46 S; circumpolar 65-35 S).
Checks: (1) the fixed-threshold overlaps from review_calc/c06a_drake_points.npy, (2) that the x3 density files are the raw
        (unsmoothed) fixed-threshold grids, (3) that obs18 of outputs/x3/x3_overlap_2018split.json is already the
        relative-threshold value, (4) the random 17/8 whole-year-split null (2000 draws, seed 20261105) against that json.
Output: <dir>/s6_occurrence_relthr.npz (drake_* and circ_*: lon_edges, lat_edges, all/early/late raw density and
        early_smooth/late_smooth, in % of months) and s6_occurrence_relthr.json.
Usage : python figures/s6_occurrence_relthr.py [outputs dir] [output dir]  (defaults: $OSR11_OUTPUTS, this folder)

Origin  : figures/s6_occurrence_relthr.py of the 10 Nov 2026 Supplement; moved here on 17 Nov 2026 with the paths routed through osr11.paths (no change of method).
Outputs : outputs/s6/s6_occurrence_relthr.npz, s6_occurrence_relthr.json (distributed; read by figures/figS6_occurrence_adt.py)
Needs   : outputs/x7/x7_points_all.npy (written by 02p, not distributed)
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import json, sys
import numpy as np
from scipy.ndimage import gaussian_filter

from osr11 import paths as PA  # noqa: E402
OUT = PA.OUT; DST = PA.out("s6")
C06A = OUT / "c06" / "c06a_drake_points.npy"   # fixed-threshold Drake points of the 4 Nov review (not distributed); only used by the checks
GRID = {"drake": (np.arange(-86, -33.99, 0.5), np.arange(-64, -45.99, 0.25), 4),      # c06 lon edges, lat edges, ncol
        "circ": (np.arange(-180, 180.01, 1.0), np.arange(-65, -35.01, 0.5), 8)}      # c07
LEV, SIG = 10.0, 1.0
SPLITS = {"split2018": (np.arange(0, 17), np.arange(17, 25)), "split2016": (np.arange(0, 15), np.arange(15, 25))}


def malvinas(p):
    return (p[:, 1] >= -62) & (p[:, 1] < -56) & (p[:, 2] > -52)


class Occ:
    def __init__(self, pts, dom):
        self.le, self.ae, self.ncol = GRID[dom]; k, lon, lat = pts[:, 0].astype(int), pts[:, 1], pts[:, 2]; yr = k // 12
        self.H = np.stack([np.histogram2d(lat[yr == y], lon[yr == y], bins=[self.ae, self.le])[0] for y in range(25)])
        self.w = np.cos(np.radians(0.5 * (self.ae[1:] + self.ae[:-1])))[:, None]

    def dens(self, years):
        return 100 * self.H[years].sum(0) / (12 * len(years) * self.ncol)

    def iou(self, ye, yl, lev=LEV):
        a, b = gaussian_filter(self.dens(ye), SIG) >= lev, gaussian_filter(self.dens(yl), SIG) >= lev
        return float((self.w * (a & b)).sum() / (self.w * (a | b)).sum())


def iou_grid(e, l, w, lev=LEV):
    a, b = gaussian_filter(e, SIG) >= lev, gaussian_filter(l, SIG) >= lev
    return float((w * (a & b)).sum() / (w * (a | b)).sum())


X3 = json.load(open(OUT / "x3" / "x3_overlap_2018split.json"))
x7 = np.load(OUT / "x7" / "x7_points_all.npy"); res = {"method": __doc__.split("Output:")[0].strip(), "relative_threshold": {}, "checks": {}}
npz = {}
for dom, x3key in (("drake", "drake"), ("circ", "circ_P1")):
    O = Occ(x7[~malvinas(x7)], dom); O0 = Occ(x7, dom)                  # with / without the Malvinas exclusion
    e, l = SPLITS["split2018"]
    r = {s: round(O.iou(*SPLITS[s]), 4) for s in SPLITS}
    r.update({f"split2018_level{v}": round(O.iou(e, l, v), 4) for v in (5, 15)})
    r["split2018_without_malvinas_exclusion"] = round(O0.iou(e, l), 4)
    rng = np.random.default_rng(20261105); null = np.empty(2000)
    for i in range(2000):
        yl = np.sort(rng.choice(25, 8, replace=False)); null[i] = O.iou(np.setdiff1d(np.arange(25), yl), yl)
    r["null_17_8_median"] = round(float(np.median(null)), 4); r["null_q05_q95"] = [round(float(np.percentile(null, q)), 4) for q in (5, 95)]
    r["P_null_le_obs"] = round(float(np.mean(null <= r["split2018"])), 4)
    r["n_points_in_grid"] = int(O.H.sum()); r["x3_json_obs18"] = round(X3[x3key]["obs18"], 4); r["x3_json_null_med"] = round(X3[x3key]["null_med"], 4)
    res["relative_threshold"][dom] = r
    d = {"lon_edges": O.le, "lat_edges": O.ae, "all": O.dens(np.arange(25)), "early": O.dens(e), "late": O.dens(l)}
    d["early_smooth"], d["late_smooth"] = gaussian_filter(d["early"], SIG), gaussian_filter(d["late"], SIG)
    npz.update({f"{dom}_{k}": v for k, v in d.items()})
    print(dom, json.dumps(r))
# checks against the fixed-threshold (0.02 m2 s-2) products -------------------------------------------------------------
ch = res["checks"]
if not C06A.exists():
    print("checks against the fixed-threshold points skipped (", C06A, "not present)")
    np.savez_compressed(DST / "s6_occurrence_relthr.npz", **npz); json.dump(res, open(DST / "s6_occurrence_relthr.json", "w"), indent=1)
    print("saved", DST / "s6_occurrence_relthr.npz"); sys.exit(0)
c06a = Occ(np.load(C06A), "drake"); XD = np.load(OUT / "x3" / "x3_density_drake.npz")
XC = np.load(OUT / "x3" / "x3_density_circ_P1.npz"); e, l = SPLITS["split2018"]
ch["fixed_drake_c06a_split2016"] = round(c06a.iou(*SPLITS["split2016"]), 4)          # c06_log: 0.829 (main Fig. 2b before relthr)
ch["fixed_drake_c06a_split2018"] = round(c06a.iou(e, l), 4)
ch["x3_drake_density_equals_raw_c06a"] = bool(all(np.allclose(XD[k], c06a.dens(y)) for k, y in (("all", np.arange(25)), ("early", e), ("late", l))))
ch["fixed_circ_from_x3_grids_split2018"] = round(iou_grid(XC["early"], XC["late"], np.cos(np.radians(0.5 * (XC["lat_edges"][1:] + XC["lat_edges"][:-1])))[:, None]), 4)
ch["x3_circ_all_equals_17_8_weighted_early_late"] = bool(np.allclose(XC["all"], (17 * XC["early"] + 8 * XC["late"]) / 25))
ch["note"] = ("x3_density_*.npz are the RAW (unsmoothed) fixed-threshold densities (Drake identical to c06a; circumpolar 'all' identical to "
              "drake/dk14 'all'); obs18 of x3_overlap_2018split.json for drake and circ_P1 is the relative-threshold value (reproduced here "
              "from x7 without the Malvinas exclusion), as its own 'note' field says; x3/README.txt is out of date on this point. The fixed-threshold "
              "circumpolar points exist here only for months 0-29 (c07a_points_000.npy), so its split-2018 IoU is taken from the x3 grids.")
print("checks", json.dumps(ch, indent=1))
np.savez_compressed(DST / "s6_occurrence_relthr.npz", **npz)
json.dump(res, open(DST / "s6_occurrence_relthr.json", "w"), indent=1)
print("saved", DST / "s6_occurrence_relthr.npz", DST / "s6_occurrence_relthr.json")
