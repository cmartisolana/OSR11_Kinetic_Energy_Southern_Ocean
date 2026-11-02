"""DK25: trend of the altimetric cross-track surface geostrophic velocity across Drake Passage by 0.25 deg latitude
along the mean LM Gould track (dk24 section fields, Product 1, same bins as the dk3 mean-profile comparison).
Points with ice fraction > 0.10 masked per month. 2001-2025: monthly test where the bin series is complete, otherwise
annual means (>= 6 months); Oct 2005-Apr 2019: annual means (>= 3 months), as dk5. Also mean profiles for 2001-2012 and
2013-2025. Usage: python dk25_profile_trends.py <revision dir> <dk24 dir> <out dir>

Purpose : trend of the altimetric cross-track surface velocity by 0.25-deg latitude along the track (Fig. S4).
Inputs  : outputs/drake/dk24/dk24_section_fields_allsat.npz (03d)
Outputs : outputs/drake/dk25/dk25_profile_trends.csv
Run     : python scripts/03_drake/03e_profile_trends.py          (no raw data needed)
Seeds   : rstats.trend_bcz default seed (0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd
from osr11 import rstats as R, paths as PA
DK24, OUT = PA.dk(24), PA.dk(25)
Z = np.load(DK24 / "dk24_section_fields_allsat.npz"); lo, la = Z["lon"], Z["lat"]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
un = (Z["u"] * nx + Z["v"] * ny) * 100; un[Z["ice"] > 0.10] = np.nan
T = pd.DataFrame({"year": np.repeat(range(2001, 2026), 12), "month": np.tile(range(1, 13), 25)}); T["x"] = T.year + (T.month - .5) / 12
ov = ((T.x >= 2005 + 9 / 12) & (T.x <= 2019 + 4 / 12)).values
def annual(years, vals, minn):
    df = pd.DataFrame(dict(year=years, v=vals)).dropna(); g = df.groupby("year").v; a = g.mean()[g.count() >= minn]
    return a.index.values.astype(float), a.values
rows = []
for c in np.arange(-63.5, -54.99, 0.25):
    s = np.abs(la - c) <= 0.125
    if not s.any(): continue
    y = np.nanmean(un[:, s], 1)
    out = dict(latb=c, npts=int(s.sum()), mean=np.nanmean(y), mean_early=np.nanmean(y[T.year.values <= 2012]),
               mean_late=np.nanmean(y[T.year.values >= 2013]), months=int(np.isfinite(y).sum()))
    if np.isfinite(y).all(): o = R.trend_bcz(T.x.values, y, month=T.month.values)
    else:
        xa, ya = annual(T.year.values, y, 6); o = R.trend_bcz(xa, ya, monthly=False) if len(xa) >= 8 else None
    out.update({k + "_0125": (np.nan if o is None else 10 * o[k] if k != "p" else o[k]) for k in ("slope", "ci_lo", "ci_hi", "p")})
    xa, ya = annual(T.year[ov].values, y[ov], 3); o = R.trend_bcz(xa, ya, monthly=False) if len(xa) >= 8 else None
    out.update({k + "_0519": (np.nan if o is None else 10 * o[k] if k != "p" else o[k]) for k in ("slope", "ci_lo", "ci_hi", "p")})
    rows.append(out)
D = pd.DataFrame(rows); D.to_csv(OUT / "dk25_profile_trends.csv", index=False); print(D.round(2).to_string())
