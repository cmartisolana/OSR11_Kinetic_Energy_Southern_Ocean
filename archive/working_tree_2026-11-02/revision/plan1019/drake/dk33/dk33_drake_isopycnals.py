"""dk33: isopycnals along the LM Gould line in Drake Passage from the RG Argo climatology (dk32 Drake-box sigma0), 2004-2025.
Annual sigma0 (TEOS-10) at 1 deg interpolated bilinearly onto the mean LM Gould track (dk24 section, 121 points).
Outputs: mean sections 2004-2014 and 2015-2025; per-year pressure of sigma0 = 27.3, 27.5, 27.7 along the track; the
cross-passage isopycnal slope (least squares of pressure on latitude over the points where the surface exists in all years;
dbar per degree, positive = deeper to the north) and its trend (% per decade; 95 % with n_eff from lag-1 autocorrelation).
Usage: python dk33_drake_isopycnals.py <dk32 out> <dk24 section npz> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from scipy.interpolate import RegularGridInterpolator as RGI
D32, SEC, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
Z = np.load(D32 / "dk32_drake_sigma0.npz"); yrs, P, la, lo, S0 = Z["years"], Z["P"], Z["lat"], Z["lon"], Z["s0"]
T = np.load(SEC); tlon, tlat = np.mod(T["lon"], 360), T["lat"]
sec = np.array([[RGI((la, lo), S0[y, k], bounds_error=False, fill_value=np.nan)(np.c_[tlat, tlon]) for k in range(P.size)] for y in range(yrs.size)])
SIG = [27.3, 27.5, 27.7]
def iso(s0):        # s0: (P, n) -> pressure of each sigma0 surface
    out = {}
    for s in SIG:
        below = s0 >= s; first = np.argmax(below, 0); has = below.any(0) & (first > 0) & ~below[0]
        k1 = np.clip(first, 1, P.size - 1); k0 = k1 - 1
        g0 = np.take_along_axis(s0, k0[None], 0)[0]; g1 = np.take_along_axis(s0, k1[None], 0)[0]
        p = P[k0] + (s - g0) / (g1 - g0) * (P[k1] - P[k0]); out[s] = np.where(has & np.isfinite(p), p, np.nan)
    return out
IP = {s: np.array([iso(sec[y])[s] for y in range(yrs.size)]) for s in SIG}
t = yrs - yrs.mean()
def trend(y):
    b, a = np.polyfit(t, y, 1); r = y - a - b * t; r1 = np.corrcoef(r[1:], r[:-1])[0, 1]; ne = max(3, len(y) * (1 - r1) / (1 + r1))
    se = np.sqrt((r ** 2).sum() / (ne - 2) / (t ** 2).sum()); return 10 * b, 10 * stats.t.ppf(0.975, ne - 2) * se, 2 * stats.t.sf(abs(b / se), ne - 2)
rows, SL = [], {}
for s in SIG:
    ok = np.isfinite(IP[s]).all(0); x = tlat[ok] - tlat[ok].mean()
    sl = (IP[s][:, ok] - IP[s][:, ok].mean(1, keepdims=True)) @ x / (x @ x)   # dbar per degree latitude: positive = deeper to the north
    SL[s] = sl; b, ci, p = trend(sl); bp, cip, pp = trend(100 * (sl / sl.mean() - 1))
    zN = np.nanmean(IP[s][:, ok][:, tlat[ok] > np.median(tlat[ok])], 1); zS = np.nanmean(IP[s][:, ok][:, tlat[ok] <= np.median(tlat[ok])], 1)
    bz, ciz, pz = trend(zN - zS)
    rows.append(dict(sigma0=s, npts=int(ok.sum()), lat_range=f"{tlat[ok].min():.1f} to {tlat[ok].max():.1f}", slope_mean=sl.mean(), slope_trend=b, slope_ci=ci, slope_p=p,
                     slope_trend_pct=bp, slope_ci_pct=cip, dz_NS_mean=(zN - zS).mean(), dz_trend=bz, dz_ci=ciz, dz_p=pz,
                     depth_trend_N=trend(zN)[0], depth_trend_S=trend(zS)[0]))
D = pd.DataFrame(rows); D.to_csv(OUT / "dk33_drake_isopycnals.csv", index=False); pd.set_option("display.width", 220); print(D.round(2).to_string(index=False))
e, l_ = yrs <= 2014, yrs >= 2015
np.savez(OUT / "dk33_drake_section.npz", years=yrs, P=P, tlat=tlat, tlon=tlon, sec_early=np.nanmean(sec[e], 0), sec_late=np.nanmean(sec[l_], 0), sec_mean=np.nanmean(sec, 0),
         **{f"p{s}": v for s, v in IP.items()}, **{f"sl{s}": v for s, v in SL.items()})
