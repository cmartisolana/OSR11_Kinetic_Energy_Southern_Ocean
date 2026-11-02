"""dk6 (19 Oct 2026): five half-power tests in Drake Passage (70-55W, envelopes with climatological half-power south of 52S,
in the ACC band), Products 1 (allsat) and 4 (twosat). Baseline half-power: 200 cm2 s-2, 25 %, 100 km, monthly KE of the
monthly-mean velocity; monthly envelopes matched to the 2001-2025 mean envelope at each longitude (as a2_run.py).
 T1  passage-mean half-power latitude (coverage >= 50 % and >= 90 %; with and without monthly ice mask)
 T2  streamline coordinates: relative ADT at the half-power latitude (psi_hp) and latitude of the climatological psi_hp streamline
 T3  EKE-weighted median latitude inside the monthly KE envelope, relative to the KE half-power latitude
 T4  fraction of KE and EKE north of the climatological half-power latitude inside the climatological envelope (+ synthetic sensitivities)
 T5  annual passage half-power vs wind stress (lags 0-3 yr): sensitivity and wind-expected shift

Purpose : half-power tests in Drake Passage (70-55 W), both products, with and without the monthly ice mask (SI S3).
Inputs  : outputs/drake/dk6/fields_<product>.npz (03l)
Outputs : outputs/drake/dk6/dk6_trends_<product>_ice<0|1>.csv, dk6_series_*.csv, dk6_tracks_*.csv (distributed)
Run     : for a in "allsat 0" "allsat 1" "twosat 0" "twosat 1"; do python scripts/03_drake/03m_section_halfpower_tests.py $a; done
          (replaces dk6_all.sh, which cd'ed into an absolute session path)
Seeds   : rstats.trend_bcz default seed (0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, numpy as np, pandas as pd
from osr11.drake_halfpower import *  # noqa: F401,F403
prod, icem = sys.argv[1], bool(int(sys.argv[2])); tag = f"{prod}_ice{int(icem)}"
D = load(prod, icem); lat, lon = D["lat"], D["lon"]; KM = 111.195
dists = F.meridional_dists(lat)
C, TB = detect_all(D["ke"], lat)
kec = np.nanmean(D["ke"], 0); ekec = np.nanmean(D["eke"], 0); adtc = np.nanmean(D["adt"], 0)
cc, ctb = F.detect_field(np.ascontiguousarray(kec), lat, dists, P["thr"], P["min_w"], P["frac"], MAXE)
meta, arr = match(C, TB, cc, ctb, lat)
M = pd.DataFrame(meta, columns=["j", "c", "lo", "hi", "hpc", "ls", "ln"]); M["j"] = M.j.astype(int); M["lon"] = lon[M.j]
ia = np.clip(np.searchsorted(lat, M.hpc), 0, lat.size - 1); M["adt_hpc"] = adtc[ia, M.j]
M["cov"] = np.isfinite(arr[:, :, FI["lat_hp"]]).mean(0)
sel0 = (M.lon >= -70) & (M.lon < -55) & (M.hpc < -52) & (M.adt_hpc >= -0.6) & (M.adt_hpc <= 0.2)
# relative ADT: monthly ADT minus its mean over the ACC band (mean ADT -0.6..0.2) of 70-55W
band = (adtc >= -0.6) & (adtc <= 0.2) & ((lon >= -70) & (lon < -55))[None, :]
adt_rel = D["adt"] - np.array([np.nanmean(a[band]) for a in D["adt"]])[:, None, None]
adtc_rel = np.nanmean(adt_rel, 0)

def wmedian(prof, lo_, hi_):
    """distance-weighted median latitude of prof between lo_ and hi_ (trapezoid, as the half-power)"""
    m = (lat >= lo_) & (lat <= hi_) & np.isfinite(prof)
    if m.sum() < 3: return np.nan
    x, y, dd = lat[m], np.clip(prof[m], 0, None), dists[m]
    a = 0.5 * (y[1:] + y[:-1]) * np.diff(dd); cum = np.r_[0, np.cumsum(a)]
    if cum[-1] <= 0: return np.nan
    return float(np.interp(0.5 * cum[-1], cum, x))

def nfrac(prof, lo_, mid, hi_):
    m = (lat >= lo_) & (lat <= hi_) & np.isfinite(prof)
    if m.sum() < 3 or not np.isfinite(prof[(lat >= lo_) & (lat <= hi_)]).all(): return np.nan
    x, y = lat[m], np.clip(prof[m], 0, None)
    tot = np.trapezoid(y, x); n = np.trapezoid(np.interp(np.r_[mid, x[x > mid]], x, y), np.r_[mid, x[x > mid]])
    return n / tot if tot > 0 else np.nan

rows, series = [], {}
for cmin in (0.5, 0.9):
    S_ = M[sel0 & (M["cov"] >= cmin)].index.values
    hp = arr[:, S_, FI["lat_hp"]]; ls_m, ln_m = arr[:, S_, FI["lat_s"]], arr[:, S_, FI["lat_n"]]
    nt, nq = hp.shape
    psi = np.full_like(hp, np.nan); phipsi = np.full_like(hp, np.nan); phie = np.full_like(hp, np.nan)
    fk = np.full_like(hp, np.nan); fe = np.full_like(hp, np.nan); G = np.full(nq, np.nan)
    for q, iq in enumerate(S_):
        j = M.j[iq]; hpc, ls_c, ln_c = M.hpc[iq], M.ls[iq], M.ln[iq]
        # climatological relative-ADT streamline at the half-power latitude, and the mean meridional ADT gradient (m per km)
        pc = adtc_rel[:, j]; psic = np.interp(hpc, lat, pc); w = (lat >= ls_c) & (lat <= ln_c) & np.isfinite(pc)
        G[q] = np.polyfit(lat[w] * KM, pc[w], 1)[0] if w.sum() > 3 else np.nan
        for k in range(nt):
            pk = adt_rel[k, :, j]
            if np.isfinite(hp[k, q]): psi[k, q] = np.interp(hp[k, q], lat, pk)
            # latitude of the climatological streamline: crossing nearest to the climatological half-power latitude
            ok = np.isfinite(pk); s = np.sign(pk[ok] - psic); idx = np.where(np.diff(s) != 0)[0]
            if idx.size:
                la = lat[ok]; pp = pk[ok]; cr = la[idx] + (psic - pp[idx]) * (la[idx + 1] - la[idx]) / (pp[idx + 1] - pp[idx])
                phipsi[k, q] = cr[np.argmin(abs(cr - hpc))]
            if np.isfinite(hp[k, q]): phie[k, q] = wmedian(D["eke"][k, :, j], ls_m[k, q], ln_m[k, q])
            fk[k, q] = nfrac(D["ke"][k, :, j], ls_c, hpc, ln_c); fe[k, q] = nfrac(D["eke"][k, :, j], ls_c, hpc, ln_c)
    def ser(a, scale=1.0):
        return np.nanmean((a - np.nanmean(a, 0, keepdims=True)) * scale, 1)
    Gm = np.nanmean(G)
    ss = {"T1 half-power latitude (km, +N)": ser(hp, KM),
          "T2 psi at half-power (relative ADT, cm)": ser(psi, 100.0),
          "T2 psi at half-power, km-equivalent (+N across streamlines)": ser(psi / G[None, :], 1.0),
          "T2 latitude of the climatological psi_hp streamline (km, +N)": ser(phipsi, KM),
          "T3 EKE median minus KE half-power (km, +N)": ser(phie - hp, KM),
          "T3 EKE-weighted median latitude (km, +N)": ser(phie, KM),
          "T4 KE fraction north of clim half-power (pp)": ser(fk, 100.0),
          "T4 EKE fraction north of clim half-power (pp)": ser(fe, 100.0)}
    for name, y in ss.items():
        r = tr(y)
        rows.append(dict(product=prod, icemask=icem, coverage_min=cmin, n_lon=nq, series=name, mean_valid=float(np.isfinite(hp).mean()),
                         std=float(np.nanstd(y)), trend_dec=r["slope"] * 10, lo=r["ci_lo"] * 10, hi=r["ci_hi"] * 10, p=r["p"]))
        series[(cmin, name)] = y
    if cmin == 0.5:
        # T2 decomposition: how much of the monthly half-power variability moves with the streamlines
        y1, yb = series[(0.5, "T1 half-power latitude (km, +N)")], series[(0.5, "T2 latitude of the climatological psi_hp streamline (km, +N)")]
        ok = np.isfinite(y1) & np.isfinite(yb); rr = np.corrcoef(y1[ok], yb[ok])[0, 1]
        rows.append(dict(product=prod, icemask=icem, coverage_min=cmin, n_lon=nq, series="T2 r(half-power latitude, streamline latitude), monthly", trend_dec=rr))
        rows.append(dict(product=prod, icemask=icem, coverage_min=cmin, n_lon=nq, series="T2 mean |dADT/dy| at the envelope (cm per 100 km)", trend_dec=Gm * 1e4))
        # T4 synthetic sensitivities on the climatological profiles (no noise)
        for X, nm in ((kec, "KE"), (ekec, "EKE")):
            dn, ds_ = [], []
            for iq in S_:
                j = M.j[iq]; hpc, ls_c, ln_c = M.hpc[iq], M.ls[iq], M.ln[iq]; p0 = X[:, j]
                f0 = nfrac(p0, ls_c, hpc, ln_c)
                p1 = np.where(lat > hpc, 1.2 * p0, p0); dn.append(nfrac(p1, ls_c, hpc, ln_c) - f0)
                p2 = np.interp(lat - 25 / KM, lat, p0); ds_.append(nfrac(p2, ls_c, hpc, ln_c) - f0)
            rows.append(dict(product=prod, icemask=icem, coverage_min=cmin, n_lon=nq, series=f"T4 synthetic: 20 % north-flank {nm} intensification (pp)", trend_dec=100 * np.nanmean(dn)))
            rows.append(dict(product=prod, icemask=icem, coverage_min=cmin, n_lon=nq, series=f"T4 synthetic: 25 km northward shift, {nm} (pp)", trend_dec=100 * np.nanmean(ds_)))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/dk6_trends_{tag}.csv", index=False)
pd.DataFrame({f"{c}|{k}": v for (c, k), v in series.items()}).assign(x=T).to_csv(f"{OUT}/dk6_series_{tag}.csv", index=False)
M[sel0].to_csv(f"{OUT}/dk6_tracks_{tag}.csv", index=False)
pd.set_option("display.width", 250); print(R.round(3).to_string())
