"""H1: local half-power meander geometry with synthetic validation (plan section 4).

SPECIFICATION (written 2026-10-19 before any geometry trend or relation was computed)
--------------------------------------------------------------------------------------
Sites and windows are taken unchanged from earlier choices made on time-mean fields only:
  SEIR150 meander 140-150E   (b1_meander_run.py; upstream 125-138E, downstream 150-158E)
  KERG    meander 76-84E     (upstream 70-77E, downstream 84-94E)
  PAR     meander 135.13-109.38W (Liu et al. 2026 sections; upstream 150-135.13W, downstream 109.38-80W)
  DRAKE   passage 70-55W     (plan: "begin in Drake Passage"; upstream 100-80W, downstream 55-35W = A5d boxes)
Detection: fronts.detect_field on monthly KE = 0.5|U_m|^2 (all-sat), D1 baseline (200 cm2 s-2, 100 km, 25 %).
Climatological envelopes from the 2001-2025 mean KE; a monthly envelope is matched to a climatological one at
the same longitude if its half-power latitude lies inside the climatological segment (strongest if several).
Branch: start at the window longitude whose climatological envelope has the largest peak KE; step west and east
choosing the climatological envelope whose half-power latitude is within 0.25 deg of the last accepted one
(gaps of up to 4 columns = 0.5 deg bridged). A branch longitude is usable if its matched coverage >= 0.9.
GATE G1 (coverage): usable branch longitudes >= 90 % of window columns. Sites failing G1 are reported, not analysed.
Monthly path: half-power latitude at usable branch longitudes; months with > 10 % of them missing -> NaN;
otherwise linear interpolation across missing longitudes. Metrics (b1_meander.path_metrics, metric coordinates):
mean latitude (displacement, km north of the mean), meridional amplitude (km), sinuosity (length / end-to-end
distance), RMS curvature and integral |kappa| ds; Gaussian smoothing 150 km primary, 100 and 200 km sensitivity.
GATE G2 (synthetic, 40 realisations per scenario, 2001-2025 monthly): synthetic KE = climatology deformed by a
linear 25-yr ramp + observed monthly KE anomalies of the same site resampled in 12-month blocks. Scenarios:
  N  no change;  T  uniform northward translation 25 km;  F  flexing: displacement 25 km x sin(2 pi s/W) along the
  window (W = window width), i.e. meander amplitude change;  U  uniform intensification x1.20;
  A  asymmetric intensification: x1.20 on the northern flank tapering to x1.00 on the southern flank (tanh over 50 km).
Pass if (i) T recovers >= 50 % of the imposed displacement, (ii) |spurious displacement| under A and U is < 50 % of
the displacement recovered under T, (iii) amplitude response to F exceeds 2 sd of its N spread.
Metrics failing (ii) or (iii) are reported as not interpretable at that site.
ANALYSIS (passing sites only): trend_bcz of each metric (150 km). Relations on monthly anomalies (deseasonalised,
Theil-Sen detrended, 2-month running mean, as b1_covar.py): upstream jet-core speed -> geometry (lags 0-3),
geometry -> meander-section EKE and -> downstream EKE (lags 0-3), downstream jet-core speed vs geometry (lags -3..3);
significance from 1000 phase-randomised surrogates of the geometry series, maximum over lags; BH-FDR (10 %) over
all sites x metrics x links. The ADT-path results (B1, B2, B3) are kept as they are; this is a separate measurement.
Usage: python h1_halfpower_meanders.py <revision dir> <out dir> [phase]   phase = gate | analyse
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy.ndimage import uniform_filter1d

REV, OUT = Path(sys.argv[1]), Path(sys.argv[2]); PHASE = sys.argv[3] if len(sys.argv) > 3 else "gate"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REV))
import fronts as F, rstats as R, b1_meander as B
FI = {k: i for i, k in enumerate(F.FIELDS)}
THR, MINW, FRAC, MAXENV = 0.02, 100.0, 0.25, 8
SITES = {
    "SEIR150": dict(reg=(120.0, 162.0, -62.0, -44.0), mea=(140.0, 150.0), up=(125.0, 138.0), down=(150.0, 158.0), band=(-60.0, -47.0)),
    "KERG": dict(reg=(66.0, 98.0, -62.0, -44.0), mea=(76.0, 84.0), up=(70.0, 77.0), down=(84.0, 94.0), band=(-57.0, -48.0)),
    "PAR": dict(reg=(-152.0, -78.0, -64.0, -45.0), mea=(-135.13, -109.38), up=(-150.0, -135.13), down=(-109.38, -80.0), band=(-60.0, -48.0)),
    "DRAKE": dict(reg=(-102.0, -33.0, -64.0, -50.0), mea=(-70.0, -55.0), up=(-100.0, -80.0), down=(-55.0, -35.0), band=(-62.0, -52.0)),
}
YRS = range(2001, 2026)
FILES = [REV/"data"/"monthly_allsat"/f"{y}"/f"{m:02d}.nc" for y in YRS for m in range(1, 13)]
NT = len(FILES)
T = pd.DataFrame({"year": np.repeat(list(YRS), 12), "month": np.tile(np.arange(1, 13), 25)}); T["x"] = T.year + (T.month-0.5)/12


def load_site(c):
    lo0, lo1, la0, la1 = c["reg"]
    U, V, S2, ICE = [], [], [], []
    for f in FILES:
        with xr.open_dataset(f) as ds:
            s = ds.sel(longitude=slice(lo0, lo1), latitude=slice(la0, la1))
            U.append(s.u_mean.values.astype(np.float32)); V.append(s.v_mean.values.astype(np.float32))
            S2.append(s.var_intra.values.astype(np.float32)); ICE.append(np.nan_to_num(s.ice_frac.values).astype(np.float32))
        lat = s.latitude.values.astype(float); lon = s.longitude.values.astype(float)
    return lat, lon, np.array(U), np.array(V), np.array(S2), np.array(ICE)


def detect_stack(KE, lat):
    d = F.meridional_dists(lat)
    cnt = np.zeros((KE.shape[0], KE.shape[2]), np.int16); tab = np.full((KE.shape[0], KE.shape[2], MAXENV, F.NF), np.nan, np.float32)
    for k in range(KE.shape[0]):
        c, tb = F.detect_field(np.asarray(KE[k], float), lat, d, THR, MINW, FRAC, MAXENV); cnt[k] = c; tab[k] = tb
    return cnt, tab


def match(cnt, tab, ccnt, ctab, lat, cols):
    """For each window column and climatological envelope: monthly matched half-power latitude (nt, ncol, MAXENV)."""
    out = np.full((cnt.shape[0], len(cols), MAXENV), np.nan)
    for q, j in enumerate(cols):
        for c in range(ccnt[j]):
            lo, hi = lat[int(ctab[j, c, FI["i_start"]])], lat[int(ctab[j, c, FI["i_end"]])]
            lh = tab[:, j, :, FI["lat_hp"]]; pk = tab[:, j, :, FI["peak_val"]]
            inside = np.isfinite(lh) & (lh >= lo) & (lh <= hi)
            best = np.argmax(np.where(inside, pk, -np.inf), 1); ok = inside.any(1)
            out[ok, q, c] = lh[ok, best[ok]]
    return out


def branch(ccnt, ctab, cols, cov):
    """Follow one climatological branch across the window; returns list of (column index q, envelope c) or None."""
    best, start = -np.inf, None
    for q, j in enumerate(cols):
        for c in range(ccnt[j]):
            if ctab[j, c, FI["peak_val"]] > best and cov[q, c] >= 0.9:
                best, start = ctab[j, c, FI["peak_val"]], (q, c)
    if start is None:
        return []
    path = {start[0]: start[1]}
    for step in (-1, 1):
        q, last = start[0], ctab[cols[start[0]], start[1], FI["lat_hp"]]; gap = 0
        while 0 <= q + step < len(cols):
            q += step; j = cols[q]
            cand = [(abs(ctab[j, c, FI["lat_hp"]] - last), c) for c in range(ccnt[j]) if abs(ctab[j, c, FI["lat_hp"]] - last) <= 0.25]
            if cand:
                c = min(cand)[1]; path[q] = c; last = ctab[j, c, FI["lat_hp"]]; gap = 0
            else:
                gap += 1
                if gap > 4:
                    break
    return sorted(path.items())


def monthly_paths(M, br, cov, lonw):
    qs = np.array([q for q, c in br if cov[q, c] >= 0.9]); cs = np.array([c for q, c in br if cov[q, c] >= 0.9])
    P = M[:, qs, cs] if qs.size else np.full((M.shape[0], 0), np.nan)
    return lonw[qs] if qs.size else np.array([]), P


def metrics_from_path(lon_b, lat_m, smooth):
    ok = np.isfinite(lat_m)
    if lon_b.size < 10 or ok.mean() < 0.9:
        return dict(meanlat_km=np.nan, amp_km=np.nan, sinuosity=np.nan, rms_kappa=np.nan, int_abs_kappa=np.nan)
    la = np.interp(lon_b, lon_b[ok], lat_m[ok])
    m = B.path_metrics(np.c_[lon_b, la], smooth=smooth)
    if m is None:
        return dict(meanlat_km=np.nan, amp_km=np.nan, sinuosity=np.nan, rms_kappa=np.nan, int_abs_kappa=np.nan)
    return dict(meanlat_km=float(np.mean(la))*111.2, amp_km=m["amp_km"], sinuosity=m["sinuosity"],
                rms_kappa=m["rms_kappa"], int_abs_kappa=m["int_abs_kappa"])


def site_geometry(KE, lat, lon, c, ccnt=None, ctab=None, smooths=(150.0,)):
    cols = np.where((lon >= c["mea"][0]) & (lon <= c["mea"][1]))[0]
    if ccnt is None:
        ccnt, ctab = F.detect_field(np.nanmean(np.asarray(KE, float), 0), lat, F.meridional_dists(lat), THR, MINW, FRAC, MAXENV)
    cnt, tab = detect_stack(KE, lat)
    M = match(cnt, tab, ccnt, ctab, lat, cols)
    cov = np.isfinite(M).mean(0)
    br = branch(ccnt, ctab, cols, cov)
    lon_b, P = monthly_paths(M, br, cov, lon[cols])
    rows = []
    for k in range(P.shape[0]):
        r = {"i": k}
        for sm in smooths:
            for kk, v in metrics_from_path(lon_b, P[k], sm).items():
                r[f"{kk}_s{int(sm)}"] = v
        rows.append(r)
    info = dict(n_cols=int(cols.size), n_branch=len(br), n_usable=int(lon_b.size), frac_usable=float(lon_b.size/cols.size),
                lon_usable_min=float(lon_b.min()) if lon_b.size else np.nan, lon_usable_max=float(lon_b.max()) if lon_b.size else np.nan,
                clim_lat_mean=float(np.nanmean(np.nanmean(P, 0))) if lon_b.size else np.nan)
    return pd.DataFrame(rows), info, (ccnt, ctab), lon_b, P


def deform(KEc, lat, lon, cols_mask, dy_km, gain_n, gain_s, jet_lat):
    """Shift climatological KE meridionally by dy_km(lon) and apply a flank-dependent gain."""
    out = np.empty_like(KEc)
    for j in range(lon.size):
        dlat = dy_km[j]/111.2
        out[:, j] = np.interp(lat - dlat, lat, KEc[:, j], left=np.nan, right=np.nan)
    if gain_n != 1.0 or gain_s != 1.0:
        yk = (lat[:, None] - jet_lat[None, :])*111.2
        wgt = 0.5*(1 + np.tanh(yk/50.0))          # 0 south of the jet, 1 north
        out = out*(gain_s + (gain_n-gain_s)*wgt)
    return out


def synthetic(KE, lat, lon, c, cc, nreal=40, seed=0):
    rng = np.random.default_rng(seed)
    KEc = np.nanmean(np.asarray(KE, float), 0); anom = np.asarray(KE, float) - KEc[None]
    cols = (lon >= c["mea"][0]) & (lon <= c["mea"][1])
    W = c["mea"][1]-c["mea"][0]; s = np.clip((lon - c["mea"][0])/W, 0, 1)
    # jet latitude per column: climatological peak KE latitude in the band
    band = (lat >= c["band"][0]) & (lat <= c["band"][1])
    jet_lat = lat[band][np.nanargmax(np.where(np.isfinite(KEc[band]), KEc[band], -1), 0)]
    ramp = (np.arange(NT) - (NT-1)/2)/(NT-1)       # -0.5 .. 0.5 over 25 yr
    SC = {"N": (0, 0, 1.0, 1.0), "T": (25, 0, 1.0, 1.0), "F": (0, 25, 1.0, 1.0), "U": (0, 0, 1.2, 1.2), "A": (0, 0, 1.2, 1.0)}
    res = []
    nblk = NT//12
    for name, (dT, dF, gn, gs) in SC.items():
        for r in range(nreal):
            blocks = rng.integers(0, nblk, nblk)
            idx = np.concatenate([np.arange(b*12, b*12+12) for b in blocks])
            SY = np.empty_like(anom, dtype=np.float32)
            for k in range(NT):
                dy = ramp[k]*(dT + dF*np.sin(2*np.pi*s))*np.ones(lon.size)
                SY[k] = deform(KEc, lat, lon, cols, dy, 1 + (gn-1)*(ramp[k]+0.5), 1 + (gs-1)*(ramp[k]+0.5), jet_lat) + anom[idx[k]]
            geo, info, *_ = site_geometry(SY, lat, lon, c, *cc)
            rec = dict(scenario=name, real=r, frac_usable=info["frac_usable"])
            for col in [q for q in geo.columns if q != "i"]:
                y = geo[col].values; ok = np.isfinite(y)
                rec[col] = R.theil_sen(T.x.values[ok], y[ok])[0]*25 if ok.sum() > 100 else np.nan   # change over 25 yr
            res.append(rec)
        print("  scenario", name, "done", flush=True)
    return pd.DataFrame(res)


def jet_speed(Usp, lat, lon, c, sec):
    m = (lon >= c[sec][0]) & (lon <= c[sec][1]); b = (lat >= c["band"][0]) & (lat <= c["band"][1])
    return np.nanmean(np.nanmax(Usp[:, b][:, :, m], 1), 1)


def box_eke(EKE, ICE, lat, lon, c, sec):
    m = (lon >= c[sec][0]) & (lon <= c[sec][1]); b = (lat >= c["band"][0]) & (lat <= c["band"][1])
    msk = b[:, None] & m[None, :] & ~(ICE.max(0) > 0.10)
    w = np.cos(np.radians(lat))[:, None]*np.ones(lon.size)[None, :]
    return np.array([np.nansum(EKE[k][msk]*w[msk])/np.sum(w[msk]) for k in range(EKE.shape[0])])


def prep(y):
    y = R.deseasonalise(np.asarray(y, float), T.month.values); ok = np.isfinite(y)
    b, a = R.theil_sen(T.x.values[ok], y[ok]); r = y - (a + b*T.x.values)
    r = pd.Series(r).interpolate(limit=2).values
    return uniform_filter1d(np.nan_to_num(r - np.nanmean(r)), 2, mode="nearest") * np.where(np.isfinite(r), 1, np.nan)


def lagcorr(a, b, lags):
    out = []
    for L in lags:   # b lags a by L months
        if L >= 0: x, y = a[:NT-L], b[L:]
        else: x, y = a[-L:], b[:NT+L]
        ok = np.isfinite(x) & np.isfinite(y); out.append(np.corrcoef(x[ok], y[ok])[0, 1])
    return np.array(out)


def surrogate_p(a, b, lags, n=1000, seed=0):
    rng = np.random.default_rng(seed); r0 = lagcorr(a, b, lags); stat = np.nanmax(np.abs(r0))
    aa = np.nan_to_num(a); Fa = np.fft.rfft(aa); cnt = 0
    for _ in range(n):
        ph = np.exp(1j*rng.uniform(0, 2*np.pi, Fa.size)); ph[0] = 1
        s = np.fft.irfft(Fa*ph, n=aa.size)*np.where(np.isfinite(a), 1, np.nan)
        cnt += np.nanmax(np.abs(lagcorr(s, b, lags))) >= stat
    return r0, (1+cnt)/(1+n)


if __name__ == "__main__":
    gate_rows, synth_all = [], []
    for name, c in SITES.items():
        print(name, flush=True)
        lat, lon, U, V, S2, ICE = load_site(c)
        KE = 0.5*(U.astype(float)**2 + V.astype(float)**2)
        geo, info, cc, lon_b, P = site_geometry(KE, lat, lon, c, smooths=(100.0, 150.0, 200.0))
        info["site"] = name; info["G1_pass"] = bool(info["frac_usable"] >= 0.9)
        geo.to_csv(OUT/f"h1_geometry_{name}.csv", index=False)
        np.savez_compressed(OUT/f"h1_paths_{name}.npz", lon=lon_b, lat_hp=P)
        print(" ", info, flush=True)
        if PHASE == "gate":
            if info["G1_pass"]:
                sy = synthetic(KE, lat, lon, c, cc); sy["site"] = name; synth_all.append(sy)
            gate_rows.append(info)
            continue
        # analyse phase
        g1 = pd.read_csv(OUT/"h1_gate.csv").set_index("site")
        if not bool(g1.loc[name, "G1_pass"]) or not bool(g1.loc[name, "G2_pass_meanlat"]) and not bool(g1.loc[name, "G2_pass_amp"]):
            continue
        Usp = np.hypot(U.astype(float), V.astype(float))
        ub, vb = U.mean(0), V.mean(0)
        EKE = 0.5*((U-ub)**2 + (V-vb)**2) + 0.5*S2
        drivers = {"up_speed": jet_speed(Usp, lat, lon, c, "up"), "down_speed": jet_speed(Usp, lat, lon, c, "down"),
                   "mea_eke": box_eke(EKE, ICE, lat, lon, c, "mea")*1e4, "down_eke": box_eke(EKE, ICE, lat, lon, c, "down")*1e4}
        pd.DataFrame(drivers).assign(year=T.year, month=T.month).to_csv(OUT/f"h1_drivers_{name}.csv", index=False)
        rows = []
        for col in [q for q in geo.columns if q.endswith("_s150")]:
            y = geo[col].values
            o = R.trend_bcz(T.x.values, y, month=T.month.values)
            rows.append(dict(site=name, metric=col, kind="trend", slope_per_decade=10*o["slope"], ci_lo=10*o["ci_lo"], ci_hi=10*o["ci_hi"],
                             p=o["p"], mean=float(np.nanmean(y)), coverage=float(np.isfinite(y).mean())))
            gp = prep(y)
            for link, drv, lags in [("up_speed->geom", "up_speed", range(0, 4)), ("geom->mea_eke", "mea_eke", range(0, 4)),
                                    ("geom->down_eke", "down_eke", range(0, 4)), ("geom~down_speed", "down_speed", range(-3, 4))]:
                dp = prep(drivers[drv])
                if link.startswith("up"):
                    r0, p = surrogate_p(gp, dp, [-L for L in lags])     # geometry lags speed
                else:
                    r0, p = surrogate_p(gp, dp, list(lags))
                rows.append(dict(site=name, metric=col, kind=link, r_by_lag=json.dumps([round(float(v), 3) for v in r0]),
                                 r_maxabs=float(r0[np.nanargmax(np.abs(r0))]), p=p))
        pd.DataFrame(rows).to_csv(OUT/f"h1_results_{name}.csv", index=False)
        print(pd.DataFrame(rows).round(3).to_string(), flush=True)
    if PHASE == "gate":
        G = pd.DataFrame(gate_rows)
        if synth_all:
            S = pd.concat(synth_all); S.to_csv(OUT/"h1_synthetic.csv", index=False)
            summ = []
            for site, s in S.groupby("site"):
                d = s.groupby("scenario").agg(["mean", "std"])
                Trec = d.loc["T", ("meanlat_km_s150", "mean")] - d.loc["N", ("meanlat_km_s150", "mean")]
                spurA = d.loc["A", ("meanlat_km_s150", "mean")] - d.loc["N", ("meanlat_km_s150", "mean")]
                spurU = d.loc["U", ("meanlat_km_s150", "mean")] - d.loc["N", ("meanlat_km_s150", "mean")]
                Famp = d.loc["F", ("amp_km_s150", "mean")] - d.loc["N", ("amp_km_s150", "mean")]
                Nsd = d.loc["N", ("amp_km_s150", "std")]
                summ.append(dict(site=site, T_recovered_km=Trec, A_spurious_km=spurA, U_spurious_km=spurU,
                                 F_amp_response_km=Famp, N_amp_sd_km=Nsd,
                                 A_amp_spurious_km=d.loc["A", ("amp_km_s150", "mean")] - d.loc["N", ("amp_km_s150", "mean")],
                                 G2_pass_meanlat=bool(Trec >= 12.5 and abs(spurA) < 0.5*Trec and abs(spurU) < 0.5*Trec),
                                 G2_pass_amp=bool(Famp > 2*Nsd)))
            G = G.merge(pd.DataFrame(summ), on="site", how="left")
        G.to_csv(OUT/"h1_gate.csv", index=False); print(G.round(3).to_string())
