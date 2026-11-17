"""Step 7: side-by-side summary of the Product 1 and Product 4 diagnostics for the Supplement comparison figure/table.

Purpose : collect, for both products, every number of the P1-vs-P4 comparison in one long table, each with its source file.
Inputs  : outputs/retest/retest_summary.csv (trend_bcz, 300 nulls: the source of the paper's band/sector trends and of Fig. S3b),
          outputs/a5/a5_trends_{allsat,twosat}.csv (older trend_calibrated test; listed because the brief names it),
          outputs/a5/a5_series_*.csv, outputs/a2/a2h_twosat_vs_allsat.csv, outputs/x7/x7_trends.csv,
          outputs/drake/dk12/dk12_band_trends_allsat.csv, outputs/a1/a1_maps_bc.nc, outputs/f4/f4_tilt_maps.nc (1-deg band cells),
          outputs/p4/p4_band_series.csv, p4_band_periods.csv, p4_bins_trends*.csv, p4_drake_domains.csv, p4_eke_trend_map.nc,
          p4_wind_series*.csv, p4_lagcorr_summary.csv, p4_exponents.csv
Outputs : outputs/p4/p4_summary.csv  quantity, unit, P1, P1_lo, P1_hi, P1_p, P4, P4_lo, P4_hi, P4_p, source_P1, source_P4, note
          outputs/p4/p4_checks.csv   re-run checks (trend_bcz on the recomputed P4 band series vs retest_summary; domain trends 2004-2025)
Method  : no new statistics except (i) the pattern correlation of the two 1-deg EKE trend maps (Pearson r over the cells finite in
          both; also cos(lat)-weighted and restricted to the 1-deg ACC band cells of f4_tilt_maps 'band'), (ii) trend_bcz(n_outer=300,
          seed 0) of the ACC band + 100 km series over 2004-2025 (as x14 / Fig. 4c) and of the recomputed P4 ACC series (check).
          Bin counts: bins with a finite trend; 'significant' = p < 0.05 (uncorrected, as Fig. 3d).
Run     : python scripts/06_product4/p4_07_summary.py
"""
import p4_common as C
import numpy as np, pandas as pd, xarray as xr

R = C.R; O = C.OUT; H = C.HERE
rows = []


def add(q, unit, p1, p4, s1, s4, note=""):
    d = dict(quantity=q, unit=unit, source_P1=s1, source_P4=s4, note=note)
    for tag, v in (("P1", p1), ("P4", p4)):
        v = v if isinstance(v, dict) else dict(val=v)
        d[tag] = v.get("val", np.nan); d[f"{tag}_lo"] = v.get("lo", np.nan); d[f"{tag}_hi"] = v.get("hi", np.nan); d[f"{tag}_p"] = v.get("p", np.nan)
    rows.append(d)


# ---------------------------------------------------------------- 1. band and sector trends
RT = pd.read_csv(O / "retest" / "retest_summary.csv").set_index("key")
SER = {p: pd.read_csv(O / "a5" / f"a5_series_{p}.csv") for p in ("allsat", "twosat")}
A5 = {p: pd.read_csv(O / "a5" / f"a5_trends_{p}.csv").query("period == '2001-2025'").set_index("series") for p in ("allsat", "twosat")}
for reg in ("ACC", "Pacific", "Atlantic", "Indian"):
    for v in ("KE", "EKE"):
        rel, ab, old = {}, {}, {}
        for p in ("allsat", "twosat"):
            r = RT.loc[f"{p}|{reg}_{v}|2001-2025"]; mu = SER[p][f"{reg}_{v}"].mean()
            rel[p] = dict(val=1000 * r.slope / mu, lo=1000 * r.ci_lo / mu, hi=1000 * r.ci_hi / mu, p=r.p_new)
            ab[p] = dict(val=10 * r.slope, lo=10 * r.ci_lo, hi=10 * r.ci_hi, p=r.p_new)
            a = A5[p].loc[f"{reg}_{v}"]
            old[p] = dict(val=a.pct_per_decade, lo=1000 * a.ci_cal_lo / a["mean"], hi=1000 * a.ci_cal_hi / a["mean"], p=a.p_cal)
        lab = "ACC band" if reg == "ACC" else f"{reg} sector"
        add(f"{lab} {v} trend 2001-2025", "% of mean per decade", rel["allsat"], rel["twosat"], "retest_summary.csv (allsat)", "retest_summary.csv (twosat)",
            "trend_bcz, 300 nulls (paper / Fig. S3b values)")
        add(f"{lab} {v} trend 2001-2025", "cm2 s-2 per decade", ab["allsat"], ab["twosat"], "retest_summary.csv (allsat)", "retest_summary.csv (twosat)", "")
        add(f"{lab} {v} trend 2001-2025 (a5_trends file)", "% of mean per decade", old["allsat"], old["twosat"], "a5_trends_allsat.csv", "a5_trends_twosat.csv",
            "older trend_calibrated test (p_cal, ci_cal); superseded by retest_summary in the paper")
c1, c4 = RT.loc["ACC EKE-KE allsat"], RT.loc["ACC EKE-KE twosat"]
add("ACC band EKE - KE contrast 2001-2025", "percentage points per decade", dict(val=c1.slope, lo=c1.ci_lo, hi=c1.ci_hi, p=c1.p_new),
    dict(val=c4.slope, lo=c4.ci_lo, hi=c4.ci_hi, p=c4.p_new), "retest_summary.csv", "retest_summary.csv", "trend of EKE/mean - KE/mean")

# ---------------------------------------------------------------- 2. half-power latitude
HP = pd.read_csv(O / "a2" / "a2h_twosat_vs_allsat.csv")
for reg in ("ACC", "Pacific", "Atlantic", "Indian"):
    g = {pn: HP[(HP["product"] == pn) & (HP.region == reg)].iloc[0] for pn in ("all-sat", "two-sat")}
    f = lambda r: dict(val=r.pos_km_dec, lo=r.pos_lo, hi=r.pos_hi, p=r.pos_p)
    add(f"Half-power latitude trend, {'circumpolar' if reg == 'ACC' else reg}, FIXED threshold 200 cm2 s-2", "km per decade (north +)",
        f(g["all-sat"]), f(g["two-sat"]), "a2h_twosat_vs_allsat.csv", "a2h_twosat_vs_allsat.csv",
        f"baseline set (02e): {int(g['all-sat'].n_lon)} P1 tracks, {int(g['two-sat'].n_lon)} P4 tracks; this is what Fig. S3c shows")
x7 = pd.read_csv(C.OUT / "x7" / "x7_trends.csv")
r = x7[(x7.region == "circumpolar") & (x7.series == "lat_hp")].iloc[0]
add("Half-power latitude trend, circumpolar, RELATIVE threshold (20 % of local max; main text)", "km per decade (north +)",
    dict(val=r.slope, lo=r.lo, hi=r.hi, p=r.p), np.nan, "outputs/x7/x7_trends.csv", "not computed",
    "1826 P1 tracks; no Product 4 relative-threshold tracking exists")

# ---------------------------------------------------------------- 3. 2-degree bins
B1 = pd.read_csv(O / "drake" / "dk12" / "dk12_band_trends_allsat.csv").rename(columns={"slope_dec": "slope"})
for tag, f, nt in (("", "p4_bins_trends.csv", "levels on P1 mean ADT (primary)"), ("_p4adt", "p4_bins_trends_p4adt.csv", "levels on P4's own mean ADT")):
    if not (H / f).exists():
        continue
    B4 = pd.read_csv(H / f)
    for q, fn, unit in (("2-deg bins with a finite EKE trend", lambda B: int(B.slope.notna().sum()), "count"),
                        ("2-deg bins with positive EKE trend", lambda B: 100 * float((B.slope > 0).sum()) / B.slope.notna().sum(), "% of bins"),
                        ("2-deg bins with p < 0.05", lambda B: int((B.p < 0.05).sum()), "count"),
                        ("2-deg bins with p < 0.05 and positive", lambda B: int(((B.p < 0.05) & (B.slope > 0)).sum()), "count"),
                        ("2-deg bins: median EKE trend", lambda B: float(B.slope.median()), "cm2 s-2 per decade"),
                        ("2-deg bins: max EKE trend", lambda B: float(B.slope.max()), "cm2 s-2 per decade")):
        add(q + (f" [P4 {nt}]" if tag else ""), unit, fn(B1), fn(B4), "dk12_band_trends_allsat.csv", f, "")
    m = B1.merge(B4, on="lon", suffixes=("_1", "_4")).dropna(subset=["slope_1", "slope_4"])
    add("2-deg bins: correlation of P1 and P4 EKE trends" + (f" [P4 {nt}]" if tag else ""), "Pearson r", float(np.corrcoef(m.slope_1, m.slope_4)[0, 1]),
        np.nan, "dk12 vs " + f, "", "single value in column P1")
    add("2-deg bins: ratio of summed EKE trends P1/P4" + (f" [P4 {nt}]" if tag else ""), "ratio", float(m.slope_1.sum() / m.slope_4.sum()), np.nan, "", "", "single value in column P1")

# ---------------------------------------------------------------- 4. map pattern correlation
if (H / "p4_eke_trend_map.nc").exists():
    M1 = xr.open_dataset(O / "a1" / "a1_maps_bc.nc"); M4 = xr.open_dataset(H / "p4_eke_trend_map.nc")
    band1 = xr.open_dataset(O / "f4" / "f4_tilt_maps.nc").band.values > 0
    wlat = np.cos(np.radians(M1.lat.values))[:, None] * np.ones((1, M1.lon.size))
    for v in ("EKE", "KE"):
        a, b = 10 * M1[f"{v}_slope"].values, 10 * M4[f"{v}_slope"].values; ok = np.isfinite(a) & np.isfinite(b)
        for sel, nm in ((ok, "all common cells"), (ok & band1, "common cells in the 1-deg ACC band")):
            r_ = float(np.corrcoef(a[sel], b[sel])[0, 1])
            aw, bw, ww = a[sel], b[sel], wlat[sel]
            ma, mb = np.sum(ww * aw) / ww.sum(), np.sum(ww * bw) / ww.sum()
            rw = float(np.sum(ww * (aw - ma) * (bw - mb)) / np.sqrt(np.sum(ww * (aw - ma) ** 2) * np.sum(ww * (bw - mb) ** 2)))
            add(f"1-deg {v} trend maps: P1-P4 pattern correlation, {nm} (n = {int(sel.sum())}), unweighted", "Pearson r", r_, np.nan,
                "a1_maps_bc.nc vs p4_eke_trend_map.nc", "", "single value for the pair (P1 column)")
            add(f"1-deg {v} trend maps: P1-P4 pattern correlation, {nm} (n = {int(sel.sum())}), cos(lat)-weighted", "Pearson r", rw, np.nan,
                "a1_maps_bc.nc vs p4_eke_trend_map.nc", "", "single value for the pair (P1 column)")
            add(f"1-deg {v} trend maps: mean trend, {nm}", "cm2 s-2 per decade", float(np.mean(a[sel])), float(np.mean(b[sel])), "a1_maps_bc.nc", "p4_eke_trend_map.nc", "")
        add(f"1-deg {v} trend maps: cells with a trend", "count", int(np.isfinite(a).sum()), int(np.isfinite(b).sum()), "a1_maps_bc.nc", "p4_eke_trend_map.nc", "")
        add(f"1-deg {v} trend maps: cells with positive trend", "% of cells", 100 * float(np.nanmean(a[np.isfinite(a)] > 0)), 100 * float(np.nanmean(b[np.isfinite(b)] > 0)),
            "a1_maps_bc.nc", "p4_eke_trend_map.nc", "")
        add(f"1-deg {v} trend maps: cells passing 10 % FDR", "count", int(M1[f"{v}_fdr"].sum()), int(M4[f"{v}_fdr"].sum()), "a1_maps_bc.nc", "p4_eke_trend_map.nc", "")

# ---------------------------------------------------------------- 5. period increases
PER = pd.read_csv(H / "p4_band_periods.csv")
for s in ("ACC_KE", "ACC_EKE"):
    g = lambda pr: float(PER[(PER["product"] == pr) & (PER.series == s)].pct_increase.iloc[0])
    add(f"ACC band {s.split('_')[1]}: 2018-2025 mean over 2001-2017 mean", "%", g("P1 (a5_series_allsat)"), g("P4 (recomputed here)"),
        "a5_series_allsat.csv", "p4_band_series.csv (= a5_series_twosat.csv)", "ratio of the monthly-value means")
for f, nt in (("p4_drake_domains.csv", ""), ("p4_drake_domains_p4adt.csv", " [P4 levels on P4's own mean ADT]")):
    if not (H / f).exists():
        continue
    D = pd.read_csv(H / f)
    for dom in ("upstream 86-68W", "passage 68-58W", "downstream 58-40W"):
        d = D[D.domain == dom]
        add(f"Drake {dom}: EKE between SAF and PF, 2018-2025 over 2001-2017{nt}", "%",
            float(d.iloc[0].pct_increase_2018_2025_vs_2001_2017), float(d.iloc[1].pct_increase_2018_2025_vs_2001_2017),
            "dk12_band_annual_allsat.csv (figS2_timing method)", f.replace("domains", "bins_annual").replace("drake_", ""), "")

# ---------------------------------------------------------------- 6. wind relations (ACC band + 100 km)
L = pd.read_csv(H / "p4_lagcorr_summary.csv")
for s in ("KE", "EKE"):
    g = lambda pr: L[(L["product"] == pr) & (L.series == s)].iloc[0]
    a, b = g("P1"), g("P4")
    add(f"{s}: maximum lagged correlation with wind stress (wind leading 0-36 months)", "r",
        dict(val=a.r_max, p=a.p_fw), dict(val=b.r_max, p=b.p_fw), "outputs/f4e/f4e_series.csv (x12 method)", "p4_wind_series.csv", "p = family-wise over 37 lags")
    add(f"{s}: lag of the maximum correlation", "months (95 % block-bootstrap range in lo/hi)",
        dict(val=a.lag_max, lo=a.lag_max_lo, hi=a.lag_max_hi), dict(val=b.lag_max, lo=b.lag_max_lo, hi=b.lag_max_hi), "p4_lagcorr_summary.csv", "p4_lagcorr_summary.csv",
        f"pointwise-significant lags P1: {a.sig_lags}; P4: {b.sig_lags}")
E = pd.read_csv(H / "p4_exponents.csv")
for s in ("KE", "EKE"):
    g = lambda pr: E[(E["product"] == pr) & (E.series == s)].iloc[0]
    a, b = g("P1"), g("P4")
    add(f"{s}: exponent b on previous-year wind stress, trends retained (2004-2025)", "dimensionless",
        dict(val=a.b, lo=a.b_lo, hi=a.b_hi, p=a.p), dict(val=b.b, lo=b.b_lo, hi=b.b_hi, p=b.p), "x9/x10 method on f4e_series.csv", "p4_wind_series.csv", "")
    add(f"{s}: exponent b, both series detrended (2004-2025)", "dimensionless",
        dict(val=a.b_detr, lo=a.bd_lo, hi=a.bd_hi, p=a.p_bd_extra), dict(val=b.b_detr, lo=b.bd_lo, hi=b.bd_hi, p=b.p_bd_extra),
        "x9/x10 method", "p4_wind_series.csv", "p: extra test (x9 surrogates, both detrended), not in x9")
    add(f"{s}: share of the trend-retaining exponent due to shared trends", "fraction", a.shared_trend_frac, b.shared_trend_frac, "", "", "(b - b_detr) / b")

# ---------------------------------------------------------------- 7. checks and domain trends 2004-2025
chk = []
S4 = pd.read_csv(H / "p4_band_series.csv")
for v in ("KE", "EKE"):
    o = R.trend_bcz(S4.x.values, S4[f"ACC_{v}"].values, month=S4.month.values, n_outer=300, seed=0)
    r = RT.loc[f"twosat|ACC_{v}|2001-2025"]
    chk.append(dict(check=f"trend_bcz on recomputed P4 ACC_{v} vs retest_summary twosat", here_slope=o["slope"], ref_slope=r.slope,
                    here_p=o["p"], ref_p=r.p_new, here_lo=o["ci_lo"], ref_lo=r.ci_lo, here_hi=o["ci_hi"], ref_hi=r.ci_hi))
x14 = pd.read_csv(C.OUT / "x14" / "x14_results.csv")
x14 = x14[(x14.domain == "buf100") & (x14.test == "trend 2004-2025")].set_index("series")
for lab, f in (("P1", O / "f4e" / "f4e_series.csv"), ("P4", H / "p4_wind_series.csv")):
    S = pd.read_csv(f); c4 = S.year.values >= 2004
    res = {}
    for v in ("wind", "KE", "EKE"):
        y = S[v].values[c4]; o = R.trend_bcz(S.x.values[c4], y, month=S.month.values[c4], n_outer=300); mu = y.mean()
        res[v] = dict(val=1000 * o["slope"] / mu, lo=1000 * o["ci_lo"] / mu, hi=1000 * o["ci_hi"] / mu, p=o["p"])
        if lab == "P1":
            chk.append(dict(check=f"P1 band+100 km {v} trend 2004-2025 vs x14 buf100", here_slope=res[v]["val"], ref_slope=x14.loc[v].pct,
                            here_p=o["p"], ref_p=x14.loc[v].p, here_lo=res[v]["lo"], ref_lo=x14.loc[v].lo, here_hi=res[v]["hi"], ref_hi=x14.loc[v].hi))
    if lab == "P1":
        P1res = res
    else:
        for v in ("KE", "EKE"):
            add(f"ACC band + 100 km {v} trend 2004-2025 (Fig. 4c domain)", "% of mean per decade", P1res[v], res[v],
                "outputs/f4e/f4e_series.csv", "p4_wind_series.csv", "trend_bcz, 300 nulls, as x14")
pd.DataFrame(chk).to_csv(H / "p4_checks.csv", index=False); print(pd.DataFrame(chk).to_string(index=False))
T = pd.DataFrame(rows)[["quantity", "unit", "P1", "P1_lo", "P1_hi", "P1_p", "P4", "P4_lo", "P4_hi", "P4_p", "source_P1", "source_P4", "note"]]
T.to_csv(H / "p4_summary.csv", index=False)
with pd.option_context("display.width", 250, "display.max_colwidth", 70):
    print(T[["quantity", "P1", "P1_lo", "P1_hi", "P1_p", "P4", "P4_lo", "P4_hi", "P4_p"]].round(3).to_string(index=False))
