"""Comparison with previous studies: SI Sect. S4 and Fig. S4 (SI v3, 3 Nov 2026).

For each earlier estimate, the trend of this study is recomputed for the region, period (overlap with our 2001-2025 record)
and units of that estimate as closely as the shipped outputs allow, and also given for 2001-2025. Published values are typed
below with the page they come from (Hogg et al. 2015: values from the repository notes only, PDF not available ->
[AUTHOR CHECK], docs/AUTHOR_TODO.md). No new analysis method: every trend is the calibrated Theil-Sen test of the paper
(rstats.trend_bcz: deseasonalised Theil-Sen, calibrated AR noise, max p of AIC and AR(12) noise models, 95 % interval).
Ported from figures_final/scripts_SI_v3/tableS4_previous_studies.py (paths through osr11.paths; the Zhang et al. 2021
Pacific-box series is split off as a raw-data step whose output is shipped).

Purpose : Fig. S4 values (SI S4).
Inputs  : outputs/m1/m1_trends.csv, m1_series_{allsat,twosat}.csv   (02l_earlier_studies: Hogg boxes, ACC band, both products)
          outputs/m2/m2_trends.csv                                    (35-65 S EKE area integral processed as Martinez-Moreno
                                                                       et al. 2021; archived m2_mm21_procedure.py, see README)
          outputs/x1/x1_tracks.npz                                    (02j: half-power and fixed-contour latitude per track)
          outputs/drake/dk12/dk12_band_series_allsat.npz, dk12_band_trends_allsat.csv  (02i: EKE between mean SAF and PF,
                                                                       2-deg bins; n_pix as area weights)
          outputs/m1/m1_zhang_pacific_62-40S_allsat.csv              (monthly EKE of the Pacific box with the 62-40 S limits
                                                                       published by Zhang et al. 2021, p. 11; Product 1 only)
Outputs : outputs/m1/m1_previous_studies.csv (read by figures/figS4_previous_studies.py)
          with --zhang-series: outputs/m1/m1_zhang_pacific_62-40S_allsat.csv (NEEDS RAW DATA: $OSR11_DATA/monthly_allsat;
          the shipped file is the output of this step)
Run     : python scripts/05_stats/05k_previous_studies.py --zhang-series   (raw data; reads the 300 monthly files twice, ~3 min)
          python scripts/05_stats/05k_previous_studies.py                  (no raw data; ~6 min)
Seeds   : rstats.trend_bcz default seed (0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys
import numpy as np, pandas as pd
from osr11 import rstats as R, paths as PA

OUT = PA.OUT
ZSER = PA.out("m1") / "m1_zhang_pacific_62-40S_allsat.csv"


def zhang_series():
    """Monthly EKE (cm2 s-2) of the Pacific box 150E-72W (288E), 62-40S, Product 1, mask as 02l (finite in all months,
    ice fraction <= 0.10 in every month), cos-latitude weighted; EKE = 0.5|U_m - ubar|^2 + 0.5 var_intra."""
    import xarray as xr
    fs = PA.monthly_files("allsat")
    assert fs[0].exists(), f"{fs[0]} not found: this step needs $OSR11_DATA/monthly_allsat"
    with xr.open_dataset(fs[0]) as d0:
        lat = d0.latitude.values.astype(float); lon = d0.longitude.values.astype(float)
    sub = (lat >= -62) & (lat <= -40); L = np.mod(lon, 360); jj = (L >= 150) & (L < 288)
    su = sv = nv = imax = 0
    for f in fs:
        with xr.open_dataset(f) as d:
            u = d.u_mean.values[sub][:, jj].astype(float); v = d.v_mean.values[sub][:, jj].astype(float)
            ic = np.nan_to_num(d.ice_frac.values[sub][:, jj])
        ok = np.isfinite(u) & np.isfinite(v); su = su + np.where(ok, u, 0); sv = sv + np.where(ok, v, 0); nv = nv + ok; imax = np.fmax(imax, ic)
    valid = (nv == len(fs)) & (imax <= 0.10); ub, vb = su / np.maximum(nv, 1), sv / np.maximum(nv, 1)
    wt = np.broadcast_to(np.cos(np.radians(lat[sub]))[:, None], valid.shape)
    y = []
    for f in fs:
        with xr.open_dataset(f) as d:
            u = d.u_mean.values[sub][:, jj].astype(float); v = d.v_mean.values[sub][:, jj].astype(float); s2 = d.var_intra.values[sub][:, jj].astype(float)
        e = 0.5 * ((u - ub) ** 2 + (v - vb) ** 2) + 0.5 * s2
        y.append(1e4 * np.sum(e[valid] * wt[valid]) / np.sum(wt[valid]))
    pd.DataFrame(dict(year=np.repeat(np.arange(2001, 2026), 12), month=np.tile(np.arange(1, 13), 25), eke_cm2s2=y,
                      n_pix=int(valid.sum()))).to_csv(ZSER, index=False)
    print("wrote", ZSER, "valid pixels:", int(valid.sum()), flush=True)


if "--zhang-series" in sys.argv[1:]:
    zhang_series(); raise SystemExit(0)

KM = 111.195
X = 2001 + (np.arange(300) + 0.5) / 12; YR = np.repeat(np.arange(2001, 2026), 12); MO = np.tile(np.arange(1, 13), 25)
PRODS = (("allsat", "Product 1"), ("twosat", "Product 4"))


def tb(y, y0, y1, scale=1.0, rel=False):
    """calibrated trend of the monthly series y over y0..y1 (inclusive): per decade x scale, or % of period mean per decade"""
    s = (YR >= y0) & (YR <= y1); yy = np.asarray(y, float)[s]
    o = R.trend_bcz(X[s], yy, month=MO[s], n_outer=300)
    f = (1000 / np.nanmean(yy)) if rel else 10 * scale
    return dict(val=o["slope"] * f, lo=o["ci_lo"] * f, hi=o["ci_hi"] * f, p=o["p"])


rows = []


def add(study, metric, region, period, pub, pub_lo, pub_hi, pub_int, our_series, our_period, our, src, note="", panel=""):
    rows.append(dict(panel=panel, study=study, metric=metric, region=region, period=period, published_value=pub,
                     published_lo=pub_lo, published_hi=pub_hi, published_interval=pub_int, our_series=our_series,
                     our_period=our_period, our_value=our["val"] if our else np.nan, our_lo=our["lo"] if our else np.nan,
                     our_hi=our["hi"] if our else np.nan, our_interval="95 %" if our else "", p=our["p"] if our else np.nan,
                     source=src, note=note))
    print(f"{study} | {region} | {our_series} {our_period}: " + (f"{our['val']:+.2f} [{our['lo']:+.2f}, {our['hi']:+.2f}] p={our['p']:.3f}" if our else "-"), flush=True)


M1 = pd.read_csv(OUT / "m1" / "m1_trends.csv").set_index("series")
def m1(prod, key, per, rel=True):
    r = M1.loc[f"{prod}|{key}|{per}"]
    if rel:
        return dict(val=r.pct_per_decade, lo=r.pct_ci_lo, hi=r.pct_ci_hi, p=r.p)
    return dict(val=r.slope_per_decade, lo=r.ci_lo, hi=r.ci_hi, p=r.p)
SER = {p: pd.read_csv(OUT / "m1" / f"m1_series_{p}.csv") for p, _ in PRODS}

# ---------------------------------------------------------------- (a) relative EKE trends, % per decade
HOGG = "Hogg et al. (2015)"; HSRC = "notes of the working repository (revision/plan1019/plan1019_results.md, sect. 0); Hogg et al. 2015 PDF not available"
for key, reg, v, e in (("Hogg_Indian_EKE", "Indian box 40-150E, 57-44S", 6.0, 1.6), ("Hogg_Pacific_EKE", "Pacific box 150E-72W, 62-48S", 6.4, 1.8)):
    for prod, pn in PRODS:
        for per in ("2001-2012", "2001-2025"):
            add(HOGG, "EKE trend (% per decade)", reg, "1993-2012", v, v - e, v + e, "90 % (relative to the 1993-2001 mean)",
                pn, per, m1(prod, key, per), HSRC, "ours: % of the period mean; crossover (theirs) vs gridded (ours) EKE", "a")

MM = "Martinez-Moreno et al. (2021)"
M2 = pd.read_csv(OUT / "m2" / "m2_trends.csv")
for prod, pn in PRODS:
    for per in ("2001-2020", "2001-2025"):
        r = M2[(M2["product"] == prod) & (M2.period == per)].iloc[0]
        lo, hi = [float(q) for q in r.bcz_ci.strip("()").split(",")]
        our = dict(val=100 * r.bcz_slope_PJm_dec / r.mean_PJm, lo=100 * lo / r.mean_PJm, hi=100 * hi / r.mean_PJm, p=r.bcz_p)
        add(MM, "EKE trend (% per decade)", "ACC region (their mask: 35-60S without boundary currents and gyres); ours: 35-65S area integral, not their ACC-region mask",
            "1993-2020", 5.1, np.nan, np.nan, "none given in text (standard error only in Fig. 5c)", pn, per, our,
            "Nature Clim. Change 11, p. 400 (text and Fig. 5c); region definition in Methods (online)",
            "ours: m2 (EKE area integral 35-65S, MM21 processing), calibrated test on monthly values", "a")

BL = "Barcelo-Llull et al. (2025)"
for prod, pn in PRODS:
    for per in ((2001, 2022), (2001, 2025)):
        our = tb(SER[prod]["ACC_EKE"].values, *per, rel=True)
        add(BL, "EKE trend (% per decade)", "high-EKE regions (global, 5 % of ocean area); ours: ACC band (not the same region)",
            "1993-2022", 5.7 if prod == "allsat" else 2.5, np.nan, np.nan, "standard error shown only in Fig. 2C",
            pn, f"{per[0]}-{per[1]}", our, "Sci. Rep. 15:21754, p. 5 (high-EKE regions: all-sat 5.7 %, two-sat 2.5 % per decade)",
            "published value is for the product in 'our_series' (all-sat = Product 1 type, two-sat = Product 4 type)", "a")

# all-sat / two-sat ratio (absolute trends), table only
for per in ((2001, 2022), (2001, 2025)):
    t1 = tb(SER["allsat"]["ACC_EKE"].values, *per); t4 = tb(SER["twosat"]["ACC_EKE"].values, *per)
    rows.append(dict(panel="", study=BL, metric="ratio of all-sat to two-sat EKE trend", region="global ocean; ours: ACC band",
                     period="2000-2022", published_value=1.8, published_lo=np.nan, published_hi=np.nan, published_interval="none",
                     our_series="Product 1 / Product 4", our_period=f"{per[0]}-{per[1]}", our_value=t1["val"] / t4["val"],
                     our_lo=np.nan, our_hi=np.nan, our_interval="", p=np.nan,
                     source="Sci. Rep. 15:21754, p. 2 (3.4 for 1993-2022, global), p. 3 (1.8 for 2000-2022, global), p. 5 (2.3, high-EKE regions, 1993-2022)",
                     note=f"ratio of absolute trends {t1['val']:.2f}/{t4['val']:.2f} cm2 s-2 per decade"))
    print("ratio", per, t1["val"] / t4["val"])

# ---------------------------------------------------------------- (b) absolute EKE trends, cm2 s-2 per decade (Zhang et al. 2021)
ZH = "Zhang et al. (2021)"
for key, reg, v, e, pan in (("Hogg_Indian_EKE", "Hogg Indian box 40-150E, 57-44S", 9.1, 3.0, "b"),
                       ("Hogg_Pacific_EKE", "Hogg Pacific box 150E-72W, 62-48S (repository m1 box, not Zhang's 62-40S)", 9.4, 5.3, "")):
    for prod, pn in PRODS:
        add(ZH, "EKE trend (cm2 s-2 per decade)", reg, "1993-2019", v, v - e, v + e, "90 %", pn, "2001-2019",
            tb(SER[prod][key].values, 2001, 2019), "JGR Oceans 126, e2020JC016973, p. 11 (crossover EKE, 27-year trends)",
            "crossover (theirs) vs gridded (ours) EKE", pan)
        add(ZH, "EKE trend (cm2 s-2 per decade)", reg, "1993-2019", v, v - e, v + e, "90 %", pn, "2001-2025",
            m1(prod, key, "2001-2025", rel=False), "JGR Oceans 126, e2020JC016973, p. 11", "", pan)
# 30-degree sector centred at 165W (180-150W), the only significant sector (Zhang p. 10)
B = np.load(OUT / "drake" / "dk12" / "dk12_band_series_allsat.npz"); TB = pd.read_csv(OUT / "drake" / "dk12" / "dk12_band_trends_allsat.csv")
sel = (B["lon"] >= -180) & (B["lon"] < -150); w = TB.set_index("lon").loc[B["lon"][sel], "n_pix"].values.astype(float)
ysec = np.nansum(B["eke"][:, sel] * w, 1) / np.nansum(np.where(np.isfinite(B["eke"][:, sel]), w, 0), 1)
for per in ((2001, 2019), (2001, 2025)):
    add(ZH, "EKE trend (cm2 s-2 per decade)", "30 deg sector 180-150W (ours: between mean SAF and PF)", "1993-2019", 18.9, 18.9 - 11.9,
        18.9 + 11.9, "90 %", "Product 1", f"{per[0]}-{per[1]}", tb(ysec, *per),
        "JGR Oceans 126, e2020JC016973, p. 10 (sector centred at 165W, Fig. 7a)", "ours: dk12 2-deg bins 179W-151W weighted by pixel count", "b")
# Zhang over Hogg's period (table only)
for key, reg, v, e in (("Hogg_Indian_EKE", "Hogg Indian box 40-150E, 57-44S", 10.9, 5.2), ("Hogg_Pacific_EKE", "Hogg Pacific box 150E-72W, 62-48S (repository m1 box, not Zhang's 62-40S)", 12.5, 7.7)):
    for prod, pn in PRODS:
        add(ZH, "EKE trend (cm2 s-2 per decade)", reg, "1993-2012", v, v - e, v + e, "90 %", pn, "2001-2012",
            m1(prod, key, "2001-2012", rel=False), "JGR Oceans 126, e2020JC016973, p. 12", "table only")
# Hogg absolute (table only, AUTHOR CHECK)
for key, reg, v, e in (("Hogg_Indian_EKE", "Indian box", 18.3, 5.1), ("Hogg_Pacific_EKE", "Pacific box", 14.9, 4.1), ("Hogg_Atlantic_EKE", "Atlantic box 35W-10E, 56-46S", 4.0, 3.7)):
    for prod, pn in PRODS:
        add(HOGG, "EKE trend (cm2 s-2 per decade)", reg, "1993-2012", v, v - e, v + e, "90 %", pn, "2001-2025",
            m1(prod, key, "2001-2025", rel=False), HSRC, "table only; crossover EKE about 1.6-1.9 x gridded EKE (repository notes)")

# Pacific box with Zhang's 62-40S limits (Product 1; series from the raw-data step --zhang-series, shipped)
if not ZSER.exists():
    zhang_series()
y = pd.read_csv(ZSER).eke_cm2s2.values
for per in ((2001, 2019), (2001, 2025)):
    add(ZH, "EKE trend (cm2 s-2 per decade)", "Hogg Pacific box 150E-72W, 62-40S (as published by Zhang et al. 2021, p. 11)", "1993-2019", 9.4, 4.1, 14.7, "90 %",
        "Product 1", f"{per[0]}-{per[1]}", tb(y, *per), "JGR Oceans 126, e2020JC016973, p. 11",
        "Zhang's published box; Product 1 only (two-satellite monthly fields not in the container); replaces the former 62-48S repository box", "b")

# ---------------------------------------------------------------- (c) front position, km per decade, positive northward
T = np.load(OUT / "x1" / "x1_tracks.npz"); lon_t = T["lon"]
def pos(field, a, b, y0, y1):
    s = (lon_t >= a) & (lon_t < b); Lt = T[field][:, s].astype(float)
    y = np.nanmean(Lt - np.nanmean(Lt, 0, keepdims=True), 1)
    return tb(y, y0, y1, scale=KM)
CH = "Chambers (2018)"
for reg, (a, b), v in (("90-170E", (90, 170), -29.0), ("160-90W (200-270E)", (-160, -90), 46.0)):
    for per in ((2001, 2016), (2001, 2025)):
        add(CH, "front-position trend (km per decade, northward positive)", reg, "1993-2016", v, np.nan, np.nan,
            "none given for the regional average", "Product 1, half-power latitude", f"{per[0]}-{per[1]}", pos("lat_hp", a, b, *per),
            "Ocean Sci. 14, 105-116, p. 113 (KE half-power point; -29 km per decade 90-170E, 46 km per decade to the north 200-270E)",
            "same half-power idea; theirs along-track CKE, 3-year estimates", "c")
GI = "Gille (2014), after Sokolov and Rintoul (2009)"
for field, nm in (("lat_c", "Product 1, fixed ADT contour"), ("lat_hp", "Product 1, half-power latitude")):
    add(GI, "front-position trend (km per decade, northward positive)", "circumpolar, SSH contours", "1993-2007 (first 15 years)",
        -40.0, np.nan, np.nan, "none (about 60 km southward in 15 years)", nm, "2001-2025", pos(field, -180, 180, 2001, 2025),
        "Phil. Trans. R. Soc. A 372: 20130273, p. 3 (60 km in 15 years) and p. 5 (consistent with 2.9 +- 0.1 mm/yr sea-level rise); "
        "transport-latitude index: no long-term trend, no value (p. 8)", "", "c")
rows.append(dict(panel="", study="Schiavon et al. (2025)", metric="jet position", region="-", period="1995-2015 vs 2080-2100 (models)",
                 published_value=np.nan, published_lo=np.nan, published_hi=np.nan, published_interval="", our_series="", our_period="",
                 our_value=np.nan, our_lo=np.nan, our_hi=np.nan, our_interval="", p=np.nan,
                 source="JGR Oceans 2025 (GFDL-CM4, HadGEM3-GC31-MM, historical vs SSP5-8.5)",
                 note="not compared: climate-model projections, no observed trend; cited only for regional, topographically anchored response"))

D = pd.DataFrame(rows)
TAB = PA.out("m1") / "m1_previous_studies.csv"
D.round(3).to_csv(TAB, index=False)
print("wrote", TAB)
