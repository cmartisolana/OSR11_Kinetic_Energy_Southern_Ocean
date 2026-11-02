# Antarctic-slope transport index after Hogg et al. (2015), 20 Oct 2026

Scripts: `plan1019/h0_prepare_bathymetry_80S_chunked.py` (GEBCO 2026 80–45°S → 1/8°, `data/bathymetry_80s.nc`),
`plan1019/h0b_hogg_pointwise.py` (index), `drake/h2_hogg_vs_insitu.py` (unchanged, in situ comparison), `h0b_diagnostics.py`.
Outputs: `h0b_allsat.csv`, `h0b_twosat.csv`, `h2_h0b_allsat.csv`, `h2_h0b_twosat.csv`.

Formula checked against the Hogg et al. (2015) PDF (Sect. 2.3): T = −3.69 (h − h_f) − 1.11 h_f [Sv, h in cm]. T is a transport **anomaly** with model-derived ("preliminary") factors, not an absolute transport.

Region: ocean 1000–3000 m south of 60°S, contiguous with Antarctica (33,701 all-sat / 8,530 two-sat points).

Why h0b and not h0: DUACS vNov2024 fills ADT under sea ice (80 % of slope points in Aug 2010 have finite ADT with ice_frac = 1), so ice points are masked (ice_frac ≤ 0.10). With h0's 50 % coverage threshold only 37/300 months survive. h0b follows Hogg et al.: per-point bias + calendar-month climatology removed, then all sampled points averaged; no coverage threshold (≥ 100 points, ≥ half of depth bins).
Coverage used: Jan–Mar 0.50, Jul–Sep 0.05 (median 4 of 12 30° bands sampled in Jul–Sep).

| | All-sat | Two-sat |
|---|---|---|
| Months with index | 296/300 | 291/300 |
| T, monthly s.d. / annual-mean s.d. | 4.7 / 2.4 Sv | 4.2 / 2.1 Sv |
| T trend 2001–2025 (trend_bcz) | −1.9 [−3.3, −0.6] Sv/dec, p = 0.02 | −1.7 [−3.3, −0.1] Sv/dec, p = 0.03 |
| of which −1.11·h_f | −1.8 | −1.6 |
| of which −3.69(h − h_f) | −0.1, p = 0.83 | −0.1, p = 0.81 |
| h trend | +1.7 [0.8, 2.6] cm/dec, p = 0.003 | +1.5 [0.2, 2.8] cm/dec, p = 0.03 |
| h trend, Jan–Mar only (best coverage) | +1.5 cm/dec, p = 0.047 | +1.2 cm/dec, p = 0.06 |
| r(T) with in situ 0–760 m total, per transect / annual 2005–2019 | 0.20 / 0.35 | 0.21 / 0.37 |
| r(T) with in situ north / south of 58.5°S, annual | 0.24 / 0.02 | 0.28 / −0.02 |
| r(T) all-sat vs two-sat, monthly | 0.96 | |

Caveats:
1. The trend is entirely the −1.11·h_f term, i.e. the slope sea-level rise. DUACS ADT includes global-mean sea-level rise; the Hughes et al. (2014) factors come from Boussinesq models without it. The slope rise (+1.7 cm/dec) is smaller than satellite GMSL rise over 2001–2025 (≈3–4 mm/yr; not computed here). Relative to GMSL, h would decline and the sign of the T trend would reverse. The trend is therefore not robust to the GMSL treatment.
2. Winter values sample ~5 % of the slope in ~4 longitude bands.
3. Agreement with in situ data is weak (annual r ≈ 0.35, n = 15, not significant); the in situ series is 0–760 m in Drake Passage only.

## Global-mean sea level removed (20 Oct 2026)

GMSL: AVISO reference series (`MSL_Serie_MERGED_Global_AVISO_GIA_Adjust_Filter2m.txt`, TOPEX/Jason, annual and semiannual signals removed), averaged to months → `data/gmsl_aviso_monthly_2001_2025.csv`. The −0.3 mm/yr GIA correction is taken back out (DUACS ADT has none): GMSL trend 3.4 mm/yr (3.7 with GIA). It is subtracted from ADT at every point and month (`H0_GMSL=... h0b_hogg_pointwise.py`), so h is slope sea level relative to the global mean.
Outputs: `h0b_{allsat,twosat}_noGMSL.csv`, `drake/out/h2_h0b_{allsat,twosat}_noGMSL.csv`. Usage: `python h0b_diagnostics.py h0b_allsat_noGMSL h0b_twosat_noGMSL`.

| | All-sat | Two-sat |
|---|---|---|
| h trend (relative to GMSL) | −1.6 [−2.4, −0.9] cm/dec, p = 0.003 | −1.8 [−2.7, −0.9] cm/dec, p = 0.007 |
| T trend 2001–2025 | **+1.8 [+0.3, +3.3] Sv/dec, p = 0.03** | **+1.9 [+0.3, +3.6] Sv/dec, p = 0.03** |
| of which −1.11·h_f / −3.69(h − h_f) | +1.8 / −0.0 (p = 0.95) | +2.0 / −0.1 (p = 0.86) |
| Jan–Mar-only h trend (×−1.11) | −2.1 cm/dec, p = 0.10 (+2.3 Sv/dec) | −2.1 cm/dec, p = 0.07 (+2.4 Sv/dec) |
| T s.d., monthly / annual | 4.8 / 2.5 Sv | 4.4 / 2.4 Sv |
| r(T) with in situ total, per transect / annual | 0.22 / 0.36 | 0.22 / 0.34 |
| r(T) with in situ north / south, annual | 0.37 / −0.11 | 0.40 / −0.16 |

Reading: removing GMSL reverses the sign of the trend (−1.9 → +1.8 Sv/dec). The trend comes entirely from the long-period term (factor −1.11, which Hogg et al. call preliminary and say may be smaller still at longer periods). The best-sampled season gives the same sign but p ≈ 0.1. Interannual agreement with in situ is essentially unchanged (weak).

## Correction to the in situ correlations above (20 Oct 2026)

`h2_hogg_vs_insitu.py` merges two per-transect tables on (year, month), duplicating rows in months with several transects (519 rows from 248 transects). The r values in both tables above are therefore wrong. Corrected by `drake/h2b_hogg_vs_insitu_fixed.py` (one row per transect) → `drake/out/h2b_*_noGMSL.csv`.
Period: 245 matched transects (241 two-sat), Feb 2005–Dec 2019, i.e. the whole in situ record; annual means over 2005–2019 (n = 15).

| T (GMSL removed) vs in situ 0–760 m | All-sat | Two-sat |
|---|---|---|
| total, per transect | 0.16 | 0.16 |
| total, annual (matched months / full-year index) | 0.16 / 0.07 | 0.13 / 0.00 |
| north of 58.5°S, per transect | 0.22 | 0.23 |
| north, annual (matched / full year) | 0.36 / 0.42 | 0.40 / 0.41 |
| south of 58.5°S, annual (matched / full year) | −0.26 / −0.40 | −0.32 / −0.45 |

No annual correlation is significant (p ≥ 0.10). Per-transect p values are not given: transects in the same month share one index value, and the samples are autocorrelated.
Other users of `d1b_insitu_split.csv` (d1b_gvtrack, d1c check, fig_v3, fig_drake) merge with monthly series (many-to-one) and are not affected.

## Trends over the in situ period, 2005–2019 (GMSL removed)

- In situ 0–760 m total, annual means (n = 15): −0.9 [−9.6, +7.8] Sv/dec, p = 0.70 (mean 75 Sv).
- T, monthly, h_f fitted on 2001–2025: all-sat +2.7 [−0.8, +6.2] Sv/dec, p = 0.14; two-sat +2.7 [−1.5, +6.8], p = 0.13.
- h: all-sat −2.2 cm/dec (p = 0.03), two-sat −2.2 (p = 0.08). Converted with −1.11 this is about +2.4 Sv/dec; with −3.69, about +8 Sv/dec.
Neither trend is significant over the overlap, and the CIs overlap. A trend difference cannot be tested with 15 years.

## Same comparison without removing GMSL (corrected h2b), 2005–2019

| T (GMSL kept) vs in situ 0–760 m | All-sat | Two-sat |
|---|---|---|
| total, per transect / annual | 0.16 / 0.21 | 0.17 / 0.22 |
| north, annual | 0.27 | 0.32 |
| south, annual | −0.13 | −0.18 |
| T trend 2005–2019 | −1.1 [−4.6, +2.3] Sv/dec, p = 0.46 | −1.1 [−5.2, +2.9], p = 0.50 |

In situ over the same period: −0.9 [−9.6, +7.8] Sv/dec, p = 0.70. Correlations barely change with the GMSL treatment, because GMSL is mostly a trend and the annual cycle is removed. Over 2005–2019 the GMSL-kept index matches the in situ trend in sign, but none of the trends is significant.
