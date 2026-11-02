# Supporting Information for "Eddies and mean flow of the Antarctic Circumpolar Current intensified together at the surface, 2001–2025"

Martí-Solana et al., State of the Planet, Ocean State Report 11, Chapter 2

<!-- v3.1 (19 Oct 2026): revised from SI_v3.md per revision_decisions_v3.md (D binding; A1, A4, A5, C-M6a, C-M8, C-M9, C-S1, C-S9, HI S4), numbercheck_v3.md items 6, 8–11, 13–17, 19, 20 and review_PO.md. Number sources (precedence brief §0 > drake_results > plan1019 > AC v2 > v23 > R2 "Revised trend test" > earlier R2): inputs/numbers_v3.md; inputs/drake_results.md; inputs/plan1019_results.md; inputs/R2_revision_results.md; inputs/AC_reply_*_v2.md; inputs/manuscript_v23.md; numbercheck_v3.md; checks/d1c_out.txt; CSVs cited per table. -->

**Conventions.** Product 1: all-satellite SEALEVEL_GLO_PHY_L4_MY_008_047 (DT2024, 1/8°), used in the main text; Product 4: two-satellite SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057 (1/4°). Unless stated, trends are Theil–Sen slopes of deseasonalised monthly series, 2001–2025, in percent of the period mean per decade, with 95 % intervals in brackets and p from the test of Sect. S2. Contrasts (trends of differences of mean-normalised series) are in percentage points ("points") per decade. Significance means p < 0.05; p = 0.050 is "at the threshold". In figures, open symbols mark p ≥ 0.05, including values that round to 0.050. Post hoc tests are labelled.

---

## S1 Two-satellite product

**Purpose.** Constellation changes over 2001–2025 may affect mapped EKE trends; a constant pair of altimeters shows which conclusions depend on the product.

**Method.** Product 4 diagnostics use the same definitions, regions, ice mask (ice fraction > 0.10 in any month), periods and test. Grid spacing is 1/8° (Product 1) and 1/4° (Product 4); effective resolution is of the order of 100–150 km in wavelength at ACC latitudes (Ballarotta et al., 2019) [CHECK: DT2024 effective resolution from product QUID]. Product differences are tested as the trend of the difference of mean-normalised series (Product 1 minus Product 4). The eddy atlas (Product 5) exists for Product 1 only.

**Results.** In Drake Passage (Table S1) upstream KE and EKE increase in both products; upstream speed is at the threshold in Product 1 (p = 0.050) and has p = 0.056 in Product 4. Box EKE is at the threshold in Product 1, about half as large and not significant in Product 4 (difference +4.7 points per decade [−1.2, 10.7], p = 0.09). Downstream EKE is significant in Product 1 only.

**Table S1.** Drake Passage trends, 2001–2025. Box 80–52°W, 63–54°S (ice-masked); upstream 100–80°W and downstream 55–35°W within the ACC band (mean ADT −0.6 to +0.2 m). Absolute trends in cm² s⁻² per decade (speed: cm s⁻¹ per decade), period mean in parentheses. "Box, no mask": the box without the ice mask (Sect. S5).
<!-- /home/claude/work/drake/out/dk1/dk1_trends.csv (period 2001-2025); relative values and p as in drake_results.md §1; no-mask rows from dk1_trends.csv *_box_nomask_* (2001-2025), as previously quoted in SI v3 S5 (EKE, KE P1; EKE P4; contrasts) -->

| Region | Quantity | Product 1, % [CI], p | Product 1, absolute (mean) | Product 4, % [CI], p | Product 4, absolute (mean) |
|---|---|---|---|---|---|
| Box | KE | +2.4 [−2.5, 7.4], 0.29 | +7.8 [−7.9, 23.4] (319) | +1.4 [−3.4, 6.2], 0.52 | +4.2 [−10.4, 18.8] (301) |
| | EKE | +8.4 [0.2, 16.6], 0.050 | +21.7 [0.6, 42.8] (258) | +3.8 [−2.5, 10.1], 0.22 | +8.3 [−5.6, 22.3] (221) |
| | Along-stream speed | −0.1 [−3.2, 3.0], 0.91 | −0.02 [−0.51, 0.47] (15.7) | 0.0 [−3.1, 3.0], 0.98 | 0.00 [−0.47, 0.47] (15.4) |
| | EKE − KE, points | +6.0 [0.7, 11.3], 0.04 | — | +2.4 [−1.7, 6.4], 0.22 | — |
| Box, no mask | KE | +2.0 [−2.8, 6.9], 0.34 | +5.3 [−7.5, 18.1] (263) | +0.9 [−3.7, 5.6], 0.63 | +2.3 [−9.4, 14.0] (251) |
| | EKE | +7.1 [−0.9, 15.1], 0.08 | +15.2 [−2.0, 32.4] (214) | +2.8 [−3.6, 9.3], 0.35 | +5.2 [−6.6, 17.0] (183) |
| | Along-stream speed | +0.2 [−1.5, 1.9], 0.74 | +0.03 [−0.21, 0.26] (14.1) | +0.2 [−1.5, 1.8], 0.76 | +0.02 [−0.20, 0.25] (13.9) |
| | EKE − KE, points | +5.0 [0.5, 9.6], 0.04 | — | +1.9 [−1.5, 5.2], 0.23 | — |
| Upstream | KE | +12.8 [6.8, 18.7], 0.007 | +19.2 [10.3, 28.2] (150) | +10.8 [5.8, 15.8], 0.007 | +16.2 [8.7, 23.8] (150) |
| | EKE | +15.3 [6.2, 24.4], 0.013 | +22.5 [9.1, 36.0] (147) | +11.5 [4.7, 18.3], 0.007 | +15.6 [6.4, 24.7] (135) |
| | Along-stream speed | +3.9 [0.3, 7.5], 0.050 | +0.33 [0.02, 0.63] (8.4) | +3.6 [0.0, 7.3], 0.056 | +0.31 [0.00, 0.61] (8.4) |
| | EKE − KE, points | +3.8 [−0.3, 7.9], 0.09 | — | +1.7 [−1.1, 4.6], 0.21 | — |
| Downstream | KE | +4.3 [−1.8, 10.3], 0.13 | +17.3 [−7.3, 41.8] (406) | +3.3 [−2.4, 9.0], 0.21 | +12.5 [−9.4, 34.5] (385) |
| | EKE | +8.0 [3.0, 13.0], 0.010 | +36.3 [13.5, 59.1] (455) | +3.8 [−1.5, 9.1], 0.13 | +14.2 [−5.6, 34.0] (375) |
| | Along-stream speed | +2.7 [−4.0, 9.3], 0.28 | +0.42 [−0.62, 1.47] (15.7) | +2.7 [−4.0, 9.3], 0.29 | +0.42 [−0.62, 1.46] (15.7) |
| | EKE − KE, points | +3.3 [−1.0, 7.7], 0.12 | — | +0.1 [−4.9, 5.1], 0.96 | — |

Eastward wind stress (WIND_GLO_PHY_L4_MY_012_006): box +8.7 [1.4, 15.9], p = 0.03; upstream +9.8 [1.3, 18.4], p = 0.03; downstream +9.4 [−1.5, 20.2], p = 0.09. Both products compare similarly with the in situ record (Table S2).

**Table S2.** Comparison with the LM Gould ADCP record (2005–2019) and section diagnostics.
<!-- drake_results.md §2–3; /home/claude/work/drake/out/dk3/dk3_summary.csv (north/south bins); dk2/dk2_vs_insitu.csv (Product 4 RMS); dk4/dk4_trends.csv (Product 4 CIs) -->
<!-- Product 4 proxy correlations 0.52 / 0.66: /home/claude/work/v3/checks/d1c_out.txt ("twosat: n=248 r_transect=0.523 r_annual=0.658"), produced by /home/claude/work/v3/checks/d1c_twosat_proxy_vs_insitu.py (inputs out/d1/d1b_insitu_split.csv, drake/out/dk4/dk4_section_series.csv), which reproduces the Product 1 values of d1b_summary.csv (0.535 / 0.670) with the same procedure -->

| Quantity | Product 1 | Product 4 |
|---|---|---|
| Mean cross-section velocity profile vs ADCP (70–118 m), r | 0.98 | 0.97 |
| Velocity per 25 km bin and transect (includes the mean cross-section structure), r (n) | 0.72 (8187) | 0.72 (8189) |
| — bins north / south of 58.5°S, r | 0.70 / 0.75 | 0.69 / 0.73 |
| Velocity anomalies (bin mean removed), r | 0.61 | 0.60 |
| RMS velocity difference | 15 cm s⁻¹ | 16 cm s⁻¹ |
| EKE north of 58.5°S, r per transect / annual | 0.67 / 0.89 | 0.65 / 0.85 |
| EKE, all ice-free points, r annual | 0.92 | 0.88 |
| Band position vs in situ half-transport latitude, r per transect / annual | 0.50 / 0.71 | 0.48 / 0.65 |
| — mean offset; RMS difference | 16 km; 114 km | 5 km; 117 km |
| Surface transport proxy (GV track) vs in situ part north of 58.55°S, r per transect / annual | 0.54 / 0.67 | 0.52 / 0.66 |
| Section half-power trend 2001–2025, km per decade | −10 [−36, +16], p = 0.43 | −8 [−44, +29], p = 0.65 |
| EKE north of 58.5°S, trend 2001–2025 | +7.4 [−5.6, 20.4], p = 0.24 | +1.8 [−8.5, 12.0], p = 0.68 |
| Surface transport proxy (GV track), trend 2001–2025 | +3.9 [−2.1, 9.9], p = 0.17 | +4.2 [−1.8, 10.2], p = 0.15 |

Over the ACC band (Table S3) all trends have the same sign in both products. Mean-flow trends agree closely; EKE trends are about half as large in Product 4, mainly in the sub-monthly part, and diverge mainly after about 2018 (Fig. S1), consistent with the time-varying constellation (Product 1 QUID; Hargous et al., 2026). No tested product difference is significant (p ≥ 0.09; p ≥ 0.12 on a common 1/4° grid). Product 4 relative KE and EKE trends are similar (+3.6 and +3.8 %); EKE − KE is not significant in either product.

**Table S3.** Circumpolar diagnostics, 2001–2025 unless stated.
<!-- v23 §3.1; plan1019 §2; R2 A4 table, A7 table and "Revised trend test"; AC-A v2 (proxy); /home/claude/work/out/m1/m1_trends.csv (EKE absolute: allsat 19.83 [8.11, 31.54], twosat 8.06 [2.82, 13.30] cm² s⁻² per decade, given here per year as numbercheck item 6; common-grid minimum p 0.12 = diff_all-two|box35_65_LF4|2001-2025, numbercheck item 10) -->

| Quantity | Product 1 | Product 4 | Difference, points per decade |
|---|---|---|---|
| ACC KE, % (absolute, cm² s⁻² yr⁻¹) | +4.9 (1.36 [0.86, 1.86]), p = 0.003 | +3.6 (0.97), p ≤ 0.013 | — |
| ACC EKE, % (absolute, cm² s⁻² yr⁻¹) | +8.0 [3.3, 12.8] (1.98 [0.81, 3.15]), p = 0.017 | +3.8 [1.3, 6.3] (0.81 [0.28, 1.33]), p = 0.020 | +4.3 [−2.2, 10.8], p = 0.12 |
| ACC EKE, 2001–2020 | +6.4 [0.9, 11.9], p = 0.04 | +3.2 [−0.1, 6.5], p = 0.07 | +3.2, p = 0.37 |
| Low-frequency EKE ½\|Uₘ − ū\|² | +5.6 | +3.2 | — |
| Low-frequency EKE, common 1/4° grid | +5.5, p = 0.003 | +3.2, p = 0.013 | +2.3, p = 0.18 |
| Sub-monthly EKE ½σ²ₘ | +12.0, p = 0.066 | +5.0, p = 0.040 | — |
| Along-stream speed | +2.0, p = 0.003 | +1.8, p ≤ 0.013 | — |
| EKE − KE, points per decade | +3.0 [−0.1, 6.1], p = 0.060 | +0.2 [−1.2, 1.5] (earlier test, p = 0.74) | — |
| Circumpolar surface transport proxy | +1.6 [0.9, 2.3], p = 0.003 | +1.5 [0.9, 2.1], p = 0.003 | — |
| Wake boxes minus rest of band, EKE, points per decade | +0.1 [−3.0, 3.2], p = 0.95 | −0.3, p = 0.80 | — |
| Box 35–65°S, EKE (ice-masked) | +6.6 [1.5, 11.6], p = 0.03 | +2.3 [0.2, 4.5], p = 0.04 | +4.5 [−1.4, 10.4], p = 0.10 |
| Hogg et al. (2015) Indian box, EKE | +9.0 [4.9, 13.0], p = 0.007 | +4.8 [1.7, 7.9], p = 0.013 | +4.5 [−5.8, 14.7], p = 0.20 |
| Hogg et al. (2015) Pacific box, EKE | +8.6 [0.9, 16.4], p = 0.04 | +4.9 [0.5, 9.4], p = 0.05 (at the threshold) | +4.3 [−0.8, 9.4], p = 0.09 |
| Hogg et al. (2015) Atlantic box, EKE | +8.8 [−0.1, 17.6], p = 0.056 | +5.2 [−0.2, 10.7], p = 0.06 | +3.9, p = 0.17 |
| PAR jet-core speed, whole meander / Lower Ridge, m s⁻¹ per decade | +0.019 / +0.026, p = 0.01 | +0.017 / +0.024, p = 0.01 | — |
| Along-stream bins (Sect. S9) and cross-bin r(EKE, speed) | main text Sect. 3.4 | [NUMBER NEEDED: Product 4 bin statistics and cross-bin r; not run] | — |
| Sector EKE, absolute (cm² s⁻² yr⁻¹), Pacific / Atlantic / Indian | +2.35 / +2.14 / +1.33 | +1.21 / +0.90 / +0.23 | — |
| Sector KE, absolute (cm² s⁻² yr⁻¹), Pacific / Atlantic / Indian | +1.87 / +1.27 / +0.92 | +1.53 / +0.89 / +0.41 | — |

**Energy-equivalent comparison.** EKE and KE scale with velocity squared, the along-stream speed and surface transport proxies with velocity, so a small relative velocity trend corresponds to about twice that trend in energy. The ratio of the relative EKE and circumpolar proxy trends is 5.0 (Product 1) and 2.5 (Product 4); in energy-equivalent terms the proxy rises at about 40 % of the EKE rate in Product 1 (3.2 against 8.0 % per decade) and at a comparable rate in Product 4 (3.0 against 3.8); the along-stream speed rises at about the KE rate (4.0 against 4.9; Product 4, 3.6 against 3.6). In Drake Passage the energy-equivalent northern proxy trend is comparable with box EKE in Product 1 (7.8 against 8.4) and larger in Product 4 (8.4 against 3.8; neither significant). This describes the surface partition between eddies and mean flow and is not a saturation test.
<!-- Energy equivalent = 2 × relative trend of the velocity-like quantity (revision_decisions_v3 A1; review_PO M2): proxy 2 × 1.6 = 3.2 (P1), 2 × 1.5 = 3.0 (P4); speed 2 × 2.0 = 4.0 (P1), 2 × 1.8 = 3.6 (P4); Drake northern proxy 2 × 3.9 = 7.8 (P1), 2 × 4.2 = 8.4 (P4). Ratios 8.0/1.6 = 5.0, 3.8/1.5 = 2.5 (numbercheck B, MS 3.6). -->

In the sectors (Table S3), Pacific KE and EKE are significant in both products (p = 0.007–0.033), Indian EKE in Product 1 only (p = 0.010) and Atlantic EKE in neither (Product 1 p = 0.053); Product 4 KE is not significant in the Atlantic (p = 0.06) or Indian (p = 0.08) sector; along-stream speed and proxy trends agree between products (proxy within 0.2 points).
<!-- R2 A4 sector table (absolute values); revised p from R2 "Revised trend test" -->

**Figure S1.** Results from the all-satellite (Product 1) and two-satellite (Product 4) products. (a) Mean KE, 2001–2025, Product 1 (light blue: sea ice in any month). (b) ACC-band KE and EKE, 12-month running means in percent of the 2001–2025 mean; solid: Product 1, dashed: Product 4. (c) ACC-band trends, 2001–2025, in percent of the mean per decade with 95 % confidence intervals: KE, EKE, its low-frequency and sub-monthly parts, along-stream speed and surface transport proxy. (d) The same for the sectors (speed: along-stream speed; transport: surface transport proxy). Filled symbols: p < 0.05; open symbols: p ≥ 0.05.

**Figure S2.** Drake Passage in both products. (a) Mean cross-track velocity along the LM Gould track from the ADCP (70–118 m, orange) and from Products 1 (green) and 4 (grey), 2005–2019. (b) Annual EKE and (c) annual transport north of 58.5°S, in percent of the 2005–2019 mean: in situ (orange, 2005–2019), Product 1 (colour) and Product 4 (grey), 2001–2025. (d) Annual position of the energetic band along the section: in situ half-transport latitude and half-power latitude in both products. (e) Trends 2001–2025 of KE, EKE, along-stream speed and eastward wind stress in the Drake box and the upstream and downstream regions, both products (filled: p < 0.05; open: p ≥ 0.05, including values that round to 0.050). Grey shading: overlap with the in situ record.

The Drake Passage diagnostics of both products are shown in Fig. S2.

**Conclusion.** Signs and mean-flow trends are product-independent; the EKE magnitude differs by about a factor of two (not significantly), so statements resting on it, including the eddy–mean-flow contrast, depend on the product.

---

## S2 Trend-test calibration

**Purpose.** Monthly series are strongly autocorrelated, and common trend tests then reject too often.

**Method.** The statistic is the Theil–Sen slope over the exact OLS-slope standard deviation under a fitted autoregressive (AR) noise model, bias-corrected by bootstrap including the trend removal (Kilian, 1998); non-stationary fits are shrunk to a maximum root of 0.99. p and the 95 % interval are calibrated against 300 null series processed identically (smallest p 0.0033, reported as 0.003). Monthly series use AIC-selected and AR(12) orders, reporting the larger p and the union of intervals; the two need not agree exactly at the threshold (Drake box EKE: interval excludes zero, p = 0.050). Annual series use AR(1) with iterated bias correction. For 1° maps the calibration is pooled over 10 persistence classes (p from |z|/κ, κ = max(q95/1.96, q99/2.576) of the class null). The false discovery rate (FDR) is controlled at 10 %, following Wilks (2016; α_FDR = 2α_global for correlated fields).

**Results.** AIC is liberal for the Drake-like AR(12) process (dropping small high-lag coefficients) and AR(12) for pure AR(1); the maximum rule takes the more conservative (Table S4). Fifteen annual values cannot be tested at the nominal level (9.5 %); no in situ trend is significant under either test. With 300 null series, p near 0.05 has a Monte Carlo uncertainty of about ±0.01, so p of about 0.03–0.07 is marginal. Of 585 recomputed trends, fewer are significant than under the earlier test in every group except the contrasts (3 → 4, including two new tests); e.g. KE/EKE terms 145 → 107. ACC KE, EKE, along-stream speed and low-frequency EKE remain significant; ACC EKE − KE (p = 0.060), Atlantic EKE (p = 0.053) and Product 1 sub-monthly EKE (p = 0.066) do not.

**Table S4.** False-positive rate at nominal 5 % (number of series in parentheses).
<!-- R2 "Revised trend test" table; earlier methods from R2 "Fixed choices" and "Statistics check" -->

| Process (n = 300 months unless stated) | Earlier calibrated sieve | New, AIC order | New, AR(12) |
|---|---|---|---|
| AR(1), φ = 0.8 | 10.7 % | 5.3 % (300) | 8 % (100) |
| AR(1), φ = 0.9 (independent check) | — | 3.5 % (200) | — |
| AR(12) fitted to Drake-box EKE | — | 9.5 % (200) | 3 % (100) |
| Annual AR(1), φ = 0.5, n = 25 | — | 5.0 % (200) | — |
| Annual AR(1), φ = 0.5, n = 15 | — | 9.5 % (200) | — |

Rejected alternatives: moving-block bootstrap 13–18 %; uncalibrated AR sieve 7–13 %; Hamed–Rao 25–33 %. An independent re-implementation matched the analytic slope standard deviation to Monte Carlo within 0.3 %.

**Pre-specified and post hoc.** Fixed before any trend was computed: half-power parameters (Sect. S3); ACC band and ice mask; wake boxes, including Drake upstream and downstream; wake test, along-stream analysis, meander gates and monthly covariability tests. Post hoc or added after results were seen: the 400 cm² s⁻² Drake split (Sect. S3); the branch-start check (Sect. S4); direct tests of EKE − KE contrasts; corrections of the box-free covariability design (Sect. S9).

**Conclusion.** Monthly trends with p ≤ 0.02 are unlikely to be artefacts of autocorrelation; p of about 0.03–0.07 is marginal; 15 annual values cannot be tested at the nominal level.

---

## S3 Half-power method: parameters, detectability and the Drake Passage section

**Purpose.** A bound on the displacement of the energetic band is meaningful only if the parameters do not drive it and the detectable shift is known.

**Method.** Following Chambers (2018), on monthly-mean KE: at each longitude, envelopes extend from peaks above a KE threshold to where KE falls to a fraction of the peak (minimum width imposed); the half-power latitude is the median of KE integrated across the envelope. Baseline: 200 cm² s⁻² (the rounded median, 195 cm² s⁻², of the 2001–2025 mean KE over the ACC band), 25 %, 100 km. As a KE-weighted median it also responds to unequal intensification of the flanks (Sect. S4), so it bounds the displacement of the energy distribution, not of a jet axis. Monthly envelopes are matched to the climatological envelope; longitudes matched in ≥ 90 % of months are kept (900, about 31 % of the band). Sensitivity: thresholds 150–400 cm² s⁻², fractions 20–50 %, widths 50–300 km. Detectability: shifts imposed on climatological profiles plus observed anomalies.

**Results.**
- Every parameter set gives p ≥ 0.16 for the circumpolar trend (revised test; Table S5, Fig. S3a).
- None of the 33 bins of 5° reaches p < 0.05 (revised test).
- Over the regions of earlier studies, trends are small and not significant (km per decade, 95 % intervals): 90–170°E −1.6 (−12.9 to +9.6) and 160–90°W +5.6 (−4.9 to +16.2), where Chambers (2018) found −29 and +46 over 1993–2016; 75–150°E −1.5 (−9.4 to +6.5), where Kim and Orsi (2014) reported a southward drift over 1992–2011; 150–80°W (Pacific–Antarctic Ridge meander) +5.9 (−5.3 to +17.1), against about +3 over 1993–2023 (Liu et al., 2026).<!-- outputs/a2/a2g_literature_regions.csv -->
- Over 25 years, a 25 km shift is detected in 100 % of realisations and a 10 km shift in 44 %; about 70 % of a shift is recovered.
- Detection frequency is 25–45 % at about 170–150°W and 115–75°W (including the Drake upstream region), which are excluded.
- The unmatched median latitude of the discussion version (at most 27 % power for a 50 km shift; −0.002 to −0.007° yr⁻¹ from amplitude growth alone) was replaced.
- Three-year blocks (8 values, OLS as in Chambers, 2018): ACC p = 0.17; Indian sector p = 0.049, a single borderline value.

**Table S5.** Circumpolar half-power trend per parameter set, 2001–2025 (revised test; km per decade, positive northward).
<!-- outputs/retest/retest_summary.csv keys <set>|ACC|lat_hp; slope in degrees per year × 1111.95 -->

| Parameter set | Trend | 95 % interval | p |
|---|---|---|---|
| Baseline: 200 cm² s⁻², 25 %, 100 km | −0.5 | −4.2 to +3.1 | 0.74 |
| Threshold 150 cm² s⁻² | +0.1 | −5.0 to +5.1 | 0.98 |
| Threshold 250 cm² s⁻² | −0.7 | −4.2 to +2.9 | 0.65 |
| Threshold 300 cm² s⁻² | −0.6 | −4.3 to +3.1 | 0.68 |
| Threshold 400 cm² s⁻² | −0.7 | −2.9 to +1.5 | 0.48 |
| Relative threshold (monthly ACC-band mean KE) | +0.9 | −6.1 to +8.0 | 0.76 |
| Fraction 20 % | 0.0 | −4.2 to +4.2 | 0.99 |
| Fraction 33 % | −0.9 | −4.6 to +2.7 | 0.56 |
| Fraction 50 % | −1.0 | −4.8 to +2.8 | 0.53 |
| Minimum width 50 km | −1.5 | −3.8 to +0.7 | 0.16 |
| Minimum width 200 km | +1.8 | −5.9 to +9.6 | 0.58 |
| Minimum width 300 km | −1.9 | −8.2 to +4.4 | 0.50 |

**Withdrawn envelope width.** The broadening reported in the discussion version is withdrawn. The threshold-defined width increased (+0.47 km yr⁻¹, earlier test p = 0.002), but an imposed amplitude increase with a fixed jet gave +0.43 to +1.6 km yr⁻¹. Amplitude-independent widths narrowed slightly (half-maximum −0.10 km yr⁻¹, p = 0.017; at 25 % −0.12 km yr⁻¹, p = 0.047), but 5 km of resolution sharpening gives −0.05 to −0.25 km yr⁻¹, and the half-maximum width has ≤ 35 % power even for large broadening. No width change can be inferred.

**Drake Passage section and the 400 cm² s⁻² split (post hoc).** Along the mean LM Gould track (55.1–62.3°S) the baseline profile has one envelope, 55.1–59.3°S (half-power latitude 57.0°S; 99.7 % of months): the Subantarctic Front (SAF) and Polar Front (PF) merge because KE stays above 200 cm² s⁻² between them, so the baseline is the energy median of the merged band, not a front. Points with ice fraction > 0.10 are masked each month, so points south of 58.5°S enter only when ice-free; the in situ half-transport latitude uses the climatological envelope, so no shared variability is built into the correlation. Mapped KE may be damped near the ice edge, so the profile was also truncated at 58.55°S (the ice-free limit). The envelope then spans 55.1–58.55°S (half-power latitude 56.7°S; 97 % of months), and its trend, −10.3 km per decade (−33.6 to +13.1; p = 0.36; Product 4 −8.4, −28.4 to +11.6; p = 0.39), is unchanged from the baseline. Against the in situ half-transport latitude within the truncated envelope, r = 0.50 per transect and 0.50 for annual means (Product 4: 0.49 and 0.64), with a mean offset of 7 km and an RMS difference of 91 km; the annual correlation is lower than with the full profile (0.71), so the band-position comparison (Sect. S5) depends partly on the seasonally ice-covered part of the section.<!-- source: drake/dk2_section_fronts.py with DK2_LATMIN=-58.55 -> drake/out/dk2_trunc5855/dk2_trends.csv, dk2_vs_insitu.csv (run 19 Oct 2026) --> The post hoc 400 cm² s⁻² run (inside the pre-defined range; Table S6) detects the PF too rarely for a trend; the SAF trend (at the threshold in Product 1) rests on 55 % of months.

**Table S6.** Section half-power at 400 cm² s⁻² (post hoc), 2001–2025.
<!-- drake_results.md §2; /home/claude/work/drake/out/dk2_thr400/dk2_trends.csv and dk2_vs_insitu.csv -->

| | SAF, Product 1 | SAF, Product 4 | PF, Product 1 | PF, Product 4 |
|---|---|---|---|---|
| Climatological half-power latitude | 56.2°S | 56.0°S | 58.3°S | 58.4°S |
| Months detected | 55 % | 61 % | 31 % | 35 % |
| Trend, km per decade | −12.6 [−25.0, −0.2], p = 0.050 | −10.3 [−21.0, +0.5], p = 0.06 | not computed | not computed |
| r vs in situ, per transect / annual | 0.47 / 0.71 | 0.45 / 0.76 | 0.48 / 0.38 | 0.31 / 0.35 |
| Mean offset from in situ | −39 km | −35 km | −13 km | −12 km |

**Figure S3.** Half-power analysis. (a) Circumpolar trend of the half-power latitude, 2001–2025, for the baseline parameters and for each alternative (threshold in cm² s⁻²; relative threshold: the monthly ACC-band mean KE; fraction: of the peak KE; minimum width of the envelope), in km per decade (positive northward). (b) Trends in 5° longitude bins (bin centre; bins with too few matched longitudes are omitted). (c) Percentage of realisations in which an injected shift, growing linearly to the value on the x axis over 2001–2025, is detected at p < 0.05 (dashed: 5 %). (d) Median recovered shift (dashed: 1:1). Open symbols: p ≥ 0.05; 95 % confidence intervals. [AUTHOR CHECK: panels (a, b) were drawn with the earlier test; the revised p-values are in Table S5 and the text.]

**Figure S4.** Half-power analysis in the all-satellite (Product 1) and two-satellite (Product 4) products, 2001–2025: trends of (a) the half-power latitude (km per decade, positive northward) and (b) the KE at the envelope peak (percent of the mean per decade), circumpolar and by sector, with 95 % confidence intervals. Filled symbols: p < 0.05.

The half-power trends of both products are compared in Fig. S4.

**Conclusion.** Within the tested ranges the circumpolar bound (no displacement of the energy distribution larger than about 25 km at the retained longitudes) does not depend on the parameters. At Drake Passage the baseline measures a merged band; its post hoc separation into fronts is too intermittent for a trend.

---

## S4 Half-power meander gates

**Purpose.** Referee C asked whether meander flexing, a step of the Thompson and Naveira Garabato (2014) mechanism, is measurable from half-power positions; two gates fixed before any geometry trend decide where.

**Method.** Four windows from bathymetry and climatology (Table S7). G1: one half-power branch continuous across ≥ 90 % of window longitudes. G2, at sites passing G1: synthetic tests (40 realisations; climatological KE plus observed anomalies resampled in 12-month blocks; linear change over 25 years; 150 km smoothing) of translation, 25 km flexing, uniform intensification and one-flank intensification (×1.2); a metric passes if translation is recovered and intensification gives < 50 % of it.

**Results.** Only SEIR150 passes G1; a post hoc change of branch start rescues no site. At SEIR150 (Table S8) the mean position passes G2 narrowly (northern-flank intensification gives +10.0 km, 47 % of the recovered translation); amplitude, sinuosity and curvature do not separate flexing from intensification. The SEIR150 position trend is −4.3 km per decade [−14.3, +5.8], p = 0.37. Upstream jet-core speed leads northward displacement (r = 0.31 at lag 3 months; surrogate p = 0.004, passing FDR over 20 relation tests), but the peak lies at the edge of the pre-specified window and, given G2, may partly reflect energy redistribution. Displacement is unrelated to meander-section (p = 0.89) and downstream EKE (p = 0.27). DRAKE fails G1 (one branch traceable across 70–55°W), a different measurement from the single monthly profile of Sect. S3.

**Table S7.** Gate G1.
<!-- /home/claude/work/out/h1/h1_gate.csv; plan1019 §3 -->

| Site (window) | Usable longitudes | G1 |
|---|---|---|
| SEIR150 (140–150°E) | 96 % | pass |
| KERG (76–84°E) | 41 % | fail |
| PAR (135.1–109.4°W) | 21 % | fail |
| DRAKE (70–55°W) | 14 % | fail |

**Table S8.** Gate G2 at SEIR150 (recovered 25-year change).
<!-- /home/claude/work/out/h1/h1_gate.csv; plan1019 §3 -->

| Scenario → metric | Response | Null spread (1 SD) |
|---|---|---|
| 25 km translation → mean position | 21.3 km (85 %) | — |
| Northern-flank intensification → mean position | +10.0 km | — |
| Uniform intensification → mean position | −1.2 km | — |
| 25 km flexing → amplitude | 8.3 km | 29.4 km |
| Translation → amplitude | +9.5 km | — |
| Uniform intensification → sinuosity | −0.076 | 0.056 |
| Flexing → sinuosity | −0.021 | 0.056 |

**Figure S5.** Half-power analysis at standing meanders. (a) Percentage of longitudes in each meander window, selected beforehand from bathymetry and climatology, at which a single half-power branch can be followed continuously in at least 90 % of months (dashed: 90 % criterion; only the Southeast Indian Ridge passes). (b–d) Synthetic tests at the Southeast Indian Ridge (40 realisations per scenario: climatological KE plus observed monthly anomalies resampled in 12-month blocks, with a linear change over 25 years): recovered 25-year change in (b) mean half-power position, (c) meander amplitude and (d) sinuosity, at 150 km smoothing, for no change, a 25 km northward shift, a 25 km change in meander amplitude (flexing), a uniform 20 % intensification and a 20 % intensification of the northern flank only. Dashes mark the imposed change.

The gate coverage and synthetic tests are shown in Fig. S5.

**Conclusion.** Half-power geometry gives an interpretable mean position at one of four meanders and cannot detect 25 km flexing; meander shape is examined on the jet-core contour (Sect. S9).

---

## S5 Drake Passage: comparison with the in situ record and sensitivity tests

**Purpose.** The main text uses altimetry alone. Here we check over 2005–2019 which altimetric quantities reproduce the variability measured from the ARSV Laurence M. Gould, and test the Drake Passage results that lie near the threshold or depend on mask, section or series choices.

**Comparison method.** The in situ record comprises 248 ADCP transects across Drake Passage over 2005–2019 (Product 3; Gutierrez-Villanueva et al., 2023): cross-track velocity at 70–118 m, EKE integrated over 0–780 m and total transport over 0–760 m [AUTHOR CHECK: Gutierrez-Villanueva et al. (2023) quote 237 transects for EKE]. Altimetric values from the monthly field of each transect's month are compared along the mean track, per transect and as annual means. As altimetry gives surface geostrophic flow and the in situ quantities include deeper and ageostrophic flow, the comparison checks variability, not magnitude. EKE and transport are compared north of about 58.5°S, ice-free in every month. The in situ band position is the latitude dividing the 0–760 m transport within the climatological altimetric envelope (55.1–59.3°S) into halves, so no shared variability is built in. Annual correlations are tested against surrogate series.

**Comparison results.** Over 2005–2019 Product 1 reproduces the mean cross-track velocity profile (r = 0.98), the velocity per 25 km bin and transect (r = 0.72 over 8187 pairs; anomalies r = 0.61; RMS difference 15 cm s⁻¹), EKE north of 58.5°S (r = 0.67 per transect, 0.89 for annual means) and transport north of 58.55°S (r = 0.54 per transect, 0.67 annual; surrogate p = 0.003) (Table S9; Fig. S2). The full-section in situ transport is reproduced neither by the northern proxy (r = 0.03 per transect and for annual means), as the in situ transports north and south of the ice-free limit vary in opposition (r = −0.54), nor by the whole-line integral including values mapped under seasonal ice (r = 0.07 per transect). This opposition is partly built into the split (a jet crossing a fixed latitude moves transport between the parts) and does not test the opposing frontal trends of Gutierrez-Villanueva et al. (2023). The half-power latitude of the merged energetic band along the section follows the in situ half-transport latitude with r = 0.50 per transect and 0.71 for annual means (0.50 with the profile truncated at the ice-free limit; Sect. S3), with a mean offset of 16 km and an RMS difference of 114 km, reflecting monthly fields compared with transects lasting about two days. Trends cannot be checked: over 15 annual values the in situ intervals (EKE −58.8 to +74.4 % per decade) are many times wider than the trends, and the overlap ends in 2019, as the products begin to diverge (Sect. S1). Over 2001–2025, EKE north of 58.5°S changes by +7.4 % per decade (−5.6 to +20.4; p = 0.24), not detectably.

**Table S9.** What the comparison of altimetry (Product 1) with the in situ record (Product 3) supports, 2005–2019.

| Quantity | r per transect / annual means | Variability reproduced? | Basis for extending the record |
|---|---|---|---|
| Mean cross-track velocity profile | 0.98 (mean profile) | n/a (mean structure) | supported |
| Velocity per 25 km bin and transect (8187 pairs) | 0.72 (anomalies 0.61) / – | transect-scale anomalies: yes | supported |
| EKE north of 58.5°S | 0.67 / 0.89 | yes | supported |
| Transport north of 58.55°S | 0.54 / 0.67 | yes | supported for the northern part |
| Position of the energetic band | 0.50 / 0.71 | yes | interannual variability only (RMS 114 km) |
| Full-section transport | 0.03 / 0.03 | no (northern vs southern in situ parts: r = −0.54) | not supported |
| Trends | – | cannot be checked (15 annual values) | not established |

**Ice mask and EKE − KE contrast.** The mask removes 37 % of box pixels. Without it, box EKE is no longer at the threshold in Product 1 (+7.1 % per decade, p = 0.08) and the Product 1 EKE − KE contrast weakens from +6.0 to +5.0 points per decade (p = 0.04 in both); in Product 4 neither is significant with or without the mask (Table S1). The contrasts were tested after being described. Ending in 2024, masked box EKE gives p = 0.057.
<!-- dk1_trends.csv (nomask rows); R2 A6 (37 %) and "Revised trend test" (2001-2024) -->

**Section geometry.** An earlier proxy that dropped ice-covered points month by month tracked the number of valid points (r = 0.89) and gave a spurious trend (+1.9 % per decade); it is withdrawn. T_perm integrates fixed points ice-free in every month, on the earlier line (65°W 55°S–60°W 63°S) and the mean LM Gould track (GV track, primary); T_all integrates all points, including values mapped under ice (Table S11).

**Table S11.** Coverage-matched comparison with the in situ transport, 2005–2019 (248 transects in 137 months; annual means).
<!-- plan1019 §1; /home/claude/work/out/d1/d1_summary.csv, d1_trends.csv (earlier-line south), d1b_summary.csv, d1b_trends.csv -->

| | Earlier line | GV track |
|---|---|---|
| Points kept / southern limit | 53 of 121 / 58.47°S | 58 / 58.55°S |
| Share of full-line mean surface proxy | 65 % | 67 % |
| Share of in situ total (0–760 m) north of the limit | 62 % | 63 % |
| T_perm vs in situ north part, r per transect / annual | 0.44 / 0.72 | 0.54 / 0.67 (surrogate p = 0.003) |
| T_perm vs in situ total, r per transect / annual | 0.05 / 0.11 | 0.03 / 0.03 |
| T_perm vs in situ south part, r per transect / annual | −0.45 / −0.69 | −0.59 / −0.72 |
| In situ north vs south part, r per transect | −0.55 | −0.54 |
| T_all vs in situ north / total / south, r per transect | — | 0.14 / 0.07 / −0.09 |
| T_perm trend 2001–2025 | +2.3 [−3.7, 8.2], p = 0.40 | +3.9 [−2.1, 9.9], p = 0.17 |
| T_all trend 2001–2025 | +0.7, p = 0.16 (earlier test) | +0.6 [−0.4, 1.6], p = 0.22 |

The northern points carry 63 % of the mean ADT drop across the earlier line. Averaged by month (n = 137) the north-part r is 0.66; deseasonalised, 0.50; independent re-coding reproduced every correlation within ±0.02. The north–south anticorrelation is partly built in (a jet crossing the fixed latitude moves transport between the parts without changing the total) and does not test the frontal compensation of Gutierrez-Villanueva et al. (2023), which concerns opposing trends among frontal regions. Seven line variants (end points moved by 0.5° or 1° of longitude, the southern limit moved by ±0.5°, and an approximation of the LM Gould track) give T_perm trends of −0.2 to +4.8 % per decade over 2001–2025 (p = 0.14–0.94, revised test) and T_all trends of +0.4 to +1.6 % (p = 0.03–0.49).<!-- outputs/retest/retest_summary.csv keys T(perm|all)_v2_*|2001-2025 and *_gould_approx -->

**Changes along the transect (main-text Sect. 3.1).** The analysis follows the two coordinate systems of Gutierrez-Villanueva et al. (2023). In passage coordinates, the altimetric surface velocity is averaged in 25 km bins of distance from the north along the mean LM Gould track; the in situ velocity is the ADCP cross-track velocity at 70–118 m, with northbound transects reordered from the north and the ADCP offset applied as in their code. In streamwise coordinates, their streamline values ψ (Niiler–Maximenko mean dynamic topography minus 0.60 m plus sea-level anomaly) were mapped onto our ADT by sampling the monthly ADT of Product 1 at the mean position of each streamline pair on each transect (4005 pairs): ADT = 0.62 + 1.04ψ (r = 0.99; median offset 0.59 m, consistent with their 0.60 m). Zones: Subantarctic Front ψ = −0.40 to −0.60 m, zone between the fronts −0.60 to −1.00 m, Polar Front −1.00 to −1.25 m (their −1.00 to −1.30 m, limited by the bin edges of the archived files). Points with ice fraction above 0.10 are masked monthly; in streamwise coordinates a month is excluded for a zone when any of its points is ice-covered. Trends over October 2005–April 2019 use annual means for both records (AR(1) test), over 2001–2025 monthly series where complete and annual means otherwise. Over 2005–2019 the bin-by-bin trend profiles of altimetry and in situ velocity correlate with r = 0.57 (Product 1) and 0.63 (Product 4); the 2001–2025 altimetric profiles correlate with the 2005–2019 in situ profile with r = 0.64 and 0.62. The zone trends (Table S10) reproduce the sign pattern of Gutierrez-Villanueva et al. (2023) in passage coordinates; in streamwise coordinates neither record shows a significant transport trend with this test, and the EKE decrease they report in the Polar Front is not seen in altimetry. The southern zones (PF flank, SACCF) are seasonally ice-covered, so their altimetric trends rest on ice-free months only.
<!-- drake/dk5_transect_vs_gv23.py, dk5b_insitu_velocity_fix.py, dk5c_zones.py -> drake/out/dk5/dk5_trends.csv, dk5_zones.csv, dk5_psi_to_adt.csv -->

**Table S10.** Trends along the LM Gould transect in the zones of Gutierrez-Villanueva et al. (2023): 95 % intervals in brackets. In situ: ADCP velocity at 70–118 m, 0–760 m transport and 0–780 m EKE; altimetry: surface geostrophic velocity, surface transport per unit depth and surface EKE.

| Zone | Quantity | In situ 2005–2019 | Product 1 2005–2019 | Product 4 2005–2019 | Product 1 2001–2025 | Product 4 2001–2025 |
|---|---|---|---|---|---|---|
| SAF–PF zone (150–300 km) | cross-track velocity (cm s⁻¹ per decade) | +4.1 [−2.8, +11.1], p = 0.23 | +3.6 [−1.8, +9.0], p = 0.16 | +3.1 [−3.3, +9.6], p = 0.24 | +3.0 [+0.6, +5.4], p = 0.02 | +2.9 [+0.5, +5.4], p = 0.02 |
| PF southern flank (350–550 km) | cross-track velocity (cm s⁻¹ per decade) | −1.9 [−9.3, +5.4], p = 0.40 | −1.9 [−4.3, +0.4], p = 0.10 | −1.4 [−4.1, +1.3], p = 0.27 | −0.4 [−1.9, +1.2], p = 0.57 | −0.5 [−2.0, +1.0], p = 0.47 |
| SACCF (750–850 km) | cross-track velocity (cm s⁻¹ per decade) | +3.1 [−2.6, +8.7], p = 0.19 | +2.9 [+1.2, +4.5], p = 0.003 | +2.5 [+1.1, +4.0], p = 0.003 | −0.3 [−1.5, +0.9], p = 0.70 | −0.1 [−1.1, +1.0], p = 0.88 |
| SAF (−0.40 to −0.60 m) | transport (% per decade) | +5.9 [−16.6, +28.4], p = 0.47 | +0.2 [−6.5, +6.9], p = 0.89 | +1.4 [−4.0, +6.7], p = 0.48 | −1.7 [−3.3, −0.2], p = 0.03 | −0.8 [−2.0, +0.4], p = 0.17 |
| SAF–PF zone (−0.60 to −1.00 m) | transport (% per decade) | −1.5 [−17.9, +14.9], p = 0.81 | +0.3 [−4.2, +4.7], p = 0.80 | +0.3 [−1.6, +2.3], p = 0.62 | +0.8 [+0.1, +1.5], p = 0.05 | +0.5 [+0.0, +1.0], p = 0.05 |
| PF (−1.00 to −1.25 m) | transport (% per decade) | +4.8 [−22.1, +31.7], p = 0.52 | +0.6 [−1.2, +2.4], p = 0.39 | −0.9 [−3.0, +1.2], p = 0.37 | −0.1 [−0.9, +0.7], p = 0.83 | +0.0 [−0.8, +0.8], p = 0.99 |
| SAF (−0.40 to −0.60 m) | EKE (% per decade) | +24.2 [−24.2, +72.6], p = 0.26 | +19.4 [−36.2, +75.0], p = 0.32 | −12.0 [−131.7, +107.6], p = 0.65 | +18.4 [−0.8, +37.6], p = 0.06 | +6.0 [−11.5, +23.6], p = 0.53 |
| SAF–PF zone (−0.60 to −1.00 m) | EKE (% per decade) | −0.4 [−75.2, +74.4], p = 0.99 | +24.3 [−27.3, +75.8], p = 0.23 | +15.2 [−40.2, +70.5], p = 0.39 | +15.0 [−4.3, +34.3], p = 0.14 | +10.9 [−6.1, +27.8], p = 0.18 |
| PF (−1.00 to −1.25 m) | EKE (% per decade) | −9.3 [−91.2, +72.7], p = 0.66 | +17.9 [−2.6, +38.3], p = 0.11 | +14.6 [−13.0, +42.2], p = 0.26 | +3.4 [−4.0, +10.8], p = 0.40 | +0.9 [−6.8, +8.5], p = 0.87 |

**In situ quantities.** Total cross-track transport over 0–760 m (mean 75.5 Sv) = depth-independent reference part 51.1 + Ekman 8.1 + geostrophic relative to 760 m 15.9 Sv. Velocities: ADCP 70–118 m bins. EKE: the archive variable, integrated over 0–780 m (calc_eke_uvprime.m, maxz = 780); Gutierrez-Villanueva et al. (2023) present it depth-averaged, a scaling that does not affect correlations. The archive files contain 248 transects for both transport and EKE [AUTHOR CHECK: Gutierrez-Villanueva et al. (2023) quote 237 transects for the EKE analysis].

**Same-period trends, 2005–2019.** No trend is significant (Table S12), and the test is liberal at n = 15 (Sect. S2). The in situ EKE trend (+7.8 % per decade, −58.8 to +74.4) lies between the products, but the intervals are far wider than the differences, favouring neither product. Gutierrez-Villanueva et al. (2023) found a significant Polar Front EKE decrease [verify against Gutierrez-Villanueva et al., 2023], beside the non-significant trend of our in situ series north of 58.5°S.

**Table S12.** Trends over 2005–2019 (% per decade; positions in km per decade).
<!-- /home/claude/work/drake/out/dk3/dk3_summary.csv (EKE, annual); /home/claude/work/out/d1/d1b_trends.csv (transport, annual); /home/claude/work/drake/out/dk4/dk4_trends.csv (T_perm monthly, 2005-2019; numbercheck D); /home/claude/work/drake/out/dk2/dk2_trends.csv and dk2_vs_insitu.csv (positions; altimetric series monthly) -->

| Quantity | In situ | Product 1 | Product 4 |
|---|---|---|---|
| EKE north of 58.5°S (annual) | +7.8 [−58.8, 74.4], p = 0.63 | +10.1 [−71.2, 91.3], p = 0.61 | +0.8 [−63.3, 65.0], p = 0.95 |
| Transport: in situ total / north / south; altimetric T_perm (GV, annual) | −1.2 [−12.7, 10.3], p = 0.70 / −0.9 [−24.4, 22.6], p = 0.87 / +0.2 [−53.9, 54.4], p = 0.98 | +4.9 [−20.9, 30.8], p = 0.49 | +5.4 [−19.1, 29.9], p = 0.47 |
| Altimetric T_perm (GV, monthly series) | — | +6.8 [−11.4, 24.9], p = 0.40 | +7.5 [−10.4, 25.3], p = 0.34 |
| Band position (in situ half-transport latitude; altimetric half-power latitude) | +4.1 [−29.4, 37.6], p = 0.80 | +19.4 [−46.3, 85.1], p = 0.51 | +17.6 [−84.3, 119.6], p = 0.72 |

**Wind lags.** Lagged wind–transport correlations were computed with the withdrawn proxy. [PENDING: wind–transport lag correlations with T_perm on the GV track.]

**Figure S6.** Coverage of the altimetric surface transport proxy on the earlier (v2) line and comparison with the in situ transport over matching latitudes. (a) Mean in situ total cross-track transport in the upper 760 m per 25 km bin, LM Gould ADCP transects 2005–2019 (orange bars), and mean surface geostrophic velocity normal to the altimetric section, 2001–2025 (blue line); dots mark the 53 points ice-free in every month that form the proxy; shading marks the part of the section where the ice fraction exceeds 0.10 in at least one month (south of about 58.5°S). (b) Annual means of the in situ transport north of 58.5°S (orange) and of the proxy (blue), in standard-deviation units, 2005–2019. (c) As (b), for the in situ transport over the full section. (d) In situ transport north (orange) and south (grey) of 58.5°S for each transect against the proxy of the same month, in standard-deviation units. The main text uses the mean LM Gould track (58 points north of 58.55°S); values for both lines are in Table S11.

The coverage of the proxy on the earlier line and its comparison with the in situ transport are shown in Fig. S6.

**Conclusion.** Altimetry reproduces the interannual variability of the ship-measured velocity anomalies, EKE, northern transport and band position, but not full-section transport, and the 15-year overlap cannot check trends. Box EKE and the EKE − KE contrast weaken without the mask or in Product 4, and 2005–2019 trends do not discriminate between products.

---

## S6 Comparison with earlier estimates

**Purpose.** Reading our EKE trends as an update of Martínez-Moreno et al. (2021; MM21) and Hogg et al. (2015) requires separating period, product version and method from ocean change.

**Method.** *MM21 procedure* (from their code repository): EKE = ½(ugosa² + vgosa²) × ρ0, area-integrated (PJ m⁻¹), 365-day rolling mean, Theil–Sen with a modified Mann–Kendall test; their result south of 35°S is +5.1 % per decade (+0.0368 PJ m⁻¹ per decade), 1993 to March 2020, all-satellite vDT2018 (anomalies relative to 1993–2012, 1° grid). We applied it over 35–65°S (Table S13; pixels valid in every month, no ice mask, 12-month running mean of monthly values), and our calibrated test to the unsmoothed series. *Hogg et al. (2015) boxes* (ice-masked): Indian 40–150°E, 57–44°S; Pacific 150°E–72°W, 62–48°S; Atlantic 35°W–10°E, 56–46°S. Their crossover along-track EKE from 1993, relative to the 1993–2001 mean (90 % intervals), gives +6.0 ± 1.6 % (Indian) and +6.4 ± 1.8 % (Pacific) per decade, and 18.3 ± 5.1, 14.9 ± 4.1 and 4.0 ± 3.7 cm² s⁻² per decade (Atlantic without a relative value). Crossover EKE is 1.6–1.9 times gridded EKE, so only relative trends are compared. *Scale matching*: Product 1 block-averaged to 1/4° (Table S3).

**Table S13.** Our data processed with the MM21 procedure, 35–65°S, no ice mask (% per decade).
<!-- plan1019 §2b; /home/claude/work/out/m2/m2_trends.csv -->

| Product, period | MM21 procedure | Calibrated test | Mean (PJ m⁻¹) |
|---|---|---|---|
| Product 1, 2001–2020 | +5.2 | +4.5 [−2.6, 11.6], p = 0.13 | 1.44 |
| Product 1, 2001–2025 | +7.5 | +6.6 [1.2, 12.0], p = 0.04 | 1.48 |
| Product 4, 2001–2020 | +1.8 | +1.2 [−1.5, 3.9], p = 0.36 | 1.28 |
| Product 4, 2001–2025 | +2.7 | +2.3 [0.2, 4.5], p = 0.04 | 1.29 |

**Ice-masked variants.** With the ice mask, the calibrated Product 1 trend over 35–65°S is +4.6 % per decade [−1.9, 11.1] (p = 0.11) for 2001–2020 and +6.6 % [1.5, 11.6] (p = 0.03) for 2001–2025 (Table S14), close to the unmasked values; the masked integral rises by +0.092 PJ m⁻¹ per decade [0.022, 0.162] on a mean of 1.40 PJ m⁻¹ (Product 1) and +0.029 (Product 4), against MM21's +0.0368.

**Level difference.** Our integrated EKE level is about twice MM21's: 1.44 against about 0.72 PJ m⁻¹, the level implied by their trends (+0.0368 PJ m⁻¹ per decade ÷ 5.1 % per decade). The cause has not been investigated. Candidates, unranked: (i) product version (DT2024 vs vDT2018); (ii) grid and order of operations (1/8° velocities vs coarsening to 1°, which removes sub-1° variance if velocities are averaged before squaring [AUTHOR CHECK: whether MM21 coarsen velocities or EKE]); (iii) anomaly definition (their 1993–2012 reference, daily fields; our 2001–2025 reference, monthly averages of daily values); (iv) region and masking (pixel validity, land and ice); (v) normalisation (ρ0, area weighting). A constant level factor leaves relative trends unchanged; a time-varying one would not.

**Table S14.** EKE trends by region and period, all regions ice-masked (% per decade; absolute in cm² s⁻² per decade for 2001–2025).
<!-- plan1019 §2; /home/claude/work/out/m1/m1_trends.csv (Atlantic 2001-2012 and 2001-2020) -->

| Region | Period | Product 1 | Product 4 |
|---|---|---|---|
| Indian | 2001–2012 | +5.8 [−17.8, 29.3], p = 0.47 | +9.6 [1.7, 17.5], p = 0.03 |
| | 2001–2020 | +8.7 [1.7, 15.7], p = 0.03 | +5.7 [1.3, 10.2], p = 0.02 |
| | 2001–2025 | +9.0 [4.9, 13.0], p = 0.007; abs. +18.0 [9.9, 26.2] | +4.8 [1.7, 7.9], p = 0.013; abs. +8.3 |
| Pacific | 2001–2012 | +0.4 [−29.2, 30.1], p = 0.97 | +2.7 [−17.9, 23.3], p = 0.68 |
| | 2001–2020 | +5.8 [0.8, 10.8], p = 0.04 | +3.1 [−2.3, 8.4], p = 0.20 |
| | 2001–2025 | +8.6 [0.9, 16.4], p = 0.04; abs. +15.7 [1.7, 29.6] | +4.9 [0.5, 9.4], p = 0.05 (at the threshold); abs. +7.5 |
| Atlantic | 2001–2012 | −5.1 [−24.1, 13.9], p = 0.50 | −3.1 [−21.0, 14.7], p = 0.65 |
| | 2001–2020 | +5.4 [−9.8, 20.5], p = 0.32 | +3.5 [−5.2, 12.2], p = 0.38 |
| | 2001–2025 | +8.8 [−0.1, 17.6], p = 0.056; abs. +13.5 | +5.2 [−0.2, 10.7], p = 0.06; abs. +6.9 |
| 35–65°S | 2001–2020 | +4.6 [−1.9, 11.1], p = 0.11 | +1.2 [−1.5, 4.0], p = 0.36 |
| | 2001–2025 | +6.6 [1.5, 11.6], p = 0.03 | +2.3 [0.2, 4.5], p = 0.04 |
| ACC band | 2001–2020 | +6.4 [0.9, 11.9], p = 0.04 | +3.2 [−0.1, 6.5], p = 0.07 |
| | 2001–2025 | +8.0 [3.3, 12.8], p = 0.017 | +3.8 [1.3, 6.3], p = 0.020 |

The 2001–2012 windows, overlapping Hogg et al.'s record, are uninformative (intervals about ±8 to ±30 % per decade; one of six sector trends significant, Indian in Product 4). For 2001–2025, Hogg et al.'s relative trends lie between our products and inside both intervals in the Indian and Pacific boxes: consistency, not reproduction, as periods, sampling and reference means differ. Extending 2001–2020 to 2001–2025 raises every trend except Product 4 Indian and narrows every interval except Product 1 Pacific.

**Product versions.** MM21 used all-satellite vDT2018, no longer distributed; Barceló-Llull et al. (2025) found a large version sensitivity of EKE trends. Two-satellite vDT2018 remains available from C3S. [PENDING: two-satellite vDT2018 comparison, if extracted.] [AUTHOR DECISION: add a 1993–2025 two-satellite comparison only if the 1993–2000 extraction is completed.]

**Figure S7.** (a) EKE trends (95 % intervals) for 2001–2020 and 2001–2025, both products, in the Hogg et al. (2015) boxes, 35–65°S and the ACC band, with published estimates (diamonds; different period, version and method; Hogg et al. at 90 %); filled, p < 0.05; open, p ≥ 0.05. (b) Product difference (all-satellite minus two-satellite, percentage points per decade). (c) EKE integrated over 35–65°S in percent of its 2001–2025 mean (MM21 procedure, no ice mask), monthly and 12-month running mean; shading, years after the MM21 period. [AUTHOR CHECK: the 0.0498 p-value of the Product 4 Pacific box plots as filled; redraw with the rounding rule before submission.]

The comparison with earlier estimates is summarised in Fig. S7.

**Conclusion.** Processed as in the earlier studies, our trends have the same sign and order. Over 35–65°S the calibrated trend is not significant for 2001–2020 and reaches p = 0.03–0.04 for 2001–2025 with or without the ice mask, a marginal result (Sect. S2); the twofold level difference with MM21 is unexplained.

---

## S7 Antarctic-slope sea-level index (after Hughes et al., 2014, and Hogg et al., 2015)

**Purpose.** The surface proxy misses the southern Drake Passage section, and the in situ record ends in 2019. Sea level over the Antarctic slope responds to the bottom-pressure (barotropic-type) part of the circumpolar transport. The index is the third transport quantity of main-text Sect. 2.4, distinct from the surface proxy and the in situ 0–760 m transport.

**Background.** In NEMO simulations (ORCA025, ORCA12), Hughes et al. (2014) found circumpolar transport correlated with bottom pressure over the Antarctic slope (model correlation below −0.5; 100 m depth bins over 0–5000 m, strongest at 1000–3000 m), with a factor of −3.69 Sv mbar⁻¹; sea level is amplified relative to bottom pressure beyond about 5-year periods. Hogg et al. (2015) applied this to AVISO sea level over shelf and slope (T = −3.69(h − h_f) − 1.11 h_f Sv cm⁻¹ [CHECK signs and units against Hogg et al., 2015]), called the factors preliminary and gave no numerical transport trend.

**Our construction.**
1. Product 1 monthly sea level south of 60°S, 2001–2025 [PENDING: extraction south of 65°S].
2. Region: the 1000–3000 m depth range of the continental slope contiguous with the Antarctic shelf, excluding isolated rises, plateaus and ridges (e.g. Maud Rise, the southern Kerguelen Plateau, the Pacific–Antarctic Ridge), by a contiguity mask on bathymetry to 80°S [CHECK: bathymetry product]; it replaces the model-correlation mask. The unrestricted selection (all 1000–3000 m points south of 60°S) is a sensitivity test (Table S15).
3. Each month, only pixels with ice fraction ≤ 0.10; the ice-free fraction is recorded.
4. Mean per 100 m depth bin, then over bins; also in 12 longitude bands of 30° to monitor zonal coverage.
5. Bias and annual and semi-annual harmonics removed (h); low-frequency part h_f from 10, 20 and 40-year sinusoid fits, as in Hogg et al. (2015).
6. Conversion to T_hogg (Sv) with the published factors, labelled model-calibrated; h is reported first.
7. Trends of h and T_hogg (Sect. S2); correlations with the in situ total, north and south parts (matched transects, annual means), the depth-independent reference part [PENDING: if the archive allows] and T_perm.

**Results.** [PENDING: h (cm) and T_hogg (Sv, model-calibrated) series.] [PENDING: trends of h and T_hogg with 95 % intervals and p.] [PENDING: correlations with in situ total, north, south and reference parts, and with T_perm.] [PENDING: monthly ice-free fraction, mean and distribution by longitude band.] [PENDING: comparison with Hogg et al. (2015).]

*Sampling caveat (not a result).* A smoke test on the slope already extracted (60–65°S, 1000–3000 m, before the contiguity mask was introduced) found 47 % of slope pixels ice-free in a typical month; sea-ice sampling is expected to dominate further south. Nothing from that test is used.

**Table S15.** Planned sensitivity tests.

| Test | What it checks | Result |
|---|---|---|
| Contiguous slope vs unrestricted 1000–3000 m selection south of 60°S | Inclusion of isolated rises, plateaus and ridges off the Antarctic margin | [PENDING] |
| Correlation of h with monthly ice-free fraction; fixed always-ice-free pixel set vs monthly set | Spurious variability or trend from varying coverage (cf. the withdrawn Drake proxy, Sect. S5) | [PENDING] |
| Minimum coverage per month and per longitude band | Dependence on poorly sampled months | [PENDING] |
| Ice-fraction threshold | Dependence on the mask | [PENDING] |
| Depth range (1000–3000 m vs wider bands within 0–5000 m) | Approximation of the slope region | [PENDING] |
| h_f construction (10/20/40-year fits vs other low-pass filters) | Low-frequency part, weighted differently by the conversion | [PENDING] |
| Conversion factors (ORCA025- vs ORCA12-based; no conversion) | Calibration uncertainty, largest for trends | [PENDING] |
| Removal of a large-scale sea-level signal (method to be specified) | Aliasing of steric or mass-driven margin sea-level change into "transport" | [PENDING] |
| Product 4 | Constellation dependence | [PENDING] |

**Interpretation (conditional).** If the index is well sampled and follows the in situ depth-independent or southern part, it adds, to 2025, a barotropic-type indicator of circumpolar transport variability that neither the surface proxy nor the 0–760 m record contains; if sea-ice sampling or large-scale sea-level change dominates, it will be reported with its coverage and not interpreted. Either way the Sv values remain model-calibrated, with the largest uncertainty for trends.

**Figure S8.** [PENDING] (a) Contiguous slope region with mean monthly ice-free fraction. (b) h and T_hogg (model-calibrated) with the in situ total, south and reference parts and T_perm. (c) Monthly ice-free fraction.

Fig. S8 will show the index and its sampling.

**Conclusion.** [PENDING: one sentence on whether the index is usable and what it shows.]

---

## S8 Tracked eddies

**Purpose.** Where eddies form and whether their nonlinearity changes, as context for the main-text reading of the EKE increase as more intense rather than more numerous eddies (main text Sect. 3.3).

**Method.** META4.0 DT all-satellite atlas (Product 5), south of 35°S, January 2001 to August 2023 (16.8 million observations, 524 000 tracks). Births per 1° cell per year. U/c = rotational speed over 7-day centre-displacement speed, for tracks of at least seven days.

**Results.** Births are enhanced along the whole energetic band, particularly at Kerguelen and Drake Passage, and the meander sections are not singular generation sites at 1° (Table S16). No wake-minus-rest trend difference is significant, monthly births and path length show no consistent link, and U/c has no detectable trend.

**Table S16.** Tracked-eddy diagnostics (Product 5, January 2001–August 2023).
<!-- R2 E (as in SI v3 S8); AC-A v2 -->

| Quantity | Value |
|---|---|
| Births per year per 1° cell: ACC band; wake boxes | 1.7; 1.8 (PAR) to 2.6 (Kerguelen) |
| Births, upstream / meander / downstream section: PAR; SEIR150; Kerguelen | 1.5 / 2.0 / 1.6; 1.5 / 1.1 / 2.1; 4.0 / 3.6 / 2.7 |
| Wake minus rest, normalised trend (points per decade): eddies present; births; amplitude; Σ amplitude² | +0.3 (p = 0.71); −1.9 (p = 0.21); −0.03 (p = 0.97); −1.9 (p = 0.55) |
| Monthly births vs path length | PAR p = 0.06 (r ≈ −0.2); SEIR150 p = 0.16; KERG p = 0.74 |
| U/c: median; share of observations > 1; trend | 6.9; 99 %; p = 0.97 |
| U and c, trend | about +3.5 and +3.2 % per decade |

**Caveats.** Only eddies above the mapping resolution are tracked (speed-radius median about 37 km, 5th percentile 24 km), so an unchanged radius is weakly informative; the atlas shares the constellation dependence of Product 1 (Sect. S1).

**Figure S9.** Tracked eddies (Product 5, January 2001 to August 2023). (a) Eddy births per year per 1° cell, with the mean boundaries of the ACC band (grey). (b) Median nonlinearity ratio U/c of the eddies in the ACC band, monthly (thin) and 12-month running mean (thick).

Eddy births and nonlinearity are shown in Fig. S9.

**Conclusion.** Eddy generation is spread along the energetic band rather than concentrated at the meanders, and tracked-eddy nonlinearity shows no detectable trend.

---

## S9 Along-stream analysis: details

**Purpose.** Choices behind main-text Sect. 3.4, monthly covariation tests and our reproduction of Liu et al. (2026).

**Path and bins.** Monthly contour of ADT relative to its ACC-band mean at +0.18 m, the level of highest climatological speed (Subantarctic Front side); the speed–level curve is broad (0.16–0.22 m within 0.5 %), so the level is weakly determined. Smoothing 150 km; 72 bins of 5°. Bins at 60–50°W are excluded because the contour follows the Malvinas Current (sinuosity coverage 0.69). Per-bin tests use the Sect. S2 statistic with 300 null series (single bootstrap for intervals) and FDR at 10 %. The EKE–sinuosity cross-bin correlation (r = 0.29, p < 0.001) is not used in the main text because EKE and monthly-path sinuosity both contain the eddies.

**Monthly covariation.** A pre-specified box test at PAR, SEIR150 and KERG found weak, episodic covariation, with path lengthening leading downstream EKE by 1–3 months (pooled p = 0.007), carried at SEIR150 by a few episodes, and no upstream-speed → path-length link. A box-free version (36 overlapping 30° windows; downstream EKE 15–35° east), corrected after its first pass for window overlap, birth normalisation and contour level, shows simultaneous path-length–EKE covariation in 22 of 36 windows (lag 0 strongest) and a downstream association in 1 of 36 (Polar Front side, −0.2 m: 14 and 3 of 36, mixed signs). What holds is local, simultaneous covariation, as expected if eddies form the path distortions; no delayed downstream response is detected.

**Wake test details.** Wake boxes cover 23 % of the band area and 43 % of its mean EKE. In the pre-specified sensitivity period 2001–2017 the wake-minus-rest EKE difference is +4.1 points per decade (p = 0.13; Product 4 +4.3, p = 0.06; earlier test); no secondary localisation test passes FDR. At Drake, upstream KE rises more than downstream (Scotia) KE (p = 0.010), with EKE at the threshold (p = 0.047–0.053).

**Reproduction of Liu et al. (2026).** The only analysis outside 2001–2025: to match Liu et al., it uses daily Product 1 fields from 1993 (1993–2023 and 1993–2025, block-averaged to 1/4°). Following their Methods, we computed daily |∇ADT|, its 4-month occurrence above 20 % of the daily maximum, and axis, speed and section positions. Their jet-core acceleration is reproduced in every section to their two-decimal rounding, survives the calibrated test except at the Upper Ridge, and appears in the 2001–2025 monthly fields of both products. Their position trends are reproduced in sign and approximately in magnitude, but only the Lower-Ridge northward shift approaches significance with an autocorrelation-aware test (1993–2023 p = 0.047; 1993–2025 p = 0.10; 2001–2025 p = 0.24, earlier test). The difference lies in the test: on trend-free AR(1) series (n = 372) their modified Mann–Kendall test rejects at nominal 5 % in 63 % (φ = 0.5) to 98 % (φ = 0.8, 4-month running mean) of cases, because its variance uses n_eff from the undetrended series; this concerns significance, not their trend magnitudes. We could not reproduce their width (about 650 km zero-to-zero against their 100 km; we may have misread the definition) and make no width comparison.

**Figure S10.** Month-to-month covariation of the length of the jet-core path and EKE, in 30° windows (2001–2025, Product 1): correlation at the lag (0–3 months, EKE lagging) with the largest magnitude, between the path length and (a) EKE within 150 km of the path in the same window and (b) EKE within 150 km of the path 15–35° downstream of the window centre. Black: path at the Subantarctic Front level (+0.18 m, as in the main text); grey: Polar Front level (−0.20 m). Filled symbols: p < 0.05 against phase-randomised surrogates.

The month-to-month covariation of path length and EKE is shown in Fig. S10.

**Conclusion.** The contour level is weakly determined but stable; monthly tests support only local, simultaneous path–EKE covariation; the reproducible signal at the Pacific–Antarctic Ridge is jet-core acceleration.

---

## S10 Southern Annular Mode

**Purpose.** Referee A asked whether KE and EKE variability follows the SAM.

**Method.** Deseasonalised, detrended monthly ACC-band KE and EKE (Product 1) correlated with the SAM index of Marshall (2003), SAM leading by 0–12 months; significance from phase-randomised surrogates of the index, accounting for the lags examined.

**Results.** Maximum r = 0.15 for KE (lag 0, p = 0.12) and 0.09 for EKE (p = 0.24); 12-month running means r = 0.11 and 0.08; annual means 0.04 and 0.05.

**Conclusion.** Circumpolar KE and EKE variability shows no detectable relation to the SAM index over 2001–2025; a lagged regional response to individual SAM events (Meredith and Hogg, 2006) would not be expected to appear in this correlation.
