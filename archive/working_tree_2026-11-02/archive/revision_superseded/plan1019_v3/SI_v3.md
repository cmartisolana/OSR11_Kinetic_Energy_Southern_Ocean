# Supporting Information for "Kinetic Energy Variability and Eddy Saturation in the Southern Ocean"

Martí-Solana et al., State of the Planet, Ocean State Report 11, Chapter 2

<!-- Number sources (precedence brief §0 > drake_results > plan1019 > AC v2 > v23 > R2 "Revised trend test" > earlier R2): inputs/numbers_v3.md; inputs/drake_results.md; inputs/plan1019_results.md; inputs/R2_revision_results.md; inputs/AC_reply_*_v2.md; inputs/manuscript_v23.md; CSVs cited per table. -->

**Conventions.** Product 1: all-satellite SEALEVEL_GLO_PHY_L4_MY_008_047 (DT2024, 1/8°), used in the main text. Product 4: two-satellite SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057 (1/4°). Unless stated, trends are Theil–Sen slopes of deseasonalised monthly series, 2001–2025, in percent of the period mean per decade, with 95 % intervals in brackets and p from the test of Sect. S2. Significance means p < 0.05; p = 0.050 is "at the threshold". Tests chosen after results had been seen are labelled post hoc.

---

## S1 Two-satellite product

**Purpose.** The all-satellite constellation changes over 2001–2025, so part of a trend in mapped EKE may come from sampling. Repeating each diagnostic with a constant pair of altimeters shows which conclusions depend on the product.

**Method.** All diagnostics were recomputed with Product 4 using the same definitions, regions, ice mask (pixels with ice fraction > 0.10 in any month excluded), periods and test. Grid spacing is 1/8° (Product 1) and 1/4° (Product 4); the effective resolution is of the order of 100–150 km in wavelength at ACC latitudes (Ballarotta et al., 2019) [CHECK: DT2024 effective resolution from product QUID]. Product differences are tested as the trend of the difference of the mean-normalised series (Product 1 minus Product 4, points per decade). The eddy atlas (Product 5) exists for Product 1 only.

**Results.** In Drake Passage (Table S1), upstream KE and EKE increase in both products and upstream speed is at the threshold in both. Inside the box, Product 1 EKE is at the threshold and Product 4 EKE is about half as large and not significant (difference +4.7 points per decade [−1.2, 10.7], p = 0.09). Downstream EKE is significant in Product 1 only.

**Table S1.** Drake Passage trends, 2001–2025. Box 80–52°W, 63–54°S (ice-masked); upstream 100–80°W and downstream 55–35°W within the ACC band (mean ADT −0.6 to +0.2 m). Relative trends in % per decade; absolute in cm² s⁻² per decade (speed: cm s⁻¹ per decade), period mean in parentheses. EKE − KE: trend of the difference of mean-normalised series.
<!-- /home/claude/work/drake/out/dk1/dk1_trends.csv (period 2001-2025); relative values and p as in drake_results.md §1 -->

| Region | Quantity | Product 1, % [CI], p | Product 1, absolute (mean) | Product 4, % [CI], p | Product 4, absolute (mean) |
|---|---|---|---|---|---|
| Box | KE | +2.4 [−2.5, 7.4], 0.29 | +7.8 [−7.9, 23.4] (319) | +1.4 [−3.4, 6.2], 0.52 | +4.2 [−10.4, 18.8] (301) |
| | EKE | +8.4 [0.2, 16.6], 0.050 | +21.7 [0.6, 42.8] (258) | +3.8 [−2.5, 10.1], 0.22 | +8.3 [−5.6, 22.3] (221) |
| | Along-stream speed | −0.1 [−3.2, 3.0], 0.91 | −0.02 [−0.51, 0.47] (15.7) | 0.0 [−3.1, 3.0], 0.98 | 0.00 [−0.47, 0.47] (15.4) |
| | EKE − KE | +6.0 [0.7, 11.3], 0.04 | — | +2.4 [−1.7, 6.4], 0.22 | — |
| Upstream | KE | +12.8 [6.8, 18.7], 0.007 | +19.2 [10.3, 28.2] (150) | +10.8 [5.8, 15.8], 0.007 | +16.2 [8.7, 23.8] (150) |
| | EKE | +15.3 [6.2, 24.4], 0.013 | +22.5 [9.1, 36.0] (147) | +11.5 [4.7, 18.3], 0.007 | +15.6 [6.4, 24.7] (135) |
| | Along-stream speed | +3.9 [0.3, 7.5], 0.050 | +0.33 [0.02, 0.63] (8.4) | +3.6 [0.0, 7.3], 0.056 | +0.31 [0.00, 0.61] (8.4) |
| | EKE − KE | +3.8 [−0.3, 7.9], 0.09 | — | +1.7 [−1.1, 4.6], 0.21 | — |
| Downstream | KE | +4.3 [−1.8, 10.3], 0.13 | +17.3 [−7.3, 41.8] (406) | +3.3 [−2.4, 9.0], 0.21 | +12.5 [−9.4, 34.5] (385) |
| | EKE | +8.0 [3.0, 13.0], 0.010 | +36.3 [13.5, 59.1] (455) | +3.8 [−1.5, 9.1], 0.13 | +14.2 [−5.6, 34.0] (375) |
| | Along-stream speed | +2.7 [−4.0, 9.3], 0.28 | +0.42 [−0.62, 1.47] (15.7) | +2.7 [−4.0, 9.3], 0.29 | +0.42 [−0.62, 1.46] (15.7) |
| | EKE − KE | +3.3 [−1.0, 7.7], 0.12 | — | +0.1 [−4.9, 5.1], 0.96 | — |

Eastward wind stress (WIND_GLO_PHY_L4_MY_012_006): box +8.7 [1.4, 15.9], p = 0.03; upstream +9.8 [1.3, 18.4], p = 0.03; downstream +9.4 [−1.5, 20.2], p = 0.09.

Along the LM Gould section (Table S2) the two products compare equally well with the in situ record.

**Table S2.** Comparison with the LM Gould ADCP record (2005–2019) and section diagnostics.
<!-- drake_results.md §2–3; /home/claude/work/drake/out/dk3/dk3_summary.csv (north/south bins); dk2/dk2_vs_insitu.csv (Product 4 RMS); dk4/dk4_trends.csv (Product 4 CIs) -->

| Quantity | Product 1 | Product 4 |
|---|---|---|
| Mean cross-section velocity profile vs ADCP (70–118 m), r | 0.98 | 0.97 |
| Velocity per 25 km bin and transect, r (n) | 0.72 (8187) | 0.72 (8189) |
| — bins north / south of 58.5°S, r | 0.70 / 0.75 | 0.69 / 0.73 |
| Velocity anomalies (bin mean removed), r | 0.61 | 0.60 |
| RMS velocity difference | 15 cm s⁻¹ | 16 cm s⁻¹ |
| EKE north of 58.5°S, r per transect / annual | 0.67 / 0.89 | 0.65 / 0.85 |
| EKE, all ice-free points, r annual | 0.92 | 0.88 |
| Band position vs in situ half-transport latitude, r per transect / annual | 0.50 / 0.71 | 0.48 / 0.65 |
| — mean offset; RMS difference | 16 km; 114 km | 5 km; 117 km |
| Surface transport proxy (GV track) vs in situ north part, r per transect / annual | 0.54 / 0.67 | [NUMBER NEEDED: Product 4 proxy correlations] |
| Section half-power trend 2001–2025, km per decade | −10 [−36, +16], p = 0.43 | −8 [−44, +29], p = 0.65 |
| EKE north of 58.5°S, trend 2001–2025 | +7.4 [−5.6, 20.4], p = 0.24 | +1.8 [−8.5, 12.0], p = 0.68 |
| Surface transport proxy (GV track), trend 2001–2025 | +3.9 [−2.1, 9.9], p = 0.17 | +4.2 [−1.8, 10.2], p = 0.15 |

Over the ACC band (Table S3) every trend has the same sign in both products. The mean-flow measures agree closely; the EKE trends are about half as large in Product 4, mainly through the sub-monthly part, and the products diverge mainly after about 2018 (Fig. S1), consistent with the time-varying constellation documented in the Product 1 QUID and by Hargous et al. (2026). The product difference is not significant in any region or period tested (p ≥ 0.09; p ≥ 0.18 on a common 1/4° grid). In Product 4 the relative KE and EKE trends are of similar size (+3.6 and +3.8 %), and the EKE − KE difference is not significant in either product.

**Table S3.** Circumpolar diagnostics, 2001–2025 unless stated.
<!-- v23 §3.1; plan1019 §2; R2 A4 table, A7 table and "Revised trend test"; AC-A v2 (proxy); /home/claude/work/out/m1/m1_trends.csv -->

| Quantity | Product 1 | Product 4 | Difference, points per decade |
|---|---|---|---|
| ACC KE, % (absolute, cm² s⁻² yr⁻¹) | +4.9 (1.36 [0.86, 1.86]), p = 0.003 | +3.6 (0.97), p ≤ 0.013 | — |
| ACC EKE, % (absolute, cm² s⁻² per decade) | +8.0 [3.3, 12.8] (+19.8), p = 0.017 | +3.8 [1.3, 6.3] (+8.1), p = 0.020 | +4.3 [−2.2, 10.8], p = 0.12 |
| ACC EKE, 2001–2020 | +6.4 [0.9, 11.9], p = 0.04 | +3.2 [−0.1, 6.5], p = 0.07 | +3.2, p = 0.37 |
| Low-frequency EKE ½\|Uₘ − ū\|² | +5.6 | +3.2 | — |
| Low-frequency EKE, common 1/4° grid | +5.5, p = 0.003 | +3.2, p = 0.013 | +2.3, p = 0.18 |
| Sub-monthly EKE ½σ²ₘ | +12.0, p = 0.066 | +5.0, p = 0.040 | — |
| Along-stream speed | +2.0, p = 0.003 | +1.8, p ≤ 0.013 | — |
| EKE − KE | +3.0 [−0.1, 6.1], p = 0.060 | +0.2 [−1.2, 1.5] (earlier test, p = 0.74) | — |
| Circumpolar surface transport proxy | +1.6 [0.9, 2.3], p = 0.003 | +1.5 [0.9, 2.1], p = 0.003 | — |
| Wake boxes minus rest of band, EKE | +0.1 [−3.0, 3.2], p = 0.95 | −0.3, p = 0.80 | — |
| Box 35–65°S, EKE (ice-masked) | +6.6 [1.5, 11.6], p = 0.03 | +2.3 [0.2, 4.5], p = 0.04 | +4.5 [−1.4, 10.4], p = 0.10 |
| Hogg et al. (2015) Indian box, EKE | +9.0 [4.9, 13.0], p = 0.007 | +4.8 [1.7, 7.9], p = 0.013 | +4.5 [−5.8, 14.7], p = 0.20 |
| Hogg et al. (2015) Pacific box, EKE | +8.6 [0.9, 16.4], p = 0.04 | +4.9 [0.5, 9.4], p = 0.05 (at the threshold) | +4.3 [−0.8, 9.4], p = 0.09 |
| Hogg et al. (2015) Atlantic box, EKE | +8.8 [−0.1, 17.6], p = 0.056 | +5.2 [−0.2, 10.7], p = 0.06 | +3.9, p = 0.17 |
| PAR jet-core speed, whole meander / Lower Ridge, m s⁻¹ per decade | +0.019 / +0.026, p = 0.01 | +0.017 / +0.024, p = 0.01 | — |
| Along-stream bins (Sect. S9) and cross-bin r(EKE, speed) | main text Sect. 3.7 | [NUMBER NEEDED: Product 4 bin statistics and cross-bin r] | — |

By sector (absolute, cm² s⁻² yr⁻¹, Product 1 / Product 4), EKE rises by +2.35 / +1.21 (Pacific), +2.14 / +0.90 (Atlantic) and +1.33 / +0.23 (Indian), and KE by +1.87 / +1.53, +1.27 / +0.89 and +0.92 / +0.41. Pacific KE and EKE are significant in both products (p = 0.007–0.033), Indian EKE in Product 1 only (p = 0.010), Atlantic EKE in neither (Product 1 p = 0.053) [NUMBER NEEDED: revised p for Product 4 Atlantic and Indian KE]. The sector along-stream speed and surface transport proxy trends agree between products (proxy within 0.2 points).
<!-- R2 A4 sector table (absolute values); revised p from R2 "Revised trend test" -->

**Figure S1.** ACC-band (a) KE of the monthly-mean flow, (b) EKE and (c) along-stream speed, 12-month running means in percent of the mean (Product 1 green, Product 4 grey); (d) difference of the mean-normalised EKE series with its trend and 95 % interval.

**Conclusion.** Signs and mean-flow trends do not depend on the product, whereas the EKE magnitude differs by about a factor of two; that difference is not significant over 2001–2025 but leaves any statement resting on the EKE magnitude product-dependent.

---

## S2 Trend-test calibration

**Purpose.** The monthly series are strongly autocorrelated, and common trend tests reject too often on such series; the main-text distinction between "significant" and "not detectable" requires a test that holds its nominal level.

**Method.** The statistic is the Theil–Sen slope divided by the exact OLS-slope standard deviation under a fitted autoregressive (AR) noise model. AR coefficients are bias-corrected by bootstrap (Kilian, 1998), including the trend removal; non-stationary fits are shrunk to a maximum root of 0.99. p and the 95 % interval are calibrated against 300 null series from the corrected model, processed identically (smallest p 0.0033, reported as 0.003). Monthly series are tested with an AIC-selected order and with AR(12); the larger p and the union of intervals are reported. Annual series use AR(1) with iterated bias correction. For the 1° maps the calibration is pooled over 10 persistence classes (p from |z|/κ, κ = max(q95/1.96, q99/2.576) of the class null), with false-discovery-rate (FDR) control at 10 %.

**Results.** Rates are in Table S4. AIC is liberal for the Drake-like AR(12) process (it drops small high-lag coefficients) and AR(12) for pure AR(1); the maximum rule takes the more conservative. A 15-value annual series cannot be tested at the nominal level (9.5 %); no in situ trend is significant under either test. Of 585 trends recomputed with this test, fewer are significant in every group than with the earlier test (e.g. KE/EKE terms 145 → 107). The ACC KE, EKE, along-stream speed and low-frequency EKE trends remain significant; the ACC EKE − KE contrast (p = 0.060), Atlantic EKE (p = 0.053) and Product 1 sub-monthly EKE (p = 0.066) do not.

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

**Pre-specified and post hoc.** Fixed before any trend was computed: half-power parameters (Sect. S3); ACC band and ice mask; wake boxes, including Drake upstream and downstream; specifications of the wake test, along-stream analysis, meander gates and monthly covariability tests. Post hoc or added after results were seen: the 400 cm² s⁻² Drake split (Sect. S3); the branch-start check (Sect. S4); the direct tests of the EKE − KE contrasts; the corrections of the box-free covariability design (Sect. S9).

**Conclusion.** The test holds its level within a few percent for monthly series, so main-text trends with p ≤ 0.02 are unlikely to be artefacts of autocorrelation, whereas results with p near 0.05 and trends from 15 annual values are undetermined.

---

## S3 Half-power method: parameters, detectability and the Drake Passage section

**Purpose.** The half-power position is used in the main text as a bound on the displacement of the energetic band; a bound is meaningful only if the parameters do not drive it and the detectable shift is known.

**Method.** Following Chambers (2018), applied to the KE of the monthly-mean flow: at each longitude, segments with KE above a threshold are found; around each peak an envelope extends to where KE falls to a fraction of the peak (minimum width imposed); the half-power latitude is the median of the KE integrated across the envelope. Baseline: 200 cm² s⁻² (the median of 2001–2025 mean KE over the ACC band, 195 cm² s⁻²), 25 %, 100 km. Monthly envelopes are matched to the climatological envelope at the same longitude, and only longitudes matched in ≥ 90 % of months are kept (900 longitudes, about 31 % of the band). Sensitivity: thresholds 150–400 cm² s⁻², fractions 20–50 %, widths 50–300 km. Detectability: known shifts imposed on the climatological profiles with observed anomalies, processed identically.

**Results.**
- Every parameter set gives p > 0.07 for the circumpolar trend (Table S5).
- None of the 33 bins of 5° reaches p < 0.05 [AUTHOR CHECK: R2_revision_results.md reports 4 of 33 bins at p < 0.05 under the earlier test (none after FDR) and a revised count of 5 "A2 longitude bins" with unstated binning].
- A 25 km shift over 25 years is detected in 100 % of realisations and a 10 km shift in 44 %, with about 70 % of the shift recovered.
- Detection frequency is 25–45 % at about 170–150°W and 115–75°W (including the Drake upstream region); these longitudes are excluded.
- The unmatched median latitude of the discussion version had at most 27 % power for a 50 km shift and responded to amplitude growth alone (−0.002 to −0.007° yr⁻¹); it was replaced.
- Three-year blocks (8 values, OLS as in Chambers, 2018): ACC p = 0.17; Indian sector p = 0.049, a single borderline value.

**Table S5.** Circumpolar half-power trend per parameter set. [NUMBER NEEDED: trend (km per decade), 95 % interval and p for each threshold, fraction and minimum-width combination; the files give only "all p > 0.07".]

**Withdrawn envelope width.** The broadening reported in the discussion version is withdrawn. The threshold-defined width increased (+0.47 km yr⁻¹, earlier test p = 0.002), but an imposed amplitude increase with a fixed jet produced +0.43 to +1.6 km yr⁻¹, so the metric tracks amplitude. Amplitude-independent widths narrowed slightly (half-maximum width −0.10 km yr⁻¹, p = 0.017; width at 25 % −0.12 km yr⁻¹, p = 0.047), but 5 km of resolution sharpening produces −0.05 to −0.25 km yr⁻¹, and the half-maximum width has ≤ 35 % power even for large broadening. No width change can be inferred.

**Drake Passage section and the 400 cm² s⁻² split (post hoc).** Along the mean LM Gould track (55.1–62.3°S) the baseline profile has one envelope, 55.1–59.3°S, half-power latitude 57.0°S, detected in 99.7 % of months: the Subantarctic Front (SAF) and Polar Front (PF) merge because KE stays above 200 cm² s⁻² between them. The baseline is therefore the energy centroid of the merged band, not a front. [AUTHOR CHECK: whether the section profile uses points south of the 58.5°S ice-free limit; the envelope reaches 59.3°S.] After seeing the profile we repeated the analysis at 400 cm² s⁻², inside the pre-defined range but chosen post hoc (Table S6). The PF is detected too rarely for a trend, and the SAF trend (at the threshold in Product 1) rests on 55 % of months.

**Table S6.** Section half-power at 400 cm² s⁻² (post hoc), 2001–2025.
<!-- drake_results.md §2; /home/claude/work/drake/out/dk2_thr400/dk2_trends.csv and dk2_vs_insitu.csv -->

| | SAF, Product 1 | SAF, Product 4 | PF, Product 1 | PF, Product 4 |
|---|---|---|---|---|
| Climatological half-power latitude | 56.2°S | 56.0°S | 58.3°S | 58.4°S |
| Months detected | 55 % | 61 % | 31 % | 35 % |
| Trend, km per decade | −12.6 [−25.0, −0.2], p = 0.050 | −10.3 [−21.0, +0.5], p = 0.06 | not computed | not computed |
| r vs in situ, per transect / annual | 0.47 / 0.71 | 0.45 / 0.76 | 0.48 / 0.38 | 0.31 / 0.35 |
| Mean offset from in situ | −39 km | −35 km | −13 km | −12 km |

**Figure S2.** (a) Half-power method on a meridional KE profile. (b) Detection frequency of the matched envelope by longitude (90 % criterion dashed). (c) Recovered vs imposed shift in the synthetic experiments. (d) Climatological KE along the LM Gould track with the 200 and 400 cm² s⁻² thresholds and envelopes.

**Conclusion.** The circumpolar bound does not depend on the parameters within the tested ranges and means "no displacement larger than about 25 km over the record"; at Drake Passage the baseline measures a merged band, and its separation into fronts is post hoc and too intermittent for a trend.

---

## S4 Half-power meander gates

**Purpose.** Referee C asked whether meander flexing, a step of the Thompson and Naveira Garabato (2014) mechanism, can be measured from half-power positions. Two gates, written before any geometry trend was computed, decide where this is interpretable.

**Method.** Four windows chosen from bathymetry and climatology (Table S7). G1: one half-power branch followed continuously across ≥ 90 % of window longitudes. G2, at sites passing G1: synthetic tests (40 realisations; climatological KE plus observed anomalies resampled in 12-month blocks; linear change over 25 years; 150 km smoothing) of translation, 25 km flexing, uniform intensification and intensification of one flank (×1.2). A metric passes if translation is recovered and intensification produces < 50 % of the recovered translation.

**Results.** Only SEIR150 passes G1; a post hoc change of branch start rescues no failed site. At SEIR150 (Table S8) the mean position passes G2 narrowly (northern-flank intensification alone gives +10.0 km, 47 % of the recovered translation); amplitude, sinuosity and curvature do not separate flexing from intensification. The SEIR150 position trend is −4.3 km per decade [−14.3, +5.8], p = 0.37. Upstream jet-core speed leads northward displacement (r = 0.31 at lag 3 months; surrogate p = 0.004, passing FDR over 20 relation tests), but the peak lies at the edge of the pre-specified window and, given G2, may partly reflect redistribution of energy. Displacement is unrelated to meander-section (p = 0.89) and downstream EKE (p = 0.27). The DRAKE failure concerns a different measurement from Sect. S3: G1 needs one branch traceable across 70–55°W, whereas the section analysis uses one meridional profile per month along the LM Gould track.

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

**Figure S3.** (a) Percentage of window longitudes with a continuous branch in ≥ 90 % of months. (b–d) SEIR150 synthetic tests: recovered change in (b) mean position, (c) amplitude, (d) sinuosity, per scenario; dashes mark the imposed change.

**Conclusion.** Half-power geometry gives an interpretable mean position at one of four meanders and cannot detect a 25 km flexing, so meander shape is examined on the jet-core contour (Sect. S9).

---

## S5 Drake Passage sensitivity tests

**Purpose.** The Drake Passage results anchor the main text, and several lie near the threshold or depend on mask, section or series choices.

**Ice mask.** The mask removes 37 % of box pixels. Without it, box EKE rises by +7.1 % [−0.9, 15.1], p = 0.08 (+15.2 [−2.0, 32.4] cm² s⁻² per decade, mean 214) and KE by +2.0 % [−2.8, 6.9], p = 0.34; Product 4 EKE +2.8 % [−3.6, 9.3], p = 0.35. Ending in 2024, masked box EKE gives p = 0.057.
<!-- dk1_trends.csv (nomask rows); R2 A6 (37 %) and "Revised trend test" (2001-2024) -->

**EKE − KE contrast** (tested after it had been described): +6.0 points per decade [0.7, 11.3], p = 0.04 (Product 1, masked); +5.0 [0.5, 9.6], p = 0.04 (Product 1, no mask); +2.4 [−1.7, 6.4], p = 0.22 (Product 4, masked); +1.9 [−1.5, 5.2], p = 0.23 (Product 4, no mask).

**Section geometry.** An earlier proxy that dropped ice-covered points month by month correlated with the number of valid points (r = 0.89) and gave a spurious trend (+1.9 % per decade); it is withdrawn. The proxy is now integrated over fixed points ice-free in every month (T_perm), on the earlier line (65°W 55°S–60°W 63°S) and along the mean LM Gould track (GV track, primary). T_all integrates all points, including values mapped under ice (Table S9).

**Table S9.** Coverage-matched comparison with the in situ transport, 2005–2019 (248 transects in 137 months; annual means).
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

The northern points carry 63 % of the mean ADT drop across the earlier line. Averaged by month (n = 137) the north-part correlation is 0.66; deseasonalised, 0.50. Independent re-coding reproduced every correlation within ±0.02. The north–south anticorrelation is partly built in (a jet crossing the fixed latitude moves transport between the parts with no change in the total); it is consistent with the frontal compensation of Gutierrez-Villanueva et al. (2023) but does not test it. Seven line variants (end points moved by up to 1°) give T_perm trends of −0.2 to +4.8 % per decade (p = 0.13–0.96) and T_all trends of +0.4 to +1.6 % (p = 0.03–0.41), both under the earlier test [NUMBER NEEDED: revised p for the seven line variants].

**In situ quantities.** Total cross-track transport over 0–760 m (mean 75.5 Sv) = depth-independent reference part 51.1 + Ekman 8.1 + geostrophic relative to 760 m 15.9 Sv. Velocities: ADCP 70–118 m bins. EKE: depth-integrated 0–780 m [CHECK: 0–780 vs 0–760 m in the archive].

**Same-period trends, 2005–2019** (Table S10). No trend is significant, and the test is liberal at n = 15 (Sect. S2). The in situ EKE trend lies between the two products, but the intervals are several times wider than the differences, so the comparison does not favour either product.

**Table S10.** Trends over 2005–2019 (% per decade; positions in km per decade).
<!-- /home/claude/work/drake/out/dk3/dk3_summary.csv (EKE, annual); /home/claude/work/out/d1/d1b_trends.csv (transport, annual); /home/claude/work/drake/out/dk2/dk2_trends.csv and dk2_vs_insitu.csv (positions; altimetric series monthly) -->

| Quantity | In situ | Product 1 | Product 4 |
|---|---|---|---|
| EKE north of 58.5°S (annual) | +7.8 [−58.8, 74.4], p = 0.63 | +10.1 [−71.2, 91.3], p = 0.61 | +0.8 [−63.3, 65.0], p = 0.95 |
| Transport: in situ total / north / south; altimetric T_perm (GV, annual) | −1.2 [−12.7, 10.3], p = 0.70 / −0.9 [−24.4, 22.6], p = 0.87 / +0.2 [−53.9, 54.4], p = 0.98 | +4.9 [−20.9, 30.8], p = 0.49 | [NUMBER NEEDED: Product 4 T_perm annual 2005–2019] |
| Band position (in situ half-transport latitude; altimetric half-power latitude) | +4.1 [−29.4, 37.6], p = 0.80 | +19.4 [−46.3, 85.1], p = 0.51 | +17.6 [−84.3, 119.6], p = 0.72 |

**Wind lags.** Lagged wind–transport correlations were computed with the withdrawn proxy. [PENDING: wind–transport lag correlations with T_perm on the GV track.]

**Figure S4.** (a) Mean in situ 0–760 m total transport per 25 km bin, 2005–2019 (bars), and mean altimetric velocity normal to the GV track (line); dots, the ice-free points of the proxy; shading, ice fraction > 0.10 in at least one month. (b, c) Annual in situ transport north of 58.55°S (b) and over the full section (c) with the proxy, in standard-deviation units. (d) In situ north and south parts per transect against the proxy of the same month.

**Conclusion.** The box EKE trend and the EKE − KE contrast weaken without the mask or in Product 4; the surface proxy follows the in situ transport only over the latitudes it covers; and the 2005–2019 trends cannot discriminate between products.

---

## S6 Comparison with earlier estimates

**Purpose.** Reading our EKE trends as an update of Martínez-Moreno et al. (2021; MM21) and Hogg et al. (2015) requires separating period, product version and method from changes in the ocean.

**Method.** *MM21 procedure* (from their code repository): EKE = ½(ugosa² + vgosa²), ×ρ0, area-integrated (PJ m⁻¹), 365-day rolling mean, Theil–Sen with a modified Mann–Kendall test; their result south of 35°S is +5.1 % per decade (+0.0368 PJ m⁻¹ per decade), 1993 to March 2020, all-satellite vDT2018 (anomalies relative to 1993–2012, 1° grid). We applied it to our EKE over 35–65°S (pixels valid in every month, no ice mask, 12-month running mean of monthly values), and also our calibrated test to the unsmoothed series. *Hogg et al. (2015) boxes*: Indian 40–150°E, 57–44°S; Pacific 150°E–72°W, 62–48°S; Atlantic 35°W–10°E, 56–46°S (ice-masked). Theirs use crossover along-track EKE from 1993, relative to the 1993–2001 mean, 90 % intervals: Indian +6.0 ± 1.6 %, Pacific +6.4 ± 1.8 % per decade; absolute 18.3 ± 5.1, 14.9 ± 4.1 and 4.0 ± 3.7 cm² s⁻² per decade (Atlantic without a relative value). Crossover EKE is 1.6–1.9 times gridded EKE, so only relative trends are compared. *Scale matching*: Product 1 block-averaged to 1/4° (Table S3).

**Table S11.** Our data processed with the MM21 procedure, 35–65°S (% per decade).
<!-- plan1019 §2b; /home/claude/work/out/m2/m2_trends.csv -->

| Product, period | MM21 procedure | Calibrated test | Mean (PJ m⁻¹) |
|---|---|---|---|
| Product 1, 2001–2020 | +5.2 | +4.5 [−2.6, 11.6], p = 0.13 | 1.44 |
| Product 1, 2001–2025 | +7.5 | +6.6 [1.2, 12.0], p = 0.04 | 1.48 |
| Product 4, 2001–2020 | +1.8 | +1.2 [−1.5, 3.9], p = 0.36 | 1.28 |
| Product 4, 2001–2025 | +2.7 | +2.3 [0.2, 4.5], p = 0.04 | 1.29 |

With the ice mask, the 35–65°S integral rises by +0.092 PJ m⁻¹ per decade [0.022, 0.162] on a mean of 1.40 PJ m⁻¹ (Product 1) and +0.029 (Product 4), against MM21's +0.0368.

**Level difference.** Our integrated EKE level is about twice that of MM21 (1.44 against about 0.72 PJ m⁻¹, the value implied by their absolute and relative trends [AUTHOR CHECK: derivation of 0.72]). We have not investigated the cause. Candidate causes, unranked: (i) product version (DT2024 vs vDT2018); (ii) grid and order of operations (1/8° velocities vs coarsening to 1°, which removes sub-1° variance if velocities are averaged before squaring [AUTHOR CHECK: whether MM21 coarsen velocities or EKE]); (iii) anomaly definition (their 1993–2012 reference with daily fields; our 2001–2025 reference with monthly averages of daily values); (iv) region and masking (pixel validity, land and ice treatment); (v) normalisation (ρ0, area weighting). A constant level factor does not affect relative trends, but a time-varying one would.

**Table S12.** EKE trends by region and period (% per decade; absolute in cm² s⁻² per decade for 2001–2025).
<!-- plan1019 §2; /home/claude/work/out/m1/m1_trends.csv (Atlantic 2001-2012 and 2001-2020) -->

| Region | Period | Product 1 | Product 4 |
|---|---|---|---|
| Indian | 2001–2012 | +5.8 [−17.8, 29.3], p = 0.47 | +9.6 [1.7, 17.5], p = 0.03 |
| | 2001–2020 | +8.7 [1.7, 15.7], p = 0.03 | +5.7 [1.3, 10.2], p = 0.02 |
| | 2001–2025 | +9.0 [4.9, 13.0], p = 0.007; abs. +18.0 [9.9, 26.2] | +4.8 [1.7, 7.9], p = 0.013; abs. +8.3 |
| Pacific | 2001–2012 | +0.4 [−29.2, 30.1], p = 0.97 | +2.7 [−17.9, 23.3], p = 0.68 |
| | 2001–2020 | +5.8 [0.8, 10.8], p = 0.04 | +3.1 [−2.3, 8.4], p = 0.20 |
| | 2001–2025 | +8.6 [0.9, 16.4], p = 0.04; abs. +15.7 [1.7, 29.6] | +4.9 [0.5, 9.4], p = 0.05; abs. +7.5 |
| Atlantic | 2001–2012 | −5.1 [−24.1, 13.9], p = 0.50 | −3.1 [−21.0, 14.7], p = 0.65 |
| | 2001–2020 | +5.4 [−9.8, 20.5], p = 0.32 | +3.5 [−5.2, 12.2], p = 0.38 |
| | 2001–2025 | +8.8 [−0.1, 17.6], p = 0.056; abs. +13.5 | +5.2 [−0.2, 10.7], p = 0.06; abs. +6.9 |
| 35–65°S | 2001–2020 | +4.6 [−1.9, 11.1], p = 0.11 | +1.2 [−1.5, 4.0], p = 0.36 |
| | 2001–2025 | +6.6 [1.5, 11.6], p = 0.03 | +2.3 [0.2, 4.5], p = 0.04 |
| ACC band | 2001–2020 | +6.4 [0.9, 11.9], p = 0.04 | +3.2 [−0.1, 6.5], p = 0.07 |
| | 2001–2025 | +8.0 [3.3, 12.8], p = 0.017 | +3.8 [1.3, 6.3], p = 0.020 |

The 2001–2012 windows, which overlap Hogg et al.'s record, are uninformative (intervals about ±8 to ±30 % per decade; one of six sector trends significant, Indian in Product 4). For 2001–2025, Hogg et al.'s relative trends lie between our products and inside both intervals in the Indian and Pacific boxes: consistency, not reproduction, since periods, sampling and reference means differ. Extending 2001–2020 to 2001–2025 raises every trend and narrows every interval.

**Product versions.** MM21 used all-satellite vDT2018, no longer distributed; Barceló-Llull et al. (2025) found a large version sensitivity of EKE trends. Two-satellite vDT2018 remains available from C3S. [PENDING: two-satellite vDT2018 comparison, if extracted.] [AUTHOR DECISION: add a 1993–2025 two-satellite comparison only if the 1993–2000 extraction is completed.]

**Figure S5.** (a) EKE trends (95 % intervals) for 2001–2020 and 2001–2025, both products, in the Hogg et al. (2015) boxes, 35–65°S and the ACC band, with published estimates (diamonds; different period, version and method; Hogg et al. at 90 %); filled, p < 0.05. (b) Product difference. (c) EKE integrated over 35–65°S in percent of its 2001–2025 mean (MM21 procedure, no ice mask), monthly and 12-month running mean; shading, years after the MM21 period.

**Conclusion.** Processed as in the earlier studies, our trends have the same sign and order; over 35–65°S the 2001–2020 increase is not significant with the calibrated test and becomes so with the extension to 2025, and the twofold level difference with MM21 remains unexplained.

---

## S7 Antarctic-slope sea-level index (after Hughes et al., 2014, and Hogg et al., 2015)

**Purpose.** The surface proxy misses the southern Drake Passage section and the in situ record ends in 2019. Sea level over the Antarctic slope responds to the bottom-pressure (barotropic-type) part of the circumpolar transport, a quantity of a different kind from the three in main-text Sect. 3.4.

**Background.** In NEMO simulations (ORCA025, ORCA12), Hughes et al. (2014) found circumpolar transport correlated with bottom pressure over the Antarctic slope (region where the model correlation is below −0.5), averaged in 100 m depth bins over 0–5000 m with the strongest correlations at 1000–3000 m; their factor is −3.69 Sv mbar⁻¹ for bottom pressure, and sea level is amplified relative to bottom pressure at periods longer than about 5 years. Hogg et al. (2015) applied this to AVISO sea level over the shelf and slope (bias and annual/semi-annual cycles removed; low-frequency part h_f from 10, 20 and 40-year sinusoid fits; T = −3.69(h − h_f) − 1.11 h_f Sv cm⁻¹ [CHECK signs and units against Hogg et al., 2015]), called the factors preliminary and gave no numerical transport trend, only indications of a weak decline.

**Our construction.**
1. Product 1 monthly sea level south of 60°S, 2001–2025 [PENDING: extraction south of 65°S].
2. Slope region: 1000–3000 m depth south of 60°S, from bathymetry to 80°S [CHECK: bathymetry product]; this depth band replaces the model-correlation mask.
3. Each month, only pixels with ice fraction ≤ 0.10; the ice-free fraction of slope pixels is recorded.
4. Monthly mean per 100 m depth bin, then over bins; also in 12 longitude bands of 30° to monitor zonal coverage.
5. Removal of annual and semi-annual harmonics (h); h_f as in Hogg et al. (2015).
6. Conversion to T_hogg (Sv) with the published factors, labelled model-calibrated; h is reported before T_hogg.
7. Trends of h and T_hogg (Sect. S2); correlations with the in situ total, north and south parts (matched transects, annual means), with the depth-independent reference part [PENDING: if the archive allows], and with T_perm.

**Results.** [PENDING: h (cm) and T_hogg (Sv, model-calibrated) series.] [PENDING: trends of h and T_hogg with 95 % intervals and p.] [PENDING: correlations with in situ total, north, south and reference parts, and with T_perm.] [PENDING: monthly ice-free fraction, mean and distribution by longitude band.]

*Sampling caveat (not a result).* A smoke test of the scripts on the slope already extracted (60–65°S, 1000–3000 m) found only 47 % of slope pixels ice-free in a typical month, so sea-ice sampling is expected to dominate the full index further south. No value, trend or correlation from that test is used.

**Table S13.** Planned sensitivity tests.

| Test | What it checks | Result |
|---|---|---|
| Correlation of h with monthly ice-free fraction; fixed always-ice-free pixel set vs monthly set | Spurious variability or trend from varying coverage (cf. the withdrawn Drake proxy, Sect. S5) | [PENDING] |
| Minimum coverage per month and per longitude band | Dependence on poorly sampled months | [PENDING] |
| Ice-fraction threshold | Dependence on the mask | [PENDING] |
| Depth range (1000–3000 m vs wider bands within 0–5000 m) | Approximation of the slope region | [PENDING] |
| h_f construction (10/20/40-year fits vs other low-pass filters) | Low-frequency part, weighted differently by the conversion | [PENDING] |
| Conversion factors (ORCA025- vs ORCA12-based; no conversion) | Calibration uncertainty, largest for trends | [PENDING] |
| Removal of a large-scale sea-level signal (method to be specified) | Aliasing of steric or mass-driven margin sea-level change into "transport" | [PENDING] |
| Product 4 | Constellation dependence | [PENDING] |

**Interpretation (conditional).** If the index is well sampled and follows the in situ transport (particularly its south or reference parts), it extends a barotropic-type indicator to 2025; if sea-ice sampling or large-scale sea-level change dominates, it will be reported with its coverage and not interpreted. The Sv values remain model-calibrated in either case, with the largest uncertainty for trends.

**Figure S6.** [PENDING] (a) Slope region with mean monthly ice-free fraction. (b) h and T_hogg (model-calibrated) with the in situ total, south and reference parts and T_perm. (c) Monthly ice-free fraction.

**Conclusion.** [PENDING: one sentence on whether the index is usable and what it shows.]

---

## S8 Tracked eddies

**Purpose.** To document where eddies are generated and whether their nonlinearity changes, supporting the main-text reading of the EKE increase as more intense rather than more numerous eddies (Sect. 3.6).

**Method.** META4.0 DT all-satellite atlas (Product 5), south of 35°S, January 2001 to August 2023 (16.8 million observations, 524 000 tracks). Births per 1° cell per year. U/c = rotational speed over 7-day centre-displacement speed, for tracks of at least seven days.

**Results.** Births per year per 1° cell average 1.7 over the ACC band and 1.8 (PAR) to 2.6 (Kerguelen) in the wake boxes; upstream/meander/downstream sections give 1.5/2.0/1.6 (PAR), 1.5/1.1/2.1 (SEIR150) and 4.0/3.6/2.7 (Kerguelen). Births are enhanced along the whole energetic band, particularly at Kerguelen and Drake Passage, and the meander sections are not singular generation sites at 1°. Wake-minus-rest differences in normalised trends (points per decade) are +0.3 for eddies present (p = 0.71), −1.9 for births (p = 0.21), −0.03 for amplitude (p = 0.97) and −1.9 for Σ amplitude² (p = 0.55). Monthly births and path length show no consistent link (PAR p = 0.06, r ≈ −0.2; SEIR150 p = 0.16; KERG p = 0.74). U/c has median 6.9, 99 % of observations exceed 1, and no detectable trend (p = 0.97); U and c increase by about 3.5 and 3.2 % per decade.

**Caveats.** Only eddies above the mapping resolution are tracked (speed-radius median about 37 km, 5th percentile 24 km), so an unchanged radius is weakly informative; the atlas shares the constellation dependence of Product 1 (Sect. S1).

**Figure S7.** (a) Births per year per 1° cell, 2001–2023, with ACC band, wake boxes and 3000 m isobath. (b) U/c distribution and monthly median series.

**Conclusion.** Eddy generation is spread along the energetic band rather than concentrated at the meanders, and the nonlinearity of tracked eddies shows no detectable trend.

---

## S9 Along-stream analysis: details

**Purpose.** To document the choices behind main-text Sect. 3.7, the monthly covariation tests and our reproduction of Liu et al. (2026).

**Path and bins.** Monthly contour of ADT relative to its ACC-band mean at +0.18 m, the level of highest climatological speed (Subantarctic Front side); the speed–level curve is broad (0.16–0.22 m within 0.5 %), so the level is weakly determined. Smoothing 150 km; 72 bins of 5°. The bins at 60–50°W are excluded because the contour follows the Malvinas Current (sinuosity coverage 0.69). Per-bin tests use the Sect. S2 statistic with 300 null series (single bootstrap for intervals) and FDR at 10 %. The EKE–sinuosity cross-bin correlation (r = 0.29, p < 0.001) is not used in the main text because EKE and monthly-path sinuosity both contain the eddies.

**Monthly covariation.** A pre-specified box-based test at PAR, SEIR150 and KERG found weak, episodic covariation, with path lengthening leading downstream EKE by 1–3 months (pooled p = 0.007), carried at SEIR150 by a few episodes; no upstream-speed → path-length link was found. A box-free version (36 overlapping 30° windows; downstream EKE 15–35° east), corrected after its initial pass for window overlap, birth normalisation and contour level, shows simultaneous path-length–EKE covariation in 22 of 36 windows (lag 0 strongest) and a downstream association in 1 of 36 (Polar Front side, −0.2 m: 14 and 3 of 36, mixed signs). What holds is local, simultaneous covariation, as expected if eddies form the path distortions; a delayed downstream response is not detected.

**Wake test details.** Wake boxes cover 23 % of the band area and 43 % of its mean EKE. In the pre-specified sensitivity period 2001–2017 the wake-minus-rest EKE difference is +4.1 points per decade (p = 0.13; Product 4 +4.3, p = 0.06); no secondary localisation test passes FDR. At Drake, upstream KE rises more than downstream (Scotia) KE (p = 0.010), with EKE at the threshold (p = 0.047–0.053).

**Reproduction of Liu et al. (2026).** Following their Methods, we computed daily |∇ADT| from Product 1 (1993–2025, 1/4° block-averaged), its 4-month occurrence above 20 % of the daily maximum, and the axis, speed and section positions. Their jet-core acceleration is reproduced in every section to their two-decimal rounding and survives the calibrated test except at the Upper Ridge; it also appears in the 2001–2025 monthly fields of both products. Their position trends are reproduced in sign and approximately in magnitude, but only the Lower-Ridge northward shift approaches significance with an autocorrelation-aware test (1993–2023 p = 0.047; 1993–2025 p = 0.10; 2001–2025 p = 0.24). The difference lies in the test: on trend-free AR(1) series (n = 372) their modified Mann–Kendall test rejects at nominal 5 % in 63 % (φ = 0.5) to 98 % (φ = 0.8, 4-month running mean) of cases, because its variance uses n_eff from the undetrended series. This concerns significance only, not their trend magnitudes. We could not reproduce their width (about 650 km zero-to-zero against their 100 km; we may have misread the definition), so no width comparison is made.

**Figure S8.** (a) Mean jet-core path with bins, obstacles and excluded Malvinas bins. (b) Same-window and downstream path-length–EKE correlations by window (SAF and PF levels). (c) PAR jet-core speed by section, this study vs Liu et al. (2026).

**Conclusion.** The along-stream diagnostics rest on a weakly determined but stable contour level, the monthly tests support only local, simultaneous path–EKE covariation, and the reproducible signal at the Pacific–Antarctic Ridge is jet-core acceleration.

---

## S10 Southern Annular Mode

**Purpose.** Referee A asked whether KE and EKE variability follows the SAM.

**Method.** Deseasonalised, detrended monthly ACC-band KE and EKE (Product 1) correlated with the SAM index of Marshall (2003), SAM leading by 0–12 months; significance from phase-randomised surrogates of the index, accounting for the lags examined.

**Results.** Maximum r = 0.15 for KE (lag 0, p = 0.12) and 0.09 for EKE (p = 0.24); 12-month running means r = 0.11 and 0.08; annual means 0.04 and 0.05.

**Figure S9.** Lagged correlation of ACC-band KE and EKE with the SAM index, with the 95 % surrogate envelope.

**Conclusion.** Circumpolar KE and EKE variability shows no detectable relation to the SAM index over 2001–2025; a lagged regional response to individual SAM events (Meredith and Hogg, 2006) is not expected to appear in this correlation.
