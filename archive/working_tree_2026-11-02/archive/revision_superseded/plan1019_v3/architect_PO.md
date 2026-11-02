# v3 scientific architecture: a physical oceanographer's view

Sources: [DR] drake_results.md · [P19] plan1019_results.md · [R2] R2_revision_results.md ("Revised trend test" section) · [MS] manuscript_v23.md · [AC-A/B/C] replies v2. Where sources disagree, the most recent is used and the conflict is flagged (⚑).

---

## 1 Central question and defensible findings

**Physical question.** Over 2001–2025, as the westerlies over and upstream of Drake Passage strengthened, did the surface eddy field intensify while the ACC's surface transport and pathway changed little (the surface expression expected under eddy saturation)? And is altimetry good enough to extend the one in situ transport record to answer this question circumpolarly?

**What the data show (defensible, with caveats)**

1. *Altimetry reproduces in situ interannual variability, for some quantities only.* Over 2005–2019: annual EKE north of 58.5°S r = 0.89; surface velocity per bin r = 0.72; northern-section transport r = 0.67–0.72; band position r = 0.71. The full-section transport is **not** reproduced (r ≈ 0.05) [DR §3]. This supports extending EKE and the northern-section transport. It does not support extending full-section transport or checking trends (15 years).
2. *Around Drake Passage the robust change is upstream.* At 100–80°W, KE rises +12.8 %/dec (p = 0.007), EKE +15.3 % (p = 0.013) and wind stress +9.8 % (p = 0.03), in both products. The along-stream speed rises +3.9 % (p = 0.050; two-sat p = 0.056), at the limit. Inside the box, KE (+2.4 %, p = 0.29) and along-stream speed (−0.1 %) are flat. EKE (+8.4 %, p = 0.050) is at the threshold in all-sat only (two-sat +3.8 %, p = 0.22) [DR §1].
3. *No detectable trend in the Drake transport indicators or the frontal band.* Northern surface proxy: +3.9 %/dec [−2.1, 9.9], p = 0.17 [DR §3]. In situ total: no trend over 15 annual values ⚑. Section half-power point: −10 km/dec [−36, +16], p = 0.43 [DR §2]. These are wide intervals, not evidence of no change.
4. *Circumpolar update.* In the ACC band, KE rises +4.9 %/dec (p = 0.003) and EKE +8.0 % (p = 0.017). In two-sat, EKE rises by about half as much (3.8 %), but the difference between products is not significant (p ≥ 0.09) [P19 §2]. Along-stream speed rises +2.0 % and the surface transport proxy +1.6 %, both robust to the product [MS, AC-A]. These results agree in sign and order with MM21 and Hogg et al. (2015). Extending the record to 2025 is what makes the 35–65°S increase detectable with the conservative test (2001–2020 p = 0.13 → 2001–2025 p = 0.04) [P19 §2b].
5. *Pathway geometry is stable within bounds.* Circumpolar half-power position −0.5 km/dec (−4.2, +3.1), with 25 km over the record detectable. Sinuosity and latitude along the jet-core path show no tendency [MS].

**What the data merely allow (keep as interpretation, never as result)**
- Eddy saturation as the explanation. The circumpolar forcing trend is **not detected** (45–60°S: +2.5 %, p = 0.08), so a flat transport cannot be read as insensitivity to wind [R2 A6].
- Frontal compensation (GV2023) as the cause of the flat in situ total.
- The TNG14 chain. Only the first link (faster jet) is seen, and only partially.
- Any product preference for the EKE trend magnitude.
- Equilibrium versus slow adjustment (Böning 2008).

---

## 2 Section-by-section outline (binding Drake-first order)

Suggested word budget (≤ 5000): Intro 600 · Methods 1300 · Results 2100 · Discussion 900.

### 1 Introduction
- State the eddy-saturation prediction and its observable surface expression.
- State the OSR purpose: an update that is informative whether it confirms, weakens or changes earlier results [plan §1].
- Cite earlier evidence (Hogg 2015; MM21; Zhang 2021). Name what v3 adds: the Drake anchoring with an explicit altimetry–in situ consistency check, the Antarctic-slope transport index, the extension to 2025 with a calibrated test, and the mean flow alongside EKE.
- No submesoscale.
- Introduce the mechanisms neutrally: baroclinic, wakes/barotropic, standing meanders.

### 2 Data and methods
- 2.1 Products; grid vs effective resolution; KE/EKE/along-stream speed definitions (EKE includes low-frequency meander variability).
- 2.2 Regions: Drake box, upstream, downstream, LM Gould track, ACC band.
- 2.3 Half-power point.
- 2.4 In situ record and consistency-check metrics.
- 2.5 Transport quantities: surface proxy (m² s⁻¹), in situ 0–760 m, Hogg/Hughes index (Sv, model-calibrated).
- 2.6 Along-stream path and tracked eddies.
- 2.7 Statistics.

Each subsection opens with the question it serves.

### 3.1 Drake Passage: KE, EKE, along-stream speed, wind (box/upstream/downstream)
- **Question:** has the eddy field intensified relative to the mean flow where local forcing increased?
- **Consistent with ES:** EKE up with flat mean-flow speed.
- **Inconsistent:** along-stream speed/KE rising in proportion to the forcing with no EKE change.
- **Report:** the table in [DR §1], all-sat only, in % and **absolute** units (⚑ the absolute upstream values are not in the files and should be computed, because the upstream mean is low). Add one sentence that the in-box EKE trend is not confirmed by two-sat. Report the wind for the three regions.
- **Honest reading:** inside the passage the pattern is weakly ES-like. Upstream, both mean flow and eddies intensify. This is the first TNG14 link *or* simply part of the Pacific-wide intensification (the circumpolar surface proxy's significant bins cluster at 145–90°W [MS 3.3]). It cannot be resolved here.
- **SI:** two-sat; no-mask Drake EKE (p = 0.083 [R2]); EKE−KE contrast; 2001–2024 sensitivity.

### 3.2 Front position along the LM Gould section (half-power)
- **Question:** has the energetic band in the passage moved while energy changed? This is not a transport test (B1: the f-change effect is small).
- **Consistent with ES:** a stationary band.
- **Inconsistent:** none per se. A shift would mean the box trends are relocation, not intensification.
- **Report:**
  - At baseline there is one merged SAF+PF envelope: 55.1–59.3°S, half-power point 57.0°S, detected in 99.7 % of months.
  - Trend −10 km/dec [−36, +16].
  - Against in situ: r = 0.50 per transect and 0.71 annual; offset 16 km; RMS 114 km [DR §2].
- **Must be stated:** this is the energy centroid of the merged band, not a front.
- **SI:** the 400 cm² s⁻² split (post hoc; SAF detected in 55 % of months, PF in 31 %); two-sat; the sensitivity range.
- ⚑ State whether the profile uses points south of the 58.5°S ice-free limit. The envelope reaches 59.3°S.

### 3.3 Altimetry vs in situ (purpose: extend the record)
- **Question:** does altimetry reproduce the in situ variability well enough to extend it to 2001–2025 and beyond Drake?
- **No ES test here.**
- **Report:** a compact Table 2, "what the check supports" [DR §3]:
  - Mean velocity profile r = 0.98.
  - Per-bin velocity r = 0.72 (8,187 pairs).
  - Annual EKE r = 0.89.
  - Northern transport r ≈ 0.7.
  - Position r = 0.71.
  - Full-section transport r ≈ 0.05; north–south r = −0.55.
  - Trends: not checkable.
- **Then the extension, 2001–2025:** EKE north of 58.5°S +7.4 % [−5.6, 20.4], p = 0.24; northern transport +3.9 %, p = 0.17.
- **Wording:** "comparison/consistency check", since the quantities differ (surface vs 0–760 m; ADCP 70–118 m bins).
- **Drop from the main text:** the same-period trends (+7.8 / +10.1 / +0.8 %, none significant) cannot adjudicate products. SI only, without the "closer to all-sat" reading.
- **SI:** line variants, T_perm vs T_all, monthly-averaged correlations (n = 137, r = 0.66) [P19 §1], deseasonalised check.

### 3.4 Transport
- **Question:** has the ACC transport changed as eddy energy rose?
- **Consistent with ES:** a trend in transport indicators that is small relative to the forcing change. This is weakly posed, because the circumpolar forcing trend is not detected.
- **Inconsistent:** a significant transport increase scaling with the wind.
- **Keep three quantities distinct** [plan §3]:
  - (a) The surface proxy, northern partial section (53/121 points, north of 58.47°S, 65 % of the full-line proxy) [P19 §1].
  - (b) The in situ 0–760 m **total**: mean 75.5 Sv = reference 51.1 + Ekman 8.1 + relative 15.9 [P19 §0]. Also its north and south parts (trends −0.9 % and +0.2 %, both n.s.).
  - (c) T_hogg: methods in full (slope region and depth bins after Hughes; h_f fit; 12 longitude bands; monthly ice coverage). **[PENDING: index series, T_hogg (model-calibrated, Hughes et al. 2014 factors), its trend with CI, its correlation with in situ total/north/south on matched transects and annual means, and ice-free fraction]**.
- **Physical pairing:** the slope index is a bottom-pressure/barotropic-type quantity. Its closest in situ counterpart is the depth-independent (reference) part. If the archive allows, add T_hogg vs the reference part [PENDING]. Also compare with the south part, which lies on the same (Antarctic) side.
- **SI:** conversion-factor sensitivity, filtering, averaging, ice sampling.

### 3.5 Circumpolar extension
1. **KE/EKE/along-stream speed and surface transport proxy in the ACC band.**
   - Question: is the Drake pattern regional or basin-wide?
   - Report the ACC values in §1.4 and the circumpolar proxy +1.6 % (Pacific 2.6, Atlantic 1.8, Indian 0.4 n.s.).
   - Do **not** rank sectors (ranking not supported, p = 0.05–0.64 [R2]).
2. **Comparison with earlier estimates.**
   - Question: is this an update consistent with earlier work?
   - Report:
     - MM21 procedure on our data: +5.2 % (2001–2020) and +7.5 % (2001–2025), vs MM21 +5.1 % (1993–2020).
     - Calibrated test: +4.5 % (p = 0.13) → +6.6 % (p = 0.04).
     - Hogg sector boxes 2001–2025: Indian +9.0 and Pacific +8.6 % (all-sat), vs Hogg 6.0 ± 1.6 and 6.4 ± 1.8 [P19 §2, §2b].
   - Give absolute trends too. Disclose that our integrated level is twice MM21's (1.44 vs about 0.72 PJ m⁻¹), not investigated ⚑. Check whether MM21 coarsened velocities before squaring.
3. **Fronts.**
   - Question: is the band intensifying in place?
   - Report the bound −0.5 (−4.2, +3.1) km/dec and the Chambers/Sokolov–Rintoul/Gille/Liu comparison.
   - State coverage: 900 longitudes (about 31 %). The gaps include 115–75°W, i.e. the Drake upstream box [R2 A2].
4. **Along-stream/standing meanders.**
   - Question: is the eddy increase localised at standing meanders, and do they flex?
   - Consistent with TNG14: localisation, flexing, upstream acceleration.
   - Report EKE positive in 58/70 bins (0 after FDR), with similar values over/downstream/elsewhere (7.0/7.3/6.9 %); sinuosity and latitude show no tendency; wake test +0.1 %, p = 0.95; PAR core speed +0.019 m s⁻¹/dec (p = 0.01), upstream p = 0.06.
   - Report the cross-bin r = 0.40 (speed) with its caveat.
5. **Tracked eddies (short).**
   - Question: number or intensity?
   - Report amplitude +4.6 % and rotation +3.8 %; number (+0.8 %) and radius unchanged.
   - Caveat: same maps, so the same constellation sensitivity.

### 4 Discussion
1. What Drake shows and what the check licenses.
2. The forcing problem (circumpolar stress not detected).
3. Surface vs depth-integrated transport; Böning/adjustment time scales.
4. Mechanisms: which TNG14 links are observable, and the upstream-vs-Scotia contrast.
5. Limitations and outlook.

---

## 3 Figures (5 main; 4 if Hogg stays pending)

| Fig. | Content |
|---|---|
| 1 Drake setting & trends | (a) Map: mean EKE, 3000 m isobath, box/upstream/downstream, LM Gould track (ice-free part solid), mean fronts. (b) Forest of % trends: KE, EKE, speed, wind × 3 regions (all-sat). (c) Box/upstream EKE & KE 12-month series |
| 2 Altimetry vs in situ | (a) Mean cross-track velocity profile, alt vs ADCP. (b) Annual EKE north of 58.5°S, in situ 2005–19 over altimetry 2001–25. (c) Half-power point vs in situ half-transport latitude, extended. (d) Transport: in situ total/north/south vs northern proxy, showing the north–south opposition |
| 3 Transport | (a) Slope-index map with monthly ice-free fraction. (b) Index and T_hogg (labelled model-calibrated) vs in situ total/south/reference parts and the northern proxy. [PENDING]. If pending at submission, merge into Fig. 2d |
| 4 Circumpolar energy & earlier estimates | (a) EKE trend map with FDR cells, isobath, Hogg boxes. (b) Comparison forest: this study (2001–2020/2025) vs MM21 and Hogg, relative and absolute. (c) ACC band KE/EKE/speed/surface proxy trends; eddy amplitude/number |
| 5 Pathway & meanders | (a) Half-power position trends vs earlier studies. (b–e) Along-stream EKE, jet speed, sinuosity, latitude per 5° bin with obstacle shading |

**SI:**
- S1: two-satellite product, all diagnostics.
- S2: trend-test calibration.
- S3: half-power method, parameters 150–400/20–50 %/50–300 km, detectability, longitude bins, the 400 split at Drake, gaps, withdrawn width.
- S4: half-power meander gates G1/G2 and SEIR150.
- S5: Drake sensitivity (masks, 7 line variants, T_all, wind lags, coverage matching).
- S6: MM21 procedure, scale-matched (1/4°) comparison, 2001–2012 windows.
- S7: Hogg index construction and sensitivity.
- S8: eddy births map, U/c.
- S9: along-stream details (Malvinas exclusion, PF-level run, monthly covariation B2, Liu reproduction).
- S10: SAM correlation.

---

## 4 Referee check

| Referee point | Answered by | Status |
|---|---|---|
| C: novelty vs MM21; Fig. 1 similar | 3.5(2) matched comparison; Intro framing as OSR update | Answered |
| C: half-power and Drake comparison weakest | 3.2–3.3: section half-power checked against in situ (r = 0.71 annual); coverage-matched transport | Answered, though the half-power point remains a merged-band centroid |
| C: abstract lacks results; "direct observational evidence" | Rewrite with numbers; phrase removed | Answered |
| C: no transitions; what supports/rejects ES | Question/consistent/inconsistent at each subsection (above) | Answered by design |
| C: TNG14 chain untested (upstream acceleration, flexing, EKE, slowing) | Drake upstream (3.1), PAR speed, along-stream sinuosity (3.5(4)); half-power curvature tried → SI S4 | Partly. Flexing is not detectable by half-power and not seen in sinuosity. **"Jet slowing" is not addressed explicitly**; say it is not observed |
| C: curvature at specific meanders from half-power | Gates: only SEIR150 passes; flexing undetectable (8.3 km vs null 29.4) | Answered negatively (SI) |
| C: 200 threshold rationale and sensitivity | Methods + S3 | Answered |
| C: averaging hides meander variability; Pacific basin size; upstream/inside meanders | Along-stream bins, wake test, Drake up/down | Answered |
| C: altimetry the right dataset; resolution; grid spacing | Grid vs effective resolution; **consistency check with in situ is the new core argument**; reanalysis declined (binding) | Answered, but the declined suggestion needs one firm sentence |
| C: width ~3 km undetectable; small latitude trends | Width withdrawn; positions as detectability bounds | Answered |
| C minors (submeso, early ES papers, p < 0.05, bathymetry contours, gaps, Schiavon, wind variable, refs) | Text/figure fixes | Answered |
| C: **report cut off** | — | **Unanswered. Obtain the full report before submission** |
| A: start from Drake, extend; sector surface transport | Structure; circumpolar proxy | Answered |
| A: vertical structure | Limitation; Hogg index adds a barotropic-type indicator [PENDING] | Partly (pending) |
| A: resolution / range of eddies | Methods; eddy atlas | Answered |
| A: fronts–EKE link; instability indicators | Along-stream; no subsurface data, stated | Answered within limits |
| A: SAM, KE definition, abstract | S10; 2.1 | Answered |
| B1: fronts vs transport | Stated as not a transport test | Answered |
| B2/B5: proxies oversold; lines 80–81 | Three distinct quantities; coverage-matched correlation | Answered. **⚑ AC-B5 and MS 3.3 must be updated** ("uncorrelated r = 0.01/−0.04" is superseded) |
| B3: Böning, adjustment | Discussion | Answered |
| B4: wakes/barotropic vs baroclinic | Wake test; Drake upstream > Scotia (KE p = 0.010 [R2 A7]); "cannot distinguish" | Answered |

The AC letters (especially the closing of AC-C, with its old four-figure order) must be rewritten to match v3.

---

## 5 Physics pitfalls

1. **Surface vs depth-integrated.** Eddy saturation concerns the full-depth (baroclinic) transport. Surface speed (+2.0 %) and proxy (+1.6 %) increases neither confirm nor refute it. The ratio "EKE five times faster than transport" is not a saturation metric [plan]. Describe it, don't test with it.
2. **The forcing must be detected first.** The Drake-box stress (+8.7 %) is local. ACC transport responds to circumpolar-integrated stress, which is not significant (p = 0.08).
3. **EKE includes low-frequency meander/jet variability** (low-frequency part +5.6 %/dec [R2]). "EKE up, KE flat" can partly reflect path variability, not eddies. Prefer along-stream speed as the mean-flow measure.
4. **Upstream intensification ≠ the TNG14 first link.** EKE rises with KE upstream, so it is not mean-flow-only. It sits within the Pacific-wide signal. Relative trends are inflated where the mean is low (report absolute values), and half-power detection fails there.
5. **The in situ "0–760 m" is a total transport** (reference 51.1 + Ekman 8.1 Sv). It is not geostrophic relative to 760 m. ⚑ In situ EKE is labelled 0–780 m in [DR]; check.
6. **Hughes factors are model-derived** (NEMO), "preliminary", and frequency-dependent. Sea level is amplified at periods longer than 5 years, and −3.69 is per mbar of bottom pressure [P19 §0]. A trend is the lowest frequency, so the conversion there is the least constrained.
   - Antarctic-margin steric/freshening or mass-change sea-level trends alias into "transport".
   - Ice sampling dominates: 47 % of slope points are ice-free in a typical month [P19 §4].
   - Report the index first and Sv second. Check the sign convention against the PDF.
7. **15-year overlap:** a consistency check for variability only. The annual in situ trend test has a 9.5 % false-positive rate [R2].
8. **North/south anticorrelation** (r = −0.55) is partly built in: a jet crossing a fixed latitude moves transport between the parts. It is consistent with GV2023's compensation but does not test it.
9. **Merged envelope:** the baseline section half-power merges SAF+PF. The 400 split is post hoc with 31 % PF detection. The "Drake" meander window fails G1 (14 %) for along-longitude tracking. This is a different test from the along-section profile; explain the difference to avoid apparent contradiction.
10. **p = 0.050 is not p < 0.05** under the single-level rule. Drake EKE, upstream speed and the SAF-400 trend must be described as "at the threshold", not significant.
11. **Common mapping changes** (post-2018 constellation) can raise EKE and resolved jet speed together. The cross-bin r = 0.40 may partly reflect this; show the two-sat r in SI.
12. **The eddy radius is mapping-truncated** (5th percentile 24 km), so "unchanged size" is weakly informative.

---

## 6 Reconsider or cut

- **Cut** the Drake EKE−KE contrast as a headline. It is all-sat only, two-sat p = 0.22, and mask-dependent.
- **Cut** "EKE grows most where the jet accelerates most" and the sinuosity–EKE r = 0.29 from the abstract. Sinuosity of the monthly path contains the eddies. Keep both in Results, hedged.
- **Shrink** the circumpolar half-power to one paragraph: the bound plus the earlier-study comparison. It covers 31 % of longitudes and is no longer central.
- **Shrink** the tracked eddies to 3–4 sentences, or move them to SI if words are short. They are the least necessary for the Drake-first story.
- **Move** SAM, U/c, births, B1/B2 lead–lag and the Liu-test critique to SI. Keep the critique collegial.
- **Do not use** the 2005–2019 trend comparison to favour all-sat.
- **Reconcile before writing** ⚑:
  - In situ trend 2005–2019: −2.8 % (p = 0.47) [MS] vs −1.2 % (p = 0.70) [P19].
  - 244 vs 248 transects.
  - Mean 75.5 vs 75.7 Sv.
  - Full-section r: 0.03 vs 0.05 vs 0.11.
  - North–south r: −0.54 vs −0.55.
  - Drake EKE CI upper bound: 16.6 [DR] vs 16.7 [MS].
- **Contingency for Hogg:** if ice sampling makes the index unreliable, keep the methods and the index with its coverage in SI, and give one main-text sentence. Do not force a Sv trend into the abstract.
- **Do not promise 1993–2025** unless the extraction is done.
- **Consider a title** that signals the Drake anchoring and the consistency framing. Optional.
