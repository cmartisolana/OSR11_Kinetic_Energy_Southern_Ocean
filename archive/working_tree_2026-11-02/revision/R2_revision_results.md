# Reviewer 2 revision: results so far (Option 2)

Status on 2026-10-13. **The trend test was replaced on 2026-10-13** (section "Revised trend test" near the end). p-values quoted in the sections below are from the earlier calibrated test unless marked "revised"; the revised value of every reported trend is in `outputs/retest/retest_summary.csv`, and the conclusions that change are listed in that section.

All numbers come from the scripts in `revision/` and their outputs in `revision/outputs/`. The ocean data are satellite altimetry only:
- all-sat: SEALEVEL 008_047, DT2024, 1/8°;
- two-sat: 008_057, 1/4°.

Other inputs are WIND 012_006 (monthly `eastward_stress`), GEBCO 2026 bathymetry (used as context only) and the in-situ LM Gould ADCP transport.

## Fixed choices

These were set before any trend was computed.
- **D1 baseline:** 200 cm² s⁻², minimum width 100 km, 25 % fraction. The detection field is the KE of the monthly-mean velocity.
- **ACC band:** 2001–2025 mean ADT between −0.6 and +0.2 m.
- **Ice:** a pixel is dropped if its monthly ice fraction exceeds 0.10 in any month.
- **Trends:** deseasonalised Theil–Sen slope. The CI and p-value come from an AR-sieve bootstrap (AIC, p ≤ 12) with a double-bootstrap calibration (`p_cal`). Hamed–Rao is reported for comparison only. Maps use BH-FDR at α = 0.10.
- **Synthetic tests of the nominal 5 % level:**

  | Method | False-positive rate |
  |---|---|
  | Moving-block bootstrap | 13–18 % (rejected) |
  | Sieve | 7 % (φ = 0.5), 13 % (φ = 0.9) |
  | Hamed–Rao | 25–32 % |

  The calibrated version was the reference until 2026-10-13; see "Revised trend test" (it was still 6–11 % at nominal 5 %).
- **3-year blocks (8 values):** OLS t-test on independent blocks, as in Chambers (2018). An AR sieve cannot be fitted on 8 points; the first run of a5 gave a meaningless p = 1.0, which is now fixed.
- **A5d boxes:** fixed from bathymetry and climatological EKE before any box trend was computed.
- **B1:** sections from Liu et al. (2026, Table 1). The gate was set in advance.

## A2/A3: frontal position and width (D3)

**Position (matched fixed sample: 900 longitudes, coverage ≥ 90 %, 97 % on average)**
- ACC half-power latitude: −0.0005 ° yr⁻¹, 95 % CI [−0.0032, 0.0022], p = 0.72.
- No sector trend:

  | Sector | Trend (° yr⁻¹) | p |
  |---|---|---|
  | Atlantic | +0.0018 | 0.46 |
  | Indian | −0.0020 | 0.37 |
  | Pacific | −0.0011 | 0.62 |

- Every parameter set gives p > 0.07.
- The v2-style unmatched median latitude gives Indian −0.014 (p = 0.013) and Pacific −0.011 (p = 0.003). This metric has ≤ 27 % power even for a 50 km shift. Amplitude growth alone (S3) produces −0.002 to −0.007 ° yr⁻¹ in it.
- **By longitude:** 4 of 33 bins have p < 0.05, and none survive FDR.
- **Against Chambers' regional values:**
  - 90–170°E: −1.3 km decade⁻¹ (Chambers −29).
  - 200–270°E: +3.7 (Chambers +46).
- **Detectability (matched metric):** a 10 km shift over 25 yr is detected 44 % of the time; a 25 km shift 100 % of the time. About 70 % of an injected shift is recovered.
- **3-year blocks (8 values, OLS):** ACC p = 0.17. Indian −0.0068 ° yr⁻¹, p = 0.049. This single borderline value does not survive any multiplicity consideration.

**Width: remove the broadening claim (D3)**

| Metric | Trend (km yr⁻¹) | p |
|---|---|---|
| Segment width, matched | +0.47 | 0.002 |
| Segment width, unmatched | +0.14 | — |
| FWHM | −0.10 | 0.005 |
| W25 | −0.12 | — |
| σ | −0.03 | — |

- S3 (amplitude-only growth) produces segment-width trends of +0.43 to +1.6 km yr⁻¹. The threshold-defined width therefore tracks amplitude.
- Resolution sharpening of 5 km gives FWHM trends of −0.05 to −0.25, so the narrowing cannot be interpreted.
- FWHM has ≤ 35 % power even for large broadening.

**Detection frequency** is low (25–45 %) at about 170–150°W and 115–75°W. This explains the gaps in Fig. 2a.

## A4/A5: KE, EKE and their parts, both products (D4, D5)

ACC band, 2001–2025, calibrated trend [95 % CI], % decade⁻¹:

| Term | All-sat | Two-sat | Ratio all/two |
|---|---|---|---|
| KE (monthly-mean flow) | +1.36 [0.87, 1.85], +4.9 % | +0.97 [0.56, 1.39], +3.6 % | 1.4 |
| EKE (fixed reference) | +1.98 [1.07, 2.89], +8.0 % | +0.81 [0.30, 1.31], +3.8 % | 2.5 |
| EKE, low-frequency part ½\|Uₘ−ū\|² | +0.85, +5.6 % | +0.46, +3.2 % | 1.8 |
| EKE, sub-monthly part ½σ²ₘ | +1.14, +12.0 % | +0.34, +5.0 % | 3.4 |
| EKE w.r.t. 36-month running mean | +1.97, +8.3 % | +0.82, +4.1 % | 2.4 |
| Along-stream speed Uₘ·ū/\|ū\| (cm s⁻¹ yr⁻¹) | +0.025 [0.018, 0.033], +2.0 % | +0.024 [0.016, 0.032], +1.8 % | 1.05 |
| KE of annual-mean flow | +0.52, +3.5 % | +0.50, +3.4 % | 1.0 |
| KE of 3-yr-mean flow (OLS, 8 blocks) | +0.53, +4.0 %, p = 0.008 | +0.48, +3.6 %, p = 0.012 | 1.1 |

All terms are significant (p_cal ≤ 0.01) in both products.

**What this means**
- The surface **mean flow strengthens** by about 2 % decade⁻¹ in along-stream speed and about 3.5 % decade⁻¹ in the KE of the annual and 3-year mean flow. The two products agree to within 5–10 %.
- The **EKE trend is product-dependent.** All-sat is 2.5 × two-sat, mostly through the sub-monthly part (3.4 ×). The two products diverge mainly after about 2018 (Fig. S-a5). This is consistent with the time-varying constellation documented in the 008_047 QUID and by Hargous et al. (2026).
- The **sign** of the EKE increase holds in both products.
- **"EKE trend > KE trend" holds only in all-sat.** In two-sat, KE ≥ EKE.
- v2's "decoupling robust to the product" and "EKE ≈ 1.5 × KE" therefore do not hold.

**Sectors (ACC band)**

| Sector | All-sat EKE | Two-sat EKE (p) | All-sat KE | Two-sat KE (p) |
|---|---|---|---|---|
| Pacific | +2.35 | +1.21 (0.02) | +1.87 | +1.53 (0.01) |
| Atlantic | +2.14 | +0.90 (0.11) | +1.27 | +0.89 (0.05) |
| Indian | +1.33 | +0.23 (0.41) | +0.92 | +0.41 (0.09) |

- Along-stream speed: Pacific +3.5 % decade⁻¹, Atlantic +1.8 % (both products). The Indian sector is not significant.
- In the ACC band, the Pacific has the largest trend in both products and for both KE and EKE.

**Comparison with Martínez-Moreno et al. (2021): box south of 35°S, 2001–2020**
- All-sat EKE +0.91 cm² s⁻² yr⁻¹ (+4.6 % decade⁻¹, p_cal = 0.06).
- Two-sat +0.21 (p = 0.23).
- Over 2001–2025: all-sat +1.33 (+6.6 %), two-sat +0.42 (+2.3 %, p = 0.03).

**A5d boxes (all-sat / two-sat; sieve p only, 8 boxes × 7 terms, so read as descriptive)**
- **SE Indian Ridge–Macquarie:** KE and EKE rise both upstream and downstream, in both products (p ≤ 0.02).
- **Pacific–Antarctic Ridge:**
  - Upstream along-stream speed +3.5 % decade⁻¹ (p = 0.002, both products).
  - Downstream KE +6.9 / +5.9 % decade⁻¹ (p ≤ 0.004).
  - Downstream EKE +9.1 % (p = 0.006) in all-sat, but +4.9 % (p = 0.15) in two-sat.
- **Drake:**
  - Upstream (100–80°W): KE and EKE +11–15 % decade⁻¹ in both products.
  - Downstream (Scotia, 55–35°W): EKE +7.6 % (p = 0.001) in all-sat vs +3.1 % (p = 0.09) in two-sat.
- **Kerguelen:**
  - Upstream: nothing.
  - Downstream: EKE +6.4 % (p = 0.02) in all-sat only. Two-sat +1.0 % (p = 0.71).
- **Pattern:** mean-flow increases are product-robust. Downstream EKE increases are often present in all-sat only.

## A6: Drake Passage (D6)

**Altimetric transport proxy (corrected; `a6c_transport_fixed.py`)**
- **Correction.** `a6_drake.py` excluded section points month by month where ice fraction > 0.10. That series (T_v2) correlates with the fraction of valid points (r = 0.89; e.g. 2011-08, half the points valid, is a −6σ outlier). Its earlier trend (+1.9 % decade⁻¹, p_cal = 0.03) is therefore an artefact of varying section length and is withdrawn.
- **Fixed geometry.** Two fixed point sets:
  - T_perm: only points that are ice-free (≤ 0.10) in every month. This keeps 53 of 121 points on the v2 line, north of 58.5°S, and 49–59 points on the variants.
  - T_all: all points, including mapped under-ice values.

| Series | 2001–2025 | 2005–2019 |
|---|---|---|
| T_perm, v2 line | +2.3 % decade⁻¹, p_cal = 0.31 | +5.2 %, p = 0.41 |
| T_perm, 7 line variants | −0.2 to +4.8 %, p 0.13–0.96 | p 0.12–0.97 |
| T_all, v2 line | +0.7 %, p = 0.16 | +1.7 %, p = 0.03 |
| T_all, 7 line variants | +0.4 to +1.6 %, p 0.03–0.41 | +1.8 to +2.7 %, p 0.01–0.03 |

- **2001–2025:** no robust trend in the surface proxy.
- **2005–2019:** only T_all, which includes under-ice values, gives a small positive trend (p ≤ 0.03). It has the opposite sign to the in-situ estimate below, which is itself non-significant.
- **ΔADT across the passage:** +0.0009 m yr⁻¹ (p = 0.09).

**In situ**
- 0–760 m total ADCP transport, annual means 2005–2019: −0.21 Sv yr⁻¹, 95 % CI [−0.66, +0.24], p = 0.33.

**Agreement between the two (independently verified)**
- Same-month correlation: r = 0.05 (n = 244 transects, of which 111 share a month with another transect).
- Annual means: r = −0.09 (15 years).
- Fixed geometry: T_perm r = 0.01 (transects), −0.04 (annual); T_all r = 0.05, 0.03; Gould-line approximation r ≤ 0.06.
- **The surface proxy does not track the in-situ transport at these time scales.** Any transport statement must rest on each record separately.

**Drake box (80–52°W, 63–54°S) energy**

| Mask and period | EKE trend | p_cal | KE trend | p_cal |
|---|---|---|---|---|
| a5 ice mask (≈ v2 cube), 2001–2025 | +2.18 | 0.04 | +0.78 | 0.25 |
| a5 ice mask, 2001–2024 | +2.49 | 0.01 | +0.95 | 0.18 |
| No mask, 2001–2025 | +1.52 | 0.09 | +0.53 | 0.28 |
| No mask, 2001–2024 | +1.73 | 0.02 | — | — |
| No mask, 2005–2019 | +3.57 | 0.03 | n.s. | — |

- The ice mask removes 37 % of the box pixels.
- An independent check with an AR(1) effective-n test gives smaller p-values throughout. KE becomes p ≈ 0.03–0.04 with the mask.
- The EKE/KE contrast is therefore **test- and mask-dependent**. It should be described as borderline, not as "the clearest signature".

**Wind (`eastward_stress`, N m⁻²)**

| Domain | 2001–2025 | 2008–2025 |
|---|---|---|
| Circumpolar 45–60°S | +2.5 % decade⁻¹, p = 0.07 | p = 0.30 |
| Pacific 45–60°S | p = 0.12 | p = 0.74 |
| Drake box | +8.7 % decade⁻¹, p = 0.03 | p = 0.31 |

- **Lags of 0–3 yr** (detrended annual means): no association with Drake EKE or the transport proxy (all p > 0.05; largest |r| = 0.40 at lag 0, circumpolar stress vs transport proxy, p = 0.06). These lags were computed with T_v2 and have not been recomputed with T_perm.
- **Consequence:** the circumpolar stress increase is not robustly detected over this record, so the flat in-situ transport cannot be read as insensitivity to wind.

## B1 pilot: Pacific–Antarctic Ridge (Gate G)

- **Pathway detection:** the ADT contour at the climatological core level (relative ADT −0.18 m, mean speed 0.28 m s⁻¹) spans the Lower-Ridge section in **100 % of months**, at all three levels. Gate part 1 passes.
- **Curvature (∫|κ|ds, 100 km smoothing):**
  - Trend +0.057 rad yr⁻¹ [−0.12, +0.24], p = 0.42.
  - Interannual SD 2.63; minimum detectable 25-yr change 4.15.
  - **Gate part 2 fails.** It also fails for every curvature, length, sinuosity and amplitude metric at the core level and all smoothing scales.
  - An independent check (AR(1)-adjusted OLS) gives MDC 3.87, so it fails as well.
- **Isolated hit at core − 0.10 m:** that level shows a decrease in length and curvature (p 0.02–0.05). It is not present at the core or core + 0.10 m levels, and is one of 45 geometry tests, so it should not be interpreted.
- **Chain variables, all-sat:**

  | Variable | Trend | p_cal |
  |---|---|---|
  | Upstream jet-core speed | +0.0015 m s⁻¹ yr⁻¹ | 0.03 |
  | Meander-section core speed | +0.0026 | 0.01 |
  | Downstream EKE | +1.12 cm² s⁻² yr⁻¹ | 0.01 |

- **Lagged associations** (0–24 months; n_eff ≈ 7–12): upstream speed → curvature, curvature → downstream EKE, and upstream speed → downstream EKE are all non-significant (|r| ≤ 0.47, p ≥ 0.24).
- **Original decision:** do not extend B1. **Superseded.**

### Why the gate was revised

Two flaws in the pilot design (mine):
1. **The power criterion could not be met.** With 25 annual values, even white noise gives MDC = 1.39 × SD. Upstream speed and downstream EKE "failed" it despite significant trends.
2. **Wrong time scale.** TNG14's curvature–EKE link is a month-to-month covariation (8-yr OFES run, 2-month smoothing, r = 0.34–0.46). The pilot's 12-month low-pass removed it.

What stands from the pilot: the PAR curvature **trend** is only bounded (25-yr change between about −23 % and +45 % of the mean).

### B1 revised: covariability test at three meanders

**Specification**, fixed before any cross-correlation was computed (`b1_covar.py`):
- **Series:** monthly; deseasonalised; Theil–Sen detrended; 2-month running mean.
- **Curvature:** core-level contour length at 150 km smoothing (primary, TNG14's measure); ∫|κ|ds (secondary).
- **Links** at lags of 0–3 months:
  - L1: upstream speed → curvature;
  - L2: curvature → downstream EKE (primary);
  - L3: curvature → meander-section EKE (secondary, co-located).
- **Significance:** n_eff t-test; phase-randomised surrogates, including the maximum over lags.
- **Pooling:** Stouffer z across meanders.

**Meanders.** Sections were fixed from time-mean fields only (`b1_regions.py`, `b1_meander_run.py`):
- **PAR:** Liu et al. (2026) sections.
- **SEIR150:** SE Indian Ridge crossing, TNG14's 150–155°E meander. Upstream 125–138°E, meander 140–150°E, downstream 150–158°E.
- **KERG:** Fawn Trough branch. Upstream 70–77°E, meander 76–84°E, downstream 84–94°E.

In all three, the EKE boxes exclude pixels with ice fraction > 0.10 in any month, and the pathway is tracked in 100 % of months.

**Results**

| Link | PAR | SEIR150 | KERG | Pooled |
|---|---|---|---|---|
| L1 speed → length | r ≤ 0.15, n.s. | n.s. | n.s. | max-lag p = 0.09 |
| **L2 length → downstream EKE** | r 0.07–0.11, p = 0.45 | **r 0.20–0.28 (lags 0–2), p = 0.005** | r up to 0.21 (lag 3), p = 0.054 | **z 2.1–2.6 at lags 1–3, same sign in all three; max-lag p = 0.007** |
| L2 ∫\|κ\| → downstream EKE | n.s. | r 0.21–0.27, significant | negative at lags 0–1 | p = 0.20 |
| L3 length / ∫\|κ\| → meander EKE | — | — | — | p = 0.001–0.002 (co-located; not independent evidence) |

**Lag asymmetry (exploratory, outside the pre-specified window).** At SEIR150 and KERG the correlation sits at positive lags (length leads); with EKE leading, r ≈ 0 or negative. At KERG the peak lies at lags 4–6 months.

**Independent check (separate code).**
- Every r, n_eff, z and p is reproduced exactly; surrogate p-values agree within Monte-Carlo noise.
- Robust to no 2-month smoothing (pooled surrogate p = 0.010) and to OLS detrending (p = 0.012).
- **Weaker** with Spearman (p = 0.039; lag 1 n.s.).
- **SEIR150 is carried by a few episodes:** trimming the 5 % most extreme months lowers lag-1 r from 0.28 to 0.15. The 10 largest cross-products give 54 % of the covariance (2009-12/2010-01, 2010-11, 2012-06, 2024-03/04).

**Interpretation for the paper**
- Weak (r ≈ 0.1–0.3), episodic evidence, consistent in sign across three meanders, that pathway lengthening precedes downstream EKE by 1–3 months. This is compatible with TNG14's "second effect", but weaker than their model correlations, and not supported by the ∫|κ| metric.
- **No surface evidence for the first link** (upstream acceleration → curvature).
- Lead–lag ordering does not establish causation.

Figure: `outputs/figures/figS_b1_covariability.png`.


### B1 revisited: order of the meander–EKE sequence (2026-10-13)

Checks on the pre-specified L2 result (same prepared series; surrogate p over lags 0–3):

| Meander | length → meander EKE | meander EKE → length (lags 1–3) | meander EKE → downstream EKE | length → downstream EKE given meander EKE (partial r, lags 1–3) | amplitude → downstream EKE | circumpolar stress → length |
|---|---|---|---|---|---|---|
| PAR | r 0.19, p 0.32 | r 0.23, p 0.10 | r 0.37, p 0.03 | 0.06, 0.06, 0.02 | p 0.83 | r 0.23, p 0.03 |
| SEIR150 | r 0.23, p 0.014 | r 0.12, p 0.38 | r 0.62, p 0.001 | 0.20, 0.22, 0.18 | p 0.90 | p 0.97 |
| KERG | r 0.33, p 0.002 | r 0.15, p 0.41 | r 0.58, p 0.001 | −0.18, −0.05, 0.09 | p 0.13 | p 0.10 |

- At SEIR150 and KERG the order is path lengthening → EKE in the meander section (1–2 months) → EKE downstream (1–3 months more). The reverse order shows nothing, so the lengthening is not the contour wrapping around eddies already present.
- The pooled L2 result is mostly mediated by the meander-section EKE; a direct length → downstream relation survives only at SEIR150.
- The crest-to-trough amplitude carries no signal: the lengthening is sub-meander-scale distortion (growth of lobes before shedding), not flexing of the meander envelope.
- No forcing → geometry link: upstream speed (L1) and circumpolar stress do not precede the lengthening, except weakly at PAR for the stress.
- Fig. 4 (v2.1) now shows the meander-section EKE composite (dashed) next to the downstream one; Discussion paragraph rewritten accordingly.

### B1 paths over time (descriptive; `b1_maps.py`, `b1_paths_stats.py`)

**Path definition.** ADT minus its regional mean each month, which removes uniform sea-level rise (Gille 2014), contoured at the climatological jet-core level.

**Figures** (`figS_b1_paths_{PAR,SEIR150,KERG}.png`):
- annual-mean paths coloured by year;
- monthly path occurrence with the 2001–08 / 2009–16 / 2017–25 mean paths;
- Hovmöller of path-latitude anomaly;
- within-year SD of path latitude, an OSR8-style variability diagnostic (Sánchez-Román et al. 2024, after Andres 2016).

**Section statistics** (25 annual values, Theil–Sen with sieve p, BH-FDR over 18 tests):
- No section of any meander shows an FDR-significant change in either mean path latitude or month-to-month path variability.
- Path-latitude trends are between −0.04 and +0.05 ° decade⁻¹. The 95 % CIs lie within about ±0.15 ° decade⁻¹ (±17 km decade⁻¹).
- The epoch-mean paths overlap to within about 0.1° in the meander sections, i.e. they are topographically locked.
- This is consistent with Liu et al. (2026): no significant amplitude trend at PAR.

**Caveat (Kerguelen).** Near 77–80°E the path runs north–south, so path latitude per longitude is ill-defined there.

### A7: is the energy increase concentrated downstream of the meanders? (`a7_localisation.py`)

**Specification**, written before computing the statistics:
- **Regions:** WAKE is the union of the four A5d downstream boxes (Kerguelen, SEIR–Macquarie, PAR, Drake/Scotia); REST is the remainder of the ACC band. Ice-affected pixels are excluded.
  - WAKE has 23 % of the band area and 42–43 % of its mean EKE.
- **Statistic:** the trend of n_WAKE − n_REST, where each series is normalised by its own mean. This is the difference in % decade⁻¹. Null hypothesis: the same proportional change everywhere. Tested with the calibrated test.
- **Primary test:** EKE, 2001–2025, each product.
- **Secondary tests:** KE, SPD, SUB, LF, and per-feature downstream vs upstream and downstream vs REST; BH-FDR at 0.10 within each product and period.
- **Sensitivity:** 2001–2017.
- **Caveat:** the A5d box trends had already been seen descriptively.

**Results, 2001–2025**

| Variable | All-sat WAKE / REST (% decade⁻¹) | Difference [95 % CI], p | Two-sat WAKE / REST | Difference, p |
|---|---|---|---|---|
| **EKE (primary)** | +8.4 / +8.0 | +0.1 [−3.0, +3.2], p = 0.94 | +3.9 / +3.9 | −0.3, p = 0.80 |
| KE | +4.5 / +5.1 | −0.7, p = 0.49 | +3.3 / +3.6 | −0.5, p = 0.56 |
| Along-stream speed | +1.4 / +2.2 | −0.9, p = 0.21 | +1.3 / +2.0 | −0.7, p = 0.29 |
| SUB / LF | — | n.s. | — | n.s. |

- The wakes hold 45 % (all-sat) and 43 % (two-sat) of the ACC-integrated EKE trend. This is in proportion to their share of mean EKE (43 % and 42 %).
- **Per feature:**
  - PAR and Kerguelen: downstream EKE rises more than upstream, but not significantly (PAR +7.4 points, p = 0.12 / 0.20).
  - Drake: the reverse. Upstream (100–80°W) EKE and KE rise more than in Scotia (EKE −10 points, p = 0.02–0.03; KE −9 to −10, p = 0.01, both products). These do not pass FDR: with 20 tests and the calibrated-p floor of about 0.01, a single p ≈ 0.01 cannot pass.

**Sensitivity, 2001–2017**
- EKE: WAKE +6.7 vs REST +2.3 % decade⁻¹ in all-sat (difference +4.1, p = 0.13); two-sat +5.3 vs +0.5 (+4.3, p = 0.06).
- Low-frequency EKE: concentrated in WAKE (+5.0 and +5.1 points, p = 0.02 and 0.04; FDR-significant in both products).
- Along-stream speed: increased outside the wakes but not in them (−3.2 and −3.1 points, p = 0.01; FDR-significant in both).
- **Reading:** over 2001–2017, mean-flow speed increased away from the wakes, while the increase in low-frequency EKE was concentrated in them. After 2017 the REST region catches up, and over the full record the increase is proportionally uniform. This is a sensitivity period (pre-specified for product reasons), not the primary test.

**Interpretation**
- Over 2001–2025, the ACC EKE increase is proportionally uniform. It is not preferentially located downstream of the main topographic features, in either product.
- The surface record therefore does not show the localisation of the eddy response at standing meanders expected under eddy saturation over the full period. There is a suggestion of it in 2001–2017.
- Together with B1: meander geometry is stable; path lengthening precedes downstream EKE on monthly scales; the long-term EKE increase is basin-wide in relative terms.

### Reproducing Liu et al. (2026), PAR

- **Scripts.** `extract_liu.py` (to be run on the Mac: daily 008_047 fields, 1993–2025) computes daily |∇ADT| over 150–60°W, 70–30°S, applies the 20 % threshold of the daily domain maximum, and writes monthly detection counts, summed gradient and summed speed. It does this on the native 1/8° grid and on a 2 × 2 block-averaged 1/4° grid (Liu used 1/4°).
- **Reproduction.** `liu_reproduce.py` builds the 4-month occurrence, dominant axis, zero-to-zero width, axis speed, Table 1 peaks/troughs and the Tr2–Pk3 amplitude. It computes trends with Liu's own tests (ported from their MATLAB), with Hamed–Rao, and with our calibrated test.
- **Liu et al.'s MATLAB code** (Zenodo 17810036; the copy in `Downloads/19211657` has 8 scripts) contains no detection or axis script, so those steps follow the paper's Methods.
- **Their "modified Mann–Kendall" (06-…m).** It computes S from all n values but its variance from n_eff = n/τ, and the ACF is taken from the raw, not detrended, series. On trend-free AR(1) series (n = 372) it rejects at the nominal 5 % level in:

  | Series | False-positive rate (plain MK in brackets) |
  |---|---|
  | φ = 0.5 | 63 % (23 %) |
  | φ = 0.5, 4-month running mean | 89 % (39 %) |
  | φ = 0.8, 4-month running mean | 98 % (60 %) |

  Check: `liu_mmk_check.py`. Their pMKmod < 0.001 values therefore do not establish significance. This statement concerns the test only, not their trend magnitudes.

**Monthly preview, 2001–2025** (`liu_monthly_preview.py`; not a reproduction)

Method: 4-month windows of monthly-mean fields. The axis is the latitude of the maximum summed |∇ADT|; speed is |U| at the axis; the "band width" is the extent above 20 % of the 4-month domain maximum. Values below are deseasonalised Theil–Sen per decade, calibrated p in brackets.

| Metric | All-sat | Two-sat | Liu et al. (1993–2023) |
|---|---|---|---|
| Speed, whole | +0.019 m s⁻¹ (0.01) | +0.017 (0.01) | +0.01 |
| Speed, Lower-Ridge | +0.026 (0.01) | +0.024 (0.01) | +0.02 |
| Speed, Upper-Ridge / Downstream / Flat | all p ≤ 0.05 | all p ≤ 0.05 | +0.01 each |
| Position, Lower-Ridge | +0.06° (0.24) | +0.04° (0.44) | +0.09° (significant) |
| Position, Flat | −0.18° (0.12) | −0.16° (0.11) | −0.12° (significant) |
| Tr2–Pk3 amplitude | −0.08° (0.25) | −0.09° (0.06) | +0.03° (n.s.) |
| Peaks / troughs | none significant | none significant | — |
| Threshold band width, whole | +4.4 km (0.04) | +4.2 km (0.27) | — |

**Reading the table**
- **Speed:** acceleration is robust in both products and in every section, and its magnitude is consistent with Liu.
- **Position:** no position shift survives the calibrated test.
- **Band width:** product-dependent, and not Liu's width (which measures day-to-day wandering). A threshold band also grows with gradient strength (cf. D3).
- **Liu's own tests on these same series:** their modified MK gives p < 0.001 for nearly every series, including Pk1, whose Sen slope is exactly 0. Plain MK gives p < 0.05 for the Lower-Ridge position.

**Daily reproduction, 1993–2025** (`extract_liu.py` daily 008_047 fields → `liu_reproduce.py`; 1/4° block-averaged grid; trends were computed on the 1/4° grid only)

Method: 4-month occurrence of |∇ADT| > 20 % of the daily maximum, as in Liu et al. Our trends are deseasonalised Theil–Sen with the calibrated p (p_cal floor ≈ 0.01). Values are per decade, period 1993-01 to 2023-01 (n = 360), unless stated.

| Metric | Liu Table 2 | This reproduction | p (Liu's MK / Liu's modified MK / Hamed–Rao / calibrated) |
|---|---|---|---|
| Position, whole | +0.03° | +0.036° | 0.003 / <0.001 / 0.12 / **0.30** |
| Position, Upper-Ridge | +0.04° | +0.035° | 0.046 / <0.001 / 0.17 / 0.42 |
| Position, Lower-Ridge | +0.09° | +0.11° | <0.001 / <0.001 / 0.003 / **0.02** |
| Position, Downstream | −0.03° | −0.025° | 0.25 / <0.001 / 0.65 / 0.75 |
| Position, Flat | −0.12° | −0.04° | 0.20 / <0.001 / 0.31 / 0.60 |
| Speed, whole | +0.01 m s⁻¹ | +0.016 | all ≤ 0.01 |
| Speed, Lower-Ridge | +0.02 | +0.028 | all ≤ 0.01 |
| Speed, Upper-Ridge | +0.01 | +0.008 | 0.001 / <0.001 / 0.05 / 0.10 |
| Speed, Downstream / Flat | +0.01 / +0.01 | +0.011 / +0.011 | p_cal 0.02 / 0.01 |
| Width (zero-to-zero), whole | 1.44 km (mean 100 km) | +26 km (mean 652 km) | p_cal 0.02 |
| Half-max width, whole (our alternative) | — | +5.2 km (mean 190 km) | p_cal 0.02 |
| Half-max width, Lower-Ridge | (3.52 km, Liu's width) | +10.4 km (mean 205 km) | p_cal 0.01 |
| Amplitude Tr2–Pk3 | +0.03° (n.s.) | −0.05° (mean 0.51°) | p_cal 0.48 |
| Peaks / troughs | Pk2, Tr2 significant | none significant | sieve p ≥ 0.29 (Liu's MK: Pk2 p < 0.001, Pk3 p < 0.001) |

**Extended to 1993–2025.**
- Speed: still robust in every section except Upper-Ridge (p_cal 0.08).
- Lower-Ridge position: +0.08° decade⁻¹, p_cal = 0.04.
- Whole-meander position: +0.008°, p = 0.75.

**Reading**
- **Reproduced:**
  - **Speed:** acceleration in every section, with magnitudes matching Liu's to their two-decimal rounding. Robust to the calibrated test, and also seen in the 2001–2025 monthly preview, in both products.
  - **Position:** the section position trends reproduce in sign and approximately in magnitude (Flat Region smaller).
- **Survives a calibrated test, position:** only the Lower-Ridge northward shift (+0.08 to +0.11° decade⁻¹). It depends on including 1993–2000: in 2001–2025 it is +0.06°, p = 0.24. Liu's whole-meander shift and their Flat-Region shift are not significant under autocorrelation-aware tests.
- **Width is not reproducible as described.**
  - The zero-to-zero width of the 4-month occurrence is about 650 km, not 100 km. Liu's Fig. 2c description (zeros near 60°S and 50°S) implies a width of the same order as ours, so the paper is internally inconsistent on this point.
  - Their width trend therefore cannot be compared.
  - Our half-maximum width (about 190 km) widens by 5 km decade⁻¹ (Lower-Ridge 10 km decade⁻¹), p_cal ≤ 0.02. Because the threshold is set relative to the daily maximum gradient, this widening may still partly reflect changes in jet strength (Liu's own caveat); it has not been tested here.
- **Peak and trough latitudes** are quantised on the grid (Sen slope 0). Liu's MK gives small p-values for them that the sieve test does not support.
- **Implication for our manuscript.** The robust, reproducible PAR signal is acceleration of the jet core. It agrees with our ACC-wide along-stream speed increase (A5). Meander geometry (amplitude, paths) shows no detectable change; this agrees with B1 paths.


## E: tracked eddies, AVISO META4.0 DT all-sat atlas (2001-01 to 2023-08; `eddies_subset_aviso.py`, `eddies_analysis.py`, `fig_eddies.py`)

Atlas subset south of 35°S: 16.8 M observations, 524 k tracks (track id 0 = untracked, excluded from track statistics). ACC band = mean ADT −0.6..+0.2 m, ice-free. Only eddies larger than the mapping resolution are tracked (median speed radius 36–38 km, 5th percentile 24 km).

**E2 — what the EKE trend is made of** (ACC band, revised test, % decade⁻¹; gridded EKE/KE over the same months)

| Quantity | Trend | p |
|---|---|---|
| Gridded EKE (all-sat) | +8.6 | 0.027 |
| Gridded KE (all-sat) | +5.6 | 0.003 |
| Σ amplitude² of tracked eddies per day (eddy-energy proxy) | +10.9 | 0.007 |
| Mean eddy amplitude | +4.6 | 0.007 |
| Mean rotational speed | +3.8 | 0.020 |
| Eddies present per day | +0.8 | 0.37 |
| Tracks born per month | −0.4 | 0.63 |
| Mean speed radius | −0.1 | 0.60 |

- The number and size of tracked eddies do not change; their amplitude and rotational speed increase. The EKE increase in the gridded product is therefore an increase in eddy *intensity*, not in eddy *number*. The Σ amplitude² proxy follows the gridded EKE closely (Fig. E b), including the rise after 2018.
- Caveat: the atlas is built from the same all-satellite maps, so the amplitude trend shares the constellation sensitivity (S1); the eddy count does not depend on it in the same way.

**E3 — wakes vs rest of the band** (difference of normalised trends, points decade⁻¹): eddies present +0.3 (p 0.71), births −1.9 (p 0.21), amplitude −0.03 (p 0.97), Σ amplitude² −1.9 (p 0.55). Same conclusion as A7: no concentration of the change downstream of the features.

**E1 — generation sites.** Births per year per 1° cell: ACC-band mean 1.7; wake boxes 1.8 (PAR) to 2.6 (Kerguelen); upstream / meander / downstream sections: PAR 1.5 / 2.0 / 1.6, SEIR150 1.5 / 1.1 / 2.1, Kerguelen 4.0 / 3.6 / 2.7. Births are enhanced along the whole energetic band and particularly at Kerguelen and Drake; the meander sections are not singular generation sites at 1° resolution.

**E4 — shedding vs path distortion.** Monthly births in the meander section + downstream box against the path length (b1 preparation; lags 0–3, surrogate p): PAR p = 0.06 (r ≈ −0.2, i.e. fewer births when the path is longer), SEIR150 p = 0.16 (r ≈ +0.15), KERG p = 0.74. No consistent link. The month-to-month lengthening → EKE sequence of B1 is therefore not reflected in the number of tracked eddies born; it is more likely an intensification of the existing eddies or of the meander lobes themselves than an increase in detachments at the scales the atlas resolves.

Figures: `python fig_eddies.py main` → `outputs/figures_v21/Figure4` (v2.1 main text, 2 panels: ACC-band series and forest of trends); `python fig_eddies.py` → `FigureS_eddies` (SI, with the birth-density map). `FigureE_eddies` is the earlier 3-panel version.

**In v2.1 (2026-10-14).** Atlas added as Product 5 in Table 1 (META4.0 DT allsat, https://doi.org/10.24400/527896/a01-2026.001); new Sect. 2.5 "Tracked eddies" (the old 2.4 "Time series and linear trends" becomes 2.6); new paragraph at the end of Sect. 3.1 with the E2 numbers and Figure 4; the meander figure becomes Figure 5; one clause in the abstract, one sentence in the Discussion (where the joint intensification of mean flow and eddies is argued, with the constellation caveat); Pegliasco et al. (2022) added to the references. The SI figure is cited as "Supporting Information S4" — S3 (half-power sensitivity, longitude-resolved trends) and S4 are not written yet.

## B2: box-free versions of the meander → EKE test (`b2_pathfollowing.py`, `b2_tests.py`, `b2_figure.py`; 2026-10-14)

Pre-specified before results: 36 overlapping 30° windows every 10°; per window the jet-core path (relative-ADT contour at the climatological core level, as b1) and its 150 km-smoothed length; "path band" = grid points within 150 km of the monthly path; EKE band mean; eddy births (META4.0 tracks) within the band; same preparation and lags 0–3 as b1_covar; 1000 phase surrogates. After the first pass two corrections were made because the design was contaminated by construction: downstream EKE taken from the window centred at c+30 restricted to (c+15, c+35] (the c+10 window overlaps by two thirds), and births expressed per 1000 km of path. Outputs: `outputs/b2/b2_series.csv`, `b2_eke_bins.npy`, `b2_window_tests_v2.csv`, `b2_pooled_v2.csv`, `b2_alt1_lonbins.csv`, `b2_quicklook.png`.

- Coverage: the path is tracked in ≥ 98 % of months in 29 of 36 windows; it fails at 30–50°E and 90–100°E (multiple jets, Crozet/Kerguelen region; the core-level search hits its ±0.6 m bound there) and at −170°E.
- Alt. 2, same-window link (length → EKE in the band of the same window): r = 0.25–0.49, p < 0.05 in 25 of 29 windows, strongest at lag 0 and decaying with lag everywhere. This is a simultaneous covariation: the mesoscale eddies are the distortions of the path. No precedence of the path over the EKE.
- Alt. 2, downstream link (length → EKE 15–35° east, non-overlapping): r = −0.1 to +0.3, **no window significant**; pooled sign positive in 6 of 6 testable non-overlapping windows (Stouffer z ≈ 3.5 at lags 0–1, but with n_eff = 290 and without accounting for autocorrelation, so overstated).
- Alt. 1, along-stream profiles for the three meander windows (KERG 80, SEIR 150, PAR −120): the response is concentrated within the window whose path is used (r up to 0.25 at lag 0–2); east of the window r drops to ≈ 0 and fluctuates ±0.1–0.2 with no consistent downstream maximum.
- Alt. 3, births per 1000 km of path vs path length: incoherent (significant positive in some Atlantic/Indian windows, significant negative at −60/−50°E and 170°E); pooled nothing.
**Revision of the same day (circumpolar core level).** The per-window core level picked non-ACC jets in several windows (the Brazil–Malvinas confluence at 60–40°W, the Agulhas Return Current at 10–50°E) and failed in 7 windows. Second pass: one circumpolar relative-ADT level (ADT minus its monthly mean over the ACC band), chosen from the climatology as the level whose contour carries the highest mean speed: **+0.18 m** (the SAF; the speed-vs-level curve is broad, 0.16–0.22 m within 0.5 %, so the choice is weakly determined and a southern-front sensitivity run is pending). First pass kept in `outputs/b2_v1/`. With this path: coverage 99–100 % in all 36 windows; same-window link r = 0.05–0.62, p < 0.05 in 22 of 36 windows, lag 0 strongest everywhere; downstream link (15–35° east, non-overlapping) p < 0.05 in 1 of 36 (−20°E, r = 0.27), 8 of 12 non-overlapping windows positive, pooled z = 2.0 at lag 0 (overstated, autocorrelation ignored); births per 1000 km p < 0.05 in 2 of 36 with opposite signs. Figure: `v21_fig5_b2.py` → `outputs/figures_v21/Figure5` (map of mean EKE with the mean path + r by window); the box version is kept as `Figure5_boxes_old`.

**Sensitivity to the front tracked (`B2_LEVEL=-0.2`, Polar-Front side; outputs in `outputs/b2_pf/`).** Coverage 100 % in all windows. Same-window link: p < 0.05 in 14 of 36 (median r 0.16; SAF level 20 of 36, median 0.28), lag 0 strongest. Downstream link: p < 0.05 in 3 of 36 with mixed signs (+0.32 at 60°W, +0.34 at 20°E, −0.30 at 70°W; median r 0.12, as for the SAF). Births per path length: 1 of 36. The conclusion does not depend on the level: local, simultaneous covariation between path distortion and EKE; no detectable downstream response.

- Reading: the b1 sequence (length → meander-section EKE → downstream EKE, with lags of 1–3 months) was partly a product of the boxes. Box-free, the robust result is a local, simultaneous path-length–EKE covariation around the ACC, and a weak, individually non-significant downstream association. The Thompson & Naveira Garabato shape → downstream EKE step is not demonstrated at monthly time scales with these data.

## v2.3 additions (2026-10-15)

**B3: along-stream structure of the trends** (`b3_alongstream.py` series, `b3_trends_fast.py` trends, `b3_relations.py`; spec fixed before results; Fig. 4 of v2.3 via `v23_fig4.py`). Path = monthly +0.18 m relative-ADT contour (SAF side), smoothed 150 km; 72 bins of 5°; per bin and month: latitude, sinuosity (length / zonal length), speed along the path, EKE within 150 km. Trend test per bin = trend_bcz statistic (both noise models, larger p) with 300 null series per bin (single bootstrap for CI; the double bootstrap of trend_bcz was ~2 min per series); FDR 10 % across bins. Malvinas bins (60–50°W) excluded (contour follows the Malvinas Current; sinuosity coverage 0.69).
- EKE: positive 58/70; 12 p<0.05; FDR 0; mean 7.0 %/dec over obstacles (depth_min < 3000 m), 7.3 % 1–3 bins downstream, 6.9 % elsewhere (no localisation).
- Jet speed: positive 51/70; median +3.9 %/dec; 4 p<0.05; FDR 0.
- Sinuosity: positive 36/70; 3 p<0.05; FDR 0. Latitude: positive 35/70; 8 p<0.05; FDR 3 (isolated: ~147°W north, ~32°W and ~157°E south).
- Spatial r across bins (2000 circular shifts): EKE trend vs speed trend +0.40 (p 0.008); vs sinuosity trend +0.29 (p 0.0005); vs latitude trend −0.02 (p 0.95); vs log mean EKE +0.25 (p 0.20).
- Reading: faster jet nearly everywhere; meanders neither flex nor move; EKE increase not localised; EKE grows more where the jet accelerates and the path gets more sinuous (decadal link; cannot separate cause and effect because EKE and monthly-path sinuosity both contain the eddies).

**A2g: half-power trends over the regions of earlier studies** (`a2g_literature_regions.py`, 2001–2025, km/dec, + north): circumpolar −0.5 (−4.2, +3.1); 90–170°E −1.6 (−12.9, +9.6) [Chambers 2018: −29, 1993–2016]; 160–90°W +5.6 (−4.9, +16.2) [Chambers: +46]; 75–150°E −1.5 (−9.4, +6.5) [Kim & Orsi 2014: southward drift 1992–2011, no number used]; 150–80°W PAR +5.9 (−5.3, +17.1) [Liu et al. 2026: +0.03°/dec ≈ +3 km/dec]. SSH contours: Sokolov & Rintoul 2009b 60 km over 1992–2007 (≈ −40 km/dec); Gille 2014 60 km over two decades (≈ −30 km/dec); Gille transport latitude and Freeman et al. 2016 zonal-mean PF (2002–2014): no trend. Sources: Chambers (2018) Sect. 4; Gille (2014) Sect. 3; Chapman (2017) intro for SR09b; Freeman et al. (2016) Sect. 3.3/4.1.

**Fig. 3 (v2.3)**: all trends in % of mean per decade: Drake altimetric proxy +2.3 (−3.7, +8.2) p 0.40; in situ −2.8 (−16.4, +10.8) p 0.47 (mean 75.7 Sv); circumpolar proxy +1.6; Pacific 2.6; Atlantic 1.8; Indian 0.4 n.s.; Drake EKE +8.4 (0.2, 16.7) p 0.050; Drake KE +2.4 n.s.; ACC EKE 8.0; KE 4.9; wind stress Drake +8.7 (1.4, 15.9) p 0.03; 45–60°S +2.5 (−0.3, 5.3) p 0.08.

**Two-sat** moved to SI S1 (`figS1_twosat.py` → `figures_v23/FigureS1_products`).

## Fig. 1 maps (A1)

**Method.** Stored v2 all-sat cube, coarsened to 1°, 2001–2025. Each cell gets a sieve null (200 surrogates). The p-value is the Gaussian tail of that null, because the empirical p has a floor of 1/201 that makes BH over about 8,000 cells unreachable. BH-FDR is applied at α = 0.10.

| Field | Cells with p < 0.05 | Cells surviving FDR |
|---|---|---|
| KE | 679 of 8,184 (8.3 %) | 58 (0.7 %) |
| EKE | 905 (11.1 %) | 133 (1.6 %) |

Nominal p < 0.05 would give about 5 % of cells by chance, and the sieve is somewhat liberal, so few local trends are field-significant. The FDR-significant EKE cells cluster:
- south-east of Kerguelen (80–90°E, 52–56°S);
- south of Tasmania and around Macquarie (150–160°E, 53–58°S);
- along the Pacific ACC (130–100°W);
- in parts of the Argentine Basin.

Figure: `outputs/figures/fig1_trend_maps_fdr.png` (3,000 m isobath and ACC-band contours). The v2 maps' "significant" stippling should be replaced by this.

## Statistics check (`stats_size_check.py`, `stats_contrast_tests.py`)

**1. False-positive rate of the trend tests.** Trend-free AR(1) series with a seasonal cycle, n = 300 months, 150 replicates each, nominal α = 0.05. The calibrated test ran at reduced settings (n_outer 40, n_boot 150).

| Process | Sieve p | Calibrated p | Hamed–Rao |
|---|---|---|---|
| φ = 0.8 | 14.7 % | 10.7 % (95 % CI ≈ 6–17 %) | 30 % |
| φ = 0.95 | 18.7 % | 6.0 % (≈ 3–11 %) | 33 % |

- The calibrated test is still somewhat liberal at φ = 0.8: p_cal < 0.024 corresponds to a true 5 % level.
- **Headline trends are unaffected.** ACC KE, EKE, SUB, LF and SPD have z = 3.5–7.2 in both products.
- **At risk are results with p_cal between 0.02 and 0.05:**
  - Drake-box EKE (0.04);
  - Atlantic EKE, all-sat (0.02);
  - T_all, 2005–2019 (0.03);
  - A7 Drake downstream–upstream (0.02–0.03);
  - A7 2001–2017 low-frequency EKE (0.02 / 0.04);
  - Liu reproduction, Lower-Ridge position (0.02 / 0.04).
- **Other notes:**
  - The plain sieve p (used for the A5d boxes, the maps' null and the B1 path statistics) is liberal by about 3×.
  - Hamed–Rao is liberal by about 6×.
  - The p_cal floor (≈ 0.0099 with n_outer = 100) prevents BH-FDR from passing a single strong result among ≥ 11 tests.

**2. Contrast claims tested directly.** Trend of the difference of mean-normalised series, in % decade⁻¹:

| Contrast | Difference [95 % CI], p_cal |
|---|---|
| ACC, EKE − KE, all-sat | +3.0 [1.3, 4.7], p ≈ 0.01 |
| ACC, EKE − KE, two-sat | +0.2 [−1.2, 1.5], p = 0.74 |
| Drake box, EKE − KE, ice mask (all-sat) | +6.0 [2.0, 10.1], p ≈ 0.01 |
| Drake box, EKE − KE, no mask (all-sat) | +5.0 [2.0, 8.0], p ≈ 0.01 |

- The Drake contrast has not been computed in two-sat.
- "EKE grows faster than KE" therefore holds at the ACC scale in all-sat only. In the Drake box it is significant in all-sat with either mask.

**3. Pooling across meanders (Stouffer).** Lag-1 cross-products of path length and downstream EKE are uncorrelated or negatively correlated between meanders (r = −0.16 to 0.00), so the independence assumption is not anticonservative. Downstream EKE itself co-varies between meanders after detrending (PAR–SEIR150 r = 0.29), possibly through common constellation effects.

## Revised trend test (option b, 2026-10-13)

**Why.** The 2026-10-08 check found the calibrated test still liberal (about 11 % false positives at nominal 5 % for AR(1) φ = 0.8). Cause: the AR model fitted to detrended residuals of a short persistent record underestimates persistence, and the calibration simulated its nulls from that underestimated model.

**Test now used for every reported trend** (`rstats.trend_bcz`; earlier functions kept unchanged for reproducibility).
- Statistic: z = deseasonalised Theil–Sen slope / sd, where sd is the exact OLS-slope sd under the fitted AR noise model (Toeplitz autocovariance).
- Noise model: AR coefficients bias-corrected by bootstrap (Kilian 1998), with the trend removal included in the bias simulation; non-stationary fits are shrunk to a maximum root of 0.99.
- p and 95 % CI: calibrated against 300 null series simulated from the bias-corrected model and processed identically (p floor 0.0033).
- Monthly series are tested with two noise models, AIC-selected order and fixed AR(12). **Reported p = the larger of the two; CI = the union.** Neither alone holds its level for every persistence structure (table below).
- Annual series: AR(1) with iterated bias correction.
- Maps (`a1_maps_bc.py`): same statistic per 1° cell, with the calibration pooled over 10 persistence classes; p = Gaussian tail of |z| / κ, κ = max(q95 / 1.96, q99 / 2.576) of the class null; BH-FDR at 0.10.

**False-positive rate at nominal 5 %** (`stats_size_check_bc.py`; trend-free series with a seasonal cycle)

| Process (n = 300 months unless stated) | Old calibrated | Bias-corrected, not calibrated | New, AIC order | New, AR(12) |
|---|---|---|---|---|
| AR(1) φ = 0.8 | 10.7 % | 9 % | 5.3 % (300 series) | 8 % (100) |
| AR(1) φ = 0.9 (independent check) | — | — | 3.5 % (200) | — |
| AR(12) fitted to Drake-box EKE | — | — | 9.5 % (200) | 3 % (100) |
| Annual AR(1) 0.5, n = 25 | — | — | 5.0 % (200; AR(1), iterated correction) | — |
| Annual AR(1) 0.5, n = 15 | — | — | 9.5 % (200) | — |

- With the max rule, the rate is at most the smaller of the two columns in each row.
- AIC fails on the Drake-like process because it selects AR(2) and drops the small high-lag coefficients that carry part of the low-frequency persistence; the fixed AR(12) is liberal for a pure AR(1).
- n = 15 annual series (in-situ transport) cannot be tested at the nominal level; that trend is not significant under either test.

**Independent check** (separate agent, own code): analytic slope sd matches Monte Carlo within 0.3 % (AR(1) 0.8, 0.95, AR(2), with gaps); bias correction reduces the AR-sum bias from −0.027 to −0.008 (AR(1) 0.8); gaps and month alignment handled correctly. It found a crash for non-stationary fits (fixed) and the liberal annual size (led to the iterated correction).

**Results** (`retest_all.py`, `retest_compare.py`; 585 tests). Significant at 5 %, old → revised: A5 145 → 107; A2 22 → 13; A6 14 → 9; A7 8 → 5; B1 14 → 9; Liu 35 → 29; contrasts 3 → 4 (one lost, two new tests significant). Every change is a loss except one geometry metric in the B1 pilot (core −0.10 m s50 sinuosity, p 0.089 → 0.033; one of 45 geometry tests, not interpreted).

| Result | Old p | Revised p | Status |
|---|---|---|---|
| ACC KE, along-stream speed, low-frequency EKE (both products) | ≈ 0.01 | 0.003–0.013 | **stands** |
| ACC EKE, all-sat / two-sat | ≈ 0.01 | 0.017 / 0.020 | stands |
| ACC sub-monthly EKE (SUB), all-sat / two-sat | ≈ 0.01 | **0.066** / 0.040 | lost in all-sat |
| KE of annual-mean flow, all-sat / two-sat | 0.001 | 0.017 / 0.020 | stands |
| **ACC EKE − KE (relative), all-sat** | ≈ 0.01 | **0.060** | **lost: "EKE grows faster than KE" is not supported in either product** |
| Pacific KE, EKE (both products) | ≈ 0.01–0.02 | 0.007–0.033 | stands |
| Atlantic EKE, all-sat | 0.02 | 0.053 | lost |
| Indian EKE, all-sat | 0.01 | 0.010 | stands (two-sat n.s., as before) |
| Sector ranking (Pacific vs Atlantic / Indian, relative trends) | not tested | 0.05–0.64 | not supported, except Pacific > Indian KE in two-sat (p = 0.0498) |
| Box south of 35°S, EKE 2001–2025, all-sat / two-sat | 0.01 / 0.03 | 0.027 / 0.043 | stands |
| A2 frontal latitude (ACC, sectors) | n.s. | n.s. | unchanged |
| A2 envelope width w25 / FWHM (ACC) | 0.013 / 0.005 | 0.047 / 0.017 | borderline / stands |
| A2 longitude bins with p < 0.05 | 9 | 5 | 4 lost |
| Drake-box EKE, ice mask, 2001–2025 / 2001–2024 | 0.04 / 0.01 | **0.0498 / 0.057** | borderline |
| Drake-box EKE, no mask | 0.09 | 0.083 | n.s. |
| Drake-box EKE − KE, ice mask / no mask | ≈ 0.01 | 0.037 / 0.043 | stands, weaker |
| Drake transport proxy T_perm, T_all (2001–2025) | n.s. | n.s. | unchanged |
| T_all 2005–2019; ΔADT 2005–2019 | 0.03; 0.01 | 0.060; 0.060 | lost |
| Wind stress, Drake box / circumpolar 45–60°S | 0.03 / 0.07 | 0.030 / 0.076 | unchanged |
| A7 WAKE − REST, EKE (primary) | 0.94 / 0.80 | 0.95 / 0.80 | unchanged (no concentration) |
| A7 2001–2017 low-frequency EKE WAKE − REST | 0.02 / 0.04 (FDR) | 0.070 / 0.086 | **lost** |
| A7 2001–2017 speed WAKE − REST | 0.01 / 0.01 (FDR) | 0.040 / 0.047 | nominal only; no A7 secondary passes FDR |
| A7 Drake downstream − upstream, KE / EKE | 0.01 / 0.02–0.03 | 0.010 / 0.047–0.053 | KE stands; EKE borderline |
| B1 pilot: meander core speed, downstream EKE | 0.01 | 0.020 / 0.023 | stands |
| B1 pilot: upstream core speed | 0.03 | 0.060 | lost |
| B1 path statistics (18 tests) | none FDR (2 nominal) | none at p < 0.05 | unchanged; path-latitude 95 % CIs now up to −0.19 to +0.26 ° decade⁻¹ (PAR meander section), wider than the ±0.15 quoted above |
| A5d boxes downstream of meanders, EKE (all-sat / two-sat) | Kerguelen 0.02 / 0.71; SEIR–Macq. 0.001 / 0.004; PAR 0.006 / 0.15; Scotia 0.001 / 0.09 | 0.073 / 0.69; 0.013 / 0.066; 0.056 / 0.25; 0.027 / 0.18 | only SEIR–Macquarie (all-sat) and Scotia (all-sat) remain; none in both products |
| Liu reproduction: speed (all sections except Upper-Ridge), half-max width | ≤ 0.02 | 0.003–0.02 | stands |
| Liu reproduction: Lower-Ridge position 1993–2023 / 1993–2025 | 0.02 / 0.04 | 0.047 / 0.10 | borderline / lost |
| Liu monthly preview: Upper-Ridge speed, all-sat / two-sat | 0.05 | 0.073 / 0.027 | product-dependent |

Not affected: B1 covariability (phase-randomised surrogates, not a trend test) and the 3-year block OLS tests.

**Maps** (revised, all-sat, 1° cells, `outputs/a1_maps_bc.nc`; the null |z| tails are heavier than Gaussian in the most persistent classes, q99/2.576 up to 1.39, so κ uses the larger of the q95 and q99 scalings):

| Field | Cells p < 0.05 (old) | FDR-significant (old) |
|---|---|---|
| KE | 236 of 8,184, 2.9 % (679) | 5 (58) |
| EKE | 328, 4.0 % (905) | 18 (133) |

- Fewer than 5 % of cells reach nominal p < 0.05, so local trends are not field-significant beyond a few cells.
- EKE FDR cells: south-east of Kerguelen (83–87°E, 52–54°S, 5 cells), south of Macquarie (156–162°E, 56.5°S, 4), near Campbell Plateau / New Zealand (165–172°E, 4), and 3 near the South American shelf edge (75–66°W, 48–52°S; possibly coastal artefacts, not checked). KE: 5 cells, 4 of them near New Zealand at 40–46°S.
- Figs 1 and 5 now use these maps; their stippling shows only these cells.

## Open items

- **Main figures, rebuilt in the v1 style** (`figstyle.py`, `fig1_main.py`–`fig4_main.py` → `outputs/figures/Fig1–4`):
  - Fig. 1: KE/EKE maps and trends;
  - Fig. 2: frontal position;
  - Fig. 3: standing meanders;
  - Fig. 4: Drake Passage, with the fixed-geometry proxy.
- Recompute the wind–transport lags with T_perm.
- Main-text p-values and stars in Figs 1, 2, 4 and 5 now come from the revised test (`figstyle.retest_override`); the statements in the sections above still quote the old p-values (see "Revised trend test").
- D2 list of v2 statements to change.
- Manuscript text: v2.1 built (`v21_text_changes.py`, `v21_build_docx.py`); SI sections S3 and S4 still to write.
- B2 (along-track) is optional and not started.
