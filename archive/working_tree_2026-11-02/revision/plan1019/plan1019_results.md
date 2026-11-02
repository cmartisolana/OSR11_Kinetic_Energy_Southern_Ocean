# Revision plan of 19 Oct 2026: results of the first pass

Status 19 Oct 2026. Everything uses the existing monthly fields (2001–2025, 65–35°S) and `rstats.trend_bcz` for every trend (deseasonalised Theil–Sen; calibrated AR noise; larger p of the AIC and AR(12) models). Scripts and outputs are in `revision/plan1019/`. No existing file was changed. Values come from first runs: verify them before transcription.

## 0 Corrections to premises of the plan

- **Hogg et al. (2015) is mainly an EKE paper.** It uses crossover along-track EKE in three sector boxes:
  - Indian 40–150°E, 57–44°S;
  - Pacific 150–288°E, 62–48°S;
  - Atlantic 325–10°E, 56–46°S.

  Trends: Indian 18.3 ± 5.1, Pacific 14.9 ± 4.1, Atlantic 4.0 ± 3.7 cm² s⁻² per decade (90 % level). Relative to the 1993–2001 mean: Indian +6.0 ± 1.6 % and Pacific +6.4 ± 1.8 % per decade. Crossover EKE is 1.6–1.9 times the gridded EKE, so only the relative trends are comparable.

  The transport part is short. It uses AVISO sea level over the shelf and slope with T = −3.69(h − h_f) − 1.11 h_f Sv cm⁻¹ (signs as extracted; check the PDF). The factors are model-derived (Hughes et al. 2014, NEMO ORCA025/ORCA12) and called "preliminary". **No numerical transport trend is given**, only "indications of a weak decline".
- **Hughes et al. (2014)** define the slope region from model correlations (< −0.5). They average 50 depth bins of 100 m (0–5000 m); the strongest correlations are at 1–3 km. Their −3.69 Sv mbar⁻¹ is for bottom pressure. Sea level becomes amplified at periods longer than 5 years.
- **Martínez-Moreno et al. (2021)** used all-sat vDT2018, which is no longer distributed. Barceló-Llull et al. (2025, co-authors) found large version sensitivity and could not test all-sat; two-sat vDT2018 is still on C3S. MM21's procedure, taken from their code repository (github.com/josuemtzmo/EKE_SST_trends):
  - EKE = ½(ugosa² + vgosa²) from DUACS v2019, i.e. anomalies relative to the 1993–2012 mean;
  - coarsened to 1°, 365-day rolling mean, ×ρ0, area-integrated (PJ m⁻¹);
  - Theil–Sen slope on the smoothed daily series; significance from a modified Mann–Kendall test.

  Their values (1993 to Mar 2020): "ACC" region (all ocean south of 35°S, within a 65°S–65°N domain) +0.0368 PJ m⁻¹ per decade, **+5.1 %**; Southern Ocean basin +0.05 PJ m⁻¹, **+4.0 %**; global +0.09 PJ m⁻¹, +1.2 %.
- **In situ series:** GV2023 **total** cross-track velocity, 0–760 m (mean 75.5 Sv). This includes the reference (depth-independent, 51.1 Sv) and Ekman (8.1 Sv) parts. The geostrophic part relative to 760 m is 15.9 Sv. The series used throughout `revision/` is the total one.

## 1 Drake Passage: coverage-matched comparison (`d1_drake_coverage.py`, `d1b_gvtrack.py`; outputs `d1/`)

Inputs: the GV2023 per-bin transports (25 km bins, cloned from manuelogtzv/TrendsDrakePassageTransp) and a surface proxy along the **actual mean LM Gould track**, keeping only the points that are ice-free in every month.

| | v2 line | GV track |
|---|---|---|
| Points kept / 121 | 53 | 58 |
| Southern limit of kept points | 58.47°S | 58.55°S |
| Share of full-line mean surface proxy carried by kept points | 65 % | 67 % |
| Share of mean ADT drop across the line | 63 % | — |
| Share of in situ total (0–760 m) north of the limit | 62 % | 63 % |

Correlations, GV track (248 transects; annual means 2005–2019):

| In situ part | vs T_perm, transects | vs T_perm, annual | vs T_all, transects |
|---|---|---|---|
| Total | 0.03 | 0.03 | 0.07 |
| North of 58.55°S | **0.54** | **0.67** (surrogate p = 0.003) | 0.14 |
| South of 58.55°S | −0.59 | −0.72 | −0.09 |

- The northern and southern parts of the in situ transport are anticorrelated (r = −0.54).
- The v2 line gives the same picture: 0.44 for transects, 0.72 for annual means.

Trends (2005–2019 annual means): in situ total −1.2 % per decade (p = 0.70), north part −0.9 % (p = 0.87), south part +0.2 % (p = 0.98). T_perm on the GV track is +3.9 % per decade over 2001–2025 (p = 0.17).

**Reading.** The "no agreement" between the two proxies comes from mismatched coverage. Over the same latitudes the surface proxy tracks the in situ total transport moderately well. The total is flat partly because transport is redistributed across about 58.5°S. This fits GV2023's frontal compensation but does not test it.

**Independent check** (separate code, written without seeing these scripts): every number reproduced within ±0.02. Deseasonalising changes nothing (north part vs T_perm r = 0.50). On the v2 line the transect correlation is 0.44 and the annual one 0.72.

**Caveats**
- The north/south anticorrelation is partly built into the method: a jet crossing a fixed latitude moves transport from one side to the other even when the total is unchanged.
- The 248 transects fall in 137 months and are serially correlated; averaged by month (n = 137) the north-part correlation is 0.66.
- The annual correlation rests on 15 values.

Note: the "gould_approx" line in `a6c` ends near 58.5°W; the real track ends near 62.4°W, 62.3°S.

**Text affected:** AC B5 ("the two proxies are not correlated … at monthly or annual time scales") and the agreement paragraph of the results file (A6).

## 2 Matched comparison with earlier estimates (`m1_comparison.py`; outputs `m1/`)

EKE trend in % of the period mean per decade, [95 % CI], p. Mask: pixels ice-free in every month (as in `a5`). The 35–65°S box is the domain of the monthly fields, not MM21's definition.

| Region | Period | All-sat | Two-sat | All − two (points per decade) |
|---|---|---|---|---|
| Hogg Indian | 2001–2012 | +5.8 [−17.8, +29.3] 0.47 | +9.6 [+1.7, +17.5] 0.03 | −1.3 [−91, +89] 0.95 |
| | 2001–2020 | +8.7 [+1.7, +15.7] 0.03 | +5.7 [+1.3, +10.2] 0.02 | +3.3 [−10.2, +16.9] 0.48 |
| | 2001–2025 | +9.0 [+4.9, +13.0] 0.007 | +4.8 [+1.7, +7.9] 0.013 | +4.5 [−5.8, +14.7] 0.20 |
| Hogg Pacific | 2001–2012 | +0.4 [−29.2, +30.1] 0.97 | +2.7 [−17.9, +23.3] 0.68 | +0.2 0.97 |
| | 2001–2020 | +5.8 [+0.8, +10.8] 0.04 | +3.1 [−2.3, +8.4] 0.20 | +3.1 0.31 |
| | 2001–2025 | +8.6 [+0.9, +16.4] 0.04 | +4.9 [+0.5, +9.4] 0.05 | +4.3 [−0.8, +9.4] 0.09 |
| Hogg Atlantic | 2001–2025 | +8.8 [−0.1, +17.6] 0.06 | +5.2 [−0.2, +10.7] 0.06 | +3.9 0.17 |
| Box 35–65°S | 2001–2020 | +4.6 [−1.9, +11.1] 0.11 | +1.2 [−1.5, +4.0] 0.36 | +3.5 0.35 |
| | 2001–2025 | +6.6 [+1.5, +11.6] 0.03 | +2.3 [+0.2, +4.5] 0.04 | +4.5 [−1.4, +10.4] 0.10 |
| ACC band | 2001–2020 | +6.4 [+0.9, +11.9] 0.04 | +3.2 [−0.1, +6.5] 0.07 | +3.2 0.37 |
| | 2001–2025 | +8.0 [+3.3, +12.8] 0.017 | +3.8 [+1.3, +6.3] 0.020 | +4.3 [−2.2, +10.8] 0.12 |

Absolute trends, 2001–2025 (cm² s⁻² per decade):

| Region | All-sat | Two-sat |
|---|---|---|
| Hogg Indian | +18.0 [9.9, 26.2] | +8.3 |
| Hogg Pacific | +15.7 [1.7, 29.6] | +7.5 |
| Hogg Atlantic | +13.5 | +6.9 |
| ACC band | +19.8 | +8.1 |

Area integral of the 35–65°S box (ρ0 = 1025 kg m⁻³), 2001–2025: all-sat +0.092 PJ m⁻¹ per decade [0.022, 0.162] on a mean of 1.40 PJ m⁻¹; two-sat +0.029.

Low-frequency part on a common 1/4° grid (all-sat block-averaged): ACC all-sat +5.5 % (p = 0.003), two-sat +3.2 % (p = 0.013), difference +2.3 points (p = 0.18).

**Readings**
- 12-year windows (2001–2012, the overlap with Hogg's record) are uninformative: CIs of ±20–30 % per decade.
- Hogg's relative trends (+6.0 and +6.4 % per decade, 1993–~2012) lie between our two products for 2001–2025, inside the CIs of both products in both sectors. This is consistency, not a reproduction: the periods, sampling and reference means all differ.
- Extending 2001–2020 to 2001–2025 increases every trend and narrows every CI. The 35–65°S box all-sat trend goes from p = 0.11 to p = 0.03.
- **The all-sat − two-sat difference is not significant in any region or period** (p ≥ 0.09; with scale matching p ≥ 0.18). The products differ by about a factor of two in magnitude, but over 25 years that difference is within the trend uncertainty. "About half as large in the two-satellite product" should be stated with this caveat.

### 2b Our data processed with the MM21 procedure (`m2_mm21_procedure.py`; outputs `m2/`)

Setup: EKE ×ρ0 integrated over 35–65°S (all pixels valid in every month, no ice mask), 12-month running mean, Theil–Sen. The calibrated test (`trend_bcz`) is applied to the unsmoothed monthly series.

| Product | Period | MM21 procedure, % per decade | Calibrated test, % per decade [95 % CI], p | Mean (PJ m⁻¹) |
|---|---|---|---|---|
| All-sat | 2001–2020 | +5.2 | +4.5 [−2.6, +11.6], p = 0.13 | 1.44 |
| All-sat | 2001–2025 | +7.5 | +6.6 [+1.2, +12.0], p = 0.04 | 1.48 |
| Two-sat | 2001–2020 | +1.8 | +1.2 [−1.5, +3.9], p = 0.36 | 1.28 |
| Two-sat | 2001–2025 | +2.7 | +2.3 [+0.2, +4.5], p = 0.04 | 1.29 |

- Same order as MM21's +5.1 % (1993–2020, vDT2018 all-sat), but the periods and versions differ.
- Our integrated EKE level is about twice theirs (1.44 vs about 0.72 PJ m⁻¹). Possible causes are the grid and version (1/8° DT2024 vs 1/4° v2019) and our EKE definition; this was not investigated.
- 2001–2020 is not significant with the calibrated test, so it is the extension to 2025 that makes the increase detectable.

## 3 Half-power meanders (`h1_halfpower_meanders.py`; specification in the docstring, written before any geometry trend; outputs `h1/`)

**Gate G1: continuous branch at ≥ 90 % of window longitudes**

| Site (window) | Usable longitudes | G1 |
|---|---|---|
| SEIR150 (140–150°E) | 96 % | **pass** |
| KERG (76–84°E) | 41 % | fail |
| PAR (135.1–109.4°W) | 21 % | fail |
| DRAKE (70–55°W) | 14 % | fail |

Post-hoc check, disclosed and not part of the gate: starting the branch inside each site's ACC band does not rescue any failed site. Even the fraction of longitudes with *any* well-covered envelope is 69 % (KERG), 24 % (PAR) and 70 % (DRAKE). **The plan's "begin in Drake Passage" is not feasible with half-power positions.**

**Gate G2: synthetic tests at SEIR150** (40 realisations per scenario, observed anomalies resampled in 12-month blocks, imposed change over 25 years)
- 25 km northward translation: 21.3 km recovered (85 %).
- Asymmetric intensification (north flank ×1.2): spurious displacement **+10.0 km**, 47 % of the recovered translation. This passes the < 50 % criterion only narrowly.
- Uniform intensification: −1.2 km displacement.
- Flexing (25 km meander-amplitude change): amplitude response 8.3 km, against a null spread (1 sd) of 29.4 km, so it is **not detectable**. A pure translation also produces +9.5 km of spurious amplitude.
- Uniform intensification changes sinuosity by −0.076 (null sd 0.056) and ∫|κ| by −3.7 (null sd 3.8), while flexing changes neither (−0.021, +0.2).

So at SEIR150 only the mean position is interpretable, and even that partly responds to asymmetric intensification. Amplitude, sinuosity and curvature from half-power positions **cannot separate flexing from intensification**.

**SEIR150 results**
- Position trend: −4.3 km per decade [−14.3, +5.8], p = 0.37.
- Upstream jet-core speed leads northward displacement: r rises from 0.07 at lag 0 to 0.31 at lag 3 months (surrogate p = 0.004; passes BH-FDR over the 20 relation tests). The peak sits at the edge of the pre-specified window. Given the G2 result, part of this could be energy redistribution rather than displacement.
- Displacement shows no link to meander-section EKE (p = 0.89) or downstream EKE (p = 0.27).
- Amplitude–downstream EKE (r = −0.29, lag 3) and amplitude–downstream speed (r = −0.31) also pass FDR, but amplitude failed G2 and is not interpretable.

**Reading.** The half-power measurement does not add a usable meander-geometry test beyond one site. It supports the ACs' choice of the ADT jet-core contour (B1–B3). If anything enters the paper, a short SI note on the gates is the honest version.

## 4 Antarctic-slope index (Hogg/Hughes): tooling only (`h0_*.py`)

The script is ready (depth-bin averaging as in Hughes et al., 12 longitude bands, ice coverage per month, h_f fit, and conversion labelled "model-calibrated"). It needs new extractions south of 65°S and bathymetry to 80°S; the commands are in the docstring.

Smoke test on the part of the slope already on disk (60–65°S, 1000–3000 m): in a typical month only **47 %** of the slope points are ice-free, and 164 of 300 months have < 50 % coverage. Sea-ice sampling will dominate the real index further south.

## 5 Items left to the authors
See `plan_tasklist.md`: 1993–2000 extraction, the MM21 Methods, two-sat vDT2018, the southern extraction for the Hogg index, and the AC text changes.
