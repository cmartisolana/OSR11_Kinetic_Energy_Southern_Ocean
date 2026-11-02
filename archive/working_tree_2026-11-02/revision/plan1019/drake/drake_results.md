# Drake Passage first: altimetry checked against in situ data, then extended (19 Oct 2026, first pass)

**Logic of the section.** Over the overlap with the LM Gould record (2005–2019), check that altimetry reproduces the in situ velocity, EKE, transport and front position. Then use altimetry for the longer record (2001–2025, later 1993–2025).

Data and methods:
- Test: `rstats.trend_bcz`.
- In situ data: Gutierrez-Villanueva et al. (2023), cloned from manuelogtzv/TrendsDrakePassageTransp.
- Scripts: `revision/plan1019/drake/`. Figure: `FigureD1_drake_validation_extension`.

## 1 KE, EKE, along-stream speed and wind, 2001–2025 (`dk1_drake_trends.py`)

Trends in % of mean per decade, [95 % CI], p.

| Region | Quantity | All-sat | Two-sat |
|---|---|---|---|
| Drake box (ice-masked) | KE | +2.4 [−2.5, 7.4] 0.29 | +1.4 [−3.4, 6.2] 0.52 |
| | EKE | +8.4 [0.2, 16.6] 0.050 | +3.8 [−2.5, 10.1] 0.22 |
| | along-stream speed | −0.1 0.91 | 0.0 0.98 |
| | EKE − KE (points per decade) | +6.0 [0.7, 11.3] 0.04 | +2.4 [−1.7, 6.4] 0.22 |
| Upstream 100–80°W (ACC band) | KE | +12.8 [6.8, 18.7] 0.007 | +10.8 [5.8, 15.8] 0.007 |
| | EKE | +15.3 [6.2, 24.4] 0.013 | +11.5 [4.7, 18.3] 0.007 |
| | along-stream speed | +3.9 [0.3, 7.5] 0.050 | +3.6 [0.0, 7.3] 0.056 |
| Downstream 55–35°W (ACC band) | KE | +4.3 0.13 | +3.3 0.21 |
| | EKE | +8.0 [3.0, 13.0] 0.010 | +3.8 0.13 |

Eastward wind stress: Drake box +8.7 % [1.4, 15.9] (p = 0.03); upstream +9.8 % (p = 0.03); downstream +9.4 % (p = 0.09).

**Reading**
- Inside the passage, only all-sat EKE rises, and only at the limit of significance; the two-satellite product does not confirm it.
- The robust signal is **upstream** (southeast Pacific approach to the passage): KE, EKE and the along-stream speed all increase in both products, alongside increased wind stress.
- Downstream EKE is robust only in all-sat.

## 2 Front position by half-power along the section (`dk2_section_fronts.py`)

Method: KE of the monthly-mean flow along the mean LM Gould track, D1 baseline (200 cm² s⁻², 100 km, 25 %).

**Baseline: one envelope.** The climatological profile gives a single envelope (55.1–59.3°S, half-power point 57.0°S), i.e. the SAF and PF merge into one energetic band because KE stays above 200 cm² s⁻² between them. The band is detected in 99.7 % of months.
- Trend 2001–2025: −10 km per decade [−36, +16], p = 0.43 (all-sat); −8 [−44, +29], p = 0.65 (two-sat).
- Against the in situ position (the latitude splitting the 0–760 m transport within the same band into halves): r = 0.50 per transect and 0.71 for annual means (all-sat; two-sat 0.48 and 0.65). Mean offset 16 km (all-sat) and 5 km (two-sat).
- RMS difference 114 km: a monthly field against a 2-day transect, with altimetric spread 130 km vs in situ 77 km.

**Sensitivity: 400 cm² s⁻²** (inside the S3 range, chosen after seeing the profile, so post hoc). The SAF and PF separate.
- SAF (climatological half-power point 56.2°S): detected in 55 % of months; −12.6 km per decade (p = 0.050) in all-sat, −10.3 (p = 0.06) in two-sat. Against in situ: r = 0.47 per transect, 0.71 annual; offset −39 km.
- PF (58.3°S): detected in only 31 % of months, so no trend is computed. Against in situ: r = 0.48 per transect, 0.38 annual.
- With these coverage gaps the SAF trend is not robust.

## 3 Altimetry against in situ over 2005–2019 (`dk3_validation.py`, `d1b_gvtrack.py`, `dk4_section_series.py`)

| Quantity | All-sat | Two-sat |
|---|---|---|
| Cross-section velocity: mean profile, altimetry vs ADCP at 70–118 m | r = 0.98 | 0.97 |
| Velocity per bin and transect (8,187 pairs) | r = 0.72 | 0.72 |
| Velocity anomalies (bin mean removed) | r = 0.61 | 0.60 |
| RMS velocity difference | 15 cm s⁻¹ | 16 cm s⁻¹ |
| EKE north of 58.5°S, per transect | r = 0.67 | 0.65 |
| EKE north of 58.5°S, annual means | **r = 0.89** | **0.85** |
| EKE over all ice-free points, annual | r = 0.92 | 0.88 |
| Transport north of 58.5°S, annual (GV track / v2 line) | r = 0.67 / 0.72 | — |
| Band position, annual | r = 0.71 | 0.65 |

Same-period EKE trends north of 58.5°S, 2005–2019 (none significant with 15 annual values):
- in situ: +7.8 % per decade (0–780 m);
- all-sat: +10.1 %;
- two-sat: +0.8 %.

So the in situ EKE change over the overlap is closer to the all-satellite product.

**What the check supports**

| Quantity | Interannual variability reproduced? | Basis for extending the record |
|---|---|---|
| Surface velocity | yes (r ≈ 0.7) | supported |
| EKE | yes (r ≈ 0.85–0.9) | supported |
| Transport north of 58.5°S | yes (r ≈ 0.7) | supported for the northern part |
| Band position | yes (r ≈ 0.7) | supported |
| Full-section transport | **no** (r ≈ 0.05) | **not supported**: the part south of 58.5°S, where the record is affected by sea ice, is anticorrelated with the north (r = −0.55) |
| Trends | cannot be checked | 15 years is too short for a trend test |

The in situ EKE is 0–780 m depth-integrated while altimetry is the surface, so the agreement concerns variability, not magnitude.

**Extended record, 2001–2025 (same quantities, along the section)**
- EKE north of 58.5°S: +7.4 % per decade [−5.6, 20.4], p = 0.24 (all-sat); +1.8 %, p = 0.68 (two-sat).
- Transport north of 58.5°S: +3.9 % [−2.1, 9.9], p = 0.17 (all-sat); +4.2 %, p = 0.15 (two-sat).

## 4 Transport after Hogg et al. (2015) (`h0_antarctic_slope_index.py` → `h2_hogg_vs_insitu.py`)

The scripts are ready and smoke-tested on the 60–65°S part only; that run is not a result. Once the southern extraction exists:
1. `python revision/plan1019/h0_antarctic_slope_index.py revision/data/monthly_allsat_south revision/data/bathymetry_80s.nc revision/plan1019/drake/out/h0_allsat.csv 2001 2025`
2. `python revision/plan1019/drake/h2_hogg_vs_insitu.py revision/plan1019/drake/out/h0_allsat.csv revision/plan1019/outputs/d1 <GV repo> revision revision/plan1019/drake/out/h2_allsat.csv`

Output: index and T_hogg (labelled model-calibrated, Hughes et al. 2014 factors) against the in situ total, north and south parts, on matched transects and annual means, plus the T_hogg trend.
