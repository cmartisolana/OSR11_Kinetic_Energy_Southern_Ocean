# Verification (2 November 2026)

Environment: Linux container, Python 3.13.18, packages as pinned in `requirements.txt`; 2 CPU cores. Available data: the
300 all-satellite monthly files (`monthly_allsat`, 2.4 GB), the distributed `outputs/`, and `preprocessing.nc` (from the archived working tree,
for `01b`). **Not available**: two-satellite monthly files, wind (`wind_taux.nc`), Roemmich–Gilson Argo, the
Gutierrez-Villanueva et al. (2023) archive and the 1.6 GB KE/EKE cache. Every run used `PYTHONPYCACHEPREFIX` outside the
repository (no `__pycache__` written).

## 1. Byte-compilation

`python -m py_compile` on all 81 Python files outside `archive/` (`src/osr11/**`, `scripts/**`, `figures/*`): **all pass**
(`make check`). The working copy had one failure (`x1_drake.py`), repaired as `scripts/03_drake/03t_contour_vs_energy_drake.py`.

## 2. Steps run here

"Reference copy" = the file as distributed, i.e. the working-repository output copied into `outputs/` before any run
(or, for files that first appeared in this repository, the output of the first run here, so the comparison then checks
run-to-run reproducibility). "values equal to x" = largest absolute difference of all numeric values (floating-point
summation order); "identical" = byte-identical. Times are wall-clock on 2 cores.

| Step | Command | Result | Time (s) | Outputs compared with the reference copy | SHA-256 (first 16) |
|---|---|---|---|---|---|
| 01b --check (vs first-submission 2001-01) | `python3 scripts/01_base/01b_halfpower_fields.py --check /mnt/user-data/uploads/OSR11/repository/outputs/monthly_half-power_points/2001/01.nc` | PASS | 0.7 | stdout: eke NaN pattern identical: True max /diff/: 6.700162746842153e-08 | – |
| 01d --check (vs working-repo clim_allsat.nc) | `python3 scripts/01_base/01d_climatology.py --check` | PASS | 18.8 | stdout: v NaN pattern identical: True max /diff/: 5.364418029785156e-07 | – |
| 02h | `python3 scripts/02_circumpolar/02h_halfpower_density.py` | PASS | 17.2 | `drake/dk14/dk14_hp_density_circ.npz`: identical | c986c13cedeb17ba |
| 02j | `python3 scripts/02_circumpolar/02j_contour_vs_energy.py` | PASS | 14.8 | `x1/x1_tracks.npz`: identical<br>`x1/x1_band_adt.npy`: identical | 2fca16ffe479ab06<br>db690b2f00d422b8 |
| 02k | `python3 scripts/02_circumpolar/02k_contour_vs_energy_stats.py` | PASS | 72.2 | `x1/x1_circumpolar_trends.csv`: identical<br>`x1/x1_track_bathy.csv`: identical<br>`x1/x1_bathy_split.csv`: identical<br>`x1/x1_band_adt_trend.csv`: identical<br>`x1/x1_implied_shift.csv`: identical | e0c54d63ffdb4944<br>fb90c8d2543e1af4<br>ea97a2a788be3043<br>51a2652726cd8083<br>175da284e355e93e |
| 03e | `python3 scripts/03_drake/03e_profile_trends.py` | PASS | 288.6 | `drake/dk25/dk25_profile_trends.csv`: values equal to 1.06e-11 | 2f15696de5cdd22b |
| 03f | `python3 scripts/03_drake/03f_transport_bins.py` | PASS | 154.5 | `drake/dk39/dk39_transport_bins.csv`: values equal to 4.58e-10 | bd63943006dcfc49 |
| 03g | `python3 scripts/03_drake/03g_net_surface_transport.py` | PASS | 8.8 | `drake/dk40/dk40_net_trends.csv`: values equal to 2.84e-12<br>`drake/dk40/dk40_net_annual.csv`: identical<br>`drake/dk40/dk40_net_monthly.csv`: identical | 401e9acae0197cfe<br>958241b6e145a5df<br>6d3a2c887bfcb5e5 |
| 03i | `python3 scripts/03_drake/03i_along_acc.py` | PASS | 134.5 | `drake/dk8/dk8_along_acc_fronts_allsat.csv`: values equal to 5.12e-13 | 5349a64594a83443 |
| 03j | `python3 scripts/03_drake/03j_along_acc_periods.py` | PASS | 168.2 | `drake/dk8/dk8_along_acc_fronts_periods_allsat.csv`: values equal to 1.11e-10 | fa4b7ffcbd40ba05 |
| 03k | `python3 scripts/03_drake/03k_drake_maps.py` | PASS | 13.9 | `drake/dk34/dk34_drake_fields.npz`: identical | ddbe9dbf089f0111 |
| 03o | `python3 scripts/03_drake/03o_halfpower_density_drake.py` | PASS | 16.3 | `drake/dk10/dk10_hp_density_p2.npz`: identical | 5d5a81f2a100ed51 |
| 03t | `python3 scripts/03_drake/03t_contour_vs_energy_drake.py` | PASS | 59.3 | `x1/x1_drake_trends.csv`: identical<br>`x1/x1_drake_period_shift.json`: identical<br>`x1/x1_drake_implied_shift.csv`: identical | f47e9223a6b28504<br>72e90cb05a5c7f52<br>d69aa76821e52393 |
| 04c | `python3 scripts/04_wind_argo/04c_slope_fields.py` | PASS | 1826.0 | `f4/f4_tilt_maps.nc`: values equal to 3.55e-15 (bytes differ)<br>`f4/f4_regional_trends.csv`: values equal to 1.03e-13<br>`f4/f4_lon5_trends.csv`: values equal to 1.90e-13 | 0fc003d6693cb8d7<br>fa446a79cd675bf6<br>f786daae6cd6fe4e |
| 04f | – | NOT RUN | – | – | – |
| 04h | `python3 scripts/04_wind_argo/04h_wind_sensitivity.py` | PASS | 0.8 | `f4s/f4s_results.csv`: values equal to 1.11e-16<br>`f4s/f4s_points.csv`: identical | b30bd8522357b09e<br>561682148a533923 |
| 05a | `python3 scripts/05_stats/05a_retest.py A5 A7 CONTRAST` | PASS | 0.4 | `retest/retest_A5.csv`: identical<br>`retest/retest_A7.csv`: identical<br>`retest/retest_CONTRAST.csv`: identical | 31b513876998ee50<br>e954f249342f781e<br>aecd767173db3268 |
| 05b | `python3 scripts/05_stats/05b_retest_summary.py` | PASS | 0.3 | `retest/retest_summary.csv`: identical | 61b7501000b221bc |
| 05c | `python3 scripts/05_stats/05c_detectability_summary.py` | PASS | 0.5 | `a3/a3_summary.csv`: identical | a28d170f1ed232a4 |
| 05d | `python3 scripts/05_stats/05d_drake_period_change.py` | PASS | 0.3 | `drake/dk8/dk8_period_change_regions.csv`: shape (3, 9) vs (3, 8) | 47d29481ab65c7eb |
| 05h | `python3 scripts/05_stats/05h_size_table.py` | PASS | 0.3 | `stats/size_summary.csv`: identical | 6daa9c0658afd32c |
| 05i | `python3 scripts/05_stats/05i_common_cells_count.py` | PASS | 32.6 | `f4e/f4e_common_cells.csv`: identical | dad864c2ce2a6e1c |
| fig1 (cache present) | `python3 figures/fig1_halfpower_method.py` | PASS | 4.0 | `fig1/clim_profile_147.9.npy`: identical | aa3a9dacc48951b8 |
| fig2 | `python3 figures/fig2_drake.py` | PASS | 2.7 | stdout: inset values: {'overlap_pct': np.float64(82.9), 'halfpower_late_minus_early_km': -2.6, 'contour_late_minus_early_km': -1 | – |
| fig3 | `python3 figures/fig3_circumpolar.py` | PASS | 4.9 | stdout: values: KE [4.940765147738973, 0.0033222591362126] EKE [8.015255528417716, 0.0166112956810631] positive bins 163 of 180 | – |
| fig4 | `python3 figures/fig4_wind_slope.py` | PASS | 4.4 | stdout: f4s: {'KE': {'elasticity': 0.491, 'el_lo': 0.064, 'el_hi': 0.793, 'p': 0.018}, 'EKE': {'elasticity': 0.912, 'el_lo': 0.1 | – |
| figS1 | `python3 figures/figS1_products.py` | PASS | 1.4 | stdout: ACC_KE {'allsat': np.float64(4.94), 'twosat': np.float64(3.63)} | – |
| figS2 | `python3 figures/figS2_drake_periods.py` | PASS | 1.5 | stdout: S2d period differences (%): {('allsat', 'up'): np.float64(25.3), ('allsat', 'box'): np.float64(16.5), ('allsat', 'down') | – |
| figS3 | `python3 figures/figS3_velocity_transect.py` | PASS | 1.2 | stdout: S3 in situ bins 34 sig 3 ; altimetry 2001-2025 sig 3 | – |
| figS4 | `python3 figures/figS4_drake_validation.py` | PASS | 1.5 | stdout: S4 correlations: mean velocity profile r = 0.98; transect bins r = 0.72; rms diff = 15.3 cm/s; annual EKE r = 0.89; tran | – |
| figS5_S6 | `python3 figures/figS5_S6_halfpower.py` | PASS | 1.4 | stdout: saved /home/claude/work/repository_publication/rendered/png/FigureS6.png 170 x 50 mm | – |
| figS7 | `python3 figures/figS7_contour_vs_energy.py` | PASS | 1.8 | stdout: S7 uncalibrated Theil-Sen check (km or cm per decade): {'lat_hp': np.float64(-0.66), 'lat_c': np.float64(-11.22), 'adt_h | – |

Notes on the table
- `04f_slope_trend_map.py` (Fig. 4b map calibration) was started and stopped by hand after a few minutes to keep the verification time bounded (8 workers on
  2 cores; the 04c calibration of the same type took 30 min; 04f runs from the distributed `f4m_p273_monthly.npz` without raw data); `outputs/f4m/f4m_map.nc` is the distributed
  file, untouched. Not run here: all steps needing the two-satellite files, wind, Argo, the in situ archive or the KE/EKE
  cache (`00*`, `01a`, `01c`, `02a`, `02b`–`02g`, `02i`, `02l`, `02m`, `03a`–`03d`, `03h`, `03l`–`03n`, `03p`–`03s`, `03u`,
  `04a`, `04b`, `04d`, `04e`, `04g`) and the long synthetic size checks (`05e`–`05g`; the distributed logs are used).
- `05d`: the reference has one column fewer (the `anomaly_difference_pct` cross-check column was added afterwards); the
  other columns are identical.
- `05a` is a no-op on the distributed family files (all rows present); `05b` regenerates `retest_summary.csv` from the four
  kept families (A2, A5, A7, CONTRAST) and reproduces every numeric value of those 375 rows of the working-repository summary
  exactly (the A6, B1 and LIU rows are no longer included; they are archived).
- `01b --check`: KE/EKE of 2001-01 rebuilt from the monthly file vs the first-submission file: identical NaN pattern,
  max |difference| 7.0e-8 (ke) and 6.7e-8 m² s⁻² (eke). `01d --check`: `clim_allsat.nc` rebuilt from the 300 monthly files:
  identical NaN pattern, max |difference| 1.4e-6 m (ADT), 8.3e-7 / 5.4e-7 m s⁻¹ (u, v).
- `fig1` was also run once with the cache removed: it recomputed the 2001–2025 profile from the 300 monthly files in 17 s,
  byte-identical to the distributed `outputs/fig1/clim_profile_147.9.npy`, and rendered the same figure.

## 3. Figures

All figure scripts were run from the repository (`make figures`) and compared with the figures delivered for the manuscript
(`figures_final/png`, same date): all eleven are identical pixel for pixel. PDF bytes differ between runs only by the embedded
creation date.

| Figure | Size (px) | Pixel comparison with figures_final/png | SHA-256 png (first 16) | SHA-256 pdf (first 16) |
|---|---|---|---|---|
| Figure1 | 4015 x 1913 | identical | 5cbbf512a26f7036 | 36eb08ab7cf9e90d |
| Figure2 | 4015 x 3448 | identical | 80d3bba968d5c386 | 1803699b81e44a23 |
| Figure3 | 4015 x 5196 | identical | 8701477273a9811b | 6ff9b8f42812d882 |
| Figure4 | 4015 x 5338 | identical | dbeee36c560632f3 | 8e103f6b6d2e89b6 |
| FigureS1 | 4015 x 3968 | identical | 3c1bd8ce240a7e8a | 9537fc1eb2fcc2b4 |
| FigureS2 | 4015 x 3307 | identical | 6afde401f6747cf5 | fd93d3110a096bca |
| FigureS3 | 1960 x 1653 | identical | 499974f4851c7a1e | a30a50fd08019b7f |
| FigureS4 | 4015 x 2976 | identical | c85b0e82b97319b4 | 58ff718c2d36eb70 |
| FigureS5 | 4015 x 2362 | identical | 45fed066b7f027c7 | afedb66e09a96019 |
| FigureS6 | 4015 x 1181 | identical | 2ecd352e1640494d | 4df97f748df63386 |
| FigureS7 | 4015 x 2220 | identical | 81629d8177b980d3 | dbe21e390e873cc0 |

## 4. Manuscript numbers recomputed from the repository outputs

Manuscript values from the revised text (`text/numbers_used.csv`) and the numbers audit (`audit/numbers_refs_report.md`).

| Sect. | Statement | Value from the repository outputs | Output (step) | Manuscript | Status |
|---|---|---|---|---|---|
| 2.2 | 900 retained longitudes | 900 | a2/tracks_baseline.npz (02c; shipped) | 900 | match |
| 2.3 | Qs: no detectable change (25 annual values, 150 surrogates) | +0.64 %/dec (-0.15 to +1.43), p = 0.146, n = 25 | drake/dk40/dk40_net_trends.csv (03g, re-run here) | p = 0.146 | match |
| 2.5 | false positives "3–5 %" monthly, "up to 9.5 %" for 15 annual values | monthly trend_bcz: AR(1) 0.8 0.053 (300 reps), AR(12) Drake 0.040 (125 reps); annual n=15 0.095 | stats/size_summary.csv (05h, parsed logs) | 3–5 %; 9.5 % | match, see note N3 |
| 3.1 | late vs early: 21–24 % upstream and passage, 12 % downstream | 23.9 / 21.2 / 11.6 % | drake/dk8/dk8_period_change_regions.csv (05d, new) | 23.9 / 21.2 / 11.6 | match |
| 3.1 | half-power latitude 70–55°W: no significant trend | -2.9 km/dec (-7.2 to +1.3), p = 0.16 | x1/x1_drake_trends.csv (03t, re-run here) | p = 0.16 | match |
| 3.1 | fixed contour −13 km/dec (−17 to −9), p = 0.003 | -12.75 (-16.9 to -8.6), p = 0.0033 | same | −13 (−17 to −9) | match (rounding) |
| 3.1 | sea-level rise at the half-power point 2.7 cm/dec | 2.65 | same | 2.7 | match |
| Fig. 2b | half-power −3 km, contour −19 km (2016–2025 minus 2001–2015) | -2.6 / -18.8 km | x1/x1_drake_period_shift.json (03t, re-run here) | −2.6 / −18.8 | match |
| S3 | implied contour shift, Drake | -8.4 km/dec (|dADT/dy| 3.2 cm per 10 km) | x1/x1_drake_implied_shift.csv (03t) | −8.4 (3.2) | match |
| S3/4 | implied contour shift, circumpolar | -9.5 km/dec (|dADT/dy| 2.8 cm per 10 km) | x1/x1_implied_shift.csv (02k) | −9.5 (2.8) | match |
| 3.1 | one-third of synthetic 25 km shifts detected | 0.33 of 30 realisations | drake/dk6/dk6_detect_S25_seed25.csv (03n; shipped) | 0.33 | match |
| 3.1 | streamline bins between PF and SAF: +3.1 and +5.2 cm/s/dec | +3.07, +5.20 (stream/allsat/u/2001-2025, x10) | drake/dk5/dk5_trends.csv (03b/03c; shipped) | 3.1, 5.2 | match (audit value) |
| 3.2 | ACC EKE 8.0 %/dec (19.8 cm²/s²/dec), p = 0.017 | 8.02 %, 19.8, p = 0.0166 | retest/retest_summary.csv (05b, re-run here) | 8.0; 19.8; 0.017 | match |
| 3.2 | ACC KE 4.9 %, p = 0.003 | 4.94 %, p = 0.0033 | same | 4.9; 0.003 | match |
| 3.2 | Atlantic EKE p = 0.053 | p = 0.0532 | same | 0.053 | match |
| 3.2 | EKE − KE +3.0 points/dec, p = 0.060 | 3.00, p = 0.0598 | same | 3.0; 0.060 | match |
| 3.2/S1 | relative wake-vs-rest EKE difference, p = 0.95 | p = 0.953 | retest/retest_A7.csv via summary | 0.95 | match |
| 3.2 | 91 % of the 2° bins positive | 163/180 = 90.6 % | drake/dk12/dk12_band_trends_allsat.csv (02i; shipped) | 91 % | match |
| 3.2 | no significant half-power trend over 900 longitudes | p = 0.74 | retest_summary (A2) | 0.74 | match |
| 3.2/4 | contour −11.5 km/dec; pathway interval excludes > ~4 km/dec | contour -11.50; half-power -0.54 (-4.17 to +3.09) | x1/x1_circumpolar_trends.csv (02k, re-run here) | −11.5; ±4 | match |
| 4 | over rough and smooth topography alike | rough -1.0 (p 0.53); smooth -0.1 (p 0.99) | x1/x1_bathy_split.csv (02k, new) | −1.0/0.53; −0.1/0.99 | match |
| 3.3 | widened-band EKE 8.9 % (3.3–14.5), p = 0.02 | 8.87 (3.29–14.46), p = 0.0166 | f4e/f4e_results.csv (04g; shipped) | 8.9 | match |
| 3.3 | slope 1.2 % (0.3–2.0), p = 0.02 | 1.18 (0.31–2.04), p = 0.0199 | same | 1.2 | match |
| 4 | wind +2.4 %/dec, not significant | 2.36, p = 0.146 | same | 2.4 | match |
| 3.3 | lag-1 r = 0.61, 0.57; p = 0.004, 0.008 | KE 0.614 (p_fw 0.0045); EKE 0.572 (p_fw 0.0080) | same | yes | match |
| 3.3 | coefficients EKE 0.91 (0.20–1.49), KE 0.49 (0.06–0.79), slope 0.09 (−0.01–0.22) | EKE 0.91 (0.20–1.49) ; KE 0.49 (0.06–0.79) ; slope 0.09 (-0.01–0.22) | f4s/f4s_results.csv (04h, re-run here) | as text | match |
| 3.3 | 13 field-significant slope cells, all positive | 13 cells, 13 positive | f4m/f4m_map.nc (04f; shipped) | 13, all positive | match |
| S5 | 2830 common cells | slope domain 2835.0; with KE/EKE every month 2835.0; wind not checkable here | f4e/f4e_common_cells.csv (05i, new) | 2830 | TODO (needs wind file) |

Notes
- N1. Values from steps not re-run here (e.g. `04g`, `02i`, `03n`) are read from the distributed outputs; their scripts were
  byte-compiled only.
- N2. The "rough and smooth topography" split, the implied contour shifts (circumpolar −9.5, Drake −8.4 km per decade), the
  ACC-band mean ADT and the late-vs-early percentages had no generating script before; the new code (`02k`, `03t`, `05d`)
  reproduces the values previously computed by hand (x1_results.md, SI tables, writer's aggregation).
- N3. False-positive rates: the parsed logs give 5.3 % (AR(1) φ = 0.8, 300 replicates) and 4.0 % (AR(12) fitted to Drake EKE,
  125 replicates) for monthly series with the final test, and 9.5 % for 15 annual values. Two logs are incomplete
  (`z_ar1_0.95`: non-stationary fit; `z_real_ACC_EKE`: no result line) and `ord12_z_ar1_0.8` stopped at 125 replicates with
  8.0 %; the "3–5 %" statement rests on the first two values and on the revision log (R2_revision_results.md, archived),
  which also lists AR(1) φ = 0.9 → 3.5 %. **AUTHOR CHECK**: whether 8.0 % (ord12, 125 replicates) should be mentioned or the
  run completed. A re-run of `05g` now uses a stable seed and will differ slightly from the logs (README).
- N4. 2830 common cells (SI S5): not reproducible here without the wind file (2835 cells pass every other condition).

## 5. Repository checks

- No `__pycache__`, `*.pyc`, `.DS_Store` or symbolic links in the repository (checked with `find` before the commit).
- `CITATION.cff` validated with cffconvert 2.0.0 ("valid according to schema version 1.2.0").
- `outputs/MANIFEST.csv`: 142 files, 72.1 MB, SHA-256 of every file, producing step, raw-data need and users.
