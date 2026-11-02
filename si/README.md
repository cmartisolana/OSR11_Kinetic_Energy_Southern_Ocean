# Supporting Information: text and tables

- `SI.md`: SI text (S1–S5) as revised on 2 Nov 2026 (copy of the SI working file; the SI .docx is built outside this repository
  from this text and the figures in `rendered/`).
- `tables/*.csv`: the SI tables, one file per table, each with a `source` column. The table files were assembled from the
  outputs listed below (values copied and rounded; the `source` paths use the working-repository names, `revision/...`).
  The table below gives, for each table, the output files in this repository and the script that writes them.

| SI table (file) | Output read (this repository) | Generator | Needs raw data to regenerate the output? |
|---|---|---|---|
| Table S1 sector trends (`tableS1_sector_trends.csv`) | `outputs/retest/retest_summary.csv` (keys `allsat\|ACC_KE\|2001-2025`, ..., `ACC EKE-KE allsat/twosat`); `outputs/m1/m1_trends.csv` (`diff_all-two\|ACC_EKE\|2001-2025`); `outputs/retest/retest_A7.csv` (`<product>\|2001-2025\|EKE\|WAKE-REST`) | `05a_retest.py A5 A7 CONTRAST` → `05b_retest_summary.py`; series from `02a_energy.py`, `02m_wake_localisation.py`; `02l_earlier_studies.py` | `05b` no; `05a` no (reads saved series; rows present are skipped); `02a`, `02l`, `02m` yes |
| Table S2 Qs trends (`tableS2_Qs.csv`) | `outputs/drake/dk40/dk40_net_trends.csv` | `03g_net_surface_transport.py` (from `outputs/drake/dk24/dk24_section_fields_allsat.npz`, `03d`) | no (`03g` runs from the shipped section fields) |
| Table S2b Qs annual sampling (`tableS2b_Qs_annual_sampling.csv`) | `outputs/drake/dk40/dk40_net_annual.csv` | `03g_net_surface_transport.py` | no |
| Table S3 contour vs energy (`tableS3_contour_vs_energy.csv`) | `outputs/x1/x1_circumpolar_trends.csv`, `x1_drake_trends.csv`, `x1_drake_period_shift.json`, `x1_implied_shift.csv`, `x1_drake_implied_shift.csv` | `02k_contour_vs_energy_stats.py` (circumpolar; from `x1_tracks.npz`, `02j`), `03t_contour_vs_energy_drake.py` (Drake) | `02k` no; `02j`, `03t` need the monthly files |
| Table S4 slope by sector (`tableS4_slope_sectors.csv`) | `outputs/f4m/f4m_results.csv` (27.3, signed, monthly test); cell and FDR counts from `outputs/f4m/f4m_map.nc` | `04e_slope_monthly.py`; `04f_slope_trend_map.py` | `04e` yes (Argo); `04f` no (runs from the shipped `f4m_p273_monthly.npz`) |
| Table S5 wind sensitivity (`tableS5_wind_sensitivity.csv`) | `outputs/f4e/f4e_results.csv` (D1, D2); `outputs/f4s/f4s_results.csv` | `04g_common_cells.py`; `04h_wind_sensitivity.py` | `04g` yes; `04h` no |
| Drake domains, late vs early (`drake_domains_period_difference.csv`) | `outputs/drake/dk8/dk8_period_change_regions.csv` (columns `pct_change`, `anomaly_difference_pct`) | `05d_drake_period_change.py` (from `03j` and `02i` outputs) | no |

Notes
- Number of common cells of Table S5 / S5 text (2830): `04g_common_cells.py` now writes it to `outputs/f4e/f4e_common_cells_04g.csv`,
  but that step has not been re-run (it needs the wind file). `05i_common_cells_count.py` reproduces the mask from shipped files and
  the monthly altimetry: 2835 slope cells, all 2835 with KE/EKE in every month; the wind condition (wind file, author's Mac) is
  expected to remove the remaining 5. TODO (author): run `05i` with `data/wind_taux.nc` present to confirm 2830.
- The SI text also quotes false-positive rates (S4) from `outputs/stats/size_summary.csv` (`05h_size_table.py`, parsed from the
  size-check logs of `05f`/`05g`).
