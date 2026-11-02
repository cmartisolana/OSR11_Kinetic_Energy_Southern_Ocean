# Working repository → publication repository

Old paths are relative to `archive/working_tree_2026-11-02/` (the working repository as of 2 Nov 2026, kept unchanged).
Logic was not changed unless stated in the new script's docstring; changes were limited to paths (`osr11.paths`), imports
and the repairs listed in README "Known limitations".

| Old | New |
|---|---|
| `revision/rstats.py`, `revision/fronts.py` | `src/osr11/rstats.py`, `src/osr11/fronts.py` (unchanged) |
| `revision/a2_run.py` (library part / main) | `src/osr11/halfpower.py` / `scripts/02_circumpolar/02b_halfpower_sets.py` |
| `revision/plan1019/drake/dk6/dk6_common.py` | `src/osr11/drake_halfpower.py` |
| `revision/style_nc.py`, `revision/figstyle.py` | `src/osr11/style.py` (= `figures_final/scripts/style.py`, final figure style) |
| `archive/first_submission/scripts/utils/{chambers_fronts,regions_and_masks}.py` | `src/osr11/fs_utils/` (unchanged; used by `01c` only) |
| `revision/extract_monthly.py`, `extract_wind.py`, `prepare_bathymetry.py` | `scripts/00_extract/00a_extract_monthly.py`, `00b_extract_wind.py`, `00c_prepare_bathymetry.py` |
| `revision/plan1019/drake/dk32/dk32_rg_extract.py` | `scripts/00_extract/00d_argo_rg_extract.py` |
| `archive/first_submission/notebooks/00_preprocessing.ipynb` | `scripts/01_base/01a_preprocessing.py` (code cells exported) |
| *(first-submission `half_power_EKE.py` output)* | `scripts/01_base/01b_halfpower_fields.py` (new generator, verified) |
| `archive/first_submission/scripts/run_all_two_sat.py` | `scripts/01_base/01c_halfpower_fields_twosat.py` |
| *(none; `clim_allsat.nc` had no script)* | `scripts/01_base/01d_climatology.py` (new, verified) |
| `revision/a5_energy.py` | `scripts/02_circumpolar/02a_energy.py` |
| `revision/a2_tracks.py` | `02c_halfpower_tracks.py` |
| `revision/a3_run.py` | `02d_halfpower_detectability.py` (+ `05_stats/05c_detectability_summary.py`, new) |
| `revision/a2h_twosat.py` | `02e_halfpower_twosat.py` |
| `revision/a1_maps.py`, `a1_maps_bc.py` | `02f_eke_maps.py`, `02g_eke_maps_bc.py` (`coarsen` moved to `osr11.halfpower`) |
| `revision/plan1019/drake/dk14/dk14_hp_density_circ.py` | `02h_halfpower_density.py` |
| `revision/plan1019/drake/dk12/dk12_circumpolar_band.py` | `02i_band_lonbins.py` |
| `revision/x1_contour_vs_energy.py` | `02j_contour_vs_energy.py` (+ `x1_band_adt.npy`) |
| `revision/x1_stats.py` | `02k_contour_vs_energy_stats.py` (+ bathymetric split, implied shift) |
| `revision/plan1019/m1_comparison.py` | `02l_earlier_studies.py` |
| `archive/revision_superseded/analyses_not_in_v2.1/a7_localisation.py` | `02m_wake_localisation.py` (back in the pipeline: Table S1 wake-vs-rest) |
| `revision/plan1019/d1b_gvtrack.py` | `scripts/03_drake/03a_gould_track.py` |
| `revision/plan1019/drake/dk5_transect_vs_gv23.py`, `dk5b_insitu_velocity_fix.py` | `03b_transect_trends.py`, `03c_insitu_velocity_fix.py` |
| `.../dk24/dk24_section_streamwise.py`, `dk25/dk25_profile_trends.py`, `dk39/dk39_transport_bins.py`, `dk40/dk40_net_transport.py` | `03d_section_streamwise.py`, `03e_profile_trends.py`, `03f_transport_bins.py`, `03g_net_surface_transport.py` |
| `.../dk3_validation.py` | `03h_insitu_validation.py` |
| `.../dk8/dk8_along_acc.py`, `dk8_along_acc_periods.py` | `03i_along_acc.py`, `03j_along_acc_periods.py` |
| `.../dk34/dk34_drake_fields.py` | `03k_drake_maps.py` |
| `.../dk6/dk6_load.py`, `dk6_run.py`, `dk6_detect.py`, `dk6_all.sh` | `03l_section_halfpower_load.py`, `03m_section_halfpower_tests.py`, `03n_section_halfpower_detect.py`, Makefile `drake` |
| `.../dk10/dk10_hp_density_p2.py`, `dk9/dk9_hp_periods.py` | `03o_halfpower_density_drake.py`, `03p_halfpower_periods.py` (import fixed) |
| `.../dk2_section_fronts.py`, `dk1_drake_trends.py`, `dk4_section_series.py` | `03q_section_fronts.py`, `03r_drake_regions.py`, `03s_section_series.py` |
| `revision/x1_drake.py` + `revision/extract_drake_cols.py` | `03t_contour_vs_energy_drake.py` (repaired; no file deletion) |
| `revision/a6b_drake_box_masks.py` | `03u_drake_box_masks.py` |
| `.../dk32/dk32_isopycnal_slope.py` | `scripts/04_wind_argo/04a_isopycnal_depth.py` |
| `revision/f3w_wind.py`, `f4_tilt.py`, `f4d_eke_wind.py`, `f4m_monthly.py`, `f4m_map.py`, `f4e_sat.py`, `f4s_sens.py` | `04b_wind.py`, `04c_slope_fields.py`, `04d_eke_wind_bins.py`, `04e_slope_monthly.py`, `04f_slope_trend_map.py`, `04g_common_cells.py`, `04h_wind_sensitivity.py` |
| `revision/retest_all.py` (families A2 A5 A7 CONTRAST), `retest_compare.py` | `scripts/05_stats/05a_retest.py`, `05b_retest_summary.py` |
| `revision/stats_size_check.py`, `stats_size_check_annual.py`, `stats_size_check_bc.py` | `05e_size_monthly.py`, `05f_size_annual.py`, `05g_size_bc.py` (stable seed) |
| *(new)* | `05d_drake_period_change.py`, `05h_size_table.py`, `05i_common_cells_count.py` |
| `figures_final/scripts/fig1_methods.py` … `figS7_contour_vs_energy.py` (2 Nov 2026) | `figures/fig1_halfpower_method.py`, `fig2_drake.py`, `fig3_circumpolar.py`, `fig4_wind_slope.py`, `figS1_products.py`, `figS2_drake_periods.py`, `figS3_velocity_transect.py`, `figS4_drake_validation.py`, `figS5_S6_halfpower.py`, `figS7_contour_vs_energy.py` (paths only) |
| `figures_final/scripts/x1_tracks.py` | not copied: its outputs (`x1_tracks.npz`, `x1_band_adt.npy`, `x1_drake_period_shift.json`) are written by `02j` and `03t` and are identical |

Output folders: `revision/outputs/X` → `outputs/X`; `revision/outputs/clim_allsat.nc`, root `outputs/{mean_adt,sea_ice_mask}.nc`,
`revision/data/bathymetry_so.nc` → `outputs/base/`; `revision/outputs/a1_maps*.nc` → `outputs/a1/`;
`revision/plan1019/drake/out/dkN` and `revision/plan1019/drake/dkN` → `outputs/drake/dkN`; `revision/plan1019/outputs/{d1,m1}` → `outputs/{d1,m1}`;
root `outputs/monthly_half-power_points*` → `$OSR11_DATA/` (not distributed).
