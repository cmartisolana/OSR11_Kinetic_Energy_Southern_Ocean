# Run order of the analysis. `make figures stats si` runs from the distributed outputs/ (no raw data; Fig. 1 needs one
# monthly file, see README). The other targets need the raw data under $OSR11_DATA (default data/) and, for `data`,
# a Copernicus Marine login. Each target lists its steps in dependency order.
PY ?= python
S := scripts
.PHONY: all data base process circumpolar drake windargo stats sizecheck figures si check clean-cache

all: data base process stats figures si

# ---- 00: extraction (Copernicus Marine login; GEBCO and Roemmich-Gilson files downloaded by hand, see README) ----
data:
	$(PY) $(S)/00_extract/00a_extract_monthly.py --product allsat
	$(PY) $(S)/00_extract/00a_extract_monthly.py --product twosat
	$(PY) $(S)/00_extract/00b_extract_wind.py
	@echo "Then, by hand: 00c_prepare_bathymetry.py <GEBCO nc>; 00d_argo_rg_extract.py <RG T nc> \$$OSR11_DATA/argo_rg/so_T and <RG S nc> .../so_S"

# ---- 01: base fields (needs raw data) ----
base:
	$(PY) $(S)/01_base/01a_preprocessing.py
	$(PY) $(S)/01_base/01b_halfpower_fields.py
	$(PY) $(S)/01_base/01c_halfpower_fields_twosat.py --years 2001 2025
	$(PY) $(S)/01_base/01d_climatology.py

# ---- 02-04: processing (needs raw data; hours) ----
process: drake circumpolar windargo

drake:            # the Drake steps first: 02i and several figures read the dk5 streamline edges
	$(PY) $(S)/03_drake/03a_gould_track.py
	$(PY) $(S)/03_drake/03b_transect_trends.py
	$(PY) $(S)/03_drake/03c_insitu_velocity_fix.py
	$(PY) $(S)/03_drake/03d_section_streamwise.py
	$(PY) $(S)/03_drake/03e_profile_trends.py
	$(PY) $(S)/03_drake/03f_transport_bins.py
	$(PY) $(S)/03_drake/03g_net_surface_transport.py
	$(PY) $(S)/03_drake/03h_insitu_validation.py
	$(PY) $(S)/03_drake/03i_along_acc.py
	$(PY) $(S)/03_drake/03j_along_acc_periods.py
	$(PY) $(S)/03_drake/03k_drake_maps.py
	$(PY) $(S)/03_drake/03l_section_halfpower_load.py
	for a in "allsat 0" "allsat 1" "twosat 0" "twosat 1"; do $(PY) $(S)/03_drake/03m_section_halfpower_tests.py $$a; done
	for s in 0 10 25; do $(PY) $(S)/03_drake/03n_section_halfpower_detect.py $$s 30 $$s; done
	$(PY) $(S)/03_drake/03o_halfpower_density_drake.py
	$(PY) $(S)/03_drake/03p_halfpower_periods.py
	$(PY) $(S)/03_drake/03q_section_fronts.py
	$(PY) $(S)/03_drake/03r_drake_regions.py
	$(PY) $(S)/03_drake/03s_section_series.py
	$(PY) $(S)/03_drake/03u_drake_box_masks.py

circumpolar:      # 02c (tracks_baseline.npz) must precede 03t
	$(PY) $(S)/02_circumpolar/02a_energy.py allsat
	$(PY) $(S)/02_circumpolar/02a_energy.py twosat
	$(PY) $(S)/02_circumpolar/02b_halfpower_sets.py
	$(PY) $(S)/02_circumpolar/02c_halfpower_tracks.py
	$(PY) $(S)/02_circumpolar/02d_halfpower_detectability.py 41
	$(PY) $(S)/02_circumpolar/02e_halfpower_twosat.py
	$(PY) $(S)/02_circumpolar/02f_eke_maps.py
	$(PY) $(S)/02_circumpolar/02g_eke_maps_bc.py
	$(PY) $(S)/02_circumpolar/02h_halfpower_density.py
	$(PY) $(S)/02_circumpolar/02i_band_lonbins.py
	$(PY) $(S)/02_circumpolar/02j_contour_vs_energy.py
	$(PY) $(S)/02_circumpolar/02k_contour_vs_energy_stats.py
	$(PY) $(S)/02_circumpolar/02l_earlier_studies.py
	$(PY) $(S)/02_circumpolar/02m_wake_localisation.py allsat
	$(PY) $(S)/02_circumpolar/02m_wake_localisation.py twosat
	$(PY) $(S)/03_drake/03t_contour_vs_energy_drake.py

windargo:         # needs 02i (dk12) for 04b
	$(PY) $(S)/04_wind_argo/04a_isopycnal_depth.py
	$(PY) $(S)/04_wind_argo/04b_wind.py map
	$(PY) $(S)/04_wind_argo/04c_slope_fields.py
	$(PY) $(S)/04_wind_argo/04d_eke_wind_bins.py
	$(PY) $(S)/04_wind_argo/04e_slope_monthly.py
	$(PY) $(S)/04_wind_argo/04f_slope_trend_map.py
	$(PY) $(S)/04_wind_argo/04g_common_cells.py
	$(PY) $(S)/04_wind_argo/04h_wind_sensitivity.py
	$(PY) $(S)/05_stats/05i_common_cells_count.py

# ---- 05: statistics on the saved series (no raw data; minutes) ----
stats:
	$(PY) $(S)/05_stats/05a_retest.py A5 A7 CONTRAST
	$(PY) $(S)/05_stats/05b_retest_summary.py
	$(PY) $(S)/05_stats/05c_detectability_summary.py
	$(PY) $(S)/05_stats/05d_drake_period_change.py
	$(PY) $(S)/05_stats/05h_size_table.py
	$(PY) $(S)/02_circumpolar/02k_contour_vs_energy_stats.py
	$(PY) $(S)/03_drake/03g_net_surface_transport.py
	$(PY) $(S)/04_wind_argo/04h_wind_sensitivity.py

# synthetic false-positive checks of SI S4 (no data; hours). Logs go to outputs/stats/size_bc/.
sizecheck:
	$(PY) $(S)/05_stats/05f_size_annual.py 15 0.5 200 > outputs/stats/size_bc/ann15.log
	$(PY) $(S)/05_stats/05f_size_annual.py 25 0.5 200 > outputs/stats/size_bc/ann25.log
	$(PY) $(S)/05_stats/05g_size_bc.py ar1_0.8 300 outputs/stats/size_bc z > outputs/stats/size_bc/z_ar1_0.8.log
	$(PY) $(S)/05_stats/05g_size_bc.py real_Drake_EKE 200 outputs/stats/size_bc z > outputs/stats/size_bc/z_real_Drake_EKE.log
	RSTATS_MONTHLY_ORDER=12 $(PY) $(S)/05_stats/05g_size_bc.py real_Drake_EKE 125 outputs/stats/size_bc z > outputs/stats/size_bc/ord12_z_real_Drake_EKE.log
	$(PY) $(S)/05_stats/05h_size_table.py

# ---- figures (from outputs/; Fig. 1 also reads $OSR11_DATA/monthly_allsat/2001/01.nc) -> rendered/png, rendered/pdf ----
figures:
	$(PY) figures/fig1_halfpower_method.py
	$(PY) figures/fig2_drake.py
	$(PY) figures/fig3_circumpolar.py
	$(PY) figures/fig4_wind_slope.py
	$(PY) figures/figS1_products.py
	$(PY) figures/figS2_drake_periods.py
	$(PY) figures/figS3_velocity_transect.py
	$(PY) figures/figS4_drake_validation.py
	$(PY) figures/figS5_S6_halfpower.py
	$(PY) figures/figS7_contour_vs_energy.py

# ---- SI: figures S1-S7 (above) and the tables of si/tables from the outputs listed in si/README.md ----
si: figures
	$(PY) $(S)/05_stats/05d_drake_period_change.py
	@echo "SI tables: see si/README.md (table -> output -> script). The SI text is si/SI.md."

check:            # byte-compile everything outside archive/ without writing __pycache__
	PYTHONPYCACHEPREFIX=$${TMPDIR:-/tmp}/osr11_pyc $(PY) -m py_compile src/osr11/*.py src/osr11/fs_utils/*.py scripts/*/*.py figures/*.py

clean-cache:
	rm -rf cache outputs/drake/dk6/fields_*.npz outputs/drake/dk5/dk5_section_fields.npz
