# Run order of the analysis. `make figures notebooks test` run from the distributed outputs/ (no raw data). The other targets need
# the raw data under $OSR11_DATA (default data/) and, for `data`, a Copernicus Marine login. Each target lists its steps in
# dependency order.
PY ?= python
S := scripts
FIGOUT ?= figures_out
.PHONY: all test lint data base process circumpolar drake windargo product4 stats sizecheck figures notebooks manifest check clean-cache

all: data base process stats figures

test:               # package, trend test, half-power detection, every figure script, manifest
	$(PY) -m pytest -q tests
lint:
	ruff check src scripts figures tests

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
process: drake circumpolar windargo product4

drake:            # the Drake steps first: 02i and several figures read the dk5 streamline edges
	$(PY) $(S)/03_drake/03a_gould_track.py
	$(PY) $(S)/03_drake/03b_transect_trends.py
	$(PY) $(S)/03_drake/03c_insitu_velocity_fix.py
	$(PY) $(S)/03_drake/03d_section_streamwise.py
	$(PY) $(S)/03_drake/03e_profile_trends.py
	$(PY) $(S)/03_drake/03f_transport_bins.py
	$(PY) $(S)/03_drake/03h_insitu_validation.py
	$(PY) $(S)/03_drake/03i_along_acc.py
	$(PY) $(S)/03_drake/03j_along_acc_periods.py
	$(PY) $(S)/03_drake/03k_drake_maps.py
	$(PY) $(S)/03_drake/03o_halfpower_density_drake.py
	$(PY) $(S)/03_drake/03u_drake_box_masks.py
	$(PY) $(S)/03_drake/03v_transect_zones.py
	for v in 0 1 2; do $(PY) $(S)/03_drake/03w_adt_bin_drift.py $$v; done; $(PY) $(S)/03_drake/03w_adt_bin_drift.py

circumpolar:      # 02c (tracks_baseline.npz) must precede 03t
	$(PY) $(S)/02_circumpolar/02a_energy.py allsat
	$(PY) $(S)/02_circumpolar/02a_energy.py twosat
	$(PY) $(S)/02_circumpolar/02b_halfpower_sets.py
	$(PY) $(S)/02_circumpolar/02c_halfpower_tracks.py
	$(PY) $(S)/02_circumpolar/02e_halfpower_twosat.py
	$(PY) $(S)/02_circumpolar/02f_eke_maps.py
	$(PY) $(S)/02_circumpolar/02g_eke_maps_bc.py
	$(PY) $(S)/02_circumpolar/02h_halfpower_density.py
	$(PY) $(S)/02_circumpolar/02i_band_lonbins.py
	$(PY) $(S)/02_circumpolar/02j_contour_vs_energy.py
	$(PY) $(S)/02_circumpolar/02k_contour_vs_energy_stats.py
	$(PY) $(S)/02_circumpolar/02l_earlier_studies.py
	$(PY) $(S)/05_stats/05k_previous_studies.py --zhang-series   # raw-data step of 05k (output shipped)
	$(PY) $(S)/02_circumpolar/02m_wake_localisation.py allsat
	$(PY) $(S)/02_circumpolar/02m_wake_localisation.py twosat
	$(PY) $(S)/03_drake/03t_contour_vs_energy_drake.py
	$(PY) $(S)/02_circumpolar/02n_pathway_tracks.py      # energetic pathway of the paper (longitude-dependent threshold)
	$(PY) $(S)/02_circumpolar/02o_pathway_bins.py
	$(PY) $(S)/02_circumpolar/02p_pathway_products.py
	$(PY) $(S)/02_circumpolar/02q_pathway_occurrence.py

windargo:         # needs 02i (dk12) for 04b
	$(PY) $(S)/04_wind_argo/04a_isopycnal_depth.py
	$(PY) $(S)/04_wind_argo/04b_wind.py map
	$(PY) $(S)/04_wind_argo/04c_slope_fields.py
	$(PY) $(S)/04_wind_argo/04d_eke_wind_bins.py
	$(PY) $(S)/04_wind_argo/04e_slope_monthly.py
	$(PY) $(S)/04_wind_argo/04f_slope_trend_map.py
	$(PY) $(S)/04_wind_argo/04g_common_cells.py
	$(PY) $(S)/04_wind_argo/04i_wind_exponents.py
	$(PY) $(S)/04_wind_argo/04j_exponent_intervals.py
	$(PY) $(S)/04_wind_argo/04k_lag_correlation.py
	$(PY) $(S)/04_wind_argo/04l_cell_series.py
	$(PY) $(S)/04_wind_argo/04m_buffer_width.py

product4:         # two-satellite product (Fig. S3); after circumpolar and windargo
	$(PY) $(S)/06_product4/p4_01_band_series.py
	$(PY) $(S)/06_product4/p4_02_lonbins.py
	$(PY) $(S)/06_product4/p4_03_trend_map.py
	$(PY) $(S)/06_product4/p4_04_wind_domain.py
	$(PY) $(S)/06_product4/p4_05_lagcorr.py
	$(PY) $(S)/06_product4/p4_06_exponents.py
	$(PY) $(S)/06_product4/p4_08_relthr_tracks.py --thr relthr --ice none
	$(PY) $(S)/06_product4/p4_09_relthr_trends.py
	$(PY) $(S)/06_product4/p4_07_summary.py

# ---- 05: statistics on the saved series (no raw data; minutes) ----
stats:
	$(PY) $(S)/05_stats/05a_retest.py A5 A7 CONTRAST
	$(PY) $(S)/05_stats/05b_retest_summary.py
	$(PY) $(S)/05_stats/05h_size_table.py
	$(PY) $(S)/05_stats/05j_retention_trends.py
	$(PY) $(S)/05_stats/05k_previous_studies.py
	$(PY) $(S)/02_circumpolar/02k_contour_vs_energy_stats.py
	$(PY) $(S)/05_stats/05l_sam_lag.py

# synthetic false-positive checks of SI S1.7 (no data; hours). Logs go to outputs/stats/size_bc/.
sizecheck:
	$(PY) $(S)/05_stats/05f_size_annual.py 15 0.5 200 > outputs/stats/size_bc/ann15.log
	$(PY) $(S)/05_stats/05f_size_annual.py 25 0.5 200 > outputs/stats/size_bc/ann25.log
	$(PY) $(S)/05_stats/05g_size_bc.py ar1_0.8 300 outputs/stats/size_bc z > outputs/stats/size_bc/z_ar1_0.8.log
	$(PY) $(S)/05_stats/05g_size_bc.py real_Drake_EKE 200 outputs/stats/size_bc z > outputs/stats/size_bc/z_real_Drake_EKE.log
	RSTATS_MONTHLY_ORDER=12 $(PY) $(S)/05_stats/05g_size_bc.py real_Drake_EKE 125 outputs/stats/size_bc z > outputs/stats/size_bc/ord12_z_real_Drake_EKE.log
	$(PY) $(S)/05_stats/05h_size_table.py

# ---- figures 1-4 and S1-S7 from outputs/ -> $(FIGOUT)/png, $(FIGOUT)/pdf (rendered/ keeps the figures of the manuscript) ----
figures:
	$(PY) figures/fig1_halfpower_method.py $(FIGOUT)
	$(PY) figures/fig2_drake.py $(FIGOUT)
	$(PY) figures/fig3_circumpolar.py $(FIGOUT)
	$(PY) figures/fig4_wind_slope.py $(FIGOUT)
	$(PY) figures/figS1_drake_section.py outputs $(FIGOUT)
	$(PY) figures/figS2_timing.py outputs $(FIGOUT)
	$(PY) figures/figS3_products.py outputs $(FIGOUT)
	$(PY) figures/figS4_previous_studies.py outputs $(FIGOUT)
	$(PY) figures/figS5_pathway_robustness.py outputs $(FIGOUT)
	$(PY) figures/figS6_occurrence_adt.py outputs $(FIGOUT)
	$(PY) figures/figS7_position_longitude.py outputs $(FIGOUT)
	$(PY) scripts/compare_figures.py $(FIGOUT)/png rendered/png

# ---- notebooks: execute in place from outputs/ ----
notebooks:
	cd notebooks && for n in 0*.ipynb; do jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=900 $$n; done

manifest:           # rewrite outputs/MANIFEST.csv after any step that changes outputs/
	$(PY) scripts/make_manifest.py

check:              # byte-compile everything without writing __pycache__
	PYTHONPYCACHEPREFIX=$${TMPDIR:-/tmp}/osr11_pyc $(PY) -m py_compile src/osr11/*.py src/osr11/fs_utils/*.py scripts/*.py scripts/*/*.py figures/*.py

clean-cache:
	rm -rf cache outputs/drake/dk6/fields_*.npz outputs/drake/dk5/dk5_section_fields.npz
