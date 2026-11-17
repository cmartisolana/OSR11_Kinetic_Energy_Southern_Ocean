"""Write outputs/MANIFEST.csv: every distributed output with size, SHA-256, producing step, raw-data need and the figures that read it.

Purpose : inventory of outputs/ for the archive (re-run after any step that rewrites outputs).
Inputs  : outputs/** (all files except MANIFEST.csv)
Outputs : outputs/MANIFEST.csv (path, bytes, sha256, producer, needs_raw_data, read_by_figures)
Run     : python scripts/make_manifest.py
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import fnmatch
import hashlib
import pandas as pd
from osr11 import paths as PA

# (glob relative to outputs/, producer, needs raw data to regenerate, used by)
RULES = [
    ("base/mean_adt.nc", "01a_preprocessing", "yes", "all maps, ACC band"),
    ("base/sea_ice_mask.nc", "01a_preprocessing", "yes", "all maps"),
    ("base/clim_allsat.nc", "01d_climatology", "yes (monthly files)", "fronts in Figs. 2-4; 02i, 03i-03k, 04b"),
    ("base/bathymetry_so.nc", "00c_prepare_bathymetry", "yes (GEBCO)", "isobaths, depths; Figs. 2-4"),
    ("fig1/monthly_allsat_2001_01_uv.nc", "00a_extract_monthly (copy of one monthly file: u, v)", "yes", ""),
    ("x6/x6_relthr_tracks.npz", "02n_pathway_tracks", "yes", ""),
    ("x6/*", "02o_pathway_bins", "no (x6 tracks)", ""),
    ("x7/*", "02p_pathway_products", "yes", ""),
    ("s6/*", "02q_pathway_occurrence", "no (x7 points, not distributed)", ""),
    ("x3/*", "review of 5 Nov 2026 (occurrence overlap; read by 02q for its checks)", "yes", ""),
    ("c04/*", "03w_adt_bin_drift", "no (dk24 npz)", ""),
    ("x9/*", "04i_wind_exponents", "no (f4e series)", ""),
    ("x10/*", "04j_exponent_intervals", "no (f4e series)", ""),
    ("x12/*", "04k_lag_correlation", "no (f4e series)", ""),
    ("x13/*", "04l_cell_series", "yes (monthly files, wind, Argo cache)", ""),
    ("x14/*", "04m_buffer_width", "no (x13 cells)", ""),
    ("x11/marshall_sam_monthly_1997_2025.csv", "Marshall (2003) SAM index, British Antarctic Survey (transcribed 9 Nov 2026)", "external", ""),
    ("x11/*", "05l_sam_lag", "no (m1 series)", ""),
    ("p4/p4_09*", "06_product4/p4_08_relthr_tracks + p4_09_relthr_trends", "yes (Product 4 monthly files)", ""),
    ("p4/p4_lagcorr*", "06_product4/p4_05_lagcorr", "no", ""),
    ("p4/p4_exponents*", "06_product4/p4_06_exponents", "no", ""),
    ("p4/p4_summary.csv", "06_product4/p4_07_summary", "no", ""),
    ("p4/p4_checks.csv", "06_product4/p4_07_summary", "no", ""),
    ("p4/*", "06_product4/p4_01-p4_04", "yes (Product 4 monthly files)", ""),
    ("a1/*", "02g_eke_maps_bc", "yes (cache)", "Fig. 3a"),
    ("a2/tracks_baseline.npz", "02c_halfpower_tracks", "yes (cache)", "02e, 02j, 03t, 05a A2"),
    ("a2/a2h_*", "02e_halfpower_twosat", "yes", "Fig. S1c"),
    ("a2/a2_retention_trends.csv", "05j_retention_trends", "no", "SI S5 (retention 75 % / 50 %), Fig. S5a"),
    ("a2/*", "02b_halfpower_sets", "yes (cache)", "05a A2 (old p column); SI S5"),
    ("a5/*", "02a_energy", "yes", "05a A5/CONTRAST, Fig. 3c, Fig. S1, 05g"),
    ("a6/*", "03u_drake_box_masks", "yes", "05a CONTRAST, 05g"),
    ("drake/dk5/dk5_zones.csv", "03v_transect_zones", "yes (+GV23, dk5 cache)", "SI S3 (passage-coordinate zones)"),
    ("a7/*", "02m_wake_localisation", "yes", "05a A7 (Sect. 3.2, SI S2)"),
    ("m1/m1_previous_studies.csv", "05k_previous_studies", "no", "Fig. S4 (SI S4)"),
    ("m1/m1_zhang_pacific_62-40S_allsat.csv", "05k_previous_studies --zhang-series", "yes (monthly files)", "05k (Fig. S4b)"),
    ("m1/*", "02l_earlier_studies", "yes", "05k (Fig. S4); SI S2 (product difference)"),
    ("m2/*", "m2_mm21_procedure (working repository; not part of this version)", "yes", "05k (Fig. S4a: Martinez-Moreno et al. 2021)"),
    ("d1/*", "03a_gould_track", "yes (+GV23)", "Fig. 2, 03d"),
    ("retest/retest_summary.csv", "05b_retest_summary", "no", "Figs. 3, S1, S5; Sect. 3.2; SI S2"),
    ("retest/*", "05a_retest", "no (saved series)", "05b"),
    ("stats/size_summary.csv", "05h_size_table", "no", "Sect. 2.5, SI S1.7"),
    ("stats/size_bc/*", "05g_size_bc / 05f_size_annual", "no (synthetic)", "05h"),
    ("stats/size_phi*", "05e_size_monthly", "no (synthetic)", "SI S1.7 (background)"),
    ("x1/x1_tracks.npz", "02j_contour_vs_energy", "yes", "02k, 05k, Fig. S6"),
    ("x1/x1_band_adt.npy", "02j_contour_vs_energy", "yes", "02k, Fig. S6b"),
    ("x1/x1_drake_*", "03t_contour_vs_energy_drake", "yes", "Figs. 2b, S6c; Sect. 3.1; SI S6"),
    ("x1/*", "02k_contour_vs_energy_stats", "no", "Figs. S4c, S5a, S6; Sect. 3.2; SI S4-S6"),
    ("fig1/*", "figures/fig1_halfpower_method (cache)", "yes", "Fig. 1"),
    ("drake/dk3/*", "03h_insitu_validation", "yes (+GV23)", "Fig. S3a, SI S3"),
    ("drake/dk5/dk5_trend_profile_r.csv", "figures/figS3_drake_section", "no", "SI S3 (r of the trend profiles)"),
    ("drake/dk5/*", "03b_transect_trends + 03c_insitu_velocity_fix", "yes (+GV23)", "Figs. 2d, S3b; 02i, 03d, 03i-03k; SI S3"),
    ("drake/dk8/dk8_along_acc_fronts_periods_allsat.csv", "03j_along_acc_periods", "yes", "05d, Fig. S2b"),
    ("drake/dk8/*", "03i_along_acc", "yes", "Fig. 2c"),
    ("drake/dk10/*", "03o_halfpower_density_drake", "yes", "Fig. 2b (83 %)"),
    ("drake/dk12/*", "02i_band_lonbins", "yes", "Figs. 3d, 4, S2a; 04b; 05d; 05k (Fig. S4); Sect. 3.2"),
    ("drake/dk14/*", "02h_halfpower_density", "yes", "Fig. 3b (79 %)"),
    ("drake/dk24/*", "03d_section_streamwise", "yes", "Figs. 2d, S3; 03e-03g"),
    ("drake/dk25/*", "03e_profile_trends", "no (dk24 npz)", "Fig. S3 (ice-covered part)"),
    ("drake/dk32/*", "04a_isopycnal_depth", "yes (Argo)", "04c"),
    ("drake/dk34/*", "03k_drake_maps", "yes", "Figs. 2a, S3 (front levels)"),
    ("drake/dk39/*", "03f_transport_bins", "no (dk24 npz)", "Fig. S3b (bin latitudes)"),
    ("f3w/*", "04b_wind", "yes (wind)", "Fig. 4a"),
    ("f4/f4d_*", "04d_eke_wind_bins", "yes", "04e"),
    ("f4/*", "04c_slope_fields", "no (dk32 npz)", "Fig. 4, 04d-04g"),
    ("f4m/f4m_map.nc", "04f_slope_trend_map", "no (f4m cache)", "Fig. 4b"),
    ("f4m/*", "04e_slope_monthly", "yes (Argo)", "04f, 04g; Sect. 3.3, SI S1.5 (sector trends)"),
    ("f4e/*", "04g_common_cells", "yes", "Fig. 4c,e; Sect. 3.3"),
]


def rule(rel):
    for g, prod, raw, use in RULES:
        if fnmatch.fnmatch(rel, g):
            return prod, raw, use
    return "?", "?", "?"


import json
FIGS = json.load(open(PA.ROOT / "docs" / "figure_inputs.json"))          # files each figure script opens (traced; docs/figure_inputs.json)
LABEL = lambda n: ("Fig. S" + n[4] if n.startswith("figS") else "Fig. " + n[3])
rows = []
for f in sorted(p for p in PA.OUT.rglob("*") if p.is_file() and p.name not in ("MANIFEST.csv", "README.md")):
    rel = f.relative_to(PA.OUT).as_posix()
    prod, raw, _ = rule(rel)
    figs = ", ".join(LABEL(n) for n, files in FIGS.items() if rel in files)
    rows.append(dict(path=rel, bytes=f.stat().st_size, sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
                     producer=prod, needs_raw_data=raw, read_by_figures=figs))
M = pd.DataFrame(rows); M.to_csv(PA.OUT / "MANIFEST.csv", index=False)
print(len(M), "files,", round(M.bytes.sum() / 1e6, 1), "MB; unmatched:", list(M.path[M.producer == "?"]))
