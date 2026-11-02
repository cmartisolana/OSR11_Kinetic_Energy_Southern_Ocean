"""Write outputs/MANIFEST.csv: every distributed output with size, SHA-256, producing step, raw-data need and users.

Purpose : inventory of outputs/ for the archive (re-run after any step that rewrites outputs).
Inputs  : outputs/** (all files except MANIFEST.csv)
Outputs : outputs/MANIFEST.csv (path, bytes, sha256, producer, needs_raw_data, used_by)
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
    ("a1/*", "02g_eke_maps_bc", "yes (cache)", "Fig. 3a"),
    ("a2/tracks_baseline.npz", "02c_halfpower_tracks", "yes (cache)", "02e, 02j, 03t, 05a A2"),
    ("a2/a2h_*", "02e_halfpower_twosat", "yes", "Fig. S6"),
    ("a2/*", "02b_halfpower_sets", "yes (cache)", "05a A2 (old p column); SI S3"),
    ("a3/a3_realisations.csv", "02d_halfpower_detectability", "yes (cache)", "05c"),
    ("a3/a3_summary.csv", "05c_detectability_summary", "no", "Fig. S5"),
    ("a5/*", "02a_energy", "yes", "05a A5/CONTRAST, Fig. 3c, Fig. S1, 05g"),
    ("a6/*", "03u_drake_box_masks", "yes", "05a CONTRAST, 05g"),
    ("a7/*", "02m_wake_localisation", "yes", "05a A7 (Table S1)"),
    ("m1/*", "02l_earlier_studies", "yes", "Fig. S1d, Table S1"),
    ("d1/*", "03a_gould_track", "yes (+GV23)", "Fig. 2, 03d"),
    ("retest/retest_summary.csv", "05b_retest_summary", "no", "Figs. 3, S1, S5; Table S1; Sect. 3.2"),
    ("retest/*", "05a_retest", "no (saved series)", "05b"),
    ("stats/size_summary.csv", "05h_size_table", "no", "Sect. 2.5, SI S4"),
    ("stats/size_bc/*", "05g_size_bc / 05f_size_annual", "no (synthetic)", "05h"),
    ("stats/size_phi*", "05e_size_monthly", "no (synthetic)", "SI S4 (background)"),
    ("x1/x1_tracks.npz", "02j_contour_vs_energy", "yes", "02k, Fig. S7"),
    ("x1/x1_band_adt.npy", "02j_contour_vs_energy", "yes", "02k, Fig. S7b"),
    ("x1/x1_drake_*", "03t_contour_vs_energy_drake", "yes", "Fig. 2b, Sect. 3.1, Table S3"),
    ("x1/*", "02k_contour_vs_energy_stats", "no", "Fig. S7, Sect. 3.2/4, Table S3"),
    ("fig1/*", "figures/fig1_halfpower_method (cache)", "yes", "Fig. 1"),
    ("drake/dk1/*", "03r_drake_regions", "yes", "Fig. S2"),
    ("drake/dk2/*", "03q_section_fronts", "yes (+GV23)", "SI S2-S3 context"),
    ("drake/dk3/*", "03h_insitu_validation", "yes (+GV23)", "Fig. S4, SI S2"),
    ("drake/dk5/*", "03b_transect_trends + 03c_insitu_velocity_fix", "yes (+GV23)", "Figs. 2d, S3, S4; 02i, 03d, 03i-03k"),
    ("drake/dk6/*", "03m_section_halfpower_tests / 03n_section_halfpower_detect", "yes", "Sect. 3.1, SI S3"),
    ("drake/dk8/dk8_period_change_regions.csv", "05d_drake_period_change", "no", "Sect. 3.1, SI Drake-domain table"),
    ("drake/dk8/dk8_along_acc_fronts_periods_allsat.csv", "03j_along_acc_periods", "yes", "05d, Fig. S2b"),
    ("drake/dk8/*", "03i_along_acc", "yes", "Fig. 2c"),
    ("drake/dk9/*", "03p_halfpower_periods", "yes", "SI S3"),
    ("drake/dk10/*", "03o_halfpower_density_drake", "yes", "Fig. 2b (83 %)"),
    ("drake/dk12/*", "02i_band_lonbins", "yes", "Figs. 3d, 4, S2; 04b; Sect. 3.2"),
    ("drake/dk14/*", "02h_halfpower_density", "yes", "Fig. 3b (79 %)"),
    ("drake/dk24/*", "03d_section_streamwise", "yes", "Figs. 2d, S4; 03e-03g"),
    ("drake/dk25/*", "03e_profile_trends", "no (dk24 npz)", "Fig. S4"),
    ("drake/dk32/*", "04a_isopycnal_depth", "yes (Argo)", "04c"),
    ("drake/dk34/*", "03k_drake_maps", "yes", "Fig. 2a, S4"),
    ("drake/dk39/*", "03f_transport_bins", "no (dk24 npz)", "Fig. S4"),
    ("drake/dk40/*", "03g_net_surface_transport", "no (dk24 npz)", "Sect. 2.3, Table S2"),
    ("f3w/*", "04b_wind", "yes (wind)", "Fig. 4a"),
    ("f4/f4d_*", "04d_eke_wind_bins", "yes", "04e"),
    ("f4/*", "04c_slope_fields", "no (dk32 npz)", "Fig. 4, 04d-04g"),
    ("f4m/f4m_map.nc", "04f_slope_trend_map", "no (f4m cache)", "Fig. 4b"),
    ("f4m/*", "04e_slope_monthly", "yes (Argo)", "04f, 04g, Table S4"),
    ("f4e/f4e_common_cells.csv", "05i_common_cells_count", "yes (monthly files)", "SI S5 cell count"),
    ("f4e/*", "04g_common_cells", "yes", "Fig. 4c,e; Sect. 3.3; Table S5"),
    ("f4s/*", "04h_wind_sensitivity", "no (f4e series)", "Fig. 4d; Sect. 3.3; Table S5"),
]


def rule(rel):
    for g, prod, raw, use in RULES:
        if fnmatch.fnmatch(rel, g):
            return prod, raw, use
    return "?", "?", "?"


rows = []
for f in sorted(p for p in PA.OUT.rglob("*") if p.is_file() and p.name != "MANIFEST.csv"):
    rel = f.relative_to(PA.OUT).as_posix()
    prod, raw, use = rule(rel)
    rows.append(dict(path=rel, bytes=f.stat().st_size, sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
                     producer=prod, needs_raw_data=raw, used_by=use))
M = pd.DataFrame(rows); M.to_csv(PA.OUT / "MANIFEST.csv", index=False)
print(len(M), "files,", round(M.bytes.sum() / 1e6, 1), "MB; unmatched:", list(M.path[M.producer == "?"]))
