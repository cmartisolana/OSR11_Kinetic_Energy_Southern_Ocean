# Archive (superseded material; nothing deleted)

`working_tree_2026-11-02/` is a complete, unchanged copy of the working repository as it stood on 2 November 2026
(`OSR11/repository/`: `revision/`, root `outputs/`, root `archive/`, README, CITATION, requirements), with three exceptions:

- the raw monthly altimetry (`revision/data/monthly_allsat/`, 2.4 GB) is not copied (raw data are not distributed);
- `__pycache__/`, `*.pyc` and `.DS_Store` are not copied;
- its `.gitignore` is renamed `_gitignore_original.txt`, because a nested `.gitignore` would otherwise exclude parts of the
  archive (its unanchored `outputs/` pattern) from version control.

The scripts and outputs of the final pipeline were **copied** out of this tree into `src/`, `scripts/`, `figures/` and
`outputs/` (mapping in `docs/script_map.md`); their originals remain here, so the archive also documents what each
publication script was derived from. Everything else here is superseded and not needed for any result of the paper:

| Where (under `working_tree_2026-11-02/`) | What |
|---|---|
| `archive/first_submission/` | discussion-paper scripts, notebooks and figures |
| `archive/revision_superseded/` | earlier revision figure scripts and manuscript builders, tracked eddies, Drake experiments (dk17–dk31, dk41–dk44, h0, h2), analyses not in v2.1 (incl. the original `a7_localisation.py`) |
| `archive/fix_paths_2026-10-26.py` | the path fix of 26 Oct 2026 |
| `outputs/` (root) | first-submission products: `trends/` (2-D trend maps, in situ series), `eof/`, `plots/`, `preprocessing.nc` (input of `01b`), `mean_adt.nc`, `sea_ice_mask.nc` (copied to `outputs/base/`) |
| `revision/a1_recompute.py`, `a2e_*`, `a2g_*`, `a6_drake.py`, `a6c_*`, `a6d_*`, `a8*`, `a9_sam.py`, `b1_*`, `b2_*`, `b3_*`, `f2s_slope.py`, `f4a*`, `f4b*`, `f4c*`, `f4j*`, `f4v*`, `f4w*`, `f4x*`, `stats_contrast_tests.py` | analyses withdrawn or superseded (transport proxies, SAM, meander paths, along-stream speed, older wind/slope variants, older contrast test) |
| `revision/extract_drake_cols.py` | **do not run**: deletes every monthly file it reads; replaced by `scripts/03_drake/03t_contour_vs_energy_drake.py` |
| `revision/x1_drake.py` | does not compile (interleaved docstrings); repaired as `03t` |
| `revision/x1_figure.py`, `revision/figures_v2.1/`, `revision/archive/` | earlier figure scripts (replaced by `figures/`) and the old SI source/builder (`SI_v2.1_lean.md`, `build_si_lean.py`) |
| `revision/plan1019/` (except the scripts mapped in `docs/script_map.md`) | d1, m2, h1, fig_plan1019, extract_c3s_twosat, run_extractions.sh, author comments on referee replies (`ac/`), dk1x–dk4x experiments, `dk6_all.sh` (absolute session path) |
| `revision/outputs/` (except the files copied to `outputs/`) | outputs of the archived analyses; `retest/retest_{A6,B1,LIU}.csv`; `stats/contrast_tests.csv` (older contrast test, p = 0.0099, superseded by the CONTRAST family, p = 0.060); `figures_final/` (stale renders of 19 Oct) |
| `revision/R2_revision_results.md`, `plan1019/plan1019_results.md`, `plan1019/drake/drake_results.md`, `dk6/dk6_results.md`, `outputs/x1/x1_results.md` | result logs of the revision (they quote p-values of older tests) |
