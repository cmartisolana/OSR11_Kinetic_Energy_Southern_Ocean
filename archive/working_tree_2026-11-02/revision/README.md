# revision/ — analyses and figures of the revised manuscript (v2.1, October 2026)

Everything behind the revised manuscript and its Supporting Information (SI) lives here. The first submission is
kept under `../archive/first_submission/`; superseded figure versions, manuscript builders and analyses that are
not used in v2.1 are under `../archive/revision_superseded/` (see `../README.md`). Results and their history are
in `R2_revision_results.md` (circumpolar and along-stream analyses) and `plan1019/plan1019_results.md`,
`plan1019/drake/drake_results.md` (Drake Passage); the p-values quoted in those files before 13 October 2026 come
from the earlier calibrated test and the revised value of every reported trend is in `outputs/retest/retest_summary.csv`.

## Figures of the manuscript and SI (`figures_v2.1/`)

All final figures share one style, `style_nc.py` (Nature-family layout: 180 mm width, 5–7 pt sans-serif,
lowercase bold panel letters, no in-panel statistics, one colour per quantity; helpers from `figstyle.py`).
Run from `revision/` (the scripts put `revision/` on the path for `style_nc`, `figstyle`, `rstats`, `fronts`).

| Figure | Script | Usage |
|---|---|---|
| Figure 1 — half-power method | `fig1_methods.py` | `python figures_v2.1/fig1_methods.py . <figure dir> [example lon, default 147.9]` |
| Figure 2 — Drake Passage | `fig2_drake.py` | `python figures_v2.1/fig2_drake.py . plan1019/drake plan1019/drake/out plan1019/drake/out/dk10/<density npz> <figure dir>` |
| Figure 3 — circumpolar | `fig3_circumpolar.py` | `python figures_v2.1/fig3_circumpolar.py . plan1019/drake/out plan1019/outputs <figure dir>` |
| Figure S1 — two-satellite product | `figS1_twosat.py` | `python figures_v2.1/figS1_twosat.py <figure dir>` |
| Figure S2 — Drake Passage, both products | no generating script archived (interactive session); panels from `plan1019/drake/out/dk{1,2,3,4,8,12}` and `plan1019/outputs/d1` | |
| Figure S3 — cross-track velocity trend | `figS3_velocity_transect.py` | `python figures_v2.1/figS3_velocity_transect.py` (writes next to the script) |
| Figure S4 — fixed-latitude profiles | `figS4_drake_validation.py` | see docstring |
| Figures S5, S6 — half-power sensitivity, both products | `figS5_S6_halfpower.py S2 S3` | `OSR11_CACHE=<cache> python figures_v2.1/figS5_S6_halfpower.py S2 S3` (writes `outputs/figures_final/FigureS2`, `FigureS3`; renumbered S5, S6 in the SI) |
| SI document | `build_si_lean.py` + `SI_v2.1_lean.md` | `python figures_v2.1/build_si_lean.py <unpacked SI docx> <unpacked main docx> figures_v2.1/SI_v2.1_lean.md <dir with new_S1..new_S6.png> <out dir>` |

`figS5_S6_halfpower.py` also contains the earlier S1/S4/S5 panels (two-satellite with the withdrawn transport
proxy, tracked eddies, path–EKE covariation); they are not used.

## Data extraction (Copernicus Marine login)

```bash
python extract_monthly.py --product allsat
python extract_monthly.py --product twosat
python extract_wind.py
python prepare_bathymetry.py /path/to/GEBCO_global.nc
```

Restartable; finished months are skipped. Output goes to `data/` (about 28 GB, not tracked).

## Analyses

| Script | Content | Outputs |
|---|---|---|
| `rstats.py` | Trend statistics; `trend_bcz` is the test used for every reported trend (deseasonalised Theil–Sen, bias-corrected AR noise, calibrated with 300 null series, AIC and AR(12) noise models, larger p) | — |
| `stats_size_check*.py`, `stats_contrast_tests.py` | False-positive rates of the tests (SI S4), contrasts between trends | `outputs/stats/` |
| `retest_all.py`, `retest_compare.py` | All reported trends recomputed with `trend_bcz` | `outputs/retest/retest_summary.csv` |
| `a1_recompute.py`, `a1_maps.py`, `a1_maps_bc.py` | KE/EKE trend numbers and 1° trend maps with FDR (Fig. 3a) | `outputs/a1_*` |
| `fronts.py`, `a2_run.py`, `a2_tracks.py`, `a2e_averaging.py`, `a2g_literature_regions.py`, `a2h_twosat.py` | Half-power detection, matched positions, sensitivity, longitude bins, regions of earlier studies, two-satellite product (SI S3) | `outputs/a2/` |
| `a3_run.py` | Injection–recovery (detectability; SI S3) | `outputs/a3/` |
| `a5_energy.py` | KE, EKE and its parts, along-stream speed, both products (Fig. 3c, SI S1) | `outputs/a5/` |
| `a6_drake.py`, `a6b_drake_box_masks.py`, `a6c_transport_fixed.py`, `a6d_*.py` | Drake box, upstream and downstream trends, masks, wind stress, section fields | `outputs/a6/` |
| `a8_surface_transport.py`, `a8b_transport_bins.py` | Surface transport proxy across the ACC band (withdrawn from the main text; kept for the record) | `outputs/a8/`, `a8b/` |
| `a9_sam.py` | Relation to the SAM index (one sentence in the Discussion) | `outputs/a9/` |
| `b1_*.py`, `b2_*.py` | Monthly covariation of path length and EKE at the standing meanders (along-path tests cited in the Discussion) | `outputs/b1/`, `b2/`, `b2_pf/` |
| `b3_alongstream.py`, `b3_trends_fast.py`, `b3_relations.py` | Along-stream trends on the jet-core path (Fig. 3b inputs) | `outputs/b3/` |
| `plan1019/d1_drake_coverage.py`, `d1b_gvtrack.py` | In situ record, mean LM Gould track, coverage-matched comparison (SI S2) | `plan1019/outputs/d1/` |
| `plan1019/m1_comparison.py`, `m2_mm21_procedure.py` | Comparison with Hogg et al. (2015) and Martínez-Moreno et al. (2021) (Fig. 3e) | `plan1019/outputs/m1/`, `m2/` |
| `plan1019/h1_halfpower_meanders.py` | Half-power meander gates (not in the lean SI; kept) | `plan1019/outputs/h1/` |
| `plan1019/drake/dk*` | Drake Passage: trends (dk1), section half-power (dk2, dk6, dk9–dk11), in situ comparison (dk3–dk5, dk7, dk24, dk25), EKE between the fronts (dk8, dk12, dk34, dk37, dk38), half-power density and periods (dk10, dk13–dk15), Argo isopycnals (dk32, dk33, dk36), transect profiles and net surface transport (dk35, dk39, dk40), Argo track (dk45) | `plan1019/drake/out/dk*` |

Environment variables: `OSR11_OUTPUTS` (default `outputs`), `OSR11_CACHE` (default `cache`, a 1.6 GB stack of the
monthly KE/EKE fields built on first use), `GV23_DIR` (unpacked in situ archive, default `data/TrendsDrakePassageTransp`).
On 26 October 2026 the absolute paths of the sessions in which 19 of these scripts were run were replaced by paths
relative to this directory (`../archive/fix_paths_2026-10-26.py`); the scripts have not been re-run since.

## Fixed decisions (made before looking at any trends)

| ID | Decision |
|----|----------|
| D1 | Half-power baseline: threshold 200 cm² s⁻², minimum width 100 km, envelope boundary 25 % of peak; detection on the KE of the monthly-mean velocity. |
| B2/B3 | Jet-core path: monthly +0.18 m relative-ADT contour (Subantarctic Front side), smoothed at 150 km; 30° windows (B2) and 5° bins (B3); Malvinas bins excluded. |
| A7 | Wake boxes fixed from bathymetry and climatological EKE before any trend was computed. |
| dk13 | Meander regions (upstream / ridge / downstream) fixed from bathymetry. |
