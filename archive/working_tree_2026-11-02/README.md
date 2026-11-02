# Surface intensification of the Antarctic Circumpolar Current from 2001 to 2025

Code for Martí-Solana, C., Ruiz, S., Combes, V., Barceló-Llull, B., Ballarotta, M. and Pascual, A.:
*Surface intensification of the Antarctic Circumpolar Current from 2001 to 2025*, State of the Planet,
Ocean State Report 11, Chapter 2 (revised manuscript v2.1, October 2026).

Satellite altimetry (Copernicus Marine, 2001–2025) is used to quantify trends of kinetic energy (KE), eddy kinetic
energy (EKE), along-stream speed and the position of the energetic band of the ACC, with Drake Passage compared
with the shipboard ADCP record of Gutierrez-Villanueva et al. (2023). Everything that appears in the manuscript and
its Supporting Information (SI) can be regenerated from `revision/`; the material of the first (discussion-paper)
submission and of intermediate revisions is kept under `archive/` for traceability only.

## Layout

```
revision/                      the analysis of the revised manuscript (see revision/README.md for every script)
  rstats.py                    trend statistics (trend_bcz: the test used for every reported trend)
  fronts.py                    half-power detection
  extract_*.py, prepare_bathymetry.py     data extraction (Copernicus Marine login)
  a1_… a9_…, b1_… b3_…         analyses (KE/EKE, half-power, Drake Passage, along-stream)
  retest_*.py, stats_*.py      recomputation of all trends with the revised test; size checks
  plan1019/                    Drake Passage analyses (drake/dk*), in situ comparison (d1*), earlier-study
                               comparisons (m1, m2), half-power meander gates (h1)
  figures_v2.1/                scripts of the final figures and the SI (Table 1 below), SI text and builder
  outputs/                     analysis outputs (CSV/NPZ/NC) read by the figure scripts
  data/, cache/                raw monthly fields and caches (not tracked; ~30 GB)
outputs/                       preprocessed fields of the first submission still used here
                               (mean_adt.nc, sea_ice_mask.nc, trends/ with the in situ series)
archive/                       superseded material: first submission (scripts, notebooks, figures),
                               figure and manuscript versions v2.1–v3.6x, analyses not used in v2.1
                               (tracked eddies, meander curvature, Liu et al. reproduction, wind-memory,
                               scale decomposition, Argo transport, Antarctic-slope index)
```

## Table 1: where each figure comes from

| Figure | Script (in `revision/figures_v2.1/`) | Main inputs |
|---|---|---|
| Fig. 1 half-power method | `fig1_methods.py` | `data/monthly_allsat/2001/01.nc`, `outputs/mean_adt.nc`, `fronts.py` |
| Fig. 2 Drake Passage | `fig2_drake.py` | `plan1019/drake/out/dk{3,5,8,24,25,32–35,40,45}`, `plan1019/outputs/d1`, `plan1019/drake/dk10` density |
| Fig. 3 circumpolar | `fig3_circumpolar.py` | `outputs/a1_maps_bc.nc`, `a2/`, `a5/`, `b3/`, `retest/`, `plan1019/drake/out/dk{2,12,14,15}`, `plan1019/outputs/{d1,m1,m2}` |
| Fig. S1 two-satellite product | `figS1_twosat.py` | `outputs/a5/`, `outputs/retest/retest_summary.csv` |
| Fig. S2 Drake Passage, both products | *no generating script archived*: the figure was produced in an interactive session; its panels are the Product 1 / Product 4 series and trends of `plan1019/drake/out/dk{1,2,3,4,8,12}` and `plan1019/outputs/d1` | |
| Fig. S3 cross-track velocity trend | `figS3_velocity_transect.py` | `plan1019/drake/out/dk5/dk5_trends.csv` |
| Fig. S4 fixed-latitude profiles | `figS4_drake_validation.py` | `plan1019/drake/out/dk{3,5,24,25,39}` |
| Fig. S5, S6 half-power sensitivity and products | `figS5_S6_halfpower.py S2 S3` (writes FigureS2/FigureS3; renumbered S5/S6 in the SI) | `outputs/a2/`, `a3/`, `retest/` |
| SI document | `build_si_lean.py` from `SI_v2.1_lean.md` and the six figures | |

Numbers quoted in the text and the SI tables are in `revision/outputs/` (circumpolar: `a1_*`, `a2/`, `a5/`,
`retest/retest_summary.csv`, `stats/`), `revision/plan1019/drake/out/` (Drake Passage, `dk*`) and
`revision/plan1019/outputs/` (`d1`, `m1`, `m2`, `h1`). `revision/R2_revision_results.md` and
`revision/plan1019/plan1019_results.md` record the results and their history.

## Environment and data

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Input data (not distributed here):

- Altimetry: Copernicus Marine `SEALEVEL_GLO_PHY_L4_MY_008_047` (all-satellite, 1/8°) and
  `SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057` (two-satellite, 1/4°), daily geostrophic velocities 2001–2025:
  `python revision/extract_monthly.py --product allsat|twosat` (Copernicus Marine login).
- Wind stress: `WIND_GLO_PHY_L4_MY_012_006`: `python revision/extract_wind.py`.
- Bathymetry: GEBCO: `python revision/prepare_bathymetry.py <GEBCO global nc>`.
- Drake Passage in situ record: https://doi.org/10.5281/zenodo.10044261 (Gutierrez-Villanueva et al., 2023); set
  `GV23_DIR` to the unpacked archive (default `revision/data/TrendsDrakePassageTransp`).
- Argo (Fig. 2c): Roemmich–Gilson climatology and monthly extensions, https://sio-argo.ucsd.edu/RG_Climatology.html,
  placed in `revision/data/argo_rg/`.
- SAM index (Marshall, 2003): British Antarctic Survey; `revision/data/sam_*.txt`.

Environment variables: `OSR11_OUTPUTS` (default `revision/outputs`), `OSR11_CACHE` (default `revision/cache`, a
1.6 GB stack of the monthly KE/EKE fields built on first use), `GV23_DIR` (see above).

## Reproducing

1. Extraction (above), then `revision/a5_energy.py`, `a1_recompute.py`, `a1_maps_bc.py`, `a2_run.py`, `a3_run.py`,
   `a6_drake.py`, `b3_alongstream.py` → `retest_all.py` (all trends with `rstats.trend_bcz`).
2. Drake Passage: `revision/plan1019/d1_drake_coverage.py`, `d1b_gvtrack.py`, then `plan1019/drake/dk*` in numerical
   order (each script's docstring states its inputs); `m1_comparison.py`, `m2_mm21_procedure.py`.
3. Figures: the scripts of Table 1, run from `revision/` (`python figures_v2.1/fig3_circumpolar.py revision
   plan1019/drake/out plan1019/outputs <figure dir>`; usage lines in each docstring).
4. SI: `python figures_v2.1/build_si_lean.py <unpacked SI docx> <unpacked main docx> figures_v2.1/SI_v2.1_lean.md
   <dir with new_S1..new_S6.png> <out dir>`, then zip the output directory as .docx.

Scripts were run with Python 3.10–3.13, numpy 2.x, pandas 2.x, xarray, scipy, statsmodels, matplotlib, cartopy.
Absolute paths of the sessions in which the analyses were run were replaced by paths relative to `revision/` on
26 October 2026 (`archive/fix_paths_2026-10-26.py`, 19 scripts); the scripts have not been re-run since, so each
figure script should be run once before the repository is archived.

## Citation and licence

MIT License (see `LICENSE`). Cite the article and the Zenodo archive of this repository (`CITATION.cff`;
https://doi.org/10.5281/zenodo.21485764).
