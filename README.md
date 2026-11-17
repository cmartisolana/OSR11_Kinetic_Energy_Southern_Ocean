# Observational constraints on eddy saturation in the Antarctic Circumpolar Current

Analysis code, processed outputs, figures and notebooks for

> Martí-Solana, C., Ruiz, S., Combes, V., Barceló-Llull, B., Ballarotta, M. and Pascual, A.: Observational constraints on
> eddy saturation in the Antarctic Circumpolar Current, State of the Planet (Copernicus Ocean State Report 11, Chapter 2),
> manuscript sp-2026-11, 2026.

Code: Cristina Martí-Solana (IMEDEA, CSIC-UIB), cmarti@imedea.uib-csic.es. Archive: https://doi.org/10.5281/zenodo.21485764
(`CITATION.cff`); source: https://github.com/cmartisolana/OSR11_Kinetic_Energy_Southern_Ocean.
Repository version 2.0.0 (17 November 2026): the version of the final manuscript and Supplement; version 1.0.0-alpha accompanied
the discussion paper (`CHANGELOG.md`).

## What the paper does

Satellite altimetry (2001–2025) and the Roemmich–Gilson Argo fields (2004–2025) are used to examine three components of the
eddy-saturation hypothesis together: the eddy kinetic energy (EKE) of the Antarctic Circumpolar Current (ACC), the position of
the current, and the slope of an isopycnal, with their relation to the wind stress. In Drake Passage the altimetric surface
velocities are set against the shipboard ADCP record of Gutierrez-Villanueva et al. (2023). The position of the current is
tracked with an energy-based diagnostic, the half-power method of Chambers (2018) with a longitude-dependent threshold, and
compared with fixed contours of absolute dynamic topography (ADT), which sea-level rise displaces.

## Quick start: figures and notebooks without the raw data

Everything the figures and the notebooks need is in `outputs/` (127 MB). No download is needed.

```
pip install -r requirements.txt          # or: conda env create -f environment.yml
make figures                             # Figs. 1-4 and S1-S7 -> figures_out/png, figures_out/pdf   (about 1 min)
make notebooks                           # executes the four notebooks in place                       (about 1 min)
make test                                # every figure script runs; the outputs match their manifest  (about 2 min)
```

`rendered/` holds the figures as they appear in the manuscript and the Supplement. `make figures` writes to `figures_out/` so
that it never overwrites them; compare the two sets with `python scripts/compare_figures.py`.

## Layout

```
notebooks/           four notebooks, one per part of the paper, executed (tables of the quoted values + the figures)
figures/             one script per figure (Figs. 1-4, S1-S7); _sifig.py is the shared set-up of the Supplement figures
rendered/            the figures of the manuscript and the Supplement (png 600 dpi, pdf)
outputs/             processed outputs read by the figures and the notebooks; MANIFEST.csv (size, SHA-256, producing script,
                     figures that read each file); README.md explains the folder names
src/osr11/           shared code: paths.py (all locations, environment variables), params.py, rstats.py (trend tests),
                     fronts.py and halfpower.py (half-power detection and tracking), style.py (figure style)
scripts/00_extract/  downloads and conversions (Copernicus Marine login; GEBCO; Roemmich-Gilson)
scripts/01_base/     base fields: ice mask, mean ADT, monthly KE/EKE fields, 2001-2025 climatology
scripts/02_circumpolar/  circumpolar KE and EKE, energetic pathway (02n-02q), maps, comparison with earlier studies
scripts/03_drake/    Drake Passage: Gould track, in situ comparison, section diagnostics, ADT-bin drift (03w)
scripts/04_wind_argo/    isopycnal slope (Argo), wind stress, domain means, exponents and lagged correlations (04i-04m)
scripts/05_stats/    statistics on the saved series: trend tests, summaries, size checks, SAM lag (05l)
scripts/06_product4/ the same diagnostics for the two-satellite product (Fig. S3)
natural_earth/       Natural Earth land shapefiles (public domain) used by the maps, so that cartopy needs no download
docs/                figure_inputs.json (files each figure script opens)
tests/               pytest: package, trend test, half-power detection, every figure script
```

Every script starts with a docstring giving its purpose, inputs, outputs and run command. Scripts find the package themselves
(they add `src/` to `sys.path`), so `python scripts/.../xx.py` works from any directory.

## Figures

| Figure | Content | Script | Notebook | Reads `outputs/…` | Produced by (scripts) |
|---|---|---|---|---|---|
| Fig. 1 | half-power method | `figures/fig1_halfpower_method.py` | `01` | `base`, `fig1`, `x6` | `00a`, `01a`, `02n` |
| Fig. 2 | Drake Passage | `figures/fig2_drake.py` | `02` | `base`, `d1`, `drake/dk10`, `drake/dk24`, `drake/dk34`, `drake/dk5`, `drake/dk8`, `x1`, `x7` | `00c`, `01d`, `02p`, `03a`, `03b`, `03c`, `03d`, `03i`, `03k`, `03o`, `03t` |
| Fig. 3 | circumpolar changes | `figures/fig3_circumpolar.py` | `03` | `a1`, `a5`, `base`, `drake/dk12`, `drake/dk14`, `retest`, `x7` | `00c`, `01a`, `01d`, `02a`, `02g`, `02h`, `02i`, `02p`, `05b` |
| Fig. 4 | wind stress, eddy energy and isopycnal slope | `figures/fig4_wind_slope.py` | `04` | `base`, `drake/dk12`, `f3w`, `f4`, `f4e`, `f4m`, `x10`, `x12`, `x9` | `00c`, `01a`, `01d`, `02i`, `04b`, `04c`, `04f`, `04g`, `04i`, `04j`, `04k` |
| Fig. S1 | Drake Passage section, in situ comparison, drift of the ADT bins | `figures/figS1_drake_section.py` | `02` | `c04`, `drake/dk24`, `drake/dk25`, `drake/dk3`, `drake/dk34`, `drake/dk39`, `drake/dk5` | `03b`, `03c`, `03d`, `03e`, `03f`, `03h`, `03k`, `03w` |
| Fig. S2 | timing of the EKE increase | `figures/figS2_timing.py` | `02` | `a5`, `drake/dk12` | `02a`, `02i` |
| Fig. S3 | all-satellite and two-satellite products | `figures/figS3_products.py` | `03` | `a1`, `a5`, `base`, `drake/dk12`, `p4`, `retest`, `x7` | `00c`, `01a`, `01d`, `02a`, `02g`, `02i`, `02p`, `05b`, `06_product4/p4_01-p4_04`, `06_product4/p4_05`, `06_product4/p4_07`, `06_product4/p4_08`, `p4_09` |
| Fig. S4 | comparison with previous studies | `figures/figS4_previous_studies.py` | `03` | `m1`, `x11` | `05k`, `05l` |
| Fig. S5 | sensitivity of the half-power method | `figures/figS5_pathway_robustness.py` | `01` | `a2`, `retest`, `x1`, `x7` | `02k`, `02p`, `05b`, `05j` |
| Fig. S6 | occurrence of the half-power points, ADT at the pathway | `figures/figS6_occurrence_adt.py` | `03` | `base`, `s6`, `x1`, `x7` | `00c`, `01a`, `02j`, `02p`, `02q` |
| Fig. S7 | position by longitude | `figures/figS7_position_longitude.py` | `03` | `base`, `drake/dk12`, `x6` | `00c`, `01d`, `02i`, `02o` |

The list of files each figure script opens was traced when the figures were rendered (`docs/figure_inputs.json`); the producing
script of every file is in `outputs/MANIFEST.csv`. Main-figure scripts take an optional figure directory as argument, Supplement
scripts take `[outputs dir] [figure dir]`.

## Notebooks

| Notebook | Part of the paper | Figures |
|---|---|---|
| `01_energetic_pathway_method.ipynb` | Sect. 2.2: half-power method, retained tracks, sensitivity of the position trend | 1, S5 |
| `02_drake_passage.ipynb` | Sect. 3.1: EKE between the fronts, position, cross-track velocity in ADT bins, in situ comparison | 2, S1, S2 |
| `03_circumpolar_energy_and_position.ipynb` | Sect. 3.2: KE and EKE trends, energetic pathway against fixed ADT contours, products, earlier studies | 3, S3, S4, S6, S7 |
| `04_wind_eddy_energy_isopycnal_slope.ipynb` | Sect. 3.3: trends, lagged correlations, power-law exponents, width of the domain | 4 |

Each notebook loads the tables behind a result, prints the values quoted in the paper and draws the figures by calling the
scripts in `figures/`. They are stored executed, so they can be read on GitHub without running anything.

## Data (Table 1 of the paper)

| Ref. | Product | Access | Used for | Placed in |
|---|---|---|---|---|
| 1 | SEALEVEL_GLO_PHY_L4_MY_008_047, all-satellite DUACS DT2024, 1/8°, daily | https://doi.org/10.48670/moi-00148 | all altimetric diagnostics | `$OSR11_DATA/monthly_allsat/` via `00a` |
| 2 | WIND_GLO_PHY_L4_MY_012_006, monthly | https://doi.org/10.48670/moi-00185 | wind stress (Fig. 4) | `$OSR11_DATA/wind_taux.nc` via `00b` |
| 3 | Drake Passage ADCP velocity, EKE and transport, LM Gould (Gutierrez-Villanueva et al., 2023) | https://doi.org/10.5281/zenodo.10044261 | in situ comparison (Fig. S1) | unpacked to `$GV23_DIR` |
| 4 | SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057, two-satellite, 1/4°, daily | https://doi.org/10.48670/moi-00145 | constellation sensitivity (Fig. S3) | `$OSR11_DATA/monthly_twosat/` via `00a` |
| 5 | Roemmich–Gilson Argo 1° temperature and salinity, 2004–2025 | https://sio-argo.ucsd.edu/RG_Climatology.html | isopycnal slope (Fig. 4) | `$ARGO_RG_DIR` via `00d` |
| – | GEBCO bathymetry | https://www.gebco.net | context (isobaths, depth along the band) | `outputs/base/bathymetry_so.nc` via `00c` (distributed) |

Terms: the derived outputs in `outputs/` come from the products above and keep their licences and attribution requirements
(Copernicus Marine Service, Zenodo record of Gutierrez-Villanueva et al. 2023, Argo/Roemmich–Gilson, GEBCO). Raw data are not
redistributed. The code is under the MIT licence (`LICENSE`).

## Environment

- Python 3.13; `requirements.txt` is pinned to the versions the figures were rendered and the checks were run with
  (numpy 2.3.5, pandas 3.0.5, xarray 2026.9.0, matplotlib 3.10.9, cartopy 0.26.0, numba 0.68.1, statsmodels 0.15.0). The notebooks
  also need `jupyter` (`nbconvert`, `ipykernel`).
- Environment variables (defaults in brackets; all read in `src/osr11/paths.py`): `OSR11_DATA` [`data/`] raw monthly fields, wind,
  Argo; `OSR11_OUTPUTS` [`outputs/`]; `OSR11_CACHE` [`cache/`] stacked monthly fields built on first use; `OSR11_FIGURES`
  [`rendered/`] default figure folder of the figure scripts; `GV23_DIR` unpacked archive of Gutierrez-Villanueva et al. (2023);
  `ARGO_RG_DIR` converted Roemmich–Gilson files; `NPROC` [2] worker processes.
- Fonts: the figure style asks for Arial and falls back to Liberation Sans (same metrics). Text is rasterised slightly differently
  from one system to another, so a re-rendered figure can differ from `rendered/` by a pixel around the labels; the data marks do not.

## Running the pipeline from the raw data

Raw data (about 30 GB) are not distributed. `make data` downloads them (Copernicus Marine login; GEBCO and Roemmich–Gilson by
hand), `make base` builds the base fields, `make process` runs steps 02–04 and 06 in dependency order (hours), `make stats` the
statistics on the saved series. The order of the steps is the order of the Makefile.

## What was checked for version 2.0.0, and what was not

- **Figures.** The eleven figure scripts were run from `outputs/` in the pinned environment. Compared with the figures of the
  manuscript and the Supplement (`rendered/`), every data mark is identical; the differences are confined to text
  (anti-aliasing and the spacing of super- and subscripts), which depends on the fonts of the system.
- **Steps that run from `outputs/`, re-run here:** `03w` (ADT-bin drift), `04i`, `04j`, `04k` (exponents, their intervals, lagged
  correlations), `04m` (width of the domain), `05l` (SAM lag), `06_product4/p4_05`, `p4_06`, `p4_07`. Each reproduces the
  stored values to rounding error (differences below 1e-13).
- **Notebooks.** Executed from `outputs/`; the values they print were compared with those quoted in the manuscript and the Supplement.
- **Not re-run for this version:** the steps that need the raw data (`00`–`01`, most of `02`–`04`, `06_product4/p4_01`–`p4_04`,
  `p4_08`, `p4_09`) and the long statistics (`02o`, `05a`). Their outputs are the ones produced during the revision
  (8 October–10 November 2026). The scripts were moved from the working area with their paths changed and nothing else; those
  that need raw data were byte-compiled only.
- `outputs/x7/x7_composite_annual_*.npz` (the curves of Figs. 2b and 3b) were reduced from the monthly composites (255 MB, not
  distributed) with the averaging the figure scripts used; `02p` now writes them.
- Published values shown in Fig. S4 were typed from the papers (`scripts/05_stats/05k_previous_studies.py`, column `source`).
- Analyses made during the revision that the final manuscript does not use were removed from this version; they remain in the
  git history.

## Citation and licence

Cite the article and the Zenodo archive (`CITATION.cff`, https://doi.org/10.5281/zenodo.21485764). Code: MIT licence (`LICENSE`).
The derived outputs come from the data products above and keep their licences and attribution requirements (Copernicus Marine
Service, Gutierrez-Villanueva et al. 2023, Argo/Roemmich–Gilson, GEBCO, British Antarctic Survey SAM index). Natural Earth
shapefiles are in the public domain.
