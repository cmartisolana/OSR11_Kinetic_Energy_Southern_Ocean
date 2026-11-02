# Observational constraints on eddy saturation in the Antarctic Circumpolar Current

Analysis code and processed outputs for

> Martí-Solana, C., Ruiz, S., Combes, V., Barceló-Llull, B., Ballarotta, M. and Pascual, A.: Observational constraints on
> eddy saturation in the Antarctic Circumpolar Current, State of the Planet (Copernicus Ocean State Report 11, Chapter 2),
> manuscript sp-2026-11, 2026.

Code: Cristina Martí-Solana (IMEDEA, CSIC-UIB), cmarti@imedea.uib-csic.es. Archive: https://doi.org/10.5281/zenodo.21485764
(`CITATION.cff`). Repository version 2.0.0 (2 November 2026; publication version of the revised manuscript; v1.0.0-alpha accompanied the discussion paper). Archived data arrays omitted from GitHub are listed in `archive/OMITTED_DATA_FILES.md`.

## What is reproduced

Satellite altimetry (2001–2025) is used to measure trends of kinetic energy (KE) and eddy kinetic energy (EKE) in the
Antarctic Circumpolar Current (ACC), first in Drake Passage (compared with the LM Gould ADCP record of
Gutierrez-Villanueva et al., 2023), then around the ACC, together with the position of the energetic pathway
(half-power method, energy-based) and of fixed ADT contours, and the response of KE, EKE and the σ₀ = 27.3 kg m⁻³
isopycnal slope (Roemmich–Gilson Argo) to the wind. The repository contains:

- every script of the final pipeline, from the data download to the figures (`scripts/`, `figures/`, shared code in `src/osr11/`);
- the small processed outputs (72 MB) from which **all main and SI figures and all quoted numbers** are regenerated without
  the raw data (`outputs/`, listed with checksums in `outputs/MANIFEST.csv`);
- the figures as rendered by these scripts (`rendered/png`, `rendered/pdf`; identical, pixel for pixel, to the submitted figures);
- the SI text and tables with a table → output → script map (`si/`);
- everything superseded, unchanged, in `archive/` (nothing was deleted; see `archive/README.md`).

Raw data (about 30 GB) are not distributed; they are downloaded by `scripts/00_extract`.

## Layout

```
src/osr11/           shared code: paths.py (all locations, env vars), params.py, rstats.py (trend tests),
                     fronts.py (half-power detection), halfpower.py, drake_halfpower.py, style.py (figure style),
                     fs_utils/ (first-submission helpers used only by 01c)
scripts/00_extract/  downloads and conversions (Copernicus Marine login; GEBCO; Roemmich-Gilson)
scripts/01_base/     base fields: ice mask, mean ADT, monthly KE/EKE fields, 2001-2025 climatology
scripts/02_circumpolar/  circumpolar KE/EKE, half-power pathway, maps, contour-vs-energy, earlier studies, wake test
scripts/03_drake/    Drake Passage: Gould track, in situ comparison, section diagnostics, Qs, half-power tests
scripts/04_wind_argo/    isopycnal slope (Argo), wind stress, common-cell diagnostics, wind sensitivity
scripts/05_stats/    statistics on the saved series: recalibrated tests, summaries, size checks, derived tables
scripts/make_manifest.py  rewrites outputs/MANIFEST.csv (sizes, SHA-256, producing step)
figures/             one script per figure (Figs. 1-4, S1-S7) -> rendered/
outputs/             processed outputs read by the figures and quoted in the text + MANIFEST.csv
rendered/            Figures 1-4 and S1-S7 (png 600 dpi, pdf)
si/                  SI text (SI.md), SI tables (tables/*.csv), README.md (table -> output -> script)
docs/script_map.md   old (working repository) -> new script names
archive/             superseded material, unchanged (working tree of 2 Nov 2026)
VERIFICATION.md      what was run here, timings, checksums, comparison with the manuscript numbers
```

Every script starts with a docstring giving its purpose, inputs, outputs, run command and random seeds. Scripts find the
package themselves (they add `src/` to `sys.path`), so `python scripts/.../xx.py` works from any directory; alternatively
`pip install -e .`.

## Data (Table 1 of the paper)

| Ref. | Product | Access | Used for | Placed in |
|---|---|---|---|---|
| 1 | SEALEVEL_GLO_PHY_L4_MY_008_047, all-satellite DUACS DT2024, 1/8°, daily | https://doi.org/10.48670/moi-00148 | all altimetric diagnostics | `$OSR11_DATA/monthly_allsat/` via `00a` |
| 2 | WIND_GLO_PHY_L4_MY_012_006, monthly | https://doi.org/10.48670/moi-00185 | wind stress (Fig. 4) | `$OSR11_DATA/wind_taux.nc` via `00b` |
| 3 | Drake Passage ADCP velocity, EKE and transport, LM Gould (Gutierrez-Villanueva et al., 2023) | https://doi.org/10.5281/zenodo.10044261 | in situ comparison (SI S2) | unpacked to `$GV23_DIR` |
| 4 | SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057, two-satellite, 1/4°, daily | https://doi.org/10.48670/moi-00145 | constellation sensitivity (SI S1) | `$OSR11_DATA/monthly_twosat/` via `00a` |
| 5 | Roemmich–Gilson Argo 1° temperature and salinity, 2004–2025 | https://sio-argo.ucsd.edu/RG_Climatology.html | isopycnal slope | `$ARGO_RG_DIR` via `00d` |
| – | GEBCO bathymetry | https://www.gebco.net | context (isobaths, depth along the band) | `outputs/base/bathymetry_so.nc` via `00c` (distributed) |

Terms: the derived outputs in `outputs/` come from the products above and keep their licences and attribution requirements
(Copernicus Marine Service, Zenodo record of Gutierrez-Villanueva et al. 2023, Argo/Roemmich–Gilson, GEBCO). Raw data are not
redistributed. The code is under the MIT licence (`LICENSE`).

## Environment

- Python 3.13 (verified here: Linux, Python 3.13.18). The analyses were run on the author's Mac with Python 3.10–3.13,
  numpy 2.x, pandas 2.x.
- `pip install -r requirements.txt` (pinned to the versions installed where this repository was verified) or
  `conda env create -f environment.yml`. Extraction only: `pip install -r requirements-extract.txt` (copernicusmarine;
  optional pymannkendall). cartopy needs the Natural Earth 110m/50m land shapefiles (downloaded automatically on first use).

Environment variables (defaults in brackets; all read in `src/osr11/paths.py`):

| Variable | Meaning |
|---|---|
| `OSR11_DATA` [`data/`] | raw monthly fields, wind, Argo, first-submission KE/EKE fields, `preprocessing.nc` |
| `OSR11_OUTPUTS` [`outputs/`] | processed outputs (read and written by every step) |
| `OSR11_CACHE` [`cache/`] | 1.6 GB stack of the monthly KE/EKE fields, built on first use by `osr11.halfpower` |
| `OSR11_FIGURES` [`rendered/`] | figure output folder |
| `GV23_DIR` [`$OSR11_DATA/TrendsDrakePassageTransp`] | unpacked Gutierrez-Villanueva et al. (2023) archive |
| `ARGO_RG_DIR` [`$OSR11_DATA/argo_rg`] | Roemmich–Gilson files converted by `00d` (+ `RG_ArgoClim_20*.nc` monthly extensions) |
| `ALTIMETRY_DIR` | local daily altimetry tree, only for `01a` when not streaming |
| `NPROC` [2] | worker processes for the map calibrations (`02f`, `02g`, `04f`) |
| `RETEST_N_OUTER` [300] | null series of the calibrated test in `05a` |
| `RSTATS_MONTHLY_ORDER` [`both`] | AR order(s) of the monthly test (`both` = larger p of AIC and AR(12); used for every reported trend) |
| `RSTATS_QUICK` | smoke tests only (5 outer, 50 bootstrap series) |
| `M1_Y0` [2001] | start year of `02l` |
| `DK2_THR`, `DK2_LATMIN` | sensitivity runs of `03q` |

## How to run

```
make figures      # Figs. 1-4, S1-S7 from outputs/ (Fig. 1 also reads $OSR11_DATA/monthly_allsat/2001/01.nc)   ~1 min
make stats        # statistics on the saved series (no raw data)                                               ~2 min
make si           # SI figures + derived SI tables (no raw data)
make check        # byte-compile everything outside archive/
# with the raw data:
make data         # 00: downloads (Copernicus Marine login); GEBCO and RG by hand (see the echo)
make base         # 01: ice mask, mean ADT, monthly KE/EKE fields, climatology
make process      # 02-04 in dependency order (drake, then circumpolar, then windargo); hours
make sizecheck    # SI S4 false-positive checks (synthetic; hours)
```

Which targets need what:

| Target | Raw data | Runs from the distributed `outputs/` |
|---|---|---|
| `figures`, `si` | only Fig. 1 (`monthly_allsat/2001/01.nc`; the 2001–2025 profile is cached in `outputs/fig1/`) | yes |
| `stats` | no | yes (`05a` skips rows already present; the A2 family also needs `$OSR11_CACHE`) |
| `sizecheck` | no | yes |
| `data`, `base`, `process` | yes (and `$GV23_DIR`, `$ARGO_RG_DIR`) | no |

Steps of `process` that run without raw data, from distributed intermediate files: `02k`, `03e`, `03f`, `03g`, `04c`,
`04f`, `04h` and all of `05`.

## Where each figure, table and number comes from

| Item | Script | Main inputs (`outputs/…`) | Producing steps |
|---|---|---|---|
| Fig. 1 half-power method | `figures/fig1_halfpower_method.py` | `fig1/clim_profile_147.9.npy`, `base/*`; one monthly file | `00a`, `01a`; profile recomputed if absent |
| Fig. 2 Drake Passage | `figures/fig2_drake.py` | `drake/dk34`, `drake/dk10`, `drake/dk8`, `drake/dk5`, `drake/dk24`, `d1/d1b_track.csv`, `x1/x1_drake_period_shift.json` | `03k`, `03o`, `03i`, `03b`+`03c`, `03d`, `03a`, `03t` |
| Fig. 3 circumpolar | `figures/fig3_circumpolar.py` | `a1/a1_maps_bc.nc`, `drake/dk14`, `drake/dk12`, `a5/a5_series_allsat.csv`, `retest/retest_summary.csv` | `02g`, `02h`, `02i`, `02a`, `05a`/`05b` |
| Fig. 4 wind and slope | `figures/fig4_wind_slope.py` | `f3w/f3w_wind_map.nc`, `f4/f4_tilt_maps.nc`, `f4m/f4m_map.nc`, `f4e/*`, `f4s/*` | `04b`, `04c`, `04f`, `04g`, `04h` |
| Fig. S1 products | `figures/figS1_products.py` | `a5/*`, `retest/retest_summary.csv`, `m1/m1_trends.csv` | `02a`, `05b`, `02l` |
| Fig. S2 Drake periods | `figures/figS2_drake_periods.py` | `drake/dk12`, `drake/dk1`, `drake/dk8/dk8_along_acc_fronts_periods_allsat.csv` | `02i`, `03r`, `03j` |
| Fig. S3 velocity trend | `figures/figS3_velocity_transect.py` | `drake/dk5/dk5_trends.csv` | `03b`, `03c` |
| Fig. S4 in situ validation | `figures/figS4_drake_validation.py` | `drake/dk3`, `dk5`, `dk24`, `dk25`, `dk34`, `dk39` | `03h`, `03b`, `03d`, `03e`, `03k`, `03f` |
| Figs. S5, S6 half-power sensitivity, products | `figures/figS5_S6_halfpower.py` | `retest/retest_summary.csv`, `a3/a3_summary.csv`, `a2/a2h_twosat_vs_allsat.csv` | `05b`, `02d`+`05c`, `02e` |
| Fig. S7 contour vs energy | `figures/figS7_contour_vs_energy.py` | `x1/x1_tracks.npz`, `x1/x1_band_adt.npy`, `x1/x1_circumpolar_trends.csv` | `02j`, `02k` |
| Table 1 | – (data products) | – | – |
| SI Tables S1–S5 and the Drake-domain table | see `si/README.md` | | |
| Sect. 2.2: 900 retained longitudes | – | `a2/tracks_baseline.npz` | `02c` |
| Sect. 2.3, SI S2: Qs, no detectable change | – | `drake/dk40/dk40_net_trends.csv` | `03g` |
| Sect. 2.5, SI S4: false-positive rates | – | `stats/size_summary.csv` (from the logs in `stats/size_bc/`) | `05f`, `05g`, `05h` |
| Sect. 3.1: 21–24 % / 12 % late-vs-early EKE | – | `drake/dk8/dk8_period_change_regions.csv` | `05d` (from `03j`) |
| Sect. 3.1: half-power vs contour, 70–55°W; −3 / −19 km | – | `x1/x1_drake_trends.csv`, `x1/x1_drake_period_shift.json` | `03t` |
| Sect. 3.1: one third of 25 km shifts detected | – | `drake/dk6/dk6_detect_S25_seed25.csv` | `03n` |
| Sect. 3.1: streamline bins +3.1, +5.2 cm s⁻¹ per decade | – | `drake/dk5/dk5_trends.csv` | `03b`, `03c` |
| Sect. 3.2: KE/EKE trends, sectors, EKE−KE contrast, wake vs rest | – | `retest/retest_summary.csv` | `02a`, `02m`, `05a`, `05b` |
| Sect. 3.2: 91 % positive bins, 103, 70–84, 40–63 | – | `drake/dk12/dk12_band_trends_allsat.csv` | `02i` |
| Sect. 3.2, 4: pathway trend, −11.5 km per decade contour, rough/smooth split | – | `x1/x1_circumpolar_trends.csv`, `x1/x1_bathy_split.csv`, `x1/x1_implied_shift.csv` | `02k` |
| Sect. 3.3: 8.9 %, 1.2 %, r = 0.61/0.57, coefficients 0.91/0.49/0.09 | – | `f4e/f4e_results.csv`, `f4s/f4s_results.csv` | `04g`, `04h` |
| Sect. 3.3: 13 field-significant slope cells, all positive | – | `f4m/f4m_map.nc` | `04f` |

The full list of main-text numbers, with the value recomputed from the distributed outputs, is in `VERIFICATION.md`.

## Known limitations (please read)

- **Raw data and long steps were not re-run for this version.** The distributed outputs were produced on the author's Mac
  between 8 and 30 October 2026 by the scripts as they stood then. For the publication repository, the scripts were moved,
  their paths routed through `osr11.paths`, and a few were repaired (below). In the verification container only the
  all-satellite monthly files were available: steps needing only them and the distributed outputs were re-run and reproduce
  the distributed files (`VERIFICATION.md`); steps needing the two-satellite product, the wind, the Argo fields, the in situ
  archive or the 1.6 GB KE/EKE cache were checked by byte-compilation only.
- **Repairs made here** (each documented in the script docstring): `x1_drake.py` did not compile and only printed the
  period displacement → `03t` (writes `outputs/x1/x1_drake_period_shift.json`, reproduces the stored values exactly);
  `extract_drake_cols.py` deleted every raw monthly file it read → replaced by direct reading in `03t`, original kept only in
  `archive/`; the `dk6_common` import of dk9/dk10 (broken by the 26 Oct path fix) → `osr11.drake_halfpower`; dk6/dk8/dk9 outputs
  and their readers now share `outputs/drake/dkN/`; the hash()-based seed of the size check → `zlib.crc32` (`05g`).
- **SI S4 false-positive rates may shift slightly on a re-run.** The shipped logs were produced with a seed taken from Python's
  per-process randomised string hash, which cannot be reproduced; `05g` now uses a stable seed, so a re-run draws different
  synthetic series (Monte Carlo error about ±1.3 percentage points at 5 % with 300 replicates).
- **Generators written for this version** (no script existed; each reproduces the stored file where one existed):
  `01b` (monthly KE/EKE fields; checked on 2001-01), `01d` (`clim_allsat.nc`), `02j`/`02k` additions (`x1_band_adt.npy`,
  bathymetric split, implied contour shift), `05c` (`a3_summary.csv`), `05d` (Sect. 3.1 percentages), `05h` (size table),
  `05i` (common-cell count).
- **`01a` and `01c` are the first-submission generators of record** (notebook 00 exported to a script; `run_all_two_sat.py`),
  with paths changed only; they were not re-run.
- **Fig. S2 was produced with a new script** (`figS2_drake_periods.py`) from existing outputs; the earlier SI Fig. S2 had no
  script. All figures were re-rendered here and match the figures delivered for the manuscript pixel for pixel; the fonts are
  Liberation Sans in this environment (metric-compatible with Arial, used on the Mac).
- **Fig. S1d** contains values of earlier studies typed from the papers (Hogg et al. 2015; Martínez-Moreno et al. 2021),
  marked "[AUTHOR CHECK]" in the script.
- **2830 common cells (SI S5)**: `05i` reproduces 2835 slope cells with KE/EKE in every month; the 5-cell reduction comes from the
  wind condition, which needs the wind file. `04g` now writes the count when it is re-run.
- The 2005–2019 Qs trend (+2.3 % per decade, p = 0.040) is in `dk40_net_trends.csv` but not quoted.
- The first-submission notebooks, the tracked-eddy analysis, the transport proxies, the Liu et al. reproduction and the
  retest families A6/B1/LIU are archived and are not needed for any result of the paper.

## Citation and licence

Cite the article and the Zenodo archive (`CITATION.cff`, https://doi.org/10.5281/zenodo.21485764). Code: MIT licence
(`LICENSE`). Data products: see "Data".
