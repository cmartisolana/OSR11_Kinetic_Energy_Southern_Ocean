# Changelog

## v2.0.0 (2026-11-17) — final manuscript and Supplement

Version of the revised, final manuscript "Observational constraints on eddy saturation in the Antarctic Circumpolar Current"
and of its Supplement. The repository was rebuilt around the final analysis:

- `figures/`: one script per figure (Figs. 1–4, S1–S7); `rendered/`: the figures of the manuscript and the Supplement.
- `scripts/`: the pipeline from the data download to the processed outputs (00 extraction, 01 base fields, 02 circumpolar energy
  and energetic pathway, 03 Drake Passage, 04 wind and isopycnal slope, 05 statistics, 06 two-satellite product).
- `outputs/`: the processed outputs the figures and notebooks read, with `MANIFEST.csv` (producing script, checksums, figures
  that read each file). Figures and notebooks run without the raw data.
- `notebooks/`: four executed notebooks, one per part of the paper.
- `src/osr11/`: shared code (paths, parameters, trend tests, half-power detection, figure style); `tests/`.

Changes of method since the discussion paper (1.0.0-alpha): 2001–2025 all-satellite record with the two-satellite product as a
sensitivity test; energetic pathway tracked with a longitude-dependent half-power threshold and compared with fixed ADT
contours; Drake Passage section compared with the shipboard ADCP record of Gutierrez-Villanueva et al. (2023); isopycnal slope
from the Roemmich–Gilson Argo fields; lagged correlations and power-law exponents with respect to the wind stress; trend tests
with autoregressive surrogates.

A tag `v2.0.0` made on 2 November 2026 for an intermediate state of the revision was replaced by this release.

## v1.0.0-alpha (2026-07-22)

Code of the discussion paper (sp-2026-11). Archived at https://doi.org/10.5281/zenodo.21485764.
