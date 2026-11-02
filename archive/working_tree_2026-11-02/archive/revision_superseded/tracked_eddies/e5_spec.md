# E5 — track-following eddies at the five standing meanders (specification)

Written 26 Oct 2026 before any E5 result was computed. Atlas: META4.0 DT all-sat subset (`revision/data/aviso_eddies/subset/`, south of 35°S, 2001-01 to 2023-09-08), tracked observations only (track > 0); both polarities, reported together and separately. Regions: the dk13 longitude ranges (upstream / ridge / downstream) of the five meanders; a track is assigned to a feature and a region by the position of its **first observation** (`observation_number = 0`, i.e. a true birth inside the subset); births in the ACC band (mean ADT −0.6 to +0.2 m, ice-free in every month, as E1–E4) outside all fifteen ranges form the control class "away". Tracks that enter the subset already alive are excluded from cohort statistics.

## Questions (fixed)

- **Q1 (descriptive):** how amplitude, speed radius and rotational speed of eddies evolve downstream of each ridge: composite of all ACC-band observations binned by longitude relative to the eastern edge of the ridge range (−10° to +25°, 1° bins), median and interquartile range, with the observation count; and amplitude against days since birth (0–120 d) for eddies born upstream/over the ridge versus "away" (survivor-weighted: stated as such).
- **Q2 (descriptive):** lifetime (days), propagation distance (along-track, km) and net zonal displacement of eddies by birth class (upstream, ridge, downstream, away). Births after 2022-08-31 are excluded so that every cohort has at least one year of possible lifetime before the end of the atlas; the few tracks still alive at the end are reported as censored.
- **Q3 (the only test):** has the character of the eddies born over and downstream of the features changed relative to the rest of the band? Annual birth-year series 2001–2022 of (i) mean lifetime, (ii) mean maximum amplitude over the track, (iii) mean propagation distance, for births in ridge + downstream of the five features pooled ("MEANDER") and for "away"; statistic = trend of MEANDER/mean − AWAY/mean (points per decade), as in A7 and E3; test = `rstats.trend_bcz` on annual series (AR(1), iterated bias correction); the three contrasts are one family (BH-FDR 10 %). Per-feature versions are secondary and descriptive. Period split 2001–2011 vs 2012–2022 shown for Q1 as a descriptive check only.

## What would count

Under a wake/standing-meander mechanism that has strengthened, eddies born over or downstream of the features would be larger-amplitude, longer-lived or farther-travelling than elsewhere, and more so in the later period. Equal properties and no contrast trend would mean the atlas does not distinguish the wakes from the rest of the band, as E3 found for counts and amplitude.

## Caveats stated in advance

The atlas is built from the same all-satellite maps as the EKE (constellation sensitivity, SI S1); only eddies larger than the mapping resolution are tracked (median speed radius 36–38 km); amplitude is relative to the local contour, not to a fixed reference; cohort statistics are survivor-dependent; the record ends in 2023-08, two years before the altimetry record; and the five ranges were fixed from bathymetry (dk13) but their 5° bin trends had been seen in B3 before dk13 was written.
