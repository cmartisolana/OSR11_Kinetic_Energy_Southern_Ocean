# E5: track-following eddies at the five standing meanders — descriptive pass (26 Oct 2026)

Specification: `e5_spec.md` (written before any result). Script: `e5_tracks.py`; quick-look: `e5_quicklook.py`; outputs in `outputs/eddies/e5_*`. Atlas META4.0 DT all-sat subset, tracked observations 2001-01 to 2023-09-08 (16.6 M observations, 524,857 tracks). Regions: dk13 longitude ranges; classes by the position of the first observation inside the ACC band (88,849 births; cohort with births up to 2022-08-31: 84,861, of which 2 are censored at the end of the atlas). All tests: `rstats.trend_bcz`, annual, AR(1).

## Q3 (the pre-specified test): no change in the contrast between meander-born and other eddies

Birth-year series 2001–2022, MEANDER (ridge + downstream, five features pooled, ~1,630 births per year) minus AWAY (~1,900 per year), each normalised by its mean; points per decade [95 %]:

| Contrast | Trend | p |
|---|---|---|
| mean lifetime | −0.6 [−5.0, +3.9] | 0.70 |
| mean maximum amplitude | +0.4 [−5.0, +5.9] | 0.81 |
| mean propagation distance | −22 [−55, +11] | 0.14 |

None passes; nothing to correct for. Levels (secondary): MEANDER maximum amplitude +2.0 % per decade (p = 0.17), AWAY +0.3 % (p = 0.85); lifetimes −1.2 / −1.3 % (p 0.20 / 0.30); births −1.5 / −1.3 % (n.s.); speed radius flat. The only p < 0.05 is the MEANDER mean propagation distance (−13.5 % per decade, p = 0.037): it is a mean dominated by a few long-travelling eddies (the series spikes to 1,700 km in single years), one of 10 secondary level tests, and is not interpreted. Per feature (secondary), no contrast with AWAY is significant except the Macquarie Ridge distance (p = 0.037, same caveat).

## Q2 (descriptive): meander-born eddies are larger but not longer-lived

Cohort statistics (both polarities; medians):

| Birth class | n | max amplitude (cm) | lifetime (d) | fraction ≥ 90 d | distance (km) | eastward fraction |
|---|---|---|---|---|---|---|
| upstream | 8,531 | 3.9 | 17 | 0.078 | 61 | 0.65 |
| ridge | 15,998 | 4.8 | 14 | 0.078 | 60 | 0.74 |
| downstream | 19,230 | 6.0 | 13 | 0.068 | 69 | 0.64 |
| away | 41,102 | 3.9 | 16 | 0.082 | 61 | 0.67 |

Eddies born over and downstream of the features reach 25–55 % larger amplitude than eddies born elsewhere in the band and have ~40 % higher rotational speed (0.16–0.17 vs 0.12 m s⁻¹), with the same speed radius (31–33 km, i.e. at the mapping limit). They are not longer-lived: the median lifetime is shorter (13–14 d vs 16 d) and the fraction surviving 90 days is the same or lower. Cyclones are slightly longer-lived and larger than anticyclones in every class. The amplitude-against-age composite (survivors only) shows the downstream-born eddies ahead of the others from the first weeks (5.9 vs 3.8 cm at 14 d; 9.1 vs 6.2 cm at 90 d).

## Q1 (descriptive): the amplitude maximum sits 5–15° east of each ridge and did not move

Median amplitude of all band observations against longitude relative to the eastern edge of the ridge range: Kerguelen peaks at +3 to +9° (10–11 cm vs 2–3 cm upstream of the plateau); Southeast Indian Ridge at +9 to +12° (16 cm vs 4 cm over the ridge); Macquarie at +9° (11 cm); Drake–Scotia at +6 to +9° (8–10 cm). At the Pacific–Antarctic Ridge the maximum lies inside the broad ridge range (−9 to −3°, 7–8 cm) and decays downstream. The 2012–2023 minus 2001–2011 difference of the medians is small and positive at nearly every bin (typically +0.2 to +1.5 cm, i.e. 5–15 %), with no concentration at the downstream maxima.

## Reading

The atlas reproduces the expected climatology — eddies downstream of the five features are more intense than elsewhere — but gives no sign that this has changed over 2001–2022: their amplitude, lifetime and propagation have not diverged from those of eddies born elsewhere in the band. Together with E3 and A7, the wake eddies have intensified in step with the rest of the ACC, not preferentially. This is a null for the "strengthened wake" reading of the standing-meander mechanism at the scales the atlas resolves; it does not test the mechanism itself.

Caveats (as pre-specified): same maps as the EKE (constellation sensitivity); only eddies above the mapping resolution; short median lifetimes mean many tracks are brief detections; cohort composites are survivor-dependent; record ends 2023-08.

Decision for v2.1: no main-text figure from E5 (consistent with the decision not to use tracked eddies in the main manuscript); the Q2 table and one quick-look panel could support one sentence in SI S8 if wanted.
