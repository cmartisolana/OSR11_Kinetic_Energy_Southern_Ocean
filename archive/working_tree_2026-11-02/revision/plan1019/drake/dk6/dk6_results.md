# Drake Passage half-power tests (dk6), 19 Oct 2026

Region: 70–55°W; envelopes whose 2001–2025 climatological half-power lies south of 52°S (Malvinas-Current envelopes excluded) and in the ACC band (mean ADT −0.6 to +0.2 m). Baseline half-power (200 cm² s⁻², 25 %, 100 km) on monthly KE of the monthly-mean velocity, recomputed from the archived monthly files; monthly envelopes matched to the climatological envelope at each longitude (as a2_run.py). Trends: Theil–Sen with the calibrated test (trend_bcz), per decade, 95 % intervals. Main selection: envelopes present in ≥ 50 % of months (Product 1: 106 longitudes; Product 4: 56). Ice mask (ice fraction > 0.10 in that month) changes no result by more than 0.2 km per decade.

**Caveat on reproduction.** KE recomputed from the archived monthly files reproduces the stored circumpolar half-power tracks (outputs/a2/tracks_baseline.npz) in 81 % of track-months in the passage (96 % median over tracks in 75–50°W); masking ice does not explain the rest, and the stored KE fields (outputs/monthly_half-power_points) are not in the staged repository. With the stored tracks, the passage mean is −3.1 km per decade (−11.3 to +5.1; 101 longitudes; monthly SD 23 km) against −6.8 (−21.9 to +8.2; monthly SD 34 km) recomputed.

| Test | Product 1 | Product 4 |
|---|---|---|
| T1 passage-mean half-power latitude (km per decade, +N) | −6.8 (−21.9 to +8.2), p = 0.30 | −6.4 (−23.1 to +10.3), p = 0.33 |
| T1, envelopes in ≥ 90 % of months | −5.7 (−22.2 to +10.7), p = 0.39 (70 lon.) | −5.2 (−21.4 to +11.1), p = 0.45 (34) |
| T1 detectability, 30 realisations | 25 km over 25 yr detected in 33 %, 10 km in 27 %, none (false positives) 3 %; median recovered 22 of 25 km | – |
| T2 latitude of the climatological half-power streamline (relative ADT; km per decade) | −0.1 (−4.5 to +4.3), p = 0.93 | +0.8 (−4.8 to +6.4), p = 0.73 |
| T2 half-power across streamlines (km-equivalent per decade) | −6.2 (−20.1 to +7.7), p = 0.27 | −8.0 (−22.9 to +6.9), p = 0.21 |
| T2 monthly r(half-power latitude, streamline latitude) | 0.52 (27 % of variance) | 0.45 |
| T3 EKE-weighted median minus KE half-power (km per decade) | +0.1 (−1.3 to +1.4), p = 0.92 | +0.2 (−1.3 to +1.6), p = 0.80 |
| T4 KE fraction north of the climatological half-power (pp per decade) | −1.0 (−3.2 to +1.1), p = 0.33 | −0.8 (−3.0 to +1.4), p = 0.39 |
| T4 EKE fraction north (pp per decade) | −0.3 (−3.9 to +3.3), p = 0.80 | −0.2 (−2.9 to +2.5), p = 0.82 |
| T4 synthetic (no noise): 20 % north-flank intensification / 25 km northward shift | KE +4.3 / +6.4 pp; EKE +4.0 / +4.6 pp | KE +4.2 / +5.1 pp |
| T5 half-power vs wind stress, annual, lags 0–3 yr (Drake box, upstream 100–80°W, 45–60°S) | no relation: |r| ≤ 0.35, all p ≥ 0.10 (phase-randomised); wind-expected shift −2.2 to +0.8 km per decade (95 % ranges within −8.6 to +6.0) | |r| ≤ 0.42, all p ≥ 0.057 |

Reading:
- T1: the passage mean narrows the single-section interval (−10, −36 to +16) by about 40 %, but a 25 km shift over 25 years is detected in only a third of realisations; the passage bound is weak, unlike the circumpolar one (25 km detected in 100 %).
- T2: the sea-surface-height streamline at the energy median has not moved (±4.5 km per decade), a tighter constraint than the half-power latitude; about a quarter of the monthly half-power variability moves with the streamlines, the rest is energy redistributed across them. No trend in either part.
- T3: eddy energy has not migrated relative to the energetic band (±1.4 km per decade). EKE and KE share the monthly velocity anomaly, so this separation is partly constrained by construction.
- T4: no change in the north–south partition of KE or EKE within the merged band; a 20 % north-flank KE intensification over 2001–2025 (+1.7 pp per decade if linear) lies just outside the KE interval, a south-flank one does not; for EKE neither is excluded. Not consistent with a surface EKE decrease in the PF relative to the SAF (Gutierrez-Villanueva et al., 2023).
- T5: the passage position does not respond to wind stress at 0–3 year lags; the shift expected from the observed wind trends (< 9 km per decade in magnitude) is below what the passage can detect, so the absence of a trend does not test a wind-driven displacement here.
