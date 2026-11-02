# x1 — contour-vs-energy test of frontal displacement (exploratory, 29 Oct 2026)

Scripts: `x1_contour_vs_energy.py` (extraction), `x1_stats.py` (tests), `x1_figure.py` (figure). Inputs: `outputs/a2/tracks_baseline.npz`
(the 900 retained half-power tracks of Sect. 3.3), monthly `adt_mean` (Product 1), `outputs/mean_adt.nc`. Test: `rstats.trend_bcz`
(monthly, both noise models), as in the main text. Circumpolar series = mean over tracks of each track's anomaly from its own mean.

| Quantity (900 tracks, 2001–2025) | Trend | 95 % | p |
|---|---|---|---|
| Half-power latitude (energy pathway) | −0.5 km per decade | −4.2 to +3.1 | 0.74 |
| Latitude of the fixed ADT contour through the mean half-power point | −11.5 km per decade | −13.7 to −9.2 | 0.003 |
| ADT at the half-power point | +2.6 cm per decade | +2.1 to +3.2 | 0.003 |
| ACC-band mean ADT (no ice mask) | +3.1 cm per decade | +2.1 to +4.0 | 0.003 |

Consistency: median |dADT/dlat| at the tracks = 2.8 cm per 10 km, so +2.6 cm per decade of sea level corresponds to a −9.5 km per decade
contour shift, close to the −11.5 measured. Per track (OLS): contour trend southward at 85 % of tracks (median −9.8 km per decade),
half-power trend southward at 56 % (median −1.0); correlation of the two per-track trends r = 0.25. Sectors: contour −14.5 / −12.7 / −8.3
km per decade (Atlantic / Indian / Pacific, all p < 0.01); half-power +2.6 / −2.2 / −0.6 (p = 0.41 / 0.43 / 0.84).

Robustness of the half-power result: coverage ≥ 75 % (1301 tracks) −0.8 (−4.1 to +2.5), p = 0.56; ≥ 50 % (1695) −0.9 (−3.9 to +2.0), p = 0.47.
Bathymetric split (std of depth in a ±2° box, median 448 m): rough half −1.0 (−4.7 to +2.7), p = 0.53; smooth half −0.1 (−6.7 to +6.6),
p = 0.99; depth ≥ 3500 m −1.2 (−5.4 to +2.9); < 3500 m (168 tracks) +3.1 (−4.4 to +10.6).

Reading: the height contours of the ACC drift poleward at the rate implied by regional sea-level rise, while the energy pathway does not
move (95 % bound about 4 km per decade, ~10 km over the record), over rough and smooth topography alike.
