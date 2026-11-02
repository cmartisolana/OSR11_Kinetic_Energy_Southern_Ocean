# Notebooks

Each notebook is self-contained and documented; publication figures are written
to `figures/manuscript/` or `figures/supplementary/` (all other diagnostics go
to the untracked `outputs/` folder).

## Run order

1. `00_preprocessing.ipynb` — sea-ice mask and preprocessing grids
2. `01_front_detection.ipynb` — front detection; saves the *front_trends* detection-example panels
3. `02_trend_analysis.ipynb` — regional trend tables
4. `03a_ke_eke_trends.ipynb` — main article figure *KE_EKE_trends*
5. `03b_fronts_envelope_trends.ipynb` — *front_trends* timeseries/bars panels
6. `05_wind_stress_trends.ipynb` — wind stress CSVs (required by 06 and 07)
7. `06_drake_passage_transport.ipynb` — transport CSVs (required by 07)
8. `07_drake_passage_timeseries.ipynb` — main article figure *drake_passage* (both KE+EKE and EKE-only variants via the `INCLUDE_KE` flag)

## Supplementary

Run the prerequisite script listed in the main README first.

- `supplementary/s01_ke_eke_trends_since2010.ipynb` — *KE_EKE_trends_since2010*
- `supplementary/s02_ke_eke_trends_since2016.ipynb` — *KE_EKE_trends_since2016*
- `supplementary/s03_ke_eke_trends_two_sat.ipynb` — *KE_EKE_trends_twosat*

## Trend policy

All trend slopes are computed with **Theil–Sen**. Significance is reported with
**Modified Mann–Kendall (Yue–Wang)**. No OLS/weighted-OLS trend outputs are part
of the workflow.
