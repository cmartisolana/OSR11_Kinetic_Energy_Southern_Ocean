# Revision plan of 19 Oct: task list mapped onto `revision/`

Status: **done** = run in this session (scripts in `revision/plan1019/`); **ready** = script written, waiting for data that needs your login; **decision** = needs the authors.

| # | Plan item | Script | Outputs | Status |
|---|---|---|---|---|
| 1 | §3 Drake: partial-section label, coverage-matched in situ comparison | `d1_drake_coverage.py`, `d1b_gvtrack.py` | `outputs/d1/*` | done; independently re-coded and reproduced |
| 2 | §3 In situ variable definition | — | text | done: series = TOTAL 0–760 m |
| 3 | §2 Matched comparison: period, product, scale; Hogg et al. (2015) boxes | `m1_comparison.py` | `outputs/m1/*` | done |
| 3b | §2 MM21 procedure applied to our data | `m2_mm21_procedure.py` | `outputs/m2/*` | done |
| 4 | §2 Full 1993–2025 comparison | `run_extractions.sh` step 1, then `M1_Y0=1993 python m1_comparison.py …` | — | ready (Copernicus Marine login) |
| 5 | §2 Earlier product version (two-sat vDT2018, C3S; available, checked via the CDS API) | `extract_c3s_twosat.py` (`run_extractions.sh` step 3) | — | ready (CDS account) |
| 6 | §2 MM21 definitions | read from their code repository | results §0, §2b | done |
| 7 | §4 Half-power gates G1/G2 | `h1_halfpower_meanders.py … gate` | `outputs/h1/h1_gate.csv`, `h1_synthetic.csv` | done |
| 8 | §4 Trends and relations (SEIR150 only) | `h1_halfpower_meanders.py … analyse` | `outputs/h1/h1_results_SEIR150.csv`, `h1_relations_fdr.csv` | done |
| 9 | §3 Antarctic-slope sea-level index | `run_extractions.sh` step 2, `h0_antarctic_slope_index.py` | `outputs/h0/*` | ready (login, ETOPO/GEBCO file) |
| 10 | Conversion to Sv with Hughes/Hogg model factors | — | — | decision |
| 11 | AC replies | `claude/AC_reply_*_v2.md` (project) | — | drafted; two [AUTHOR DECISION] brackets in C1 |
| 12 | Optional (later paper): frontal compensation, moving-path energy decomposition | — | — | not started (by design) |
