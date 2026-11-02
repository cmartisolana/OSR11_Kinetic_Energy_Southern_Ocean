# Cover note: manuscript v3.1 (main text)

File: `/home/claude/work/v3/manuscript_v3.1.md` (v3 is unchanged at `manuscript_v3.md`). The file covers the main text only. It applies `revision_decisions_v3.md` sections A, B and C, using the replacement texts in `review_PO.md`, `review_HI.md` and `numbercheck_v3.md`. Section D (SI fixes) is left to the SI revision; the main-text changes that the SI has to follow are listed in Sect. 4 below.

## 1 Word counts

Counted with `wc -w` on the body from "1 Introduction" to the end of "4 Discussion and conclusions". Headings, the five figure captions, Tables 1–2, the abstract and the back matter are excluded.

| Item | v3 | v3.1 |
|---|---|---|
| Body, including bracketed notes ([PENDING], [AUTHOR CHECK], [CHECK], [verify]) | 4956 | **4980** |
| Body, notes removed | 4878 | 4920 |
| Sect. 1 Introduction | 528 | 599 |
| Sect. 2 Data and methods | 1537 | 1442 |
| Sect. 3 Results | 2106 | 2142 |
| Sect. 4 Discussion and conclusions | 794 | 805 |
| Abstract (label excluded) | 221 (218 without [PENDING]) | **204 (200 without [PENDING])** |
| Figure captions (not counted) | 633 | 723 |

Per-section counts tokenise slightly differently from the single `wc`.

- **Within the 5000 limit, but 80 words over the 4900 target.** The decisions added about 400 words: the energy-equivalent comparison, the four-reason conclusion, the forcing/attribution sentence, the revised gap statement, the contiguous-slope definition, the KE-weighted-median caveat, the Monte Carlo and FDR sentences, and the in situ EKE definition. The Methods cuts (HI C1–C11) and the trims listed in Sect. 3 offset most of this.
- **Room for the slope index.** Replacing the [PENDING] block in Sect. 3.4 with results will add words. If that pushes the body past 5000, these cuts (about 80 words) lose the least:
  - Sect. 3.6: the eddy census sentence ("About 345 tracked eddies…"), 16 words;
  - Sect. 1 ¶1: the projections clause, 18 words;
  - Sect. 3.5: "With different periods, sampling and reference means…", 12 words;
  - Sect. 2.6: the circular-shift sentence, moved to SI S9, 10 words;
  - Sect. 3.2: "reflecting monthly fields compared with transects lasting about two days", 9 words;
  - Sect. 2.2: "the upstream box lies in the Pacific sector", 8 words;
  - Sect. 4 ¶3: the Donohue clause, which repeats 2.5(b), 8 words.
- **Abstract.** The upstream clause was dropped for space, as decision B allows; the result stays in Sect. 3.1. The abstract no longer gives the 95 % interval of the ACC EKE trend.

## 2 Referee points mapped to sections (updated for v3.1)

### Referee C (Gutiérrez-Villanueva and Tedesco)

| # | Referee C point | Where addressed in v3.1 | Status |
|---|---|---|---|
| C-a | No novelty relative to Martínez-Moreno et al. (2021), which is not cited | Sect. 1 ¶4–5, Sect. 3.5, Sect. 4 ¶5, Fig. 4b, SI S6. **v3.1:** the gap statement is rewritten. It now credits Gutierrez-Villanueva et al. (2023) for the in situ EKE increase consistent with altimetric trends, Hogg et al. (2015) for pairing EKE with the slope index, and Chambers (2018) and Sokolov and Rintoul (2009b) for the pathway. The claimed gap is the matched-coverage comparison of variability, plus estimates made separately with different periods, products and tests. Novelty point (ii) is reworded ("with their procedure and with an autocorrelation-calibrated test, in two altimetry products"). The OSR sentence now follows (iv). | Answered |
| C-b | The half-power method and the Drake comparison are the weakest parts | Sects. 2.3, 3.2, 2.4/3.3 and Table 2. **v3.1:** "energy median" replaces "centroid". As a KE-weighted median, the half-power latitude bounds displacement of the energy distribution, not of a jet axis, and recovers about 70 % of an imposed shift. The section profile masks ice-affected points monthly (the AUTHOR CHECK is resolved). The in situ half-transport latitude uses the climatological envelope. | Answered |
| C-c | The abstract repeats the Introduction, gives no results and contains "direct observational evidence" | Abstract rebuilt from HI M4 with the decision-B corrections: an approach sentence, velocity anomalies r = 0.61, 58.55°S annual r = 0.67, the 31 % pathway caveat, the wind-stress trend not detected, and the new conclusion | Answered |
| C-d | Methods appear without explanation | **v3.1:** the opening questions of 2.1–2.7 are removed (HI C2); they duplicated the Results questions. Each subsection still states its rationale: 2.3 (why an amplitude-based position; what it does not test), 2.4 (what the comparison checks), 2.5 ("Because they measure different things…"; why ice-free points), 2.1 (why altimetry only). | Answered; check that the rationale sentences are enough without the questions |
| C-e | Add transitions; state which outcome supports or rejects eddy saturation | Every Results subsection opens with its question and outcomes. **v3.1:** 3.4 and 3.6 frame the comparison in like (energy-equivalent) terms. 3.3 states that it is a prerequisite, not a test. | Answered |
| C-f | TNG14: upstream acceleration | Sects. 3.1, 3.7; Sect. 4 ¶4 now gives one verdict per link | Partly (may reflect the Pacific-wide speed increase) |
| C-g | TNG14: meander flexing | Sect. 3.7; Sect. 4 ¶4 ("not detected and, at 25 km, would not be detectable") | Answered negatively |
| C-h | TNG14: eddy fluxes/EKE | Sect. 3.7 (wake test +0.1 percentage points per decade, p = 0.95); Sect. 4 ¶4 | Answered: no localisation found |
| C-i | TNG14: the jet slows down | **v3.1 (PO S4):** Sect. 4 ¶4. Under slowly increasing forcing, the final link would appear as compensation limiting the net acceleration, not as an absolute slowing. The surface jet speed increases, so any compensation is at most partial at the surface and cannot be isolated at depth or within meanders. | Answered |
| C-j | Sector averaging removes variability within meanders | Sect. 3.7, Fig. 5b–e, SI S3 | Answered |
| C-k | Curvature at specific meanders from half-power positions | Sect. 3.7 (one sentence), SI S4 | Answered (negative) |
| C-l | Upstream acceleration at a particular meander | Sects. 3.1, 3.7 (Pacific–Antarctic Ridge) | Answered |
| C-m | Reason for the 200 cm² s⁻² threshold; sensitivity | Sect. 2.3: "set to 200 cm² s⁻², the rounded median (195)"; SI S3 | Answered |
| C-n | Is altimetry adequate? Grid spacing and effective resolution | Sect. 2.1, Table 1, Sect. 3.3/Table 2, Sect. 4 limitations | Answered; DT2024 effective resolution [CHECK] |
| C-o | Use reanalysis or models | Sect. 2.1, one sentence (HI C1 + PO S11): "Reanalysis and model output are not used, apart from the published model-derived factors of Sect. 2.5(c), because their trends depend on their forcing…" | Declined (binding decision 2) |
| C-p | An envelope-width trend cannot be detected | Withdrawn; reason in SI S3 | Answered |
| C-q | Small latitude trends | Sect. 2.3 (detectability, 70 % recovery), Sect. 3.2, Sect. 3.7. **v3.1:** the bound is stated for the persistently tracked longitudes only (31 %), and "none of the 33 bins of 5° reaches p < 0.05 (revised test)". | Answered |
| C-minor 1 | Submesoscale | Removed | Answered |
| C-minor 2 | Early eddy-saturation papers | Sect. 1 ¶2 | Answered |
| C-minor 4 | p < 0.05 only | Sect. 2.7. **v3.1:** p = 0.050 is "at the threshold"; p of about 0.03–0.07 is "marginal" (Monte Carlo ±0.01). Captions now read "open symbols: p ≥ 0.05, including values that round to 0.050". | Answered |
| C-minor 5 | Fig. 1 resembles MM21 | See C-a; Fig. 4b | Answered |
| C-minor 6 | Energetic regions downstream of bathymetry | 3000 m isobath in Figs. 1a and 4a. **v3.1:** wake boxes drawn in Fig. 4a; Sect. 3.5 now says "the largest absolute increases, mostly not significant at the cell level, lie where mean EKE is largest". | Answered |
| C-minor 7 | Is the Pacific result a basin-size effect? | **v3.1 (PO S5):** Sect. 3.5. Relative trends do not scale with area, but their significance can, because a larger sector averages more independent variability; no ranking is supported. Sects. 3.1 and 3.7 cover the meander-scale metrics. | Answered |
| C-minor 8 | Gaps in Fig. 2a | Sect. 2.3 (rewritten per HI R2) | Answered |
| C-minor 9 | Schiavon et al. (2025) | Sect. 3.7 | Answered |
| C-minor 10 | Wind variable and line | Sect. 2.1, Table 1, Sect. 2.2, Sects. 3.1 and 4 | Answered |
| C-minor 11 | Missing references; Chelton et al. (1998) | Hall and Visbeck (2002) and Saenko et al. (2005) cited; Grassi and Pujol (2025) and Ballarotta et al. (2024) in Table 1 [AUTHOR CHECK]. **v3.1:** Chelton is no longer "to be removed"; it is marked [AUTHOR CHECK: the comment is cut off in the copy received]. | Partly |
| C-data | In situ EKE definition and transect count (PO M8) | **v3.1:** Sect. 2.4. EKE is integrated over 0–780 m (the archive variable); GV23 present it depth-averaged, and the scaling does not affect correlations. The archive has 248 transects, with [AUTHOR CHECK: GV23 quote 237 transects for EKE]. | Open (author check) |
| C-cut | **The report is cut off at "References … [cut off]"** | — | **Unanswered. Obtain the rest before submission.** |

### Referee A

| Item | Where addressed | Status |
|---|---|---|
| A1: start from Drake Passage | Sect. 3 structure; roadmap at the end of Sect. 1 (shortened); Sect. 4 ¶1 | Answered |
| A1: surface transport in all sectors | Sects. 2.5(a), 3.6, Fig. 4c | Answered |
| A1: vertical structure | Sects. 2.5(b), 4 ¶3 and ¶6. **v3.1:** the slope index is described as a barotropic-type indicator of circumpolar transport variability that neither the surface proxy nor the 0–760 m record contains (if it follows the depth-independent part) [PENDING]. Sect. 4 ¶3 adds that steric rise north of the ACC can raise the proxy. | Partly (index pending) |
| A2: resolution and eddy range | Sects. 2.1, 3.6, 4; SI S8. **v3.1:** the atlas shares the constellation dependence of the maps. | Answered |
| A3: front position and EKE; baroclinic indicators | Sects. 2.3, 3.7, 4 ¶4 | Answered within limits |
| A specific: abstract; KE definition | Abstract; Sect. 2.1 | Answered |
| A specific: SAM | **v3.1:** moved from Sect. 3.5 to Sect. 4 ¶2, as part of the attribution caveat; SI S10 | Answered |

### Referee B (D. Marshall)

| Item | Where addressed | Status |
|---|---|---|
| B1: frontal position vs eddy saturation | Sect. 2.3 (does not test transport); Sects. 3.2, 3.7 | Answered |
| B2: proxies | Sect. 2.5 (three quantities). **v3.1:** energy and velocity-like quantities are compared in energy-equivalent terms (Sects. 3.4, 3.6), and "five times / several times less" is removed. | Answered |
| B3: Böning et al. (2008); adjustment | **v3.1:** Böning is cited in Sect. 1 ¶2 as observed isopycnal slopes (PO O1) and removed from Sect. 4 for space. Adjustment time scales are in Sect. 4 ¶3. | Answered |
| B4: wakes/barotropic vs baroclinic | Sects. 1 ¶3, 3.7, 4 ¶4 | Answered |
| B5: proxies (lines 80–81) | Sect. 3.3, Table 2. **v3.1:** adds the whole-line integral including values under ice (r = 0.07; SI S5). | Answered; the AC-B5 text still has to be updated |
| Forcing scaling (PO S3, decision A2) | **v3.1:** Sect. 4 ¶2. EKE trend 8.0 % > undetected stress trend 2.5 %, so the EKE increase is not attributable to the circumpolar wind on this evidence (Marshall et al., 2017 [verify]; Maddison et al., 2025). | Added |

## 3 Open items in the manuscript

**[PENDING]** (Antarctic-slope index; extraction pending with the contiguous-slope mask)
1. Abstract: the slope-index clause. Keep it to about 8 words; the abstract is at 200.
2. Sect. 3.4: the index results (ice-free fraction; index (cm) and model-calibrated (Sv) trends with intervals and p; correlations with the in situ total, southern and depth-independent parts).
3. Sect. 3.4: retain the applicable conditional clause.
4. Fig. 3: panels (a), (b) and the index entry in (c). If the index is still pending at submission, (c) merges into Fig. 2d.
5. Sect. 4 ¶2: "[PENDING: comparison with Hogg et al. (2015)]". This replaces v3's "weak decline indicated by Hogg et al. (2015)", which is not in our sources.
6. Sect. 4 limitations: the slope ice-free fraction.

**[CHECK]**
1. Sect. 2.1: DT2024 effective resolution, from the QUID.
2. Sect. 2.5, Eq. 3: signs, and the units (per mbar in Hughes et al., 2014, vs Sv cm⁻¹ in Hogg et al., 2015).

**[AUTHOR CHECK]**
1. Title. It is now the editor's choice (decision A5); alternatives are in Sect. 5.
2. Table 1: QUID citation for Product 1; Ballarotta et al. (2024) for Product 4.
3. **New:** Sect. 2.4. GV23 quote 237 transects for their EKE analysis; the archive files hold 248.
4. Code availability (Zenodo version); author contributions; acknowledgements (LM Gould funders; wording for the referees); AI statement.
5. Reference list: Freeman et al. (2016) and Hargous et al. (2026) (SI only, or drop).
6. **New:** Chelton et al. (1998). Referee C's comment is cut off.
7. **Resolved and removed:** v3 Sect. 3.2, "does the profile use points south of 58.5°S?" Points with ice fraction > 0.10 are masked monthly.

**[verify]**
1. Table 1: wind-product grid spacing (0.125°).
2. Data availability: GEBCO 2026 version and DOI.
3. Reference entries: Gutierrez-Villanueva et al. (2023), Hogg et al. (2015), Hughes et al. (2014), Martínez-Moreno et al. (2021), Schiavon et al. (2025), Liu et al. (2026), Ballarotta et al. (2024), Grassi and Pujol (2025).
4. **New:** Marshall et al. (2017). Confirm that this is the right source for the scaling of eddy energy with wind stress under eddy saturation, and give the full entry. No details were invented.

**Not bracketed but worth checking**
1. Sect. 1 ¶4: the wording of GV23's results (no significant trend in net 0–760 m transport; opposing trends among frontal regions; an EKE increase between SAF and PF consistent with altimetric trends and eddy saturation). Taken from review_PO M1; check against the paper.
2. Sect. 1 ¶2: "observed isopycnal slopes (Böning et al., 2008) … indicate a much weaker sensitivity of the circumpolar transport".
3. Sect. 3.7: "none of the 33 bins of 5° reaches p < 0.05 (revised test)". Confirm against `retest_summary.csv`, which is not on disk (numbercheck C).
4. Vessel name; monthly fields used for the comparison; surrogate type.

**For the SI revision (main-text changes the SI must match)**
- S1: add the energy-equivalent comparison for both products. The main text quotes 3.2 vs 8.0 and 3.0 vs 3.8, speed 4.0 vs KE 4.9, and Drake 7.8 vs 8.4.
- S2: explain that the interval is the union of the two model orders, so an interval can exclude zero at p = 0.050. The main text keeps only "the more conservative of two model orders". Also give the Monte Carlo resolution, and keep the post hoc list (the main text now points to it).
- S3: the 400 cm² s⁻² split is no longer mentioned in the main text.
- S6: receives the ice-masked 35–65°S variants (+4.6 %, p = 0.11; +6.6 % [1.5, 11.6], p = 0.03).
- S9: name the four wake boxes (they are named in the Fig. 4a caption), and keep the +0.18 m level rationale, which was cut from main Sect. 2.6.
- S5 and S7: in situ EKE definition, and contiguous-slope region (decision D).

## 4 Changes from v3 after internal review

1. **Title** changed to "Surface eddy and mean-flow trends of the Antarctic Circumpolar Current, 2001–2025, from altimetry compared with Drake Passage observations" (A5).
2. **Energy vs velocity.** "Five times" and "several times less" are removed everywhere. Sects. 3.4 and 3.6 give the ratios (5.0/2.5) and the energy-equivalent comparison (proxy at about 40 % of the EKE rate all-satellite, comparable in two-satellite; speed at about the KE rate; Drake proxy 7.8 vs box EKE 8.4), described as surface partition, not a saturation test (A1).
3. **Conclusion reframed.** The observations are "compatible with eddy saturation but do not discriminate it from a joint intensification of the surface flow", for four stated reasons. The EKE increase is not attributable to the undetected circumpolar stress trend (A2).
4. **Licensing narrowed.** The Drake comparison licenses interannual variability at that section; extension elsewhere is by plausibility (3.3 question, Sect. 4 ¶1, summary). Velocity agreement is reported as transect-scale anomalies (r = 0.61). Table 2 header and rows are revised (A3, C/PO M5).
5. **Gap and compensation.** The gap statement credits GV23, Hogg et al. (2015), Chambers (2018) and Sokolov and Rintoul (2009b). The north–south opposition is "partly built into the split" and does not test GV23's opposing frontal trends (A4).
6. **Abstract rebuilt** (B): 200 words without [PENDING]; no upstream clause.
7. **Drake proxy interval caveat.** The interval does not exclude an energy-equivalent change comparable to box EKE (Sects. 3.4, 4 ¶1).
8. **Antarctic-slope index.** The region is now the contiguous slope, excluding isolated rises, plateaus and ridges. The index is interpreted as a barotropic-type indicator. Fig. 3 caption updated; Methods caveats moved to SI S7.
9. **Half-power wording.** "Energy median"; KE-weighted-median caveat and about 70 % recovery; "rounded median (195)"; ice masking resolved; pathway bound limited to persistently tracked longitudes; 33 bins "(revised test)".
10. **Statistics.** Monte Carlo resolution: p of about 0.03–0.07 is marginal, including the 35–65°S gain with the extension to 2025. FDR justified (Wilks, 2016). Contrasts are in percentage points. The post hoc list points to SI S2. Captions now give the open-symbol rule.
11. **In situ details.** EKE definition (0–780 m integrated); 248 vs 237 [AUTHOR CHECK]; whole-line r = 0.07; trends not checkable, with the in situ EKE interval (−58.8 to +74.4 % per decade).
12. **Cuts and smaller fixes:** Methods cuts C1–C11 (including the opening questions); TNG14 final link reworded as partial compensation; Sect. 4 ¶5 qualifies Chambers (2018); sector-size answer; wake boxes in Fig. 4a; Kim and Orsi as a text annotation only; Hogg "weak decline" → [PENDING]; Chelton → [AUTHOR CHECK]; ice-masked variants → SI S6; SAM sentence → Sect. 4.

Not applied: PO O2 and O3 (optional, not listed in the decisions), and HI M2's suggested "five times" sentence, which decision A1/A2 supersedes. The roadmap sentence (HI S2) is shortened for the word budget.

## 5 Title alternatives

None includes "eddy saturation".
1. (Editor's decision, used) Surface eddy and mean-flow trends of the Antarctic Circumpolar Current, 2001–2025, from altimetry compared with Drake Passage observations
2. (HI S7 alternative) Antarctic Circumpolar Current surface kinetic energy trends, 2001–2025: an altimetric update compared with Drake Passage observations
3. Kinetic energy variability of the Antarctic Circumpolar Current from altimetry, 2001–2025, checked against Drake Passage observations
4. From Drake Passage to the circumpolar band: updated surface kinetic energy and transport-proxy trends of the Antarctic Circumpolar Current

v3 working title ("…anchored in Drake Passage observations") is dropped, because "anchored" implies the trends were checked (HI S7). The v3 alternative "…along an undisplaced … pathway" is dropped, because "undisplaced" overclaims (PO M7). The SI title must match (decision D).

## 6 Changes relative to v23 (carried over from cover_note_v3.md, corrected for v3.1)

- Drake-first structure (3.1–3.4 Drake; 3.5–3.7 circumpolar).
- Results open with questions and outcomes. Methods no longer do (v3.1).
- New coverage-matched in situ comparison (248 transects; 75.5 = 51.1 + 8.1 + 15.9 Sv).
- Updated Drake numbers (GV-track proxy +3.9 %, p = 0.17; in situ −1.2 %, p = 0.70); new upstream and downstream boxes; box EKE at the threshold, not a headline.
- Along-section half-power: the energy median of the merged band.
- Antarctic-slope index methods, with results [PENDING].
- Matched comparison with MM21 and Hogg et al. (2015).
- Two-satellite results and sensitivity tests in the SI.
- Overstatements corrected ("sinuosity", "no detectable trend", "larger where … larger"; v3.1 adds the energy-equivalent framing and the narrower licensing).
- Width result withdrawn.

**Follow-up outside this file:** the AC letters (AC-A A1, AC-B B2/B5, AC-C C1 and closing) still quote superseded numbers and v3 wording. They must be aligned with v3.1, notably on "five times", "several times less", the novelty claim, and "meanders neither flex nor move".
