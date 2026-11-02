# Review of v3: Referee C perspective (Gutiérrez-Villanueva & Tedesco)

**Verdict.** Referee C would now see a contribution. Martínez-Moreno et al. (2021) (MM21) is cited and matched, and the abstract gives results. Each Results subsection says what it tests, and the chain of Thompson and Naveira Garabato (2014) (TNG14) is handled link by link.

Two problems remain:
1. The gap statement overclaims.
2. The criterion the paper sets up, EKE rising more than the mean flow, is not met. The paper still concludes "consistent with eddy saturation" without saying which comparison supports that.

Both can be fixed in the text.

## MUST

**M1. Gap overclaims (Abstract, sentence 2; Sect. 1, para 4, last sentence).** Hogg et al. (2015) paired EKE with the slope transport index that this paper adopts. Chambers (2018) and Sokolov and Rintoul (2009b) examined the pathway. Only the in situ claim holds.
*Replace:* "To our knowledge, the altimetric EKE increases have not been compared with in situ observations, and the eddy, mean-flow, transport and pathway changes have been estimated in separate studies with different periods, products and tests."

**M2. Which comparison is "consistent"? (Sect. 4, para 3; Abstract; summary.)** Sect. 1 expects "EKE … to increase more than the mean flow and transport". Sect. 3.5 says comparable trends "would indicate a joint intensification". That is the result: the EKE excess over KE is +3.0 % per decade, with p = 0.060. The conclusion then rests silently on the eddy–transport contrast alone.
*Add after "…not significantly exceeding that of KE.":* "The part of the surface pattern consistent with eddy saturation is therefore the contrast between eddies and transport proxies (an EKE trend about five times the proxy trend in the all-satellite product and 2.5 times in the two-satellite product), not a contrast between eddies and mean flow, which the data do not show."

**M3. The forcing caveat cuts both ways (Sect. 4, para 2).** The undetected circumpolar stress trend (p = 0.08) is used only to qualify the flat transport.
*Add after "…is not detected.":* "For the same reason, the circumpolar EKE increase cannot be attributed to the wind on this evidence, and its interannual variability is not significantly correlated with the SAM (Supporting Information S10)."

**M4. Abstract.** Problems in the current version:
- There is no approach sentence.
- "In both products" follows all-satellite numbers only.
- "The pathway position changes by −0.5 km" reads as a detected change.
- KE is used before it is defined.
- It gives about 20 numbers.
- It repeats the M1 gap problem.

*Proposed (198 words by `wc -w`, placeholder excluded):*

> Eddy-resolving models suggest stronger westerlies raise the eddy kinetic energy (EKE) of the Antarctic Circumpolar Current (ACC) more than its transport (eddy saturation). Earlier altimetric EKE increases lacked an in situ check, and eddy, mean-flow and pathway changes were estimated with different periods, products and tests. We compare 1/8° altimetry with Drake Passage shipboard ADCP transects (2005–2019), then extend the record around the ACC to 2001–2025. Altimetry reproduces the interannual variability of in situ velocity, EKE (annual r = 0.89) and northern transport (r = 0.67), but not full-section transport. Upstream, EKE rises by 15.3 % per decade, but the northern surface transport proxy has no detectable trend [PENDING: slope-index clause]. ACC-wide, EKE rises by 8.0 % per decade (3.3–12.8; about half with two satellites), kinetic energy by 4.9 %, along-stream speed by 2.0 % and the surface transport proxy by 1.6 %, with no detectable pathway displacement. Over 35–65°S, EKE rises by 6.6 % per decade, of the order of earlier estimates. Eddies and mean flow strengthen together while transport proxies change several times less, consistent with, but not a test of, eddy saturation: the circumpolar wind-stress trend is not detected, and surface data cannot separate equilibrium from slow adjustment.

If the slope clause exceeds about 8 words, drop the upstream clause; that result stays in Sect. 3.1.

**M5. Methods: cut ≈ 205 words (1537 → ≈ 1330).** Everything removed is already in Table 1 or the SI.

| # | Location | Action | Saves |
|---|---|---|---|
| C1 | 2.1, last sentence | "Reanalyses and models are not used, because their trends depend on their forcing; the adequacy of altimetry is tested against in situ data instead (Sect. 3.3)." | 14 |
| C2 | Opening question of each of 2.1–2.7 | Delete all seven. They duplicate the Results questions. In 2.5, write "Because they measure different things, three transport quantities are kept distinct." | 52 |
| C3 | 2.3 and 2.6, "fixed before any trend was computed" | Delete both; Sect. 2.7 states it once. | 10 |
| C4 | 2.3, "Unlike sea-surface height (SSH) contours…" | "Unlike sea-surface height contours, it is not shifted by sea-level rise (Gille, 2014); it does not test the transport, which a displacement changes only through the Coriolis parameter." | 7 |
| C5 | 2.3, gaps sentence | Use R2 below. | ≈ 8 |
| C6 | 2.4, first sentence | "The in situ record consists of 248 ADCP transects of the ARSV Laurence M. Gould across Drake Passage over 2005–2019 (Product 3; Table 1; Gutierrez-Villanueva et al., 2023)." | 18 |
| C7 | 2.4, "Because altimetry gives…" | "As altimetry gives surface geostrophic flow and the in situ quantities include deeper and ageostrophic flow, the comparison checks variability, not magnitude." | 10 |
| C8 | 2.5(c) | Delete ", most strongly at 1000–3000 m depth" and "in 100 m depth bins and 12 longitude bands, recording the ice-free fraction" (both are in SI S7). Replace the last sentence with: "The conversion is least constrained for trends, and margin sea-level change and sea-ice sampling may alias into the index (Supporting Information S7)." | ≈ 45 |
| C9 | 2.6, last two sentences | "Four wake boxes downstream of major topography (Supporting Information S9) test whether the EKE increase concentrates in the wakes. Tracked eddies come from META4.0 (Product 5; Table 1)." | 23 |
| C10 | 2.3, "(sensitivity to 150–400 cm² s⁻², 20–50 % and 50–300 km in Supporting Information S3)" | "(sensitivity in Supporting Information S3)" | 9 |
| C11 | 2.7, AIC/AR(12) sentences | "Monthly series take the more conservative of two model orders and annual series use AR(1); false-positive rates at the nominal 5 % level are 3–5 % for monthly series and up to 9.5 % for 15 annual values (Supporting Information S2)." | 10 |

Keep the 2.5(b) Sv decomposition, because Referee C supplied these data.

**M6. Discussion (+≈ 100 words).**
- Spend about 75 words on M2 and M3.
- Replace the last sentence of para 5, which restates the Introduction, with: "It also qualifies earlier results: the regional pathway shifts reported for 1993–2016 (Chambers, 2018) are not found over 2001–2025, which may reflect the period or the method, and the EKE increase is accompanied by a smaller but detectable strengthening of the surface mean flow."
- Offset this by cutting para 1 to its second sentence. The first sentence repeats Table 2.

## SHOULD

**S1. Novelty paragraph (end of Sect. 1).** The opening Ocean State Report (OSR) sentence reads as defensive, and "one autocorrelation-calibrated test" in (ii) is unclear.
- Replace (ii) with: "repeats the estimates of Martínez-Moreno et al. (2021) and Hogg et al. (2015) with their procedure and with an autocorrelation-calibrated test, in two altimetry products".
- Move the OSR sentence after (iv): "As an Ocean State Report update, the result is informative whether it confirms or qualifies earlier estimates."

**S2. Linearity.** Sects. 3.1–3.2 present altimetric results before Sect. 3.3 shows they can be trusted. The order is fixed, so tell the reader why.
*Replace the roadmap sentence:* "We start in Drake Passage, where the wind and energy trends (Sect. 3.1) and the band position (Sect. 3.2) are set against the in situ record (Sect. 3.3) and the transport indicators (Sect. 3.4); the quantities that the comparison supports are then extended around the ACC (Sects. 3.5–3.7)."

**S3. Redundancy.**
- The northern proxy trend appears in both Sect. 3.3 and 3.4. Keep it in 3.4 only.
- "The comparison thus supports extending…" in Sect. 3.3 repeats Table 2.
- Sect. 4, para 4 repeats the numbers from Sect. 3.7. Give one verdict per link instead.

**S4. Sect. 3.5 loses the reader** (455 words, about 30 numbers, two questions). Move the ice-masked 35–65°S variants (+4.6/+6.6 %) to SI S6; they duplicate the unmasked values. Apply R3.

**S5. Table 2.** For the band position, replace "supported" with "interannual variability only (RMS 114 km)".

**S6. Figures.**
- Five figures is right. Merge Fig. 3 into Fig. 2d if the index is still pending.
- Draw the wake boxes in Fig. 4a.
- Fig. 5a includes Kim and Orsi (2014) with "no rate used". Remove it or explain how it is plotted.

**S7. Title.** "Anchored in" implies the trends were checked against in situ data, but Table 2 says they cannot be.
- *Recommended:* "Surface eddy and mean-flow trends of the Antarctic Circumpolar Current, 2001–2025, from altimetry compared with Drake Passage observations". This also names novelty point (iii).
- *Alternative:* "Antarctic Circumpolar Current surface kinetic energy trends, 2001–2025: an altimetric update compared with Drake Passage observations".

**S8. SI title.** It still reads "…and Eddy Saturation in the Southern Ocean". Align it with the main title.

## Readability rewrites

**R1 (Sect. 3.1, last sentence; "saturation-like" is jargon).** "Inside the passage the pattern thus matches the eddy-saturation expectation only weakly, whereas upstream eddies and mean flow increase together with the wind, possibly as part of the Pacific-wide increase in along-stream speed (Sect. 3.6)."

**R2 (Sect. 2.3, gaps).** "Around the ACC, monthly envelopes are matched to the mean envelope at each longitude. Only the 900 longitudes (about 31 % of the band) matched in at least 90 % of months are kept; the excluded gaps lie mainly in the Pacific, including 115–75°W, where the band is weak and broad."

**R3 (Sect. 3.5, MM21).** "With the procedure of Martínez-Moreno et al. (2021), our data give +5.2 % per decade over 2001–2020 and +7.5 % over 2001–2025, the same order as their +5.1 % over 1993–2020. Our integrated EKE level is about twice theirs, for reasons not yet investigated (Supporting Information S6)."

**R5 (Sect. 3.7, opening).** "A fixed-band trend could reflect a moving pathway. The TNG14 mechanism implies upstream acceleration, meander flexing and EKE increases localised downstream of topography. Detecting these at the surface would make it a candidate explanation; not detecting them would leave it untested."

**R6 (Sect. 3.7, subject of "and an order of magnitude below…" is unclear).** Split after "…1993–2016." and continue: "All are an order of magnitude below the SSH-contour shifts of Sokolov and Rintoul (2009b)."

## OPTIONAL

- **O1.** Referee C's report is truncated. Obtain the rest before resubmission.
