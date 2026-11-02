# Review of manuscript v3 (physical oceanography)

Quotes from manuscript_v3.md unless marked SI; binding decisions respected.

## MUST FIX

**M1. Novelty claim is contradicted by Gutierrez-Villanueva et al. (2023, = Referee C).** Sect. 1 ¶4: "To our knowledge, the altimetric EKE increases have not been compared with in situ data". Abstract: "neither compared with in situ data nor examined with the mean flow and pathway".
- GV23 reports a significant in situ EKE increase between the SAF and PF, "consistent with altimetry-based EKE trends" and with eddy saturation.
- Hogg et al. (2015) examined transport alongside EKE, and Chambers (2018) the pathway.
- GV23's "compensation" means opposing *trends* by frontal region, partly attributed to meandering. It does not mean interannual anticorrelation.
- Sect. 1 replacement: "In Drake Passage, repeat ADCP transects show no significant trend in net 0–760 m transport over 2005–2019, with opposing trends among frontal regions, and an EKE increase between the Subantarctic and Polar Fronts that Gutierrez-Villanueva et al. (2023) found consistent with altimetric EKE trends and with eddy saturation. To our knowledge, altimetric and in situ variability have not been compared at matched coverage to define what altimetry can extend, nor have EKE, mean flow and pathway been examined jointly around the ACC."
- Abstract replacement: "…have not, to our knowledge, been compared with in situ variability at matched coverage or examined jointly with the mean flow and pathway."

**M2. Energy is compared with velocity-like quantities.** EKE and KE scale as u², while the speed and the proxy scale as u. A relative velocity change x therefore gives about 2x in energy.
- Sect. 3.6 ("about five times … about 2.5 times"), the Abstract ("several times less") and the summary overstate the contrast as a result.
- Energy-equivalent ratios: circumpolar EKE 8.0 against proxy 2 × 1.6 = 3.2, about 2.5 (two-satellite: 3.8 against 3.0, about 1.3).
- Speed: 2 × 2.0 = 4.0, close to KE at 4.9.
- Drake: proxy 2 × 3.9 = 7.8, close to box EKE at 8.4.
- Sect. 3.6 replacement: "In energy-equivalent terms (twice the relative trend of a velocity-like quantity), the surface transport proxy rises at about 40 % of the EKE rate in the all-satellite product and at a similar rate in the two-satellite product (Supporting Information S1); this describes relative rates at the surface and does not test eddy saturation."
- Apply the same correction to Sect. 3.4 ("comparable to the EKE increase"), the Abstract and the summary.

**M3. Drake "no detectable trend" is read as "small".**
- The proxy interval [−2.1, +9.9] contains the box-EKE energy-equivalent (M2).
- The Abstract's "whereas upstream…" contrasts different regions.
- "in both products" attaches all-satellite numbers to both products.
- Abstract replacement: "Over 2001–2025 the northern surface transport proxy has no detectable trend (+3.9 % per decade, −2.1 to +9.9), an interval that does not exclude a change comparable to that of EKE in the passage; upstream, KE and EKE rise by 12.8 % and 15.3 % per decade on a low mean (all-satellite; significant in both products)."
- Add the interval caveat to Sect. 4 ¶1.

**M4. The consistency check is used to license circumpolar trends.** Sect. 4 summary: "anchored by the comparison in Drake Passage, the altimetric record shows…". Sect. 3.3: "beyond the section". One 15-year section licenses interannual variability at that section only.
- Sect. 3.3 replacement: "…to extend the record at this section to 2001–2025 and, by plausibility rather than test, to other ice-free parts of the ACC?"
- Summary replacement: "Altimetry reproduces the interannual variability of upper-ocean velocity, EKE and northern transport across Drake Passage. Around the ACC, the altimetric record shows that over 2001–2025…"

**M5. Velocity correlations are not interannual (Table 2, Abstract).**
- The per-bin, per-transect r = 0.72 is dominated by the mean cross-section structure (anomalies give 0.61), and those anomalies are transect-scale.
- The mean profile (r = 0.98) cannot test variability, yet its row reads "yes".
- Fix: change the column header to "Variability reproduced?"; mean-profile row "n/a (mean structure)"; velocity row "transect-scale anomalies (r = 0.61): yes".
- Abstract: "velocity anomalies (r = 0.61 per bin and transect)".

**M6. The Antarctic-slope index is wrongly specified and interpreted.**
- (a) Sect. 2.5(c) selects "between the 1000 and 3000 m isobaths south of 60°S", as h0_antarctic_slope_index.py does. That set includes Maud Rise, the Southern Kerguelen Plateau, the Pacific–Antarctic Ridge crest and the South Scotia Ridge, none of which is on the Antarctic margin used by Hughes et al. (2014). Fix before extraction: "…the 1000–3000 m depth range of the continental slope contiguous with the Antarctic shelf, excluding isolated plateaus, ridges and seamounts".
- (b) Sect. 3.4: "it indicates the transport of the part of the passage that the surface proxy misses" is wrong. The index is a circumpolar, bottom-pressure (barotropic-mode) indicator of full-depth transport variability. Replacement: "If the index follows the in situ depth-independent part, it adds a barotropic-type indicator of circumpolar transport variability that neither the surface proxy nor the 0–760 m record contains."

**M7. "Undisplaced pathway" (Abstract) and "circumpolar-mean half-power latitude" (Sect. 3.7) overclaim.**
- The bound covers 31 % of longitudes and excludes the Pacific gap that contains the upstream Drake box.
- The metric recovers only about 70 % of an imposed shift.
- One-flank intensification moves it by 10 km (SI S4).
- Abstract: "along a pathway with no detectable displacement where the energetic band is persistently detected (31 % of longitudes)".
- Add to Sect. 2.3: "As a KE-weighted median, it also responds to unequal intensification of the flanks, so it bounds the displacement of the energy distribution rather than of a jet axis."

**M8. In situ EKE definition.** GV23 Methods give *depth-averaged 0–760 m* EKE from 237 transects (Oct 2005–Apr 2019). Table 1 and Sect. 2.4 give "integrated over 0–780 m" and 248 transects, and SI S5 still carries [CHECK]. Verify against the archive and state the definition exactly.

**M9. SI–main conflicts.**
- The SI title is still v23's, which contains "Eddy Saturation".
- Sect. 3.7 says "no sector or 5° bin has a significant trend", while SI S3 [AUTHOR CHECK] notes 4–5 bins at p < 0.05 in R2. Resolve, or write "none survives FDR control".

## SHOULD FIX

**S1. p-value resolution (Sect. 2.7).** With 300 null series, p near 0.05 has a Monte Carlo SE of about 0.013, so p = 0.04–0.06 cannot be told apart. This affects "the extension to 2025 makes the increase detectable" (0.13 → 0.04) and the Hogg Pacific box.
- Add: "With 300 null series, p-values near 0.05 carry a Monte Carlo uncertainty of about ±0.01; results with p of about 0.03–0.07 are marginal."
- Explain why the Drake EKE interval (0.2–16.6) excludes zero at p = 0.050.
- Justify FDR at 10 % "following Wilks (2016), α_FDR = 2α_global for correlated fields".

**S2. Sect. 3.1: "upstream eddies and mean flow increase together".** Upstream speed is only at the threshold, and KE contains the low-frequency EKE. Replacement: "upstream, KE and EKE increase with the wind, while the along-stream speed increase is at the threshold".

**S3. Forcing (Sect. 4 ¶2).** Eddy-saturation energy scalings give eddy energy roughly proportional to the stress (e.g. Marshall et al., 2017 [verify]). EKE rises +8.0 % per decade against an undetected +2.5 % in the stress. Add: "The circumpolar EKE trend exceeds the circumpolar stress trend severalfold, so the eddy increase cannot be attributed to the circumpolar wind alone from these data." The Abstract should state that the stress trend is not detected.

**S4. TNG14 final link (Sect. 4 ¶4).** "not observed: the surface jet speed increases" reads as evidence against the mechanism. Replacement: "Under slowly increasing forcing the final link would appear as compensation limiting the net acceleration rather than as an absolute slowing; the surface jet speed increases, so any compensation is at most partial at the surface, and the record cannot isolate it at depth or within meanders."

**S5. Sector size (Referee C, not genuinely answered).** "they do not depend on sector size" answers only the magnitude. Replacement: "Relative trends do not scale with area, but their significance can, because a larger sector averages more independent variability; they do not support a ranking of sectors."

**S6. Sect. 3.5: "Mean EKE and its increase are largest downstream…".** Few cells are significant, and the increase is proportionally uniform. Replacement: "The largest absolute increases, mostly not significant at the cell level, lie where mean EKE is largest…".

**S7. Sect. 3.3: "Trends cannot be checked: neither record has a significant trend … before the two altimetric products diverge".** Replacement: "Trends cannot be checked: over 15 annual values the in situ intervals (EKE −59 to +74 % per decade) are many times wider than the trends to be extended, and the overlap ends in 2019, as the products begin to diverge (Supporting Information S1)."

**S8. Sect. 3.3, "full-section transport is not reproduced … because".** A northern proxy fails by construction. Add: "nor does the integral over the whole line, including values mapped under seasonal ice (r = 0.07; Supporting Information S5)".

**S9. Half-power details.**
- "energy centroid" should be "energy median".
- "The threshold is the median (195)" should be "set to 200 cm² s⁻², the rounded median (195)".
- The section profile reaches 59.3°S, inside seasonal ice, where mapped KE is damped. Test a 58.55°S-truncated profile (SI S3).
- Sect. 2.4: is the "altimetric envelope" for the in situ half-transport latitude climatological or monthly? A monthly envelope builds shared variability into r.

**S10. Sect. 3.6, tracked eddies.** Add: "The atlas is built from the all-satellite maps and shares their constellation dependence (Supporting Information S1, S8)."

**S11. Sect. 2.1: "No reanalysis or model output is used".** Add "apart from the published model-derived factors of Sect. 2.5(c)".

**S12. Sect. 2.7 post hoc list.** SI S2 also lists as post hoc the EKE − KE tests (used in Sect. 3.6), the branch-start check and the covariability corrections. Align the two lists.

**S13. Surface vs depth (Sect. 4 ¶3).** The circumpolar proxy is essentially the cross-band SSH difference. Add: "a larger steric sea-level rise north of the ACC would raise it without changing the depth-integrated transport". In Sect. 2.6, "removes sea-level rise" should read "removes its band-mean part".

**S14. SI S7 wording.** It calls the index "of a different kind from the three in main-text Sect. 3.4", while Sect. 2.5 lists it as the third quantity. Align.

## OPTIONAL

- O1. Sect. 1 ¶2: Böning et al. (2008) observed isopycnal slopes, not transport.
- O2. Sect. 2.1: Eulerian along-stream speed drops with more meandering at constant jet speed; distinguish it from jet-core speed.
- O3. Sect. 3.7, r = 0.40 caveat: "the monthly jet speed and EKE both contain eddy velocities". The two-satellite r is still [NUMBER NEEDED].
- O4. SI S9 uses 1993–2025, unlike Table 1; state the exception.
- O5. SI S5: note GV23's significant PF EKE decrease beside the n.s. in situ trend.
