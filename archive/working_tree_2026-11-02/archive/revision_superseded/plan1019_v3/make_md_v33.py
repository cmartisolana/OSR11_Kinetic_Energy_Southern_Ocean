import re
s = open('/home/claude/v30/manuscript_v3.2W.md', encoding='utf8').read()
R = []
def rep(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:80]); s = s.replace(a, b)

# ---------------- key points + abstract
old_abs = re.search(r'\*\*Abstract\.\*\*.*', s).group(0)
new_abs = ("**Key points.** [AUTHOR CHECK: not part of the State of the Planet template; for the cover letter, short summary and OSR synthesis]\n\n"
"- The surface eddy kinetic energy of the Antarctic Circumpolar Current rose by about 8 % per decade over 2001–2025, because existing eddies became stronger, not more numerous.\n"
"- The surface mean flow accelerated as well, by about 2 % per decade, while the energetic pathway stayed in place: the current intensified rather than migrating south.\n"
"- Drake Passage transport indicators show no detectable trend; with no detected trend in circumpolar wind stress, eddy saturation is compatible with the record but not tested by it.\n\n"
"**Abstract.** The Southern Hemisphere westerlies have strengthened, and theory predicts that the Antarctic Circumpolar Current (ACC) responds through its eddies more than through its transport (eddy saturation). We use 1/8° satellite altimetry over 2001–2025 to examine the surface expression of this response. The surface eddy kinetic energy (EKE) of the ACC rose by 8.0 % per decade (95 % interval 3.3–12.8; about half as large in a two-satellite product), because existing eddies became stronger, not more numerous. The surface mean flow accelerated too (along-stream speed +2.0 % per decade), and the energetic pathway did not shift detectably (−0.5 km per decade, −4.2 to +3.1): the current intensified in place rather than migrating south. The largest changes lie upstream of Drake Passage, where the eastward wind stress, kinetic energy and EKE increased by 10–15 % per decade. Within the passage, the flow accelerated between the Subantarctic and Polar Fronts while the surface transport between streamlines barely changed, consistent with more vigorous meandering, and no transport indicator shows a detectable trend [PENDING: Antarctic-slope index clause]. Because the circumpolar wind-stress trend is not detected and the data are surface-only, these changes are compatible with eddy saturation but do not distinguish it from a joint intensification of eddies and mean flow.")
rep(old_abs, new_abs)

# ---------------- methods: wake boxes no longer drawn
rep("Four wake boxes downstream of major topography (Fig. 2a) test whether the EKE increase concentrates in the wakes.",
    "Four wake boxes within the ACC band downstream of major topography (Kerguelen Plateau, 75–95°E; Southeast Indian and Macquarie Ridges, 145–175°E; Pacific–Antarctic Ridge, 150–125°W; Drake Passage, 55–35°W) test whether the EKE increase concentrates in the wakes.")

# ---------------- section headings (physical statements)
rep("## 3.1 Drake Passage: forcing, energy and changes along the transect", "## 3.1 Drake Passage: stronger winds and flow upstream, faster flow between the fronts")
rep("## 3.2 Circumpolar EKE increase and earlier estimates", "## 3.2 Stronger eddies around the ACC")
rep("## 3.3 Transport from altimetry", "## 3.3 Mean flow and transport")
rep("## 3.4 Pathway: half-power positions", "## 3.4 A stationary energetic pathway")
rep("## 3.5 Standing meanders", "## 3.5 Standing meanders: no localised response")

# ---------------- 3.1 figure references
rep("but not significantly downstream (+9.4 %, p = 0.09; Fig. 1b). The clearest oceanic change lies upstream:",
    "but not significantly downstream (+9.4 %, p = 0.09; Fig. 1b). The clearest oceanic change lies upstream, where EKE increased over most of the ACC band between the Subantarctic and Polar Fronts (Fig. 1a):")
rep("In streamwise coordinates, which assign meanders and eddies to their parent jets, the surface transport between streamlines changes little (Fig. 1d):",
    "In streamwise coordinates, which assign meanders and eddies to their parent jets, the surface transport between streamlines changes little (Fig. 1d):")

# ---------------- 3.2
rep("The largest absolute increases, mostly not significant at the cell level, lie where mean EKE is largest, downstream of the Kerguelen Plateau, Southeast Indian Ridge and Campbell Plateau (Fig. 2a).",
    "The largest absolute increases, mostly not significant at the cell level, lie where mean EKE is largest, downstream of the Kerguelen Plateau, Southeast Indian Ridge and Campbell Plateau and in the Brazil–Malvinas and Agulhas regions north of the band (Fig. 2a).")
rep("KE and EKE trends are positive in all sectors, significant for both in the Pacific and Indian and for KE in the Atlantic (EKE p = 0.053).",
    "The rise is gradual, with the highest EKE in 2019–2024 (Fig. 2c), and KE and EKE trends are positive in all sectors, significant for both in the Pacific and Indian and for KE in the Atlantic (EKE p = 0.053; Fig. 2d).")
rep("With the calibrated test the trend is +4.5 % per decade (p = 0.13) over 2001–2020 and +6.6 % (1.2–12.0; p = 0.04) over 2001–2025 (Fig. 2b);",
    "With the calibrated test the trend is +4.5 % per decade (p = 0.13) over 2001–2020 and +6.6 % (1.2–12.0; p = 0.04) over 2001–2025 (Supporting Information S6, Fig. S7);")
rep("whereas their amplitude increases by 4.6 % per decade (2.4–6.8; p = 0.007) and their rotational speed by 3.8 % (0.9–6.8; p = 0.02) (Fig. 2d).",
    "whereas their amplitude increases by 4.6 % per decade (2.4–6.8; p = 0.007) and their rotational speed by 3.8 % (0.9–6.8; p = 0.02) (Fig. 2e).")

# ---------------- 3.3
rep("2.6 % in the Pacific, 1.8 % in the Atlantic and 0.4 % (not significant) in the Indian sector (Fig. 2c).",
    "2.6 % in the Pacific, 1.8 % in the Atlantic and 0.4 % (not significant) in the Indian sector (Fig. 2b, d).")

# ---------------- 3.4 (no regions of earlier studies)
old34 = s[s.index("Averaged over the 900 persistently tracked longitudes"):s.index("## 3.5")]
new34 = ("Averaged over the 900 persistently tracked longitudes, the half-power latitude changes by −0.5 km per decade (−4.2 to +3.1; p = 0.74), although a 25 km shift would have been detected; its interannual excursions of 10–15 km reverse within a few years (Fig. 3b). No sector has a significant trend, and none of the 33 longitude bins of 5° reaches p < 0.05 (Fig. 3c; Supporting Information S3), including the southeast Indian and southeast Pacific sectors where shifts were reported for earlier periods (Kim and Orsi, 2014; Chambers, 2018). Our trends are an order of magnitude below the southward shifts of SSH contours (about 40 km per decade over 1992–2007; Sokolov and Rintoul, 2009b). Along the LM Gould section, where the Subantarctic and Polar Fronts merge into one energetic envelope, the half-power latitude shows no detectable displacement either (−10 km per decade, −36 to +16; Fig. 3c). In CMIP6 projections, Schiavon et al. (2025) likewise found the ACC jets largely anchored by topography, with the largest projected displacements in the topographically unconstrained southeast Indian sector.\n\n")
rep(old34, new34)

# ---------------- 3.5 figure letters
rep("Along the jet-core path, EKE increases in 58 of 70 bins (Fig. 3c),", "Along the jet-core path, EKE increases in 58 of 70 bins (Fig. 3d),")
rep("The jet speed increases in 51 of 70 bins (median +4.1 % per decade; Fig. 3d).", "The jet speed increases in 51 of 70 bins (median +4.1 % per decade; Fig. 3e).")
rep("shows a tendency (Fig. 3e, f).", "shows a tendency (Fig. 3c, f).")

# ---------------- discussion: key results as physical statements
old_d1 = s[s.index("Three results emerge from the 2001–2025 altimetric record."):s.index("The wind stress increased by about 9 % per decade")]
new_d1 = ("Three results emerge from the 2001–2025 altimetric record. First, the eddy field of the ACC has strengthened at the surface: EKE increased by 8 % per decade because existing eddies became more energetic, and the earlier altimetric increase (Hogg et al., 2015; Martínez-Moreno et al., 2021) persists in the reprocessed record with a conservative test. Second, the surface mean flow strengthened as well, so the record shows a joint intensification rather than eddies growing over an unchanged mean flow. Third, the current intensified in place: the energetic pathway has not moved detectably where it is persistently tracked, and no regional shift of the kind reported for 1993–2016 (Chambers, 2018) is detected, a difference that may reflect the period or the method. The magnitude of the EKE increase depends on the constellation by about a factor of two, a difference not significant over 25 years but large enough to quote both products.\n\n"
"These results bear on how the ACC responds to a changing climate and how it is monitored. A stationary pathway is consistent with topographic anchoring of the jets (Schiavon et al., 2025) and with the view that the southward drift of SSH contours largely reflects sea-level rise rather than a migration of the flow (Gille, 2014); warming or freshening observed at fixed latitudes in the ACC therefore need not imply a southward-moving current. Surface EKE has become a detectable trend indicator over the 25-year record, but its magnitude depends on the satellite constellation, so its monitoring needs a constant-constellation product alongside the all-satellite one. Because the mean flow strengthened together with the eddies, the eddy field alone would miss part of the surface response; the along-stream speed and the pathway position, bounded to a few kilometres per decade where tracked, are complementary indicators.\n\n")
rep(old_d1, new_d1)

# ---------------- figure captions
cap1_old = s[s.index("**Figure 1.**"):s.index("\n", s.index("**Figure 1.**"))]
cap1 = ("**Figure 1.** Drake Passage, 2001–2025 (Product 1). (a) EKE trend (cm² s⁻² per decade; least-squares trend of annual means on 0.5° cells; blank where the ice fraction exceeds 0.10 in any month), with the mean Subantarctic (SAF) and Polar (PF) Fronts (the central streamlines of the frontal zones of Gutierrez-Villanueva et al. (2023), mapped onto the mean ADT), the 3000 m isobath (grey), the Drake box and the upstream and downstream regions within the ACC band (dotted), and the mean LM Gould track (orange; dashed: seasonally ice-covered part). (b) Annual anomalies, in percent of the 2001–2025 mean, of eastward wind stress, KE and EKE upstream, in the Drake box and downstream, and of the northern surface transport proxy; lines show the trends (Theil–Sen on the deseasonalised monthly series), solid where p < 0.05 and dashed otherwise. (c) Trends of the cross-track velocity in 25 km bins of distance from the north along the LM Gould transect: in situ ADCP velocity at 70–118 m (orange) and altimetric surface geostrophic velocity (green) over October 2005–April 2019 (annual means), and altimetric velocity over 2001–2025 (black); grey shading: mean altimetric velocity (right axis); dotted: seasonally ice-covered part, where altimetric months with ice fraction above 0.10 are excluded; SAF, PF and SACCF mark the mean positions of the Subantarctic, Polar and Southern ACC Fronts. (d, e) Trends of the transport (d) and EKE (e) between pairs of sea-surface-height streamlines, in percent of the mean per decade, with the streamline values of Gutierrez-Villanueva et al. (2023) and the same colours as (c) (in situ: 0–760 m transport and 0–780 m EKE; altimetry: surface transport per unit depth and surface EKE); shading marks the SAF, SAF–PF and PF zones, and black bars with grey bands the zone trends over 2001–2025 with 95 % intervals. Filled symbols: p < 0.05; open symbols: p ≥ 0.05, including values that round to 0.050.")
rep(cap1_old, cap1)
cap2_old = s[s.index("**Figure 2.**"):s.index("\n", s.index("**Figure 2.**"))]
cap2 = ("**Figure 2.** Eddies and mean flow of the ACC, 2001–2025 (Product 1 unless stated). (a) EKE trend on 1° cells (cm² s⁻² per decade; dots: significant after FDR control), with the 3000 m isobath, the −0.6 and +0.2 m mean ADT contours bounding the ACC band, the sector boundaries (dashed), the sea-ice mask (light blue) and the main topographic obstacles. (b) Annual anomalies over the ACC band of the eastward wind stress over 45–60°S, the along-stream speed and the surface transport proxy, and (c) of EKE and KE, with EKE from the two-satellite product (Product 4, grey), in percent of the 2001–2025 mean; lines show the trends, solid where p < 0.05 and dashed otherwise. (d) Trends of EKE, KE, along-stream speed and surface transport proxy for the ACC band and its sectors, in percent of the mean per decade with 95 % intervals; open bars: p ≥ 0.05. (e) Annual anomalies of the mean amplitude, rotational speed, number present and speed radius of tracked eddies in the ACC band (Product 5), with their trends over January 2001–August 2023.")
rep(cap2_old, cap2)
cap3_old = s[s.index("**Figure 3.**"):s.index("\n", s.index("**Figure 3.**"))]
cap3 = ("**Figure 3.** Pathway and standing meanders, 2001–2025 (Product 1). (a) Mean surface speed, with the mean half-power latitude at the 900 persistently tracked longitudes (purple) and at the other longitudes of the ACC band (grey), the mean jet-core path (red; +0.18 m relative-ADT contour), the 3000 m isobath and the LM Gould track (orange). (b) Circumpolar mean anomaly of the half-power latitude over the tracked longitudes (km, positive northward): monthly (grey), 12-month running mean (purple) and trend with 95 % interval (dashed, shading). (c–f) Trends along the ACC in 5° longitude bins, flow from left to right, with 95 % intervals: (c) half-power latitude (purple), latitude of the jet-core path (grey diamonds) and half-power latitude of the merged Subantarctic–Polar Front band on the LM Gould section (orange star), in km per decade (positive northward); (d) EKE within 150 km of the path, (e) jet speed and (f) sinuosity, in percent per decade. (g) Depth under the mean path (dark: bin mean; light: shallowest point). Grey bands: bins where the path crosses depths shallower than 3000 m; hatching: excluded Malvinas Current bins. Filled symbols: p < 0.05; open symbols: p ≥ 0.05, including values that round to 0.050; ringed: significant after FDR control.")
rep(cap3_old, cap3)

# ---------------- reference notes
rep("(citation changed: now cited for the absence of significant jet-position trends and interannual variability over flat topography, not for climate projections)",
    "(citation corrected: CMIP6 SSP5-8.5 projections with two models, cited for the projected topographic anchoring of the ACC jets, with the largest displacements in the southeast Indian sector)")
rep("Kim and Orsi (2014) is now cited only in the Fig. 5 caption.", "Kim and Orsi (2014) is cited in Sect. 3.4 for the southeast Indian drift.")
open('/home/claude/v30/manuscript_v3.3W.md', 'w', encoding='utf8').write(s)
print("ok")
