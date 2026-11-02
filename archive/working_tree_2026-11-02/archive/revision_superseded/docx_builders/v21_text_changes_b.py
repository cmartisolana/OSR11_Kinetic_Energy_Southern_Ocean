"""v2.1 (second pass, 2026-10-14): text changes answering the three referees. Imports the first-pass changes and
modifies them. Numbers: outputs/retest, a8 (sector surface transport), a9 (SAM), b2 (box-free path-EKE), eddies.
Placeholders __TWOSAT_T__ are filled by build when outputs/a8/a8_trends_twosat.csv exists."""
import changes as C
REPLACE, INSERT_AFTER, NEW_REFS, NEW_TABLE_ROW = dict(C.REPLACE), dict(C.INSERT_AFTER), list(C.NEW_REFS), C.NEW_TABLE_ROW

# ---------------- Abstract (8)
REPLACE[8] = ("Here, we provide an observation-based assessment of long-term trends and regional contrasts in kinetic energy "
    "across the Southern Ocean over 2001–2025, using gridded satellite altimetry together with a tracked-eddy atlas, "
    "scatterometer wind stress and an in-situ Drake Passage transport record. We ask whether the surface expression of "
    "the Antarctic Circumpolar Current (ACC) over the altimetric era is consistent with eddy saturation, in which "
    "stronger winds energise the eddy field rather than accelerate the circumpolar transport. Kinetic energy (KE) of the "
    "monthly-mean flow and eddy kinetic energy (EKE) in the ACC band increase by 4.9 % and 8.0 % per decade in the "
    "all-satellite product, and by 3.6 % and 3.8 % per decade in a product built from a constant pair of altimeters; the "
    "EKE increase is carried by larger amplitudes and faster rotation of a constant population of eddies. The "
    "energetic ACC pathway has not moved: its circumpolar position changes by −0.5 ± 4 km per decade, and distortions of "
    "the jet path covary with EKE in the same month without a detectable downstream response. The surface transport "
    "proxy across the ACC band rises by 1.6 % per decade and the Drake Passage proxies show no detectable trend, while "
    "the local eddy field strengthens. These results are consistent with an eddy-saturated regime at the surface, "
    "although altimetry cannot resolve the vertical structure of the response and the record cannot separate an "
    "equilibrated state from a slowly adjusting circulation.")

# ---------------- Introduction
REPLACE[36] = ("The Southern Ocean plays a crucial role in the Earth system by regulating the uptake of heat and carbon and "
    "enabling exchanges between the surface and deep ocean (Gray, 2024). Its circulation is dominated by the Antarctic "
    "Circumpolar Current (ACC), whose variability strongly influences global ocean dynamics and climate. This is one of "
    "the most energetic regions in the global ocean, where mesoscale eddies and meanders drive the horizontal and "
    "vertical fluxes of momentum and tracers essential for global climate stability. A substantial fraction of the "
    "ocean’s kinetic energy (KE) is contained in mesoscale eddies, with eddy kinetic energy (EKE) providing a measure of "
    "mesoscale activity (Meredith and Hogg, 2006; Barceló-Llull et al., 2025; Cortés-Morales et al., 2026).")
REPLACE[38] = ("While early ocean models predicted that intensifying westerly winds would accelerate the ACC (e.g., Hall and "
    "Visbeck, 2002; Saenko et al., 2005), observations and eddy-resolving models have shown a much weaker sensitivity of "
    "circumpolar transport to wind forcing than in coarse-resolution models (Böning et al., 2008; Thompson and Naveira "
    "Garabato, 2014; Patara et al., 2016; Gutierrez-Villanueva et al., 2023; Stewart et al., 2023). This behaviour, "
    "first described in idealised and coarse eddy-permitting models (Straub, 1993; Hallberg and Gnanadesikan, 2001, "
    "2006) and confirmed in eddy-resolving experiments (Munday et al., 2013), is known as eddy saturation: enhanced wind "
    "stress produces relatively small changes in the time-mean circumpolar transport while increasing the energy of "
    "transient eddies and standing meanders.")
REPLACE[40] = REPLACE[40].replace("Thompson and Garabato, 2014", "Thompson and Naveira Garabato, 2014")
REPLACE[41] = ("Identifying these dynamics using standard Level 4 gridded altimetry products presents a technical challenge, as "
    "the optimal interpolation used to create these maps acts as a low-pass filter, underestimating EKE relative to "
    "along-track data (Chambers, 2018) and resolving only the larger part of the mesoscale band at the latitudes of the "
    "ACC (Ballarotta et al., 2019). A separate challenge concerns the tracking of frontal positions: conventional "
    "methods based on sea-surface height (SSH) contours are sensitive to sea-level rise, which can bias apparent frontal "
    "displacements (Gille, 2014). To address this, we adopt the “half-power point” method (Chambers, 2018), which tracks "
    "the meridional centre of enhanced KE envelopes rather than absolute SSH values, in order to establish whether the "
    "energetic ACC pathway has moved while its energy increased. In this study, each diagnostic is introduced with the "
    "question it answers: the KE and EKE trends and the tracked-eddy atlas ask whether and how the eddy field and the "
    "mean flow have intensified (Sect. 3.1); the frontal analysis asks whether the energetic pathway has moved "
    "(Sect. 3.2); the Drake Passage and circumpolar surface transport proxies ask whether the transport has changed as "
    "the eddy field strengthened (Sect. 3.3); and the jet-path analysis asks whether distortions of the current precede "
    "the growth of EKE, as expected if standing meanders mediate the eddy response (Sect. 3.4).")

# ---------------- Methods 2.1
REPLACE[44] = ("Given the horizontal geostrophic velocity field of Product 1 (all-satellite, 1/8° grid) or Product 4 "
    "(two-satellite, 1/4° grid), the KE of the monthly-mean flow is defined as half the squared speed of the monthly-mean "
    "velocity, and EKE as:")
REPLACE[46] = ("where  and  are the components of the long-term mean geostrophic velocity field averaged over the 2001-2025 "
    "period. The analysis period begins in 2001, when the altimetric constellation first reached a minimum of three "
    "simultaneous missions, ensuring consistent spatial coverage and a more stable effective resolution throughout the "
    "record (Cortés-Morales et al., 2026). The effective resolution of the merged maps is set by the mapping rather "
    "than by the grid, and is of the order of 100–150 km in wavelength at the latitudes of the ACC (Ballarotta et al., "
    "2019), so that KE and EKE here refer to the surface geostrophic flow at scales resolved by the gridded maps. "
    "The along-stream speed of the monthly-mean flow, the component of the velocity along the direction of the "
    "2001–2025 mean velocity, is used as a measure of the mean-flow intensity. Area means over the ACC band, defined "
    "by mean absolute dynamic topography (ADT) between −0.6 and +0.2 m and free of sea ice in every month, and over its "
    "Atlantic (70°W–20°E), Indian (20–147°E) and Pacific (147°E–70°W) sectors, are cosine-weighted.")

# ---------------- Methods 2.2: parameters fixed in advance
REPLACE[49] = ("For each longitude, contiguous latitude segments where KE exceeds 200 cm² s⁻² are identified. Within each "
    "segment, the peak KE value is located and an envelope is defined between the latitudes where KE first drops to "
    "25 % of the peak, subject to a minimum width of 100 km. These three parameters were fixed before any trend was "
    "computed; their influence is examined in Supporting Information S3. Within this envelope, the half-power latitude "
    " is the latitude at which the cumulative meridional integral of KE reaches half its total:")

# ---------------- Methods 2.3: transport proxies
REPLACE[52] = "2.3 Transport proxies"
REPLACE[53] = ("Surface transport across the Drake Passage is estimated from gridded altimetry (Product 1) as the geostrophic "
    "velocity component normal to a transect from 65°W, 55°S to 60°W, 63°S, integrated along the section over the "
    "points that are free of sea ice in every month of the record (53 of 121 points, north of 58.5°S). A fixed set of "
    "points is required because a section whose length varies with the seasonal ice cover produces a spurious trend. "
    "No depth integration is applied, so this quantity is a surface transport proxy rather than a volume transport. It "
    "is therefore expressed as a standardized anomaly (z-score) and used only to compare temporal variability and "
    "trends. Displacing the transect end points by up to 1° does not change the conclusions. The same construction is "
    "applied around the whole ACC: at each longitude the eastward geostrophic velocity is integrated over the fixed ACC "
    "band defined in Sect. 2.1, giving a surface transport proxy per unit depth that is averaged over each sector and "
    "over the circumpolar band and expressed in percent of its mean.")
REPLACE[54] = ("As an independent comparison, we use the Drake Passage transport time series of Gutierrez-Villanueva et al. "
    "(2023), the net transport over 0–760 m from 244 repeat ADCP transects between 2005 and 2019, which we average into "
    "annual means. Both quantities are proxies for the full-depth circumpolar transport, which has a substantial "
    "barotropic component that neither captures (Donohue et al., 2016); we refer to them collectively as Drake Passage "
    "transport proxies. Wind forcing is characterised by the eastward wind stress of Product 2 averaged over the Drake "
    "Passage box (80–52°W, 63–54°S) and over the circumpolar band 45–60°S, and the Southern Annular Mode by the "
    "station-based index of Marshall (2003).")

# ---------------- Methods 2.4: jet path (box-free), replaces the first-pass meander section
INSERT_AFTER[54] = [("Heading2", "2.4 Jet path and eddy kinetic energy"),
    ("Normal", "To test whether distortions of the ACC path precede the growth of EKE, as in the standing-meander "
     "mechanism of Thompson and Naveira Garabato (2014), the path of the main jet is tracked each month as a contour of "
     "ADT relative to its monthly mean over the ACC band, which removes sea-level rise (Gille, 2014). The level, +0.18 m, "
     "is the one whose climatological contour carries the highest mean geostrophic speed around the ACC; it follows the "
     "Subantarctic Front, and a second level on the Polar Front side (−0.20 m) is used as a sensitivity test. The path is "
     "analysed in 36 overlapping windows of 30° of longitude centred every 10°: in each window and month the longest "
     "contour piece spanning the window is smoothed along the arc at 150 km and its length is taken as the measure of "
     "distortion, and EKE is averaged over the grid points within 150 km of the monthly path, in the same window and in "
     "the band 15–35° further east, which does not overlap the window. Monthly series are deseasonalised, detrended and "
     "smoothed over two months; the correlation between path length and EKE is computed with EKE lagging by 0–3 months, "
     "and its significance from 1000 phase-randomised surrogate path-length series, taking the maximum over the four "
     "lags. The same test is applied to the number of eddies born within 150 km of the path (Sect. 2.5) per unit path "
     "length. An earlier version of this analysis, restricted to three meanders with sections defined by hand, is "
     "given in Supporting Information S4 for comparison."),
    ("Heading2", "2.5 Tracked eddies"),
    ("Normal", C.EDDY_METHODS.replace("Trends in these series are estimated as in Sect. 2.6",
     "The nonlinearity ratio U/c is computed for each tracked observation as the rotational speed over the translation "
     "speed obtained from the displacement of the eddy centre over seven days. Trends in these series are estimated as "
     "in Sect. 2.6"))]

# ---------------- Methods 2.6: SAM sentence
REPLACE[57] = REPLACE[57].replace("All trends are computed for the all-satellite product (Product 1) and for the two-satellite product (Product 4); where they differ, both are reported.",
    "All trends are computed for the all-satellite product (Product 1) and for the two-satellite product (Product 4); where they differ, both are reported. The relation between the ACC-band series and the SAM index is tested on deseasonalised, detrended monthly series with the SAM leading by 0–12 months, with significance from phase-randomised surrogates of the index.")

# ---------------- Results 3.1
REPLACE[62] = ("The spatial distribution of mean KE and EKE over the 2001-2025 period reveals that energy is strongly "
    "concentrated along the ACC, where the highest values closely follow the main frontal pathways (Figure 1a,b). Mean "
    "KE values exceed 1500 cm² s⁻² in the most energetic regions, with enhanced levels particularly evident downstream "
    "of major bathymetric features such as the Kerguelen Plateau, the Campbell Plateau, and the Southeast Indian Ridge "
    "(3000 m isobath in Figure 1). The spatial patterns of KE and EKE are broadly similar, although EKE highlights more "
    "clearly the regions of intense mesoscale variability.")
REPLACE[64] = REPLACE[64].replace("Both series also show substantial interannual variability, with fluctuations at periods consistent with the SAM.",
    "Both series also show substantial interannual variability. This variability is not explained by the SAM: the "
    "correlation of the deseasonalised, detrended KE and EKE series with the SAM index is at most 0.15 (KE, same "
    "month) and 0.09 (EKE), not significant once the lags examined are accounted for, and the 12-month running means "
    "are uncorrelated (r = 0.11 and 0.08).")
# eddy paragraph: add U/c
INSERT_AFTER[65] = [(s, t.replace("The EKE increase thus reflects an intensification of the existing eddies rather than a change in their number or size.",
    "The nonlinearity ratio U/c is 6.9 in the median, 99 % of the observations are nonlinear, and the ratio shows no "
    "trend (p = 0.97): rotational and translation speeds increase together, by 3.5 % and 3.2 % per decade. The EKE "
    "increase thus reflects an intensification of the existing eddies rather than a change in their number, size or "
    "character.") if s == "Normal" and t.startswith("The tracked-eddy atlas") else t) for s, t in INSERT_AFTER[65]]
INSERT_AFTER[65] = [(s, t.replace("(b) Trends over the same months in percent per decade with 95 % confidence intervals (filled symbols: p < 0.05); the gridded EKE and KE are shown for comparison.",
    "(b) Trends over the same months in percent per decade with 95 % confidence intervals (filled symbols: p < 0.05), including the nonlinearity ratio U/c; the gridded EKE and KE are shown for comparison.")) for s, t in INSERT_AFTER[65]]

# ---------------- Results 3.3: sector transport proxy
REPLACE[73] = "3.3 Transport proxies and the Drake Passage"
REPLACE[75] = REPLACE[75].replace("Figure 3. Eddy saturation diagnostics in the Drake Passage", "Figure 3. Transport proxies and eddy kinetic energy in the Drake Passage")
REPLACE[78] = ("The EKE time series (Figure 3e) shows this intensification directly, with the running mean rising most markedly "
    "from around 2016 onward. Around the whole ACC, the surface transport proxy across the fixed band increases by "
    "1.6 % per decade (0.9–2.3; p = 0.003) in the all-satellite product__TWOSAT_T__, with 2.6 % (1.5–3.8; p = 0.003) in "
    "the Pacific sector, 1.8 % (0.9–2.7; p = 0.007) in the Atlantic sector and no significant trend in the Indian "
    "sector (0.4 %; p = 0.49); resolved in 5° longitude bins, 53 of 72 bins are positive and the significant increases "
    "are concentrated in the Pacific sector between 145°W and 90°W. The surface expression of the mean flow therefore "
    "strengthens by a few percent per decade, as the along-stream speed does, while EKE increases three to five times "
    "faster in relative terms in every sector. The combination of a weakly strengthening surface flow, transport "
    "proxies without a detectable trend in the Drake Passage and a faster-increasing eddy field is the local "
    "expression of the decoupling between transport and mesoscale energy examined in the Discussion; whether the "
    "relevant, circumpolar wind forcing has changed over the record cannot be established from the wind product.")

# ---------------- Results 3.4: new section with Figure 5
INSERT_AFTER[78] = [("Heading2", "3.4 Jet-path distortion and eddy kinetic energy"),
    ("Normal", "The jet path is tracked in every month in all 36 windows (Figure 5a), and its length shows no trend in 34 "
     "of them (the two exceptions, near 175°W, are decreases of about 5 % per decade). When the path within a window "
     "lengthens, EKE within 150 km of it rises in the same month: the correlation reaches p < 0.05 in 20 of the 36 "
     "windows (median r = 0.28, up to 0.62 at 20°W), is largest at zero lag everywhere and decays with lag (Figure 5b). "
     "EKE 15–35° downstream shows no such relation (median r = 0.12; p < 0.05 in one window), and the number of "
     "eddies born along the path per unit length is unrelated to its distortion (p < 0.05 in 4 windows, with both "
     "signs). Along the Polar-Front level the local covariation is weaker (p < 0.05 in 14 windows, median r = 0.16) and "
     "the downstream relation is likewise absent (3 windows, mixed signs). The distortion of the jet path and the "
     "local EKE therefore vary together on the monthly time scale, as expected if the mesoscale eddies are themselves "
     "the distortions of the path, and no delayed response of the eddy field downstream of the meanders is detectable "
     "at monthly resolution."),
    ("Normal", "__FIGURE5__"),
    ("Normal", "Figure 5. Jet-path distortion and eddy kinetic energy around the ACC, 2001–2025 (Product 1). (a) Mean "
     "EKE (logarithmic scale) with the mean jet-core path (black; +0.18 m relative ADT contour, drawn for the central "
     "10° of each window) and the 3000 m isobath (grey). (b) For each 30° window, centred on the longitude shown, the "
     "correlation between the monthly path length and EKE within 150 km of the path in the same window (black) and "
     "15–35° further east (red), at the lag of 0–3 months with the largest |r|, sign retained; filled symbols: p < 0.05 "
     "from phase-randomised surrogates with the maximum over lags.")]
# remove the first-pass meander figure and caption after 85
INSERT_AFTER.pop(85, None)

# ---------------- Discussion: reorder (Drake first), rewrite meander paragraph, add limitations
REPLACE[80] = ("Using two and a half decades of Copernicus Marine Service altimetry, this study provides an observational "
    "assessment of kinetic energy trends, energetic frontal pathways, transport proxies and the relation between jet "
    "distortion and eddy energy in the Southern Ocean. The results relevant to eddy saturation are the contrast between "
    "an increase in EKE and the much weaker change of the surface transport proxies, in the Drake Passage and around "
    "the ACC; the geographical stability of the energetic ACC pathway while the flow along it strengthens; and the "
    "absence of a detectable downstream eddy response to distortions of the jet. This behaviour is consistent with an "
    "eddy-saturated ACC at the surface, although the observational record is too short to determine the equilibrium "
    "sensitivity of circumpolar transport to sustained changes in wind forcing, and altimetry does not resolve the "
    "vertical structure of the response.")
REPLACE[86] = REPLACE[86].replace("The Drake Passage diagnostics (Figure 3) provide the local signature of eddy saturation within the study domain.",
    "The Drake Passage diagnostics (Figure 3) provide the local signature of eddy saturation within the study domain, and the place where altimetry can be compared with an in-situ record.")
REPLACE[87] = REPLACE[87] + (" Around the ACC, the surface transport proxy does increase, by 1.6 % per decade, so the "
    "circumpolar picture is not one of an unchanged mean flow but of a mean flow whose surface signature strengthens "
    "several times more slowly than the eddy field; the Drake Passage section, a single transect, does not resolve a "
    "change of this size.")
# 81, 82, 83 are moved after 87: delete in place and re-insert
MOVED = [REPLACE[81], C_82 := ("This increase is consistent with EKE intensification trends reported in other energetic regions, "
    "including the Gulf Stream and the Kuroshio Extension (Barceló-Llull et al., 2025) and the Mediterranean Sea "
    "(Hargous et al., 2026), suggesting a global-scale eddy field response to changing atmospheric forcing."),
    REPLACE[83]]
REPLACE[81] = None; REPLACE[82] = None; REPLACE[83] = None
MOVED[0] = MOVED[0].replace("Thompson and Garabato (2014)", "Thompson and Naveira Garabato (2014)")
MOVED[2] = MOVED[2].replace(
    "The substantial interannual variability superimposed on these trends is broadly consistent with SAM forcing, in which wind anomalies are preferentially channelled into the eddy field rather than into changes in mean transport (Meredith and Hogg, 2006; Meredith et al., 2012), though a formal attribution analysis is beyond the scope of this study.",
    "The substantial interannual variability superimposed on these trends is not explained by the SAM index over 2001–2025 (Sect. 3.1); the lagged, regional response to a single strong SAM event described by Meredith and Hogg (2006) is not expected to appear as a circumpolar correlation over a 25-year record, and a formal attribution is beyond the scope of this study.")
MEANDER_DISC = ("The jet-path analysis (Figure 5) addresses the first step of the standing-meander mechanism of Thompson and "
    "Naveira Garabato (2014), in which a change in the shape of the meanders precedes the growth of eddy energy "
    "downstream. At monthly resolution the path length and the EKE within 150 km of the path rise together, in the same "
    "month, almost everywhere around the ACC, and EKE downstream shows no delayed response; eddy generation along the "
    "path is unrelated to its distortion. The simultaneous covariation is the signature of the mesoscale eddies "
    "themselves, which distort the contour as they form and pass, rather than of a flexing of the standing meander "
    "that would precede them; the mechanism may operate on time scales shorter than a month, which the monthly means "
    "do not resolve, or in the vertical structure, which altimetry does not see. An earlier version of this analysis "
    "with hand-drawn sections at three meanders suggested a lead of the path over the downstream EKE of one to three "
    "months; that lead does not survive the box-free formulation (Supporting Information S4). The path length itself "
    "shows no trend over 2001–2025 in 34 of 36 windows, so the standing meanders, like the frontal positions, have not "
    "changed shape detectably while the energy around them increased.")
MOVED += [REPLACE[84], MEANDER_DISC]; REPLACE[84] = None; REPLACE[85] = None
INSERT_AFTER[87] = [("Normal", m) for m in MOVED]
REPLACE[88] = ("Several limitations bound these conclusions. Altimetry observes the surface geostrophic flow only: the transport "
    "proxies are surface quantities, the in-situ record covers 0–760 m, and the barotropic component of the ACC "
    "transport (Donohue et al., 2016) and the vertical structure of EKE are not observed, so a stable surface transport "
    "does not establish a stable full-depth transport. The gridded maps resolve wavelengths above roughly 100–150 km "
    "at these latitudes, and the tracked eddies have radii of 40 km and above; the smaller end of the mesoscale band, "
    "where the deformation radius places many Southern Ocean eddies, is not represented. The EKE trend depends on the "
    "satellite constellation (Supporting Information S1), and the tracked-eddy atlas, derived from the same maps, "
    "shares this dependence. And a weak transport trend over 25 years is not unique to an eddy-saturated equilibrium: "
    "Böning et al. (2008) similarly reported little change in ACC isopycnal tilt despite stronger westerlies, while "
    "idealized adjustment studies demonstrate that the ACC and global pycnocline can respond to sustained forcing over "
    "multidecadal to centennial timescales (Allison et al., 2011; Jones et al., 2011; Samelson, 2011). The present "
    "observations therefore cannot determine whether the weak transport response represents a fully equilibrated "
    "eddy-saturated state or, in part, a slowly adjusting circulation. Our result is more narrowly that over the "
    "satellite era the surface transport proxies change by at most a few percent per decade and the position and shape "
    "of the energetic ACC pathway remain statistically stable, while the mean-flow speed and the mesoscale kinetic "
    "energy increase, a pattern that is consistent with eddy-saturation theory.")
REPLACE[89] = ("As the altimetric record lengthens and higher-resolution products become available, including those derived "
    "from the SWOT mission, continued monitoring of Southern Ocean kinetic energy, frontal structure and transport will "
    "be essential for determining whether the observed decoupling between mesoscale energy and surface transport "
    "persists, for resolving the sub-monthly and subsurface steps of the meander mechanism that this analysis cannot "
    "see, and for separating an eddy-saturated response from the intrinsically slow adjustment of the ACC under "
    "projected future increases in Southern Hemisphere westerly winds.")
REPLACE[93] = ("The altimetry sea level products used in this study, SEALEVEL_GLO_PHY_L4_MY_008_047 and "
    "SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057, and the wind product WIND_GLO_PHY_L4_MY_012_006, are freely available from "
    "the Copernicus Marine Service at https://doi.org/10.48670/moi-00148, https://doi.org/10.48670/moi-00145, and "
    "https://doi.org/10.48670/moi-00185, respectively. The Mesoscale Eddy Trajectory Atlas META4.0 DT is distributed by "
    "AVISO+ at https://doi.org/10.24400/527896/a01-2026.001. The in situ Drake Passage transport time series of "
    "Gutierrez-Villanueva et al. (2023) is archived on Zenodo at https://doi.org/10.5281/zenodo.10044261. The SAM index "
    "of Marshall (2003) is available from the British Antarctic Survey.")
REPLACE[95] = REPLACE[95].replace("half-power point, standing-meander and tracked-eddy analyses", "half-power point, jet-path and tracked-eddy analyses")
REPLACE[137] = ("Thompson, A. F. and Naveira Garabato, A. C.: Equilibration of the Antarctic Circumpolar Current by standing "
    "meanders, J. Phys. Oceanogr., 44, 1811–1828, https://doi.org/10.1175/JPO-D-13-0163.1, 2014.")

# ---------------- references
NEW_REFS += [
    "Ballarotta, M., Ubelmann, C., Pujol, M.-I., Taburet, G., Fournier, F., Legeais, J.-F., Faugère, Y., Delepoulle, A., Chelton, D., Dibarboure, G., and Picot, N.: On the resolutions of ocean altimetry maps, Ocean Sci., 15, 1091–1109, https://doi.org/10.5194/os-15-1091-2019, 2019.",
    "Ballarotta, M., Pujol, M.-I., and SL-TAC team: Quality information document for Sea Level TAC C3S products, Copernicus Marine Service, 2024. [check against the product page]",
    "Donohue, K. A., Tracey, K. L., Watts, D. R., Chidichimo, M. P., and Chereskin, T. K.: Mean Antarctic Circumpolar Current transport measured in Drake Passage, Geophys. Res. Lett., 43, 11760–11767, https://doi.org/10.1002/2016GL070319, 2016.",
    "Grassi, B. and Pujol, M.-I.: Product user manual for Sea Level Altimeter products, Copernicus Marine Service, 2025. [check against the product page]",
    "Hall, A. and Visbeck, M.: Synchronous variability in the Southern Hemisphere atmosphere, sea ice, and ocean resulting from the annular mode, J. Climate, 15, 3043–3057, https://doi.org/10.1175/1520-0442(2002)015<3043:SVITSH>2.0.CO;2, 2002.",
    "Hallberg, R. and Gnanadesikan, A.: An exploration of the role of transient eddies in determining the transport of a zonally reentrant current, J. Phys. Oceanogr., 31, 3312–3330, https://doi.org/10.1175/1520-0485(2001)031<3312:AEOTRO>2.0.CO;2, 2001.",
    "Hallberg, R. and Gnanadesikan, A.: The role of eddies in determining the structure and response of the wind-driven Southern Hemisphere overturning: Results from the Modeling Eddies in the Southern Ocean (MESO) project, J. Phys. Oceanogr., 36, 2232–2252, https://doi.org/10.1175/JPO2980.1, 2006.",
    "Saenko, O. A., Fyfe, J. C., and England, M. H.: On the response of the oceanic wind-driven circulation to atmospheric CO2 increase, Clim. Dynam., 25, 415–426, https://doi.org/10.1007/s00382-005-0032-5, 2005.",
    "Straub, D. N.: On the transport and angular momentum balance of channel models of the Antarctic Circumpolar Current, J. Phys. Oceanogr., 23, 776–782, https://doi.org/10.1175/1520-0485(1993)023<0776:OTTAAM>2.0.CO;2, 1993.",
]
