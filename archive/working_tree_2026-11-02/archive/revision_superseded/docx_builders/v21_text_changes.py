# v2.1 text changes. Keys: base paragraph index (see listing). Values: new text (str), None = delete,
# ("after", [texts...]) handled by INSERT_AFTER. Keep edits minimal; numbers from outputs/retest/retest_summary.csv.
REPLACE = {
8: ("Here, we provide an observation-based assessment of long-term trends and regional contrasts in kinetic energy "
    "across the Southern Ocean over 2001–2025. Using gridded altimetry products, we calculate regional KE and EKE "
    "trends associated with the ACC, together with the trend in the along-stream speed of the mean flow, and we "
    "compare an all-satellite with a two-satellite product. To assess eddy saturation, we apply the “half-power point” "
    "method to track the meridional position of the ACC’s high-energy frontal bands, and we compare Drake Passage "
    "transport proxies with the local eddy field. Finally, we test whether three standing meanders of the ACC have "
    "changed and how their shape relates to the downstream eddy field. We find that KE, EKE and the mean-flow speed "
    "all increase, that the energetic frontal band and the meander paths show no detectable displacement, that the "
    "transport proxies show no detectable trend while EKE in the Drake Passage rises, and that, on monthly time scales, "
    "distortions of the meander paths precede eddy energy growth within and downstream of the meanders. These results are consistent with an eddy-saturated "
    "regime, although the record cannot separate it from a slowly adjusting circulation."),

40: ("However, the ACC response is strongly modulated by bathymetry and increased EKE need not arise exclusively "
     "through baroclinic instability. Interaction with major topographic features generates standing meanders and "
     "energetic downstream wakes, where changes in the mean flow can modify both baroclinic and barotropic "
     "instability. Stronger forcing may enhance topographic wakes through barotropic instability (Tansley and "
     "Marshall, 2001), while changes in standing meanders can alter local baroclinic conversion and momentum "
     "transfer toward the deep ocean (Thompson and Garabato, 2014; Stewart et al., 2023). Consistent with this "
     "topographic control, along-track altimetry shows the largest EKE increase of the last decades downstream of "
     "the Campbell Plateau, although regional trends depend strongly on the period considered (Zhang et al., 2021). "
     "Circumpolar EKE increases of the order of 5 % per decade have been reported from gridded altimetry (Hogg et "
     "al., 2015; Martínez-Moreno et al., 2021); here we extend this record to 2025 with the latest reprocessing and "
     "test its dependence on the satellite constellation."),

48: ("Frontal positions are identified from meridional KE profiles using the half-power point method, first "
     "introduced by Chambers (2018) as a sea-level rise-resilient metric to identify trends in Southern Ocean "
     "frontal positions from cross-track altimetric observations. Here, we apply this method to the KE of the "
     "monthly-mean geostrophic velocity from gridded altimetry (Product 1), so that the detected envelopes "
     "correspond to the persistent ACC jets rather than to individual eddies. We focus on those fronts persistent "
     "enough to appear consistently in monthly gridded fields since these are the structures more relevant for "
     "detecting multi-decadal trends in frontal position and ACC dynamics."),

51: ("where φS and φN are the southern and northern envelope boundaries, respectively. Monthly envelopes are "
     "assigned to the envelope detected at the same longitude in the 2001–2025 mean KE field, and only the 900 "
     "longitudes of the ACC band where the matched envelope is found in at least 90 % of months are retained, so "
     "that the circumpolar and sector statistics are computed over a fixed set of fronts. Sensitivity tests with "
     "thresholds of 150–400 cm² s⁻², peak fractions of 20–50 % and minimum widths of 50–300 km give the same "
     "conclusions (Supporting Information S3). Synthetic experiments in which a known shift is imposed on the "
     "observed climatological profiles show that a 25 km displacement over 25 years is detected in every "
     "realisation and a 10 km displacement in 44 % of them, with about 70 % of the imposed shift recovered."),

53: ("Surface transport across the Drake Passage is estimated from gridded altimetry (Product 1) as the "
     "geostrophic velocity component normal to a transect from 65°W, 55°S to 60°W, 63°S, integrated along the "
     "section over the points that are free of sea ice in every month of the record (53 of 121 points, north of "
     "58.5°S). A fixed set of points is required because a section whose length varies with the seasonal ice "
     "cover produces a spurious trend. No depth integration is applied, so this quantity is a surface transport "
     "proxy rather than a full-depth volume transport. It is therefore expressed as a standardized anomaly "
     "(z-score) and used only to compare temporal variability and trends. Displacing the transect end points by "
     "up to 1° does not change the conclusions."),

54: ("As an independent comparison, we use the Drake Passage transport time series of Gutierrez-Villanueva et al. "
     "(2023), the total transport over 0–760 m from 244 repeat ADCP transects between 2005 and 2019, which we "
     "average into annual means. Throughout the manuscript we therefore refer to the two quantities collectively "
     "as Drake Passage transport proxies, while recognizing that they measure different aspects of the "
     "circulation. Wind forcing is characterised by the eastward wind stress of Product 2 averaged over the Drake "
     "Passage box (80–52°W, 63–54°S) and over the circumpolar band 45–60°S."),

57: ("Following previous studies (Barceló-Llull et al., 2025; Hargous et al., 2026), trends are estimated using "
     "the Theil-Sen estimator after removing the monthly climatology, chosen because of its robustness to "
     "outliers. The Theil-Sen slope is defined as the median of all pairwise slopes between data points, making it "
     "insensitive to the influence of individual extreme values (Theil, 1950; Sen, 1968). Because the series are "
     "strongly autocorrelated, confidence intervals and p-values are obtained from a parametric bootstrap: an "
     "autoregressive model is fitted to the detrended residuals, its coefficients are corrected for small-sample "
     "bias (Kilian, 1998), and the trend statistic is compared with its distribution over 300 synthetic series "
     "generated from the corrected model and processed identically. For monthly series the test is run with an "
     "AIC-selected order and with a fixed order of 12, and the larger p-value is reported. In simulations with "
     "trend-free autocorrelated series this procedure rejects at or below the nominal 5 % level, whereas the "
     "Mann–Kendall test with the Yue and Wang (2004) correction, used in an earlier version of this work, rejected "
     "far too often. A single significance level, p < 0.05, is used throughout; in the figures, * denotes "
     "p < 0.05 and ** denotes p < 0.01. For the trend maps, the same statistic is computed on 1° cells and the "
     "false discovery rate is controlled at 10 % (Benjamini and Hochberg, 1995; Wilks, 2016). All trends are "
     "computed for the all-satellite product (Product 1) and for the two-satellite product (Product 4); where they "
     "differ, both are reported."),

61: ("Figure 1. Southern Ocean KE and EKE over the 2001-2025 period. Mean (a) KE and (b) EKE derived from gridded "
     "altimetry (Product 1). Trends in (c) KE and (d) EKE on 1° cells; dots mark cells whose trend remains "
     "significant after false-discovery-rate control. In (a-d), grey contours indicate the −0.6 and +0.2 m mean "
     "absolute dynamic topography isolines delimiting the ACC band, dashed lines the sector boundaries, and the "
     "blue area the sea-ice mask. (e) Area-averaged time series of KE (red) and EKE (green) for the ACC band, "
     "showing monthly values (thin), 12-month running means (solid), and trends (dashed lines); trend magnitudes "
     "and significance are annotated. (f) Regional KE and EKE trends in percent of the regional mean per decade "
     "for the ACC band and its Indian, Pacific and Atlantic sectors, for the all-satellite (circles) and the "
     "two-satellite (squares) products; bars denote the 95 % confidence interval, filled symbols trends with "
     "p < 0.05; values are listed on the right (*, p < 0.05; **, p < 0.01)."),

63: ("Trend maps show a widespread increase in both KE and EKE across the Southern Ocean (Figure 1c,d). The "
     "strongest positive trends are located along the ACC core, particularly in the Indian sector downstream of "
     "the Kerguelen Plateau and along the Southeast Indian Ridge, and in the Pacific sector south of New Zealand, "
     "downstream of the Campbell Plateau. Negative or near-zero trends are comparatively limited and generally "
     "located outside the main ACC envelope. At the scale of individual 1° cells, however, few local trends can be "
     "distinguished from the interannual variability: only 2.9 % (KE) and 4.0 % (EKE) of the cells reach p < 0.05, "
     "no more than expected by chance, and 5 (KE) and 18 (EKE) cells remain significant after false-discovery-rate "
     "control. The field-wide signal therefore rests on the area averages below, not on local trends."),

64: ("The ACC-band time series of KE and EKE (Figure 1e) show positive trends over the full record, with KE "
     "increasing at 1.36 cm² s⁻² yr⁻¹ (95 % confidence interval 0.86–1.86; 4.9 % per decade; p = 0.003) and EKE "
     "at 1.98 cm² s⁻² yr⁻¹ (0.81–3.15; 8.0 % per decade; p = 0.017). The along-stream speed of the monthly-mean "
     "flow increases at the same time, by 2.0 % per decade (p = 0.003), and the KE of the annual-mean velocity by "
     "3.5 % per decade (p = 0.017), so the intensification involves both the eddy field and the mean flow. The "
     "difference between the relative EKE and KE trends, +3.0 % per decade (−0.1 to +6.1), does not reach "
     "significance (p = 0.06). Both series also show substantial interannual variability, with fluctuations at "
     "periods consistent with the SAM."),

65: ("Regional decomposition of the trends (Figure 1f) shows positive KE and EKE trends in all three sectors, "
     "significant for both quantities in the Pacific and Indian sectors and for KE in the Atlantic sector, where "
     "the EKE trend has a wide confidence interval (p = 0.053). The two-satellite product (Product 4), which "
     "merges a constant pair of altimeters, reproduces the sign of every trend and the KE trends to within their "
     "uncertainty, but its EKE trends are about half as large (ACC band: KE 3.6 % and EKE 3.8 % per decade, both "
     "p < 0.05), and its sector EKE trends are significant only in the Pacific. The divergence between the two "
     "products grows after about 2018 as the number of merged missions increases (Supporting Information S1). "
     "The mean-flow strengthening is therefore robust to the product, whereas the magnitude of the EKE increase, "
     "and in particular the excess of the EKE trend over the KE trend, is not."),

68: ("Figure 2. Half-power point frontal analysis over the 2001-2025 period. (a) Half-power points detected in "
     "January 2001 within the ACC band (red), over the monthly KE field (grey shading); the blue dashed line marks "
     "90°E. (b) Meridional KE profile at 90°E in January 2001, showing the KE envelope (shaded), peak KE value "
     "(triangle), the half-power point position (star), and the 200 cm² s⁻² detection threshold (dashed line). "
     "(c) Circumpolar-mean half-power point latitude over the fixed set of 900 longitudes, with 12-month running "
     "mean and trend. (d) Trends in half-power point latitude for the ACC band and by sector, in km per decade "
     "(positive northward), with 95 % confidence intervals; open symbols denote trends with p ≥ 0.05."),

69: ("The half-power point method identifies coherent frontal envelopes across the Southern Ocean at monthly "
     "resolution (Figure 2a). The detected positions trace the main ACC frontal pathway, closely following the "
     "topographically steered meanders; detection is least frequent, 25–45 % of months, between 170°W and 150°W "
     "and between 115°W and 75°W, where the band is weak and broad. An illustrative meridional KE profile at 90°E "
     "(Figure 2b) demonstrates the method: the KE envelope (shaded) is defined around the peak value (triangle), "
     "with the half-power point (star) marking the energy-weighted centre of the front near 49°S."),

70: ("The circumpolar-mean half-power point position (Figure 2c) fluctuates around 50.8°S without a detectable "
     "trend over the full record (−0.5 km per decade, 95 % confidence interval −4.2 to +3.1 km per decade; "
     "p = 0.74). Given the detection experiments described in Section 2.2, a sustained shift larger than about "
     "25 km over the record can be excluded. This stationarity contrasts with the southward shifts of about 60 km "
     "over 1993–2008 reported from SSH-contour methods (Sokolov and Rintoul, 2009; Kim and Orsi, 2014), and is "
     "consistent with the known poleward bias of SSH-contour methods under sea-level rise (Gille, 2014) and with "
     "the along-track half-power analysis of Chambers (2018)."),

71: None,

72: ("None of the three sectors shows a significant trend either (Figure 2d): −2.2 km per decade in the Indian, "
     "−1.2 in the Pacific and +2.0 in the Atlantic sector, all with confidence intervals of about ±7 km per decade. "
     "Trends resolved in 5° longitude bins (Supporting Information S3) are likewise within their uncertainty, "
     "including the regions between 90°E and 170°E and between 160°W and 90°W where Chambers (2018) found "
     "southward and northward shifts, respectively, over 1993–2016. The energetic ACC pathway has therefore "
     "remained geographically anchored over 2001–2025 at the circumpolar, sector and regional scales. A previous "
     "version of this analysis reported a small broadening of the frontal envelopes; this signal is reproduced by "
     "an increase of the KE amplitude alone when the envelope width is defined by a fixed threshold, and is not "
     "present in amplitude-independent width measures, so no change in envelope width can be inferred."),

75: ("Figure 3. Eddy saturation diagnostics in the Drake Passage over the 2001–2025 period. (a) Map of the region "
     "showing the transport transect (solid: ice-free points used for the proxy; dashed: seasonally ice-covered "
     "part, excluded) and the domain used for the wind stress, KE and EKE time series. (b) Standardized trends of "
     "the altimetric surface transport proxy (2001–2025) and of the in-situ transport (annual means, 2005–2019), "
     "in standard deviations per decade, with 95 % confidence intervals; open symbols denote p ≥ 0.05. (c) Trends "
     "of EKE and KE in the domain, in percent of the mean per decade. (d) Standardized transport anomalies from "
     "altimetry (monthly values shaded, 12-month running mean) and from the in-situ transects (running mean over "
     "12 transects), with the eastward wind stress over the domain (orange, 12-month running mean, right axis). "
     "(e) EKE in the domain: monthly values, 12-month running mean and trend."),

76: ("The Drake Passage analysis focuses on a transect across the passage and a surrounding domain used for the "
     "wind stress, KE and EKE diagnostics (Figure 3a). The standardized trends in the two Drake Passage transport "
     "proxies (Figure 3b) are small and statistically indistinguishable from zero (altimetry: +0.016 σ yr⁻¹, "
     "p = 0.40; in-situ: −0.078 σ yr⁻¹, p = 0.47), indicating that neither estimate shows a detectable change in "
     "transport over its record. The two proxies are, however, not correlated with each other at monthly or "
     "annual time scales (r = 0.01 and −0.04 over 2005–2019), so each must be read as an independent constraint. "
     "In contrast, EKE over the domain increases by 2.18 cm² s⁻² yr⁻¹ (8.4 % per decade; p = 0.050), while the "
     "KE trend is smaller and not significant (+0.77 cm² s⁻² yr⁻¹, 2.4 % per decade; p = 0.30); the difference "
     "between the relative EKE and KE trends, +6.0 % per decade, is significant (p = 0.04). The EKE trend is "
     "borderline, and it becomes p = 0.08 when the sea-ice mask is not applied, so it should be read as a "
     "likely, not an established, increase."),

77: ("The transport and wind time series (Figure 3d) show that the eastward wind stress over the Drake Passage "
     "exhibits substantial interannual variability around 0.13 N m⁻² and a positive trend of 8.7 % per decade "
     "(p = 0.03), whereas over the circumpolar band 45–60°S the trend is 2.5 % per decade and not significant "
     "(p = 0.08). The two transport proxies show variability of similar amplitude on interannual timescales, "
     "without a common signal."),

78: ("The EKE time series (Figure 3e) shows this intensification directly, with the running mean rising most "
     "markedly from around 2016 onward. The combination of transport proxies without a detectable trend and an "
     "increasing eddy field within the passage provides the local expression of the decoupling between transport "
     "and mesoscale energy examined in the Discussion; whether the relevant, circumpolar wind forcing has "
     "changed over the record cannot be established from the wind product."),

80: ("Using two and a half decades of Copernicus Marine Service altimetry, this study provides an observational "
     "assessment of kinetic energy trends, energetic frontal pathways, standing meanders and transport proxies in "
     "the Southern Ocean. The clearest results relevant to eddy saturation are the contrast between an increase "
     "in EKE and the absence of detectable trends in two independent Drake Passage transport proxies, and the "
     "geographical stability of the energetic ACC pathway and of its standing meanders while the flow through "
     "them strengthens. This behaviour is consistent with an eddy-saturated ACC, although the observational "
     "record is too short to determine the equilibrium sensitivity of circumpolar transport to sustained changes "
     "in wind forcing."),

81: ("The widespread increase in both KE and EKE across the Southern Ocean (Figure 1) is concentrated along the "
     "ACC, in agreement with earlier gridded and along-track analyses (Hogg et al., 2015; Martínez-Moreno et "
     "al., 2021; Zhang et al., 2021). In relative terms, however, the EKE increase is not larger downstream of the "
     "four main topographic features (Kerguelen Plateau, Southeast Indian Ridge–Macquarie Ridge, "
     "Pacific–Antarctic Ridge and Drake Passage) than over the rest of the ACC band: the boxes downstream of these "
     "features contain 43 % of the mean EKE of the band and 45 % of its EKE increase, and the difference between "
     "their relative trend and that of the rest of the band is +0.1 % per decade (p = 0.94) in the all-satellite "
     "and −0.3 % per decade (p = 0.80) in the two-satellite product. The localisation of the eddy response at "
     "standing meanders expected from the mechanism of Thompson and Garabato (2014) is thus not visible in the "
     "surface trends over 2001–2025. At the same time, intensified interaction of the current with topographic "
     "obstacles may alter turbulent downstream wakes and permit enhanced barotropic instability, as illustrated "
     "in idealized ACC-relevant configurations by Tansley and Marshall (2001). Our altimetric diagnostics do not "
     "distinguish among these mechanisms, and the observed EKE trends should therefore be interpreted as evidence "
     "of increased mesoscale activity rather than as a direct diagnosis of a particular instability pathway."),

83: ("The relative EKE trend exceeds the relative KE trend in the all-satellite product (8.0 % against 4.9 % per "
     "decade), but the difference is not significant, and it vanishes in the two-satellite product, whose EKE and "
     "KE trends are both close to 3.7 % per decade. The along-stream speed of the mean flow increases by 2 % per "
     "decade in both products. The surface record therefore shows an intensification of the mean flow and of the "
     "eddy field together, rather than an eddy increase at the expense of the mean flow, and in the tracked-eddy atlas it takes the form of larger amplitudes and faster rotation of an unchanged population of eddies (Figure 4); the larger EKE trend of "
     "the all-satellite product is at least partly a consequence of the growing constellation (Hargous et al., "
     "2026; Supporting Information S1). The substantial interannual variability superimposed on these trends is "
     "broadly consistent with SAM forcing, in which wind anomalies are preferentially channelled into the eddy "
     "field rather than into changes in mean transport (Meredith and Hogg, 2006; Meredith et al., 2012), though a "
     "formal attribution analysis is beyond the scope of this study."),

84: ("The half-power point analysis (Figure 2) provides a sea-level-rise-resilient perspective on frontal "
     "stability. The circumpolar-mean frontal position shows no detectable trend over the record, and neither do "
     "the sector or regional positions, in contrast with the southward shifts reported by SSH-contour methods "
     "over comparable periods. As demonstrated by Gille (2014), SSH contours shift poleward under sea-level rise "
     "even when the underlying current structure is stationary; the half-power point method is immune to this "
     "bias, indicating that the ACC's energetic frontal band has remained geographically anchored over the study "
     "period, with an upper bound of about 25 km on any sustained circumpolar displacement. This result agrees "
     "with the along-track analysis of Chambers (2018) and with the topographic anchoring of the ACC jets in "
     "climate projections (Schiavon et al., 2025)."),

85: ("The standing-meander analysis (Figure 5) extends this picture to the scale of individual meanders. At the "
     "Pacific–Antarctic Ridge, the Southeast Indian Ridge and the Kerguelen Plateau, the path of the jet core "
     "shows no detectable change in position or in month-to-month variability over 2001–2025, while the "
     "along-stream speed in the meander sections increases (by 5.6 % per decade at the Pacific–Antarctic Ridge, "
     "p = 0.02), in agreement with the acceleration reported by Liu et al. (2026). On monthly time scales the "
     "meanders behave as sites of eddy generation. At the Southeast Indian Ridge and at Kerguelen, a lengthening "
     "of the jet path within the meander precedes, by one to two months, an increase of EKE within the meander "
     "section (r = 0.23 and 0.33; p = 0.014 and 0.002), and EKE within the section precedes EKE in the downstream "
     "box by a further one to three months (r = 0.4–0.6); the reverse order, EKE leading the path length, shows "
     "no relation. The pre-specified test between path length and downstream EKE (pooled p = 0.007) reflects "
     "this sequence, and once the meander-section EKE is accounted for, a direct relation to downstream EKE "
     "remains only at the Southeast Indian Ridge. At the Pacific–Antarctic Ridge the relation is weak and has no "
     "preferred direction. Two features limit the interpretation. The signal is carried by the length of the "
     "path, that is by distortions of the jet at scales of a few hundred kilometres, and not by the "
     "crest-to-trough amplitude of the meander, which shows no relation to EKE; the growth of meander lobes "
     "before eddies are shed, rather than a flexing of the standing-meander envelope, is therefore the more "
     "likely reading. And neither the upstream speed nor the circumpolar wind stress precedes the lengthening "
     "(except for a weak association with the wind stress at the Pacific–Antarctic Ridge, r = 0.23, p = 0.03), "
     "so the first step of the sequence proposed by Thompson and Garabato (2014), from forcing to meander shape, "
     "is not visible at the surface. Such topographic constraint is relevant to theories in which standing "
     "meanders play an important role in the ACC momentum balance, but frontal stationarity itself is not a "
     "diagnostic requirement for eddy saturation."),

86: ("The Drake Passage diagnostics (Figure 3) provide the local signature of eddy saturation within the study "
     "domain. Neither the altimetric nor the in-situ transport estimate shows a significant trend over its "
     "record, while altimetry-derived EKE increases at the limit of significance. The absence of a transport "
     "trend provides independent evidence that the large-scale transport signal has changed little over the "
     "available record; the confidence intervals bound the altimetric surface-transport trend to within about "
     "±5 % per decade. Together with increasing EKE, this behaviour is consistent with the weak transport "
     "sensitivity expected under eddy saturation. Gutierrez-Villanueva et al. (2023) showed that net transport "
     "stability across the passage is maintained through compensating trends among the different ACC frontal "
     "jets, which are much reduced in streamwise coordinates, a physical pathway through which the ACC can absorb "
     "additional wind energy by redistributing momentum across frontal zones via enhanced mesoscale activity, "
     "without any change in total volume transport. This decoupling between stable transport and a positive EKE "
     "trend is the hallmark of an eddy-saturated regime (Constantinou and Hogg, 2019; Maddison et al., 2025)."),

87: ("Within the passage, the relative EKE trend exceeds the relative KE trend by 6 % per decade (p = 0.04), a "
     "contrast that is not found over the ACC band as a whole. As discussed in the context of the "
     "satellite-constellation sensitivity analysis (Supporting Information S1), gridded altimetry may inflate "
     "EKE trend magnitudes as the number of merged altimeters grows (Hargous et al., 2026); the qualitative "
     "decoupling between stable transport proxies and increasing EKE is nonetheless present in both products. "
     "Two caveats apply. The local wind stress over the passage increases over the record, so a flat local "
     "transport is not by itself evidence of insensitivity to wind; and the two transport proxies do not "
     "co-vary, so their agreement is limited to the absence of a trend."),

88: ("An important limitation is that a weak transport trend over 25 years is not unique to an eddy-saturated "
     "equilibrium. Böning et al. (2008) similarly reported little change in ACC isopycnal tilt despite stronger "
     "westerlies, while idealized adjustment studies demonstrate that the ACC and global pycnocline can respond to "
     "sustained forcing over multidecadal to centennial timescales (Allison et al., 2011; Jones et al., 2011; "
     "Samelson, 2011). The present observations therefore cannot determine whether the weak transport response "
     "represents a fully equilibrated eddy-saturated state or, in part, a slowly adjusting circulation. Our "
     "result is more narrowly that over the satellite era the observable transport proxies and the position of "
     "the energetic ACC pathway remain statistically stable while the mean flow and the mesoscale kinetic energy "
     "increase, a pattern that is consistent with eddy-saturation theory."),

95: ("CMS designed the study, carried out the kinetic energy, eddy kinetic energy, half-power point and "
     "standing-meander analyses, and prepared the manuscript with contributions from all co-authors. SR, VC and "
     "AP contributed to the conceptualization and methodology and supervised the study. BBL contributed to the "
     "interpretation of eddy kinetic energy trends. MB provided expertise on the altimetry product "
     "characteristics and satellite constellation effects. All authors reviewed and edited the manuscript."),

100: ("Generative AI assisted the organization of the code repository accompanying this study, including folder "
      "structure, script organization, and documentation, and the implementation and testing of the trend "
      "statistics. All scientific content, analysis and interpretation are the authors' own."),

107: None,   # Chelton (uncited)
117: None,   # Kendall
121: None,   # Mann
138: None,   # Yue and Wang -> kept? cited in Methods; keep
}
REPLACE.pop(138)
REPLACE[133] = ("Sokolov, S. and Rintoul, S. R.: Circumpolar structure and distribution of the Antarctic Circumpolar Current fronts: 2. Variability and relationship to sea surface height, J. Geophys. Res.-Oceans, 114, C11019, https://doi.org/10.1029/2008JC005248, 2009.")

# new paragraphs inserted AFTER the given base index (in order)
REPLACE[55] = "2.6 Time series and linear trends"

EDDY_METHODS = ("Mesoscale eddies are taken from the Mesoscale Eddy Trajectory Atlas (META4.0 DT allsat; Product 5), in which "
    "closed sea-level contours are detected daily in the all-satellite absolute dynamic topography and linked into "
    "tracks following Pegliasco et al. (2022). The atlas ends in September 2023, and we use its observations south of "
    "35°S from January 2001 to August 2023. For each month we compute, over the same ACC band as the KE and EKE series, "
    "the number of eddy observations per day (eddies present), the number of tracks that begin in the month (eddy "
    "generation; tracked eddies only), the mean amplitude, the mean speed radius, the mean rotational speed along the "
    "contour of maximum speed, and the sum of the squared amplitudes per day, which we use as a proxy of the energy "
    "contained in the identified eddies. The density of track origins on 1° cells is used to map the generation sites. "
    "Trends in these series are estimated as in Sect. 2.6 and compared with the trends of the gridded EKE and KE over "
    "the same months.")

INSERT_AFTER = {
65: [("Normal", "The tracked-eddy atlas (Product 5) shows which property of the eddy field carries this increase (Figure 4). "
     "In the ACC band about 345 eddies are present on an average day and about 330 tracks begin each month, with a mean "
     "amplitude of 7.7 cm and a mean speed radius of 41 km. Over January 2001 to August 2023 neither the number of eddies "
     "present (+0.8 % per decade; p = 0.37), the number of new tracks (−0.4 %; p = 0.63) nor the radius (−0.1 %; p = 0.60) "
     "shows a trend, whereas the mean amplitude increases by 4.6 % per decade (2.4–6.8; p = 0.007), the rotational speed "
     "by 3.8 % (0.9–6.8; p = 0.02) and the sum of the squared amplitudes by 10.9 % (4.8–17.0; p = 0.007); the gridded "
     "EKE increases by 8.6 % per decade (2.1–15.1; p = 0.027) over the same months. The EKE increase thus reflects an "
     "intensification of the existing eddies rather than a change in their number or size. Eddy generation is not "
     "concentrated in the wakes of the standing meanders (Supporting Information S4), and the trends of eddy number and "
     "amplitude in the wakes do not differ from those in the rest of the band (p > 0.2)."),
     ("Normal", "__FIGURE4__"),
     ("Normal", "Figure 4. Tracked mesoscale eddies in the ACC band, January 2001 to August 2023 (Product 5). (a) 12-month "
     "running means of the number of eddies present per day, their mean amplitude, the sum of their squared amplitudes "
     "per day and the gridded EKE (Product 1), as anomalies in percent of the 2001–2023 mean. (b) Trends over the same "
     "months in percent per decade with 95 % confidence intervals (filled symbols: p < 0.05); the gridded EKE and KE "
     "are shown for comparison. The map of eddy generation sites is given in Supporting Information S4.")],
54: [("Heading2", "2.4 Standing meanders"),
     ("Normal", "Three standing meanders with a persistent downstream EKE maximum are analysed: the Pacific–Antarctic "
      "Ridge meander (150–80°W), with the sections of Liu et al. (2026); the Southeast Indian Ridge meander near "
      "150°E (upstream 125–138°E, meander 140–150°E, downstream 150–158°E); and the Kerguelen Plateau meander at "
      "the Fawn Trough (upstream 70–77°E, meander 76–84°E, downstream 84–94°E). The path of the jet core is tracked "
      "each month as the contour of absolute dynamic topography relative to its regional mean, which removes the "
      "effect of sea-level rise (Gille, 2014), at the level of the climatological jet core; it is found in every "
      "month at the three meanders. From the path we compute its latitude by longitude, its within-year "
      "variability, its length within the meander section after smoothing at 150 km, and its crest-to-trough "
      "amplitude; from the velocity field, the along-stream speed upstream of and within the meander section and "
      "the EKE within the meander section and in the downstream box. Following Thompson and Garabato (2014), the "
      "relation between the path geometry and EKE is tested on deseasonalised, detrended monthly series smoothed "
      "over two months, at lags of 0–3 months fixed in advance, with significance from phase-randomised surrogates "
      "and the three meanders combined with a Stouffer test; the reverse lags (EKE leading) and the partial "
      "correlation with downstream EKE given the meander-section EKE are used to establish the order of the "
      "sequence."),
     ("Heading2", "2.5 Tracked eddies"),
     ("Normal", EDDY_METHODS)],
85: [("Normal", "__FIGURE5__"),
     ("Normal", "Figure 5. Standing meanders and eddy kinetic energy, 2001–2025 (Product 1). (a, c, e) Mean EKE "
      "(logarithmic scale) at the Pacific–Antarctic Ridge, the Southeast Indian Ridge (150°E) and the Kerguelen "
      "Plateau, with the mean jet-core path (black), the meander section (dashed) and the downstream EKE box "
      "(red); grey contour: 3000 m isobath. (b, d, f) Difference in EKE anomaly within the meander section (dashed) "
      "and in the downstream box (solid), in percent of the respective mean, between months with path length in "
      "the upper and in the lower tercile, as a function of lag (positive: EKE after path length). Grey band: "
      "95 % range from 1000 surrogate path-length series for the downstream box; light shading: the lags 0–3 "
      "months fixed before the analysis, with the p-value of the pre-specified test (path length to downstream "
      "EKE) over that window. (g) Lag correlation between path length and downstream EKE for each meander and "
      "their mean (solid), the mean for the meander-section EKE (dashed), and the 95 % surrogate range of the "
      "mean.")],
}

# new references (Copernicus style), inserted alphabetically
NEW_REFS = [
    'Pegliasco, C., Delepoulle, A., Mason, E., Morrow, R., Faugère, Y., and Dibarboure, G.: META3.1exp: a new global mesoscale eddy trajectory atlas derived from altimetry, Earth Syst. Sci. Data, 14, 1087–1107, https://doi.org/10.5194/essd-14-1087-2022, 2022.',
"Benjamini, Y. and Hochberg, Y.: Controlling the false discovery rate: a practical and powerful approach to multiple testing, J. R. Stat. Soc. B, 57, 289–300, https://doi.org/10.1111/j.2517-6161.1995.tb02031.x, 1995.",
"Hogg, A. M., Meredith, M. P., Chambers, D. P., Abrahamsen, E. P., Hughes, C. W., and Morrison, A. K.: Recent trends in the Southern Ocean eddy field, J. Geophys. Res.-Oceans, 120, 257–267, https://doi.org/10.1002/2014JC010470, 2015.",
"Kilian, L.: Small-sample confidence intervals for impulse response functions, Rev. Econ. Stat., 80, 218–230, https://doi.org/10.1162/003465398557465, 1998.",
"Martínez-Moreno, J., Hogg, A. M., England, M. H., Constantinou, N. C., Kiss, A. E., and Morrison, A. K.: Global changes in oceanic mesoscale currents over the satellite altimetry record, Nat. Clim. Change, 11, 397–403, https://doi.org/10.1038/s41558-021-01006-9, 2021.",
"Wilks, D. S.: “The stippling shows statistically significant grid points”: How research results are routinely overstated and overinterpreted, and what to do about it, B. Am. Meteorol. Soc., 97, 2263–2273, https://doi.org/10.1175/BAMS-D-15-00267.1, 2016.",
]

# Product 5 row for Table 1 (inserted after the row of Product 4, tracked)
NEW_TABLE_ROW = ["5", "META4.0 DT allsat – Mesoscale Eddy Trajectory Atlas (AVISO+)",
                 "https://doi.org/10.24400/527896/a01-2026.001", ["Pegliasco et al. (2022)"]]
