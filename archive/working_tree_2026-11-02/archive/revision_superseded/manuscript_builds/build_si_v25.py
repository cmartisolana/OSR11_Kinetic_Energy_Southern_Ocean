"""Build the Supporting Information (v2.5) docx from the v1.0 SI template (styles) and outputs/figures_final/FigureS*.png.
Usage: python build_si.py <SI_v1.docx> <out.docx>"""
import sys, copy
from pathlib import Path
import docx
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

HERE = Path(__file__).resolve().parent
FIG = HERE / "outputs" / "figures_final"
src, out = sys.argv[1], sys.argv[2]
d = docx.Document(src)
body = d.element.body
for el in list(body):
    if el.tag.endswith("}sectPr"):
        continue
    body.remove(el)


def P(text, style="normal", bold_prefix=None):
    p = d.add_paragraph(style=style)
    if bold_prefix:
        r = p.add_run(bold_prefix); r.bold = True
    p.add_run(text)
    return p


def H(text):
    return d.add_paragraph(text, style="Heading 1")


def FIGP(name, width=6.6):
    p = d.add_paragraph(style="normal"); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(FIG / f"{name}.png"), width=Inches(width))


def TABLE(rows, header, widths=None):
    t = d.add_table(rows=1, cols=len(header)); t.style = "Table Grid" if "Table Grid" in [s.name for s in d.styles] else t.style
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]; c.text = ""; r = c.paragraphs[0].add_run(h); r.bold = True; r.font.size = Pt(9)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""; r = cells[i].paragraphs[0].add_run(v); r.font.size = Pt(9)
    return t


tp = d.add_paragraph(style="normal")
_r = tp.add_run("Supporting Information for “Kinetic Energy Variability and Eddy Saturation in the Southern Ocean”"); _r.bold = True; _r.font.size = Pt(14)
P("Contents: S1 Results from the two-satellite product; S2 Calibration of the trend test; S3 Half-power analysis: "
  "sensitivity, longitude bins, detection and the two-satellite product; S4 Tracked eddies: generation sites and "
  "nonlinearity; S5 Month-to-month covariation of path length and EKE. Products are numbered as in Table 1 of the main text. "
  "All trends are for 2001–2025 unless stated, with the test of Sect. 2.6 of the main text (95 % confidence intervals; "
  "significance at p < 0.05).")

# ---------------------------------------------------------------- S1
H("S1 Results from the two-satellite product")
P("The all-satellite maps (Product 1) merge every altimeter available at each time, a number that varies over the record and "
  "reaches seven in recent years (International Altimetry Team, 2021). A denser constellation resolves more of the mesoscale "
  "variability, so part of an EKE trend computed from these maps can come from the changing sampling rather than from "
  "the ocean (Barceló-Llull et al., 2025; Hargous et al., 2026). Product 4 merges a constant pair of altimeters on a 1/4° "
  "grid and is designed for climate trends (Ballarotta et al., 2024). Every diagnostic of the main text was recomputed "
  "with Product 4 on its own grid, with the same ACC band, sectors and trend test.")
P("Mean KE has the same pattern as mean EKE (Fig. S1a; main text Fig. 1a), with maxima along the ACC downstream of the "
  "Kerguelen Plateau, the Southeast Indian Ridge and the Campbell Plateau and in the Agulhas and Brazil–Malvinas regions. "
  "Over the ACC band the KE series of the two products follow each other closely; the EKE series depart from each other "
  "in some periods (e.g. 2006–2008), and from about 2018 onwards the all-satellite EKE stays above the two-satellite EKE "
  "(Fig. S1b).")
P("Every ACC-band trend has the same sign in both products (Fig. S1c). The KE of the monthly-mean flow increases by 4.9 % "
  "per decade (all-satellite) and 3.6 % (two-satellite), the along-stream speed by 2.0 % and 1.8 %, and the surface "
  "transport proxy by 1.6 % (0.9–2.3) and 1.5 % (0.9–2.1; p = 0.003 in both): the strengthening of the mean flow does not "
  "depend on the constellation. The EKE trend does: 8.0 % per decade (3.3–12.7; p = 0.017) in the all-satellite and 3.8 % "
  "(1.3–6.3; p = 0.020) in the two-satellite product. The difference lies mainly in the sub-monthly part of EKE, the "
  "variance within each month (12.0 % per decade, p = 0.066, against 5.0 %, p = 0.040), which is the part most sensitive to "
  "the mapping; the low-frequency part, the variance of the monthly means about the 2001–2025 mean, increases by 5.6 % and "
  "3.2 % per decade. The EKE trend does not significantly exceed the KE trend in either product (difference of the "
  "relative trends +3.0 % per decade, −0.1 to +6.1, p = 0.06, all-satellite; +0.2 %, p = 0.74, two-satellite). The surface "
  "transport proxy therefore strengthens about five times more slowly than EKE in the all-satellite product and about "
  "2.5 times more slowly in the two-satellite product.")
P("By sector (Fig. S1d), the KE, along-stream speed and transport-proxy trends agree between the products within their "
  "uncertainty. The EKE trends are smaller in the two-satellite product in every sector and remain significant only in the "
  "Pacific (5.7 % per decade); in the Indian sector the two-satellite EKE and KE trends are not significant. Over the box "
  "south of 35° S of Martínez-Moreno et al. (2021) EKE increases by 6.6 % per decade (p = 0.027) in the all-satellite and "
  "2.3 % (p = 0.043) in the two-satellite product.")
P("The localisation test of the main text (Sect. 2.3) gives the same result in both products: the difference between the "
  "relative EKE trend in the four wake boxes and that in the rest of the band is +0.1 % per decade (−3.1 to +3.3; p = 0.95) "
  "in the all-satellite and −0.3 % (−3.5 to +2.8; p = 0.80) in the two-satellite product, and none of the KE, along-stream "
  "speed or EKE-part differences is significant (p > 0.25). At the Pacific–Antarctic Ridge the jet-core speed on the "
  "sections of Liu et al. (2026) increases by 0.019 m s⁻¹ per decade in the all-satellite and 0.017 m s⁻¹ per decade in the "
  "two-satellite product over the whole meander (p = 0.003 in both), and by 0.026 and 0.023 m s⁻¹ per decade at the Lower "
  "Ridge.")
FIGP("FigureS1")
P(" Results from the all-satellite (Product 1) and two-satellite (Product 4) products. (a) Mean KE, 2001–2025, Product 1 "
  "(light blue: sea ice in any month). (b) ACC-band KE and EKE, 12-month running means in percent of the 2001–2025 mean; "
  "solid: Product 1, dashed: Product 4. (c) ACC-band trends, 2001–2025, in percent of the mean per decade with 95 % "
  "confidence intervals: KE, EKE, its low-frequency and sub-monthly parts, along-stream speed and surface transport proxy. "
  "(d) The same for the sectors (speed: along-stream speed; transport: surface transport proxy). Filled symbols: p < 0.05; "
  "open symbols: p ≥ 0.05.", bold_prefix="Figure S1.")

# ---------------------------------------------------------------- S2
H("S2 Calibration of the trend test")
P("Monthly ocean series are strongly autocorrelated, and an autoregressive (AR) model fitted to the residuals of a short "
  "record underestimates their persistence, so that tests based on it reject a true null hypothesis too often. The test "
  "used here (main text Sect. 2.6) corrects the AR coefficients for small-sample bias (Kilian, 1998) and calibrates the "
  "trend statistic against 300 synthetic series simulated from the corrected model and processed identically. For "
  "monthly series it is run with an AIC-selected AR order and with a fixed order of 12; the larger p-value is reported "
  "and the confidence interval is the union of the two.")
P("The false-positive rate was measured on trend-free synthetic series with a seasonal cycle (Table S1). With the AIC "
  "order the rate is close to nominal for AR(1) processes and too high (9.5 %) for a process with the autocorrelation of the "
  "Drake Passage EKE, because AIC drops the small high-lag coefficients that carry part of the low-frequency persistence; "
  "with the fixed AR(12) model the opposite holds. Taking the larger p-value gives at most the smaller of the two rates, "
  "3–5 % for the monthly processes tested. For annual series an AR(1) model with iterated bias correction holds the level "
  "for 25 values (5.0 %) but not for 15 (9.5 %); the 15-year in situ transport series can therefore not be tested at the "
  "nominal level, and its trend is not significant under either test. With this test, fewer than 5 % of the 1° cells of the "
  "trend maps reach p < 0.05 (2.9 % for KE, 4.0 % for EKE), so the local trends are not significant as a field beyond the "
  "few cells retained after false-discovery-rate control (Benjamini and Hochberg, 1995; Wilks, 2016).")
P(" False-positive rate of the trend test at the nominal 5 % level, on trend-free series with a seasonal cycle "
  "(number of synthetic series in brackets). –: not run.", bold_prefix="Table S1.")
TABLE([["AR(1), φ = 0.8, 300 months", "5.3 % (300)", "8 % (100)"],
       ["AR(1), φ = 0.9, 300 months", "3.5 % (200)", "–"],
       ["AR(12) fitted to Drake Passage EKE, 300 months", "9.5 % (200)", "3 % (100)"],
       ["Annual AR(1), φ = 0.5, 25 years", "5.0 % (200)*", "–"],
       ["Annual AR(1), φ = 0.5, 15 years", "9.5 % (200)*", "–"]],
      ["Synthetic process", "AIC-selected order", "Fixed AR(12)"])
P("* AR(1) with iterated bias correction (annual series).")

# ---------------------------------------------------------------- S3
H("S3 Half-power analysis: sensitivity, longitude bins, detection and the two-satellite product")
P("The method is illustrated in Fig. 2b of the main text. The circumpolar trend of the half-power latitude does not depend "
  "on the detection parameters (Fig. S2a): thresholds of 150–400 cm² s⁻² or a threshold relative to the monthly band mean, "
  "peak fractions of 20–50 % and minimum widths of 50–300 km all give trends between −1.9 and +1.8 km per decade with "
  "p > 0.15. Trends in the 33 longitude bins of 5° that contain enough matched longitudes range from about −24 to +33 km per "
  "decade with wide confidence intervals, and none reaches p < 0.05 (Fig. S2b). Bins are omitted where too few longitudes "
  "have the envelope in at least 90 % of the months; the detection frequency is lowest (about 20–45 %) at about 180–155° W, "
  "115–75° W and 50–55° E.")
P("The detectability of a shift was estimated by injecting a uniform meridional displacement, growing linearly over the "
  "25 years, into the monthly KE fields and repeating the detection and the trend test (41 realisations per case). A "
  "trend-free field is flagged in 5 % of the realisations; a 5 km shift in 12 %, a 10 km shift in 44 % and a 25 km or larger "
  "shift in all of them (Fig. S2c). About 70–75 % of an injected shift of 10 km or more is recovered (Fig. S2d), so trends estimated with the method "
  "underestimate a real shift by about a quarter to a third.")
P("In the two-satellite product the half-power latitude shows no trend either (Fig. S3a): +0.2 km per decade (−3.2 to "
  "+3.6; p = 0.91) for the circumpolar mean, with no significant sector trend. The KE at the envelope peak increases by "
  "3.6 % per decade (1.8–5.5; p = 0.003), against 5.3 % (2.8–7.8; p = 0.003) in the all-satellite product, and is "
  "significant in the Pacific and Indian sectors in both products (Fig. S3b).")
FIGP("FigureS2")
P(" Half-power analysis. (a) Circumpolar trend of the half-power latitude, 2001–2025, for the baseline parameters and "
  "for each alternative (threshold in cm² s⁻²; relative threshold: the monthly ACC-band mean KE; fraction: of the peak "
  "KE; minimum width of the envelope), in km per decade (positive northward). (b) Trends in 5° longitude bins (bin "
  "centre; bins with too few matched longitudes are omitted). (c) Percentage of realisations in which an injected "
  "shift, growing linearly to the value on the x axis over 2001–2025, is detected at p < 0.05 (dashed: 5 %). (d) Median "
  "recovered shift (dashed: 1:1). Open symbols: p ≥ 0.05; 95 % confidence intervals.", bold_prefix="Figure S2.")
FIGP("FigureS3", width=6.6)
P(" Half-power analysis in the all-satellite (Product 1) and two-satellite (Product 4) products, 2001–2025: trends of "
  "(a) the half-power latitude (km per decade, positive northward) and (b) the KE at the envelope peak (percent of the "
  "mean per decade), circumpolar and by sector, with 95 % confidence intervals. Filled symbols: p < 0.05.",
  bold_prefix="Figure S3.")

# ---------------------------------------------------------------- S4
H("S4 Tracked eddies: generation sites and nonlinearity")
P("Eddy births, the first observation of each track in the META4.0 atlas (Product 5; method of Pegliasco et al., 2022), are "
  "enhanced along the whole energetic band rather than at particular sites (Fig. S4a). The mean over the ACC band is 1.7 "
  "births per year per 1° cell; it is 1.8 in the wake box of the Pacific–Antarctic Ridge and 2.6 in that of the Kerguelen "
  "Plateau, and births are particularly frequent around Kerguelen and in Drake Passage. The meander "
  "sections are not singular generation sites at this resolution. The trends of the number of eddies present, the "
  "births, the mean amplitude and the sum of squared amplitudes in the four wake boxes do not differ from those in the "
  "rest of the band (difference of the relative trends between −1.9 and +0.3 % per decade; p > 0.2).")
P("The nonlinearity ratio U/c, the rotational speed over the translation speed from the displacement of the centre over "
  "seven days, for tracks lasting at least seven days, has a monthly median of about 6.9, and 99 % of the observations have U/c > 1 "
  "(Fig. S4b). It has no trend (−0.02 % per decade, −1.5 to +1.4; p = 0.97): in this subset both U and c increase, by about "
  "3.5 and 3.2 % per decade. The tracked eddies are therefore nonlinear, coherent structures throughout the record, and "
  "their intensification does not change this.")
FIGP("FigureS4")
P(" Tracked eddies (Product 5, January 2001 to August 2023). (a) Eddy births per year per 1° cell, with the mean "
  "boundaries of the ACC band (grey). (b) Median nonlinearity ratio U/c of the eddies in the ACC band, monthly (thin) and "
  "12-month running mean (thick).",
  bold_prefix="Figure S4.")

# ---------------------------------------------------------------- S5
H("S5 Month-to-month covariation of path length and EKE")
P("The mechanism of Thompson and Naveira Garabato (2014) implies that a change in the shape of a standing meander is "
  "followed by a change in eddy energy downstream of it. We tested this on monthly time scales without fixed boxes. The "
  "jet-core path of the main text (Sect. 2.4; +0.18 m relative-ADT contour, Subantarctic Front side) was divided into 36 "
  "overlapping windows of 30° of longitude, centred every 10°. In each window and month we computed the length of the "
  "path smoothed at 150 km and the EKE within 150 km of the path, in the same window and in the band 15–35° downstream of "
  "the window centre (which does not overlap the window). The deseasonalised, detrended series, smoothed with a 2-month running mean, were correlated at lags of "
  "0–3 months (EKE lagging), and significance was assessed with 1000 phase-randomised surrogates of the path length, "
  "taking the maximum over lags. The design was fixed before the results were computed and corrected, before the final results, where it was "
  "contaminated by construction (the downstream band now excludes the window, and a single circumpolar contour level replaced "
  "per-window levels that had picked non-ACC jets); the path is tracked in 99–100 % "
  "of the months in every window.")
P("Path length and EKE in the same window covary positively almost everywhere (Fig. S5a): the correlation is significant "
  "in 20 of the 36 windows (median over all windows 0.28, maximum 0.62), and in 14 of these 20 it is largest at lag 0, i.e. in the same "
  "month. EKE 15–35° downstream is not related to the path length (Fig. S5b): the correlation is significant in 1 of 36 "
  "windows (r = 0.27 at 20° W), about what is expected by chance. With the path at the Polar Front side (−0.20 m), the "
  "same-window correlation is significant in 14 of 36 windows (median over all windows 0.16), at lag 0 in 8 of these 14, and the downstream "
  "one in 3 of 36, with both signs. The robust result is therefore a local, simultaneous covariation of path distortion "
  "and EKE, which is expected because the eddies themselves distort the monthly path, and no detectable delayed response "
  "downstream; the step from meander shape to downstream eddy energy is not resolved on monthly time scales with these data.")
FIGP("FigureS5")
P(" Month-to-month covariation of the length of the jet-core path and EKE, in 30° windows (2001–2025, Product 1): "
  "correlation at the lag (0–3 months, EKE lagging) with the largest magnitude, between the path length and (a) EKE "
  "within 150 km of the path in the same window and (b) EKE within 150 km of the path 15–35° downstream of the window "
  "centre. Black: path at the Subantarctic Front level (+0.18 m, as in the main text); grey: Polar Front level (−0.20 m). "
  "Filled symbols: p < 0.05 against phase-randomised surrogates.", bold_prefix="Figure S5.")

# ---------------------------------------------------------------- refs
H("References")
REFS = [
    "Ballarotta, M., Pujol, M.-I., and SL-TAC team: Quality information document for Sea Level TAC C3S products, Copernicus Marine Service, 2024.",
    "Barceló-Llull, B., Rosselló, P., Combes, V., Sánchez-Román, A., Pujol, M. I., and Pascual, A.: Kuroshio Extension and Gulf Stream dominate the Eddy Kinetic Energy intensification observed in the global ocean, Sci. Rep., 15, 21754, https://doi.org/10.1038/s41598-025-06149-9, 2025.",
    "Benjamini, Y. and Hochberg, Y.: Controlling the false discovery rate: a practical and powerful approach to multiple testing, J. R. Stat. Soc. B, 57, 289–300, https://doi.org/10.1111/j.2517-6161.1995.tb02031.x, 1995.",
    "Hargous, P., Combes, V., Barceló-Llull, B., and Pascual, A.: Eddy kinetic energy variability from 30 years of altimetry in the Mediterranean Sea, Ocean Sci., 22, 549–564, https://doi.org/10.5194/os-22-549-2026, 2026.",
    "International Altimetry Team: Altimetry for the future: Building on 25 years of progress, Adv. Space Res., 68, 319–363, https://doi.org/10.1016/j.asr.2021.01.022, 2021.",
    "Kilian, L.: Small-sample confidence intervals for impulse response functions, Rev. Econ. Stat., 80, 218–230, https://doi.org/10.1162/003465398557465, 1998.",
    "Liu, X., Yang, C., and Chen, Y.: Standing meanders of the Antarctic Circumpolar Current: Evidence for ridge-controlled eddy saturation, J. Geophys. Res.-Oceans, 131, e2025JC023527, https://doi.org/10.1029/2025JC023527, 2026.",
    "Martínez-Moreno, J., Hogg, A. M., England, M. H., Constantinou, N. C., Kiss, A. E., and Morrison, A. K.: Global changes in oceanic mesoscale currents over the satellite altimetry record, Nat. Clim. Change, 11, 397–403, https://doi.org/10.1038/s41558-021-01006-9, 2021.",
    "Pegliasco, C., Delepoulle, A., Mason, E., Morrow, R., Faugère, Y., and Dibarboure, G.: META3.1exp: a new global mesoscale eddy trajectory atlas derived from altimetry, Earth Syst. Sci. Data, 14, 1087–1107, https://doi.org/10.5194/essd-14-1087-2022, 2022.",
    "Thompson, A. F. and Naveira Garabato, A. C.: Equilibration of the Antarctic Circumpolar Current by standing meanders, J. Phys. Oceanogr., 44, 1811–1828, https://doi.org/10.1175/JPO-D-13-0163.1, 2014.",
    "Wilks, D. S.: “The stippling shows statistically significant grid points”: How research results are routinely overstated and overinterpreted, and what to do about it, B. Am. Meteorol. Soc., 97, 2263–2273, https://doi.org/10.1175/BAMS-D-15-00267.1, 2016.",
]
for r in REFS:
    P(r)
d.save(out); print("saved", out)
