# Processed outputs

Small files written by the pipeline and read by the figures and notebooks. `MANIFEST.csv` gives, for every file, its size,
SHA-256, the script that produced it, whether that script needs the raw data, and the figures that read it. Folder names are
those of the working repository (kept so that scripts and outputs match):

| Folder | Content | Script |
|---|---|---|
| `base` | mean ADT, ice mask, 2001–2025 climatology, bathymetry | `00c`, `01a`, `01d` |
| `a1` | EKE and KE trend maps, 1° cells | `02f`, `02g` |
| `a2` | fixed-threshold half-power sets and tracks (sensitivity set) | `02b`, `02c`, `02e`, `05j` |
| `a5`, `a6`, `a7` | KE and EKE series of the ACC band, sectors and boxes | `02a`, `02m`, `03u` |
| `x1` | fixed-threshold tracks against fixed ADT contours | `02j`, `02k`, `03t` |
| `x6`, `x7` | energetic pathway with the longitude-dependent threshold (the paper's tracks): tracks, 2° bins, trends, composites | `02n`, `02o`, `02p` |
| `s6`, `x3` | occurrence of the half-power points | `02q` |
| `m1`, `m2` | regions and periods of earlier studies | `02l`, `05k` |
| `x11` | lag behind the SAM index; SAM index of Marshall (2003) | `05l` |
| `drake/dkN`, `d1`, `c04` | Drake Passage steps | `03a`–`03w` |
| `f3w`, `f4`, `f4m`, `f4e` | wind stress, isopycnal slope, domain means | `04a`–`04g` |
| `x9`, `x10`, `x12` | exponents, their intervals, lagged correlations | `04i`, `04j`, `04k` |
| `x13`, `x14` | 1° cell series and the test of the width of the domain | `04l`, `04m` |
| `p4` | two-satellite product (Product 4) | `06_product4/` |
| `retest`, `stats` | trend tests of the saved series, size checks | `05a`–`05h` |
| `fig1` | profile and January 2001 velocity for Fig. 1 | `fig1_halfpower_method.py`, `00a` |
