"""Analysis constants fixed before any trend was computed (one definition, imported where needed).

The numbered scripts were written before this module existed and some still carry the same values
inline; those inline values are identical to the ones here (checked when the repository was assembled).
"""
# Half-power envelope detection (Sect. 2.2): minimum above-threshold segment 100 km, envelope bounds at 25 % of the peak, at most
# 8 envelopes per profile. The paper uses a longitude-dependent threshold, 20 % of the maximum of each longitude's 2001-2025 mean KE
# in the ACC band (HP_REL_THR; halfpower.SETS["relthr"]). The fixed threshold of 200 cm2 s-2 (0.02 m2 s-2) of Chambers (2018)
# (HP_THR; halfpower.SETS["baseline"]) is the sensitivity set of the Supplement (Sect. S4.1, Fig. S5).
HP_REL_THR = 0.20
HP_THR = 0.02
HP_MIN_W = 100.0
HP_FRAC = 0.25
HP_MAX_ENV = 8
HP_MIN_COVERAGE = 0.90   # a track is retained when matched in at least 90 % of the months

# ACC band: 2001-2025 mean ADT between these levels (m)
ADT_SOUTH, ADT_NORTH = -0.6, 0.2

# Sectors (deg E, 0-360): Indian 20-147, Pacific 147-290, Atlantic otherwise
SECTORS_360 = {"Indian": (20.0, 147.0), "Pacific": (147.0, 290.0)}

# Record
Y0, Y1 = 2001, 2025
NMONTHS = 12 * (Y1 - Y0 + 1)

# Two periods compared in the paper (timing of the EKE increase; composites of Figs. 2b and 3b)
PERIOD_EARLY = (2001, 2017)
PERIOD_LATE = (2018, 2025)

# km per degree of latitude
KM_PER_DEG = 111.195
