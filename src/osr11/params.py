"""Analysis constants fixed before any trend was computed (one definition, imported where needed).

The numbered scripts were written before this module existed and some still carry the same values
inline; those inline values are identical to the ones here (checked when the repository was assembled).
"""
# Half-power envelope detection, baseline D1 (Sect. 2.2): threshold 200 cm2 s-2 (0.02 m2 s-2),
# minimum above-threshold segment 100 km, envelope bounds at 25 % of the peak, at most 8 envelopes per profile.
HP_THR = 0.02
HP_MIN_W = 100.0
HP_FRAC = 0.25
HP_MAX_ENV = 8

# ACC band: 2001-2025 mean ADT between these levels (m)
ADT_SOUTH, ADT_NORTH = -0.6, 0.2

# Sectors (deg E, 0-360): Indian 20-147, Pacific 147-290, Atlantic otherwise
SECTORS_360 = {"Indian": (20.0, 147.0), "Pacific": (147.0, 290.0)}

# Record
Y0, Y1 = 2001, 2025
NMONTHS = 12 * (Y1 - Y0 + 1)

# Two periods compared in Sect. 3.1 (Drake Passage EKE between the fronts)
PERIOD_EARLY = (2001, 2015)
PERIOD_LATE = (2016, 2025)

# km per degree of latitude
KM_PER_DEG = 111.195
