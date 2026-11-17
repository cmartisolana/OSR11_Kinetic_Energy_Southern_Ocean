"""Half-power envelope detection on a synthetic KE profile with a known peak and known half-power latitude."""
import numpy as np


def _profile(lat, peak_lat=-52.0, width=1.0, amp=0.05):
    return amp * np.exp(-0.5 * ((lat - peak_lat) / width) ** 2)      # m2 s-2 (0.05 = 500 cm2 s-2)


def test_detect_profile_symmetric_gaussian():
    from osr11 import fronts as F, params as P
    lat = np.arange(-65.0, -40.0, 0.125)
    prof = _profile(lat, amp=0.1)                             # 1000 cm2 s-2: the 25 % level lies above the threshold
    out = np.empty((P.HP_MAX_ENV, F.NF))
    n = F.detect_profile(prof, lat, F.meridional_dists(lat), P.HP_THR, P.HP_MIN_W, P.HP_FRAC, out)
    assert n == 1
    row = dict(zip(F.FIELDS, out[0]))
    assert abs(row["peak_lat"] + 52.0) < 0.13                 # peak on the grid
    assert abs(row["lat_hp"] + 52.0) < 0.07                   # symmetric profile: half-power at the peak
    assert abs(row["lat_s"] - (-52.0 - np.sqrt(-2 * np.log(0.25)))) < 0.15   # 25 % bounds of a Gaussian
    assert abs(row["lat_n"] - (-52.0 + np.sqrt(-2 * np.log(0.25)))) < 0.15


def test_detect_profile_skewed_profile_shifts_half_power():
    """Adding energy on the northern flank moves the half-power latitude north of the peak."""
    from osr11 import fronts as F, params as P
    lat = np.arange(-65.0, -40.0, 0.125)
    prof = _profile(lat) + 0.6 * _profile(lat, peak_lat=-50.5, width=1.5)
    out = np.empty((P.HP_MAX_ENV, F.NF))
    n = F.detect_profile(prof, lat, F.meridional_dists(lat), P.HP_THR, P.HP_MIN_W, P.HP_FRAC, out)
    assert n == 1
    row = dict(zip(F.FIELDS, out[0]))
    assert row["lat_hp"] > row["peak_lat"]


def test_detect_field_threshold_and_width():
    """Below the threshold or narrower than min_w nothing is detected; two separated jets give two envelopes."""
    from osr11 import fronts as F, params as P
    lat = np.arange(-65.0, -40.0, 0.125)
    d = F.meridional_dists(lat)
    weak = _profile(lat, amp=0.01)[:, None]                                   # 100 cm2 s-2 < 200 threshold
    narrow = _profile(lat, width=0.15)[:, None]                                # ~40 km above threshold < 100 km
    two = (_profile(lat, peak_lat=-58.0) + _profile(lat, peak_lat=-47.0))[:, None]
    field = np.hstack([weak, narrow, two])
    counts, table = F.detect_field(field, lat, d, P.HP_THR, P.HP_MIN_W, P.HP_FRAC, P.HP_MAX_ENV)
    assert list(counts) == [0, 0, 2]
    assert table.shape == (3, P.HP_MAX_ENV, F.NF)
