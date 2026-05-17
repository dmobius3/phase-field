import numpy as np

from coherence.radii import flat_onset_radius, transition_radius


def test_transition_radius_simple_onset():
    rad = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    g_obs = np.array([1.0, 1.0, 1.3, 1.5, 1.6])
    assert transition_radius(rad, g_obs, np.ones(5), 1.2) == 3.0


def test_transition_radius_requires_all_outer_points():
    rad = np.array([1.0, 2.0, 3.0, 4.0])
    g_obs = np.array([1.3, 1.0, 1.3, 1.3])  # index 1 dips below threshold
    assert transition_radius(rad, g_obs, np.ones(4), 1.2) == 3.0


def test_transition_radius_none_when_last_point_below():
    rad = np.array([1.0, 2.0, 3.0])
    g_obs = np.array([1.3, 1.3, 1.0])
    assert transition_radius(rad, g_obs, np.ones(3), 1.2) is None


def test_flat_onset_radius_basic():
    rad = np.array([1.0, 2.0, 3.0, 4.0])
    v = np.array([80.0, 96.0, 99.0, 100.0])  # v_c = 100, within 5% from idx 1
    assert flat_onset_radius(rad, v, 100.0, 0.05) == 2.0


def test_flat_onset_radius_none_for_rising_curve():
    rad = np.array([1.0, 2.0, 3.0])
    v = np.array([50.0, 60.0, 70.0])
    assert flat_onset_radius(rad, v, 100.0, 0.05) is None
