import numpy as np

from coherence.radii import (flat_onset_radius, rflat_measurable,
                             transition_radius)


def test_transition_radius_simple_onset():
    rad = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    g_obs = np.array([1.0, 1.0, 1.3, 1.5, 1.6])
    assert transition_radius(rad, g_obs, np.ones(5), 1.2) == 3.0


def test_transition_radius_tolerates_single_noisy_tail_point():
    # Nine points qualify from index 2; one point at index 7 dips below.
    # A strict "all beyond" rule would push the onset to 9.0; the
    # persistence rule (0.8) keeps it at the true onset.
    rad = np.arange(1.0, 11.0)
    ratio = np.array([1.0, 1.0, 1.3, 1.3, 1.3, 1.3, 1.3, 1.0, 1.3, 1.3])
    assert transition_radius(rad, ratio, np.ones(10), 1.2,
                             persistence=0.8) == 3.0


def test_transition_radius_rejects_isolated_point():
    # A single isolated qualifying point is not a sustained transition.
    rad = np.arange(1.0, 9.0)
    ratio = np.array([1.0, 1.0, 1.3, 1.0, 1.0, 1.0, 1.0, 1.0])
    assert transition_radius(rad, ratio, np.ones(8), 1.2,
                             persistence=0.8) is None


def test_transition_radius_large_dip_moves_onset():
    # A dip exceeding the persistence slack pushes the onset outward.
    rad = np.array([1.0, 2.0, 3.0, 4.0])
    ratio = np.array([1.3, 1.0, 1.3, 1.3])  # 1 of 4 fails at index 0
    assert transition_radius(rad, ratio, np.ones(4), 1.2,
                             persistence=0.8) == 3.0


def test_flat_onset_radius_basic():
    rad = np.array([1.0, 2.0, 3.0, 4.0])
    v = np.array([80.0, 96.0, 99.0, 100.0])  # v_c = 100, within 5% from idx 1
    assert flat_onset_radius(rad, v, 100.0, 0.05) == 2.0


def test_flat_onset_radius_tolerates_single_noisy_point():
    rad = np.arange(1.0, 11.0)
    v = np.array([70, 100, 100, 100, 100, 100, 80, 100, 100, 100],
                 dtype=float)
    # one outlier at index 6; persistence 0.8 keeps the onset at idx 1
    assert flat_onset_radius(rad, v, 100.0, 0.05, persistence=0.8) == 2.0


def test_flat_onset_radius_none_for_rising_curve():
    rad = np.array([1.0, 2.0, 3.0])
    v = np.array([50.0, 60.0, 70.0])
    assert flat_onset_radius(rad, v, 100.0, 0.05) is None


def test_rflat_measurability_cut():
    # errV/v_c must sit below 0.03 for R_flat to be measurable.
    assert rflat_measurable(4.0, 200.0)        # 2.0% -> measurable
    assert not rflat_measurable(4.0, 100.0)    # 4.0% -> not measurable
    assert not rflat_measurable(4.0, 133.3)    # 3.0% -> at the cut, excluded
