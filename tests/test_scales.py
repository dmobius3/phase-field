import numpy as np

from coherence.scales import baryonic_velocity, coherence_length_kpc, g_observed


def test_coherence_length_milky_way():
    # v_c = 220 km/s gives L_f ~ 13 kpc (working note, section I)
    assert 12.0 < coherence_length_kpc(220.0) < 14.0


def test_coherence_length_scales_as_vc_squared():
    ratio = coherence_length_kpc(200.0) / coherence_length_kpc(100.0)
    assert np.isclose(ratio, 4.0)


def test_g_observed_declines_for_flat_curve():
    g = g_observed(np.full(3, 150.0), np.array([1.0, 2.0, 4.0]))
    assert g[0] > g[1] > g[2]


def test_baryonic_velocity_disk_only():
    # gas = bulge = 0, ml_disk = 0.5  ->  v_bar = sqrt(0.5) * v_disk
    vbar = baryonic_velocity(0.0, 100.0, 0.0, ml_disk=0.5)
    assert np.isclose(vbar, np.sqrt(0.5) * 100.0)
