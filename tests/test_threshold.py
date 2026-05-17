import numpy as np

from coherence.scales import coherence_length_kpc
from coherence.synthetic import make_galaxy
from coherence.threshold import XI, trigger_ratio


def test_constant_curve_gives_inverse_xi():
    # A genuinely flat curve (v == v_c everywhere) gives T/T_c = 1/xi.
    rad = np.linspace(0.1, 30.0, 80)
    v = np.full_like(rad, 150.0)
    L_f = coherence_length_kpc(150.0)
    ratio = trigger_ratio(rad, v, L_f, 150.0)
    assert np.isclose(ratio, 1.0 / XI, rtol=0.05)


def test_flat_galaxy_above_threshold_and_below_inverse_xi():
    g = make_galaxy("flat", v_c=150.0)
    L_f = coherence_length_kpc(150.0)
    ratio = trigger_ratio(g["rad"], g["vobs"], L_f, g["v_c"])
    # rising-then-flat curve: above the T/T_c = 1 threshold, below 1/xi
    assert 1.0 < ratio < 1.0 / XI


def test_rising_curve_ratio_below_flat():
    L_f = coherence_length_kpc(150.0)
    flat = make_galaxy("flat", v_c=150.0)
    rising = make_galaxy("rising", v_c=150.0)
    r_flat = trigger_ratio(flat["rad"], flat["vobs"], L_f, flat["v_c"])
    r_rise = trigger_ratio(rising["rad"], rising["vobs"], L_f, rising["v_c"])
    assert r_rise < 1.0 < r_flat
