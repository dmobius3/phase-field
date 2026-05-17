import numpy as np

from coherence.ratios import (bootstrap_median_ci, eta, median_within,
                              scatter_chi2)


def test_eta_none_passthrough():
    assert eta(None, 5.0) is None
    assert eta(5.0, None) is None
    assert eta(6.0, 3.0) == 2.0


def test_bootstrap_median_ci_brackets_median():
    rng = np.random.default_rng(1)
    values = rng.normal(1.0, 0.2, 500)
    median, lo, hi = bootstrap_median_ci(values, n_boot=2000, seed=1)
    assert lo < median < hi
    assert abs(median - 1.0) < 0.1


def test_scatter_chi2_detects_excess_scatter():
    rng = np.random.default_rng(2)
    tight = rng.normal(1.0, 0.10, 200)
    wide = rng.normal(1.0, 0.50, 200)
    _, _, p_tight = scatter_chi2(tight, sigma_pred=0.10)
    _, _, p_wide = scatter_chi2(wide, sigma_pred=0.10)
    assert p_tight > 0.01   # scatter consistent with sigma_pred
    assert p_wide < 0.01    # scatter far exceeds sigma_pred


def test_median_within_interval():
    passed, *_ = median_within(np.full(100, 1.0), interval=(0.5, 2.0))
    assert passed
    failed, *_ = median_within(np.full(100, 3.0), interval=(0.5, 2.0))
    assert not failed
