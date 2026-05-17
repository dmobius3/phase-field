"""Dimensionless ratios eta = radius / L_f and the statistics that score
them against the registered acceptance criteria (PREREGISTRATION.md,
prediction 3): a bootstrap confidence interval on the median and a
one-sided chi-squared test of the sample scatter against sigma_pred.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def eta(radius, L_f):
    """Ratio radius / L_f. Returns ``None`` if either input is ``None``."""
    if radius is None or L_f is None:
        return None
    return float(radius) / float(L_f)


def _clean(values):
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def bootstrap_median_ci(values, n_boot: int = 10_000, ci: float = 0.95,
                        seed: int = 0):
    """Sample median and a percentile bootstrap CI on the median."""
    values = _clean(values)
    if len(values) == 0:
        raise ValueError("no finite values")
    rng = np.random.default_rng(seed)
    n = len(values)
    boot = np.median(values[rng.integers(0, n, size=(n_boot, n))], axis=1)
    alpha = 1.0 - ci
    lo, hi = np.quantile(boot, [alpha / 2, 1.0 - alpha / 2])
    return float(np.median(values)), float(lo), float(hi)


def scatter_chi2(values, sigma_pred):
    """One-sided chi-squared test that the sample scatter does not exceed
    ``sigma_pred``.

    Returns ``(chi2, dof, p_value)``; ``p_value`` is the upper-tail
    probability, so small p means the observed scatter is significantly
    larger than sigma_pred.
    """
    values = _clean(values)
    dof = len(values) - 1
    if dof < 1:
        raise ValueError("need at least two finite values")
    sample_var = np.var(values, ddof=1)
    chi2 = dof * sample_var / sigma_pred**2
    p_value = stats.chi2.sf(chi2, dof)
    return float(chi2), int(dof), float(p_value)


def median_within(values, interval=(0.5, 2.0), **boot_kwargs):
    """Prediction-3 location criterion: the median and its bootstrap CI
    both lie inside ``interval``. Returns ``(passed, median, lo, hi)``.
    """
    median, lo, hi = bootstrap_median_ci(values, **boot_kwargs)
    low, high = interval
    passed = (low <= median <= high) and (low <= lo and hi <= high)
    return passed, median, lo, hi
