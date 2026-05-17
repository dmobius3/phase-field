"""Operational radius algorithms: the transition radius r_t and the
flat-onset radius R_flat. Definitions are frozen in PREREGISTRATION.md.

Both algorithms find the smallest radius beyond which a condition holds
at *every* measured point, so a single noisy inner point cannot trip the
detection. A galaxy with no such radius returns ``None`` and is counted
as a separate sub-population, not a prediction failure.
"""

from __future__ import annotations

import numpy as np


def _suffix_onset(rad, condition):
    """Smallest rad[i] such that ``condition`` is True for all j >= i.

    Returns ``None`` if the last point fails the condition.
    """
    rad = np.asarray(rad, dtype=float)
    condition = np.asarray(condition, dtype=bool)
    n = len(rad)
    onset = n
    for i in range(n - 1, -1, -1):
        if condition[i]:
            onset = i
        else:
            break
    return None if onset == n else float(rad[onset])


def transition_radius(rad, g_obs, g_bar, ratio_threshold: float = 1.2):
    """Baryon-to-total divergence radius.

    Smallest r_t such that g_obs / g_bar >= ``ratio_threshold`` for all
    measured points at r >= r_t. ``rad`` must be sorted ascending.
    """
    ratio = np.asarray(g_obs, dtype=float) / np.asarray(g_bar, dtype=float)
    return _suffix_onset(rad, ratio >= ratio_threshold)


def flat_onset_radius(rad, v, v_c, tol: float = 0.05):
    """Flat-onset radius.

    Smallest R_flat such that |v - v_c| / v_c <= ``tol`` for all measured
    points at r >= R_flat. ``rad`` must be sorted ascending.
    """
    v = np.asarray(v, dtype=float)
    return _suffix_onset(rad, np.abs(v - v_c) / v_c <= tol)
