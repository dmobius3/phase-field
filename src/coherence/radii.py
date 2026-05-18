"""Operational radius algorithms: the transition radius r_t and the
flat-onset radius R_flat. Definitions are frozen in PREREGISTRATION.md.

Each onset uses a persistence rule: the condition must hold at the
candidate radius and at no less than a fraction ``persistence`` of the
measured points beyond it. A strict "all points beyond" rule is a
one-sided ratchet under measurement error, since a single noisy point
trips it and pushes the onset outward but never inward. The persistence
rule removes that ratchet while keeping the conceptual definition of a
sustained transition (sustained divergence for r_t, sustained flatness
for R_flat).

A galaxy with no qualifying radius returns ``None`` and is counted as a
separate sub-population, not a prediction failure.
"""

from __future__ import annotations

import numpy as np

DEFAULT_PERSISTENCE = 0.8
RFLAT_MEASURABILITY = 0.03  # max errV/v_c for R_flat to be measurable


def _persistence_onset(rad, condition, persistence: float = DEFAULT_PERSISTENCE):
    """Smallest rad[i] such that ``condition[i]`` is True and at least
    ``persistence`` of the points at indices >= i satisfy ``condition``.

    Returns ``None`` if no such radius exists. ``rad`` must be sorted
    ascending.
    """
    rad = np.asarray(rad, dtype=float)
    condition = np.asarray(condition, dtype=bool)
    for i in range(len(condition)):
        if condition[i] and condition[i:].mean() >= persistence:
            return float(rad[i])
    return None


def transition_radius(rad, g_obs, g_bar, ratio_threshold: float = 1.2,
                      persistence: float = DEFAULT_PERSISTENCE):
    """Baryon-to-total divergence radius.

    Smallest r_t such that g_obs / g_bar >= ``ratio_threshold`` holds at
    r_t and at no less than ``persistence`` of the points at r >= r_t.
    """
    ratio = np.asarray(g_obs, dtype=float) / np.asarray(g_bar, dtype=float)
    return _persistence_onset(rad, ratio >= ratio_threshold, persistence)


def flat_onset_radius(rad, v, v_c, tol: float = 0.05,
                      persistence: float = DEFAULT_PERSISTENCE):
    """Flat-onset radius.

    Smallest R_flat such that |v - v_c| / v_c <= ``tol`` holds at R_flat
    and at no less than ``persistence`` of the points at r >= R_flat.
    """
    v = np.asarray(v, dtype=float)
    return _persistence_onset(rad, np.abs(v - v_c) / v_c <= tol, persistence)


def rflat_measurable(errV, v_c, max_fraction: float = RFLAT_MEASURABILITY):
    """Whether R_flat is measurable for a galaxy.

    The representative per-point fractional velocity error errV / v_c
    must sit below ``max_fraction`` (comfortably under the flatness
    tolerance). When per-point noise approaches the tolerance, R_flat is
    biased high and is not reported. ``errV`` and ``v_c`` share a unit.
    """
    return (errV / v_c) < max_fraction

