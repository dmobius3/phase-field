"""Synthetic SPARC-like galaxies spanning a range of rotation-curve
morphologies. Used by the test suite and the sigma_pred Monte Carlo. No
SPARC data is involved, so these generators may run in Phase 0.

Four shape families approximate the SPARC morphological range:

    sharp      fast (tanh) turnover
    gradual    slow (arctan) turnover
    bumpy      gradual turnover with damped spiral-arm wiggles
    declining  gradual turnover with a gently declining outer profile

Each is generated at a grid of turnover scales (the target flat-onset
radius as a fraction of L_f), so the sigma_pred Monte Carlo can report
whether the measurement bias is stable across curve shape or shape
dependent.
"""

from __future__ import annotations

import numpy as np

from .scales import coherence_length_kpc

ML_DISK = 0.5  # 3.6 um stellar mass-to-light ratio used in the generator
SHAPES = ("sharp", "gradual", "bumpy", "declining")
# Declining curves are removed from the flat-curve analysis upstream by
# the "exclude rising or falling rotation curves" sample cut, so the
# sigma_pred population covers only the genuinely flat shapes.
FLAT_SHAPES = ("sharp", "gradual", "bumpy")
ONSET_GRID = (0.5, 1.0, 1.7)  # target flat-onset radius in units of L_f

# Internal turnover scale: x = r / r_turn at which the profile first
# reaches ~0.95 of v_c, used to place the flat onset near a target radius.
_REACH = {"sharp": 1.832, "gradual": 12.0, "bumpy": 12.0, "declining": 12.0}


def _flat_profile(rad, v_c, r_turn, shape):
    """Velocity profile for a flat-type rotation curve of a given shape."""
    base = v_c * (2.0 / np.pi) * np.arctan(rad / r_turn)
    if shape == "gradual":
        return base
    if shape == "sharp":
        return v_c * np.tanh(rad / r_turn)
    if shape == "bumpy":
        damp = np.exp(-rad / rad[-1])
        return base * (1.0 + 0.10 * damp
                       * np.sin(2.0 * np.pi * rad / (6.0 * r_turn)))
    if shape == "declining":
        return base * (1.0 - 0.06 * rad / rad[-1])
    raise ValueError(f"unknown shape: {shape!r}")


def make_galaxy(kind: str = "flat", shape: str = "gradual",
                v_c: float = 150.0, onset: float = 1.0, span: float = 4.0,
                n_points: int = 30, baryon_outer: float = 0.55,
                distance: float = 10.0, inclination: float = 60.0,
                quality: int = 1):
    """Construct a synthetic galaxy as a dict of arrays and scalars.

    ``kind`` is ``"flat"`` (a rotation curve that plateaus at v_c, with
    the morphology set by ``shape``) or ``"rising"`` (linear, still
    climbing at the last point). ``onset`` is the target flat-onset
    radius in units of L_f; the shape's internal turnover scale is set to
    hit it. ``span`` is the outermost radius as a multiple of L_f, and
    ``baryon_outer`` is v_bar / v_obs at the last point.

    Velocity fields use only the disk component (gas and bulge are zero),
    so v_bar = sqrt(ML_DISK) * v_disk.
    """
    L_f = float(coherence_length_kpc(v_c))
    r_max = span * L_f
    rad = np.linspace(r_max / n_points, r_max, n_points)

    if kind == "rising":
        vobs = v_c * rad / r_max
        r_turn = onset * L_f
    elif kind == "flat":
        r_turn = (onset * L_f) / _REACH[shape]
        vobs = _flat_profile(rad, v_c, r_turn, shape)
    else:
        raise ValueError(f"unknown kind: {kind!r}")

    frac = 1.0 - (1.0 - baryon_outer) * (rad / r_max)
    vbar = vobs * frac
    vdisk = vbar / np.sqrt(ML_DISK)
    zeros = np.zeros_like(rad)

    return {
        "rad": rad,
        "vobs": vobs,
        "e_vobs": 0.05 * np.abs(vobs) + 2.0,
        "vgas": zeros,
        "vdisk": vdisk,
        "vbul": zeros.copy(),
        "v_c": float(v_c),
        "kind": kind,
        "shape": shape,
        "distance": float(distance),
        "inclination": float(inclination),
        "quality": int(quality),
        "R_d": float(r_turn),
        "R_eff": float(1.68 * r_turn),
        "R_last": float(rad[-1]),
        "M_b": float((v_c / 100.0) ** 4),
    }


def population(v_c_values):
    """Yield flat galaxies over FLAT_SHAPES x ONSET_GRID x ``v_c_values``."""
    for shape in FLAT_SHAPES:
        for onset in ONSET_GRID:
            for v_c in v_c_values:
                yield make_galaxy("flat", shape=shape, v_c=float(v_c),
                                  onset=onset)
