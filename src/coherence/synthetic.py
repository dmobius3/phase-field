"""Synthetic SPARC-like galaxies with known properties.

Used by the test suite and by the sigma_pred Monte Carlo. No SPARC data
is involved, so these generators may run in Phase 0.
"""

from __future__ import annotations

import numpy as np

from .scales import baryonic_velocity, coherence_length_kpc

ML_DISK = 0.5  # 3.6 um stellar mass-to-light ratio used in the generator


def make_galaxy(kind: str = "flat", v_c: float = 150.0, n_points: int = 30,
                r_turn: float | None = None, span: float = 4.0,
                baryon_outer: float = 0.55, distance: float = 10.0,
                inclination: float = 60.0, quality: int = 1):
    """Construct a synthetic galaxy as a dict of arrays and scalars.

    ``kind`` is ``"flat"`` (arctan curve that plateaus at v_c) or
    ``"rising"`` (linear, still climbing at the last point). ``r_turn``
    sets the turnover scale; by default it is chosen so the flat-onset
    radius lands near L_f. ``span`` sets the outermost radius as a
    multiple of L_f. ``baryon_outer`` is v_bar / v_obs at the last point,
    which controls the transition radius.

    Velocity fields use only the disk component (gas and bulge are zero),
    so v_bar = sqrt(ML_DISK) * v_disk.
    """
    L_f = float(coherence_length_kpc(v_c))
    if r_turn is None:
        r_turn = L_f / 12.0  # arctan reaches 5% of v_c near 12 * r_turn
    r_max = span * L_f
    rad = np.linspace(r_max / n_points, r_max, n_points)

    if kind == "flat":
        vobs = v_c * (2.0 / np.pi) * np.arctan(rad / r_turn)
    elif kind == "rising":
        vobs = v_c * rad / r_max
    else:
        raise ValueError(f"unknown kind: {kind!r}")

    frac = 1.0 - (1.0 - baryon_outer) * (rad / r_max)
    vbar = vobs * frac
    vdisk = vbar / np.sqrt(ML_DISK)
    zeros = np.zeros_like(rad)

    return {
        "rad": rad,
        "vobs": vobs,
        "e_vobs": 0.05 * vobs + 2.0,
        "vgas": zeros,
        "vdisk": vdisk,
        "vbul": zeros.copy(),
        "v_c": float(v_c),
        "distance": float(distance),
        "inclination": float(inclination),
        "quality": int(quality),
        "R_d": float(r_turn),
        "R_eff": float(1.68 * r_turn),
        "R_last": float(rad[-1]),
        "M_b": float((v_c / 100.0) ** 4),
    }


def baryonic(galaxy):
    """Baryonic rotation velocity for a synthetic galaxy."""
    return baryonic_velocity(galaxy["vgas"], galaxy["vdisk"],
                             galaxy["vbul"], ml_disk=ML_DISK, ml_bul=0.7)
