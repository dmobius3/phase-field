"""Coherence length and accelerations.

All functions are unit-consistent: pass SI in, get SI out. The two
``*_kpc`` / ``*_kms`` helpers handle the SPARC catalogue units explicitly.
"""

from __future__ import annotations

import numpy as np

C_LIGHT = 299_792_458.0                  # m / s
A0 = 1.2e-10                             # m / s^2, local edge-mode scale
KPC_M = 3.085_677_581_491_367_4e19       # m per kpc
KM_S = 1.0e3                             # m/s per km/s


def coherence_length(v_c, a0: float = A0):
    """L_f = v_c^2 / a0. ``v_c`` in m/s, returns metres."""
    return np.square(np.asarray(v_c, dtype=float)) / a0


def coherence_length_kpc(v_c_kms, a0: float = A0):
    """L_f in kpc from a flat velocity ``v_c_kms`` in km/s."""
    v_c = np.asarray(v_c_kms, dtype=float) * KM_S
    return coherence_length(v_c, a0) / KPC_M


def g_observed(v, r):
    """Centripetal acceleration v^2 / r. Units must be consistent."""
    return np.square(np.asarray(v, dtype=float)) / np.asarray(r, dtype=float)


def baryonic_velocity(v_gas, v_disk, v_bul, ml_disk: float = 0.5,
                      ml_bul: float = 0.7):
    """SPARC baryonic rotation velocity.

    ``v_bar^2 = |v_gas| v_gas + Y_d |v_disk| v_disk + Y_b |v_bul| v_bul``,
    where the component velocities are the Upsilon = 1 contributions from
    the SPARC ``*_rotmod.dat`` files and ``Y_d``, ``Y_b`` are the 3.6 um
    stellar mass-to-light ratios. The signed-square form preserves the
    sign convention of the catalogue components.
    """
    v_gas = np.asarray(v_gas, dtype=float)
    v_disk = np.asarray(v_disk, dtype=float)
    v_bul = np.asarray(v_bul, dtype=float)
    sq = (np.abs(v_gas) * v_gas
          + ml_disk * np.abs(v_disk) * v_disk
          + ml_bul * np.abs(v_bul) * v_bul)
    return np.sign(sq) * np.sqrt(np.abs(sq))


def g_baryonic(v_gas, v_disk, v_bul, r, ml_disk: float = 0.5,
               ml_bul: float = 0.7):
    """Newtonian acceleration from the baryonic mass model."""
    v_bar = baryonic_velocity(v_gas, v_disk, v_bul, ml_disk, ml_bul)
    return g_observed(v_bar, r)
