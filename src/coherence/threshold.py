"""Phase-field trigger index T / T_c.

Operational definition. The trigger index is

    T = (2 / (c^2 L_f)) * integral_0^{L_f} Phi_rel(l) dl,

with ``Phi_rel(l) = v(l)^2``, the squared circular velocity at radius l.
For a flat rotation curve v == v_c this gives T = 2 v_c^2 / c^2, and with
T_c = 2 xi v_c^2 / c^2 it reproduces the closure-identity result
T / T_c = 1 / xi. For a rising curve the inner v(l) < v_c, so the index
falls below 1 / xi: the index genuinely tracks curve shape rather than
restating it.

The choice Phi_rel(l) = v(l)^2 is the operational reading of the
phase-field potential integral. It should be checked against the MIT
framework before the pre-registration tag is cut.
"""

from __future__ import annotations

import numpy as np

from .scales import C_LIGHT, KM_S

XI = 0.46  # geometry factor; 0.44-0.47 across isothermal/NFW/Hernquist


def trigger_index(rad_kpc, v_kms, L_f_kpc, v_c_kms, n_grid: int = 2048):
    """Trigger index T (dimensionless) for one galaxy.

    Inside the innermost measured radius v is taken linear from the
    origin; beyond the outermost measured point v is held flat at v_c.
    """
    rad = np.asarray(rad_kpc, dtype=float)
    v = np.asarray(v_kms, dtype=float) * KM_S
    v_c = float(v_c_kms) * KM_S
    order = np.argsort(rad)
    rad, v = rad[order], v[order]

    grid = np.linspace(0.0, float(L_f_kpc), n_grid)
    v_grid = np.interp(grid, rad, v, left=v[0], right=v_c)
    inner = grid < rad[0]
    v_grid[inner] = v[0] * grid[inner] / rad[0]

    integral = np.trapz(v_grid**2, grid)              # (m/s)^2 * kpc
    return 2.0 * integral / (C_LIGHT**2 * float(L_f_kpc))


def trigger_critical(v_c_kms, xi: float = XI):
    """Critical trigger index T_c = 2 xi v_c^2 / c^2 (dimensionless)."""
    v_c = float(v_c_kms) * KM_S
    return 2.0 * xi * v_c**2 / C_LIGHT**2


def trigger_ratio(rad_kpc, v_kms, L_f_kpc, v_c_kms, xi: float = XI,
                  n_grid: int = 2048):
    """T / T_c. A flat curve gives 1 / xi; rising curves give less."""
    T = trigger_index(rad_kpc, v_kms, L_f_kpc, v_c_kms, n_grid)
    return T / trigger_critical(v_c_kms, xi)
