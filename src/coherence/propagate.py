"""Monte Carlo estimate of sigma_pred, the eta scatter expected from
SPARC measurement error alone.

A representative error budget (not per-galaxy SPARC values) is applied to
synthetic flat galaxies: distance, inclination, and an absolute per-point
velocity error errV. Because errV is absolute, the fractional velocity
error errV / v_c is larger for low-v_c galaxies, as in SPARC (Lelli,
McGaugh & Schombert 2016: typical errV is 2-5 km/s for quality 1-2
HI/Halpha rotation curves).

R_flat is only well-defined where per-point noise sits comfortably below
the flatness tolerance, so the eta_flat statistic is computed only over
galaxies passing the measurability cut errV / v_c < RFLAT_MEASURABILITY.
r_t is well-behaved throughout and uses the full sample.

Run in Phase 0. The result is written to registration/sigma_pred.json and
frozen at the pre-registration tag. No SPARC data is used.
"""

from __future__ import annotations

import numpy as np

from .radii import (RFLAT_MEASURABILITY, flat_onset_radius, rflat_measurable,
                    transition_radius)
from .scales import coherence_length_kpc, g_baryonic, g_observed
from .synthetic import ML_DISK, make_galaxy

# Representative SPARC uncertainties (Lelli, McGaugh & Schombert 2016).
# errV is an absolute per-point velocity error; the flat velocity is
# re-measured from the perturbed outer curve rather than perturbed
# directly.
SPARC_ERRORS = {
    "distance_frac": 0.10,    # fractional distance uncertainty (radius scale)
    "inclination_deg": 5.0,   # typical inclination uncertainty (degrees)
    "errV_kms": 4.0,          # representative absolute per-point velocity error
}


def _flat_velocity(rad, vobs, outer_frac: float = 0.34):
    """Flat velocity as the mean of the outer rotation-curve points."""
    k = max(3, int(round(outer_frac * len(rad))))
    return float(np.mean(np.asarray(vobs)[-k:]))


def _eta_pair(rad, vobs, vgas, vdisk, vbul, v_c,
              ratio_threshold: float = 1.2, flat_tol: float = 0.05,
              persistence: float = 0.8):
    L_f = float(coherence_length_kpc(v_c))
    g_obs = g_observed(vobs, rad)
    g_bar = g_baryonic(vgas, vdisk, vbul, rad, ml_disk=ML_DISK, ml_bul=0.7)
    r_t = transition_radius(rad, g_obs, g_bar, ratio_threshold, persistence)
    r_flat = flat_onset_radius(rad, vobs, v_c, flat_tol, persistence)
    eta_t = None if r_t is None else r_t / L_f
    eta_flat = None if r_flat is None else r_flat / L_f
    return eta_t, eta_flat


def _perturb(galaxy, rng, errors):
    """One perturbed realization of a synthetic galaxy."""
    i0 = galaxy["inclination"]
    i1 = float(np.clip(i0 + rng.normal(0.0, errors["inclination_deg"]),
                       15.0, 89.0))
    inc = np.sin(np.radians(i0)) / np.sin(np.radians(i1))
    dist = max(1.0 + rng.normal(0.0, errors["distance_frac"]), 0.5)
    noise = rng.normal(0.0, errors["errV_kms"], size=len(galaxy["rad"]))
    rad = galaxy["rad"] * dist
    vobs = galaxy["vobs"] * inc + noise
    vdisk = galaxy["vdisk"] * inc
    vgas = galaxy["vgas"] * inc
    vbul = galaxy["vbul"] * inc
    v_c = _flat_velocity(rad, vobs)
    return rad, vobs, vgas, vdisk, vbul, v_c


def sigma_pred(v_c_grid=None, n_mc: int = 400, seed: int = 0, errors=None,
               persistence: float = 0.8):
    """Estimate sigma_pred over a grid of synthetic galaxies.

    eta_flat statistics are restricted to galaxies passing the R_flat
    measurability cut; eta_t statistics use the full grid. Reports both
    scatter (std) and bias (mean) of the fractional eta deviation.
    """
    errors = SPARC_ERRORS if errors is None else errors
    if v_c_grid is None:
        v_c_grid = np.linspace(60.0, 280.0, 12)
    rng = np.random.default_rng(seed)
    dev_t, dev_flat = [], []
    n_measurable = 0

    for v_c in v_c_grid:
        galaxy = make_galaxy("flat", v_c=float(v_c))
        measurable = rflat_measurable(errors["errV_kms"], float(v_c))
        n_measurable += int(measurable)
        eta_t0, eta_flat0 = _eta_pair(
            galaxy["rad"], galaxy["vobs"], galaxy["vgas"],
            galaxy["vdisk"], galaxy["vbul"], galaxy["v_c"],
            persistence=persistence)
        for _ in range(n_mc):
            eta_t, eta_flat = _eta_pair(*_perturb(galaxy, rng, errors),
                                        persistence=persistence)
            if eta_t is not None and eta_t0:
                dev_t.append(eta_t / eta_t0 - 1.0)
            if measurable and eta_flat is not None and eta_flat0:
                dev_flat.append(eta_flat / eta_flat0 - 1.0)

    dev_t = np.asarray(dev_t, dtype=float)
    dev_flat = np.asarray(dev_flat, dtype=float)
    return {
        "sigma_pred": float(np.std(dev_flat, ddof=1)),
        "sigma_pred_eta_t": float(np.std(dev_t, ddof=1)),
        "sigma_pred_eta_flat": float(np.std(dev_flat, ddof=1)),
        "bias_eta_t": float(np.mean(dev_t)),
        "bias_eta_flat": float(np.mean(dev_flat)),
        "n_realizations_eta_t": int(dev_t.size),
        "n_realizations_eta_flat": int(dev_flat.size),
        "n_mc_per_galaxy": int(n_mc),
        "seed": int(seed),
        "persistence": float(persistence),
        "rflat_measurability_cut": float(RFLAT_MEASURABILITY),
        "n_galaxies_total": int(len(v_c_grid)),
        "n_galaxies_rflat_measurable": int(n_measurable),
        "v_c_grid_kms": [float(v) for v in v_c_grid],
        "error_budget": dict(errors),
        "note": ("Fractional eta scatter from a representative SPARC "
                 "error budget; equals the absolute scatter at eta = 1. "
                 "sigma_pred is the eta_flat scatter over measurable "
                 "galaxies, the binding case."),
    }
