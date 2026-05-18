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
from .synthetic import ML_DISK, FLAT_SHAPES, population

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


def _summary(dev):
    """Bias (mean) and scatter (std) of a deviation list."""
    dev = np.asarray(dev, dtype=float)
    if dev.size < 2:
        return {"bias": None, "scatter": None, "n": int(dev.size)}
    return {"bias": float(np.mean(dev)), "scatter": float(np.std(dev, ddof=1)),
            "n": int(dev.size)}


def sigma_pred(v_c_grid=None, n_mc: int = 400, seed: int = 0, errors=None,
               persistence: float = 0.8):
    """Estimate sigma_pred over the broadened synthetic population
    (FLAT_SHAPES x ONSET_GRID x v_c_grid).

    eta_flat statistics are restricted to galaxies passing the R_flat
    measurability cut; eta_t statistics use the full population. Reports
    bias and scatter both pooled and per shape, so the stability of the
    measurement bias across curve morphology can be judged.
    """
    errors = SPARC_ERRORS if errors is None else errors
    if v_c_grid is None:
        v_c_grid = np.linspace(90.0, 280.0, 8)
    rng = np.random.default_rng(seed)
    dev_t = {s: [] for s in FLAT_SHAPES}
    dev_flat = {s: [] for s in FLAT_SHAPES}
    n_measurable = 0
    n_total = 0

    for galaxy in population(v_c_grid):
        shape = galaxy["shape"]
        n_total += 1
        measurable = rflat_measurable(errors["errV_kms"], galaxy["v_c"])
        n_measurable += int(measurable)
        eta_t0, eta_flat0 = _eta_pair(
            galaxy["rad"], galaxy["vobs"], galaxy["vgas"],
            galaxy["vdisk"], galaxy["vbul"], galaxy["v_c"],
            persistence=persistence)
        for _ in range(n_mc):
            eta_t, eta_flat = _eta_pair(*_perturb(galaxy, rng, errors),
                                        persistence=persistence)
            if eta_t is not None and eta_t0:
                dev_t[shape].append(eta_t / eta_t0 - 1.0)
            if measurable and eta_flat is not None and eta_flat0:
                dev_flat[shape].append(eta_flat / eta_flat0 - 1.0)

    all_t = [d for s in FLAT_SHAPES for d in dev_t[s]]
    all_flat = [d for s in FLAT_SHAPES for d in dev_flat[s]]
    by_shape = {s: {"eta_t": _summary(dev_t[s]),
                    "eta_flat": _summary(dev_flat[s])} for s in FLAT_SHAPES}
    flat_biases = [by_shape[s]["eta_flat"]["bias"] for s in FLAT_SHAPES
                   if by_shape[s]["eta_flat"]["bias"] is not None]
    pooled_flat = _summary(all_flat)
    pooled_t = _summary(all_t)

    return {
        "sigma_pred": pooled_flat["scatter"],
        "sigma_pred_eta_flat": pooled_flat["scatter"],
        "sigma_pred_eta_t": pooled_t["scatter"],
        "bias_eta_flat": pooled_flat["bias"],
        "bias_eta_t": pooled_t["bias"],
        "bias_eta_flat_by_shape": {s: by_shape[s]["eta_flat"]["bias"]
                                   for s in FLAT_SHAPES},
        "bias_eta_flat_shape_spread": ([float(min(flat_biases)),
                                        float(max(flat_biases))]
                                       if flat_biases else None),
        "by_shape": by_shape,
        "n_realizations_eta_t": int(len(all_t)),
        "n_realizations_eta_flat": int(len(all_flat)),
        "n_mc_per_galaxy": int(n_mc),
        "seed": int(seed),
        "persistence": float(persistence),
        "rflat_measurability_cut": float(RFLAT_MEASURABILITY),
        "n_galaxies_total": int(n_total),
        "n_galaxies_rflat_measurable": int(n_measurable),
        "v_c_grid_kms": [float(v) for v in v_c_grid],
        "error_budget": dict(errors),
        "note": ("Fractional eta scatter from a representative SPARC "
                 "error budget over a broadened synthetic population; "
                 "equals the absolute scatter at eta = 1. sigma_pred is "
                 "the pooled eta_flat scatter over measurable galaxies."),
    }
