"""Coherence-scale analysis pipeline for the SPARC phase-field test.

Importable modules:
    loaders    parse the SPARC master table and rotation-curve files
    scales     L_f, observed and baryonic accelerations
    radii      transition radius r_t and flat-onset radius R_flat
    ratios     dimensionless ratios eta_t, eta_flat
    threshold  trigger index T / T_c
    nulls      null-model fits and partial-correlation controls
    propagate  Monte Carlo for the acceptance scatter bound sigma_pred
"""

__version__ = "0.1.0"
