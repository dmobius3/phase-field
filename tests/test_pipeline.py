import numpy as np

from coherence.propagate import sigma_pred
from coherence.radii import flat_onset_radius, transition_radius
from coherence.ratios import eta
from coherence.scales import (coherence_length_kpc, g_baryonic, g_observed)
from coherence.synthetic import ML_DISK, make_galaxy


def test_full_pipeline_on_synthetic_flat_galaxy():
    g = make_galaxy("flat", v_c=180.0)
    L_f = coherence_length_kpc(g["v_c"])
    g_obs = g_observed(g["vobs"], g["rad"])
    g_bar = g_baryonic(g["vgas"], g["vdisk"], g["vbul"], g["rad"],
                       ml_disk=ML_DISK, ml_bul=0.7)
    r_t = transition_radius(g["rad"], g_obs, g_bar)
    r_flat = flat_onset_radius(g["rad"], g["vobs"], g["v_c"])
    assert r_t is not None and r_flat is not None
    assert eta(r_t, L_f) > 0.0
    assert eta(r_flat, L_f) > 0.0


def test_sigma_pred_runs_and_returns_valid_structure():
    # Checks mechanics only. The magnitude of sigma_pred depends on the
    # provisional error budget and is a Phase 0 deliverable, not a unit
    # test target.
    out = sigma_pred(v_c_grid=np.array([100.0, 200.0]), n_mc=40, seed=0)
    assert np.isfinite(out["sigma_pred"]) and out["sigma_pred"] > 0.0
    assert out["n_realizations_eta_flat"] > 0
    assert out["n_realizations_eta_t"] > 0
    assert set(out) >= {"sigma_pred_eta_t", "sigma_pred_eta_flat",
                        "error_budget", "v_c_grid_kms"}
