#!/usr/bin/env python3
"""Post-hoc robustness checks for the SPARC coherence-scale test.

These are NOT part of the registered analysis (run_pipeline.py, frozen
at tag v1.0-preregistration). They characterize the registered result
for the paper's results section by varying choices the pre-registration
fixed, reusing the frozen pipeline functions with varied parameters.
Requires SPARC data in data/.

  upsilon_sweep    r_t and its L_f / M_b correlations as Y_disk varies
  coverage         R_last / L_f -- trigger-integral data coverage
  label_stability  kinematic labels recomputed across the grid
  per_galaxy_errV  the R_flat measurability cut with each galaxy's errV
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402

from coherence import loaders  # noqa: E402
from coherence.nulls import ols_fit, spearman  # noqa: E402
from coherence.radii import RFLAT_MEASURABILITY, flat_onset_radius  # noqa: E402
from coherence.radii import transition_radius  # noqa: E402
from coherence.scales import (coherence_length_kpc, g_baryonic,  # noqa: E402
                              g_observed)
import run_pipeline as rp  # noqa: E402

DATA = ROOT / "data"
RESULTS = ROOT / "results"
PRIMARY = rp.PRIMARY


def upsilon_sweep(flat, upsilons=(0.3, 0.5, 0.7)):
    """r_t and its correlations as the disk mass-to-light ratio varies.

    Y_disk feeds g_bar, which sets r_t; the registered run used Y_disk
    = 0.5. L_f = v_c^2/a_0 is M/L-independent; M_b is the registered
    baryonic-mass proxy. If r_t keeps tracking M_b more tightly than L_f
    across the Y_disk range, the baryonic-mass dependence is real, not an
    artifact of the mass-to-light choice.
    """
    rows = []
    for yd in upsilons:
        L_f, r_t, M_b = [], [], []
        for g in flat:
            rad = g["rad"]
            g_obs = g_observed(g["vobs"], rad)
            g_bar = g_baryonic(g["vgas"], g["vdisk"], g["vbul"], rad,
                               ml_disk=yd, ml_bul=0.7)
            rt = transition_radius(rad, g_obs, g_bar, PRIMARY["ratio"],
                                   PRIMARY["persistence"])
            if rt is None:
                continue
            L_f.append(float(coherence_length_kpc(g["v_c"])))
            r_t.append(rt)
            M_b.append(g["M_b"])
        L_f, r_t, M_b = np.array(L_f), np.array(r_t), np.array(M_b)
        fit = ols_fit(L_f, r_t)
        rho_lf = spearman(r_t, L_f)[0]
        rho_mb = spearman(r_t, M_b)[0]
        rows.append({"Y_disk": yd, "n": int(len(r_t)),
                     "ols_slope_rt_Lf": fit["slope"],
                     "spearman_rt_Lf": rho_lf,
                     "spearman_rt_Mb": rho_mb,
                     "L_f_beats_M_b": bool(rho_lf > rho_mb)})
    return rows


def coverage(detail):
    """R_last / L_f: how far the data extends relative to the coherence
    scale that bounds the trigger integral."""
    chi = np.array([d["R_last"] / d["L_f"] for d in detail if d["L_f"]])
    return {"n": int(len(chi)),
            "median_Rlast_over_Lf": float(np.median(chi)),
            "min": float(chi.min()), "max": float(chi.max()),
            "frac_Lf_within_data": float(np.mean(chi >= 1.0))}


def label_stability(sample):
    """Kinematic label (flat = a finite R_flat exists) recomputed over
    every (tolerance, persistence) combination in the grid. The
    registered run froze the label at the primary cell; this measures
    whether it would have moved."""
    combos = sorted({(t, p) for t in rp.TOL_GRID for p in rp.PERSISTENCE_GRID})
    flat_counts = []
    for g in sample:
        n_flat = sum(flat_onset_radius(g["rad"], g["vobs"], g["v_c"], t, p)
                     is not None for t, p in combos)
        flat_counts.append(n_flat)
    arr = np.array(flat_counts)
    n = len(combos)
    return {"n_combinations": n, "n_galaxies": len(sample),
            "flat_in_all": int(np.sum(arr == n)),
            "rising_in_all": int(np.sum(arr == 0)),
            "label_flips": int(np.sum((arr > 0) & (arr < n)))}


def per_galaxy_errV(sample):
    """The R_flat measurability cut (errV/v_c < 0.03) with each galaxy's
    own median errV, against the registered fixed errV = 4 km/s."""
    fixed = pergal = both = 0
    for g in sample:
        rc = loaders.load_rotation_curve(DATA / (g["name"] + "_rotmod.dat"))
        errv = float(np.median(rc["e_vobs"]))
        f = (4.0 / g["v_c"]) < RFLAT_MEASURABILITY
        p = (errv / g["v_c"]) < RFLAT_MEASURABILITY
        fixed += int(f)
        pergal += int(p)
        both += int(f and p)
    return {"n": len(sample), "registered_fixed_4kms": fixed,
            "per_galaxy_errV": pergal, "in_both": both}


def main():
    if not (DATA / "SPARC_Lelli2016c.mrt").exists():
        print("No SPARC data in data/. Run scripts/fetch_data.py first.")
        sys.exit(1)

    sample = rp.quality_filtered(rp.load_sample())
    flat = [g for g in sample
            if rp.measure(g, PRIMARY["ratio"], PRIMARY["tol"],
                          PRIMARY["persistence"])["r_flat"] is not None]
    detail = json.loads((RESULTS / "primary_cell.json").read_text())

    report = {
        "note": ("Post-hoc robustness checks, NOT part of the registered "
                 "analysis (tag v1.0-preregistration)."),
        "n_quality_filtered": len(sample),
        "n_flat": len(flat),
        "upsilon_sweep": upsilon_sweep(flat),
        "coverage": coverage(detail),
        "label_stability": label_stability(sample),
        "per_galaxy_errV": per_galaxy_errV(sample),
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "robustness.json").write_text(json.dumps(report, indent=2))

    print("post-hoc robustness (%d quality-filtered, %d flat)"
          % (len(sample), len(flat)))
    print("\nUpsilon_disk sweep -- r_t correlations:")
    print("  Y_disk   n    OLS slope   rho(r_t,L_f)   rho(r_t,M_b)   L_f>M_b")
    for r in report["upsilon_sweep"]:
        print("   %.1f    %3d     %6.3f       %6.3f         %6.3f        %s"
              % (r["Y_disk"], r["n"], r["ols_slope_rt_Lf"],
                 r["spearman_rt_Lf"], r["spearman_rt_Mb"], r["L_f_beats_M_b"]))
    c = report["coverage"]
    print("\ncoverage: R_last/L_f median %.2f; %.0f%% of galaxies have L_f "
          "within the data" % (c["median_Rlast_over_Lf"],
                               100 * c["frac_Lf_within_data"]))
    ls = report["label_stability"]
    print("label stability: of %d galaxies, %d flat in all %d "
          "(tol,persistence) combos, %d rising in all, %d flip"
          % (ls["n_galaxies"], ls["flat_in_all"], ls["n_combinations"],
             ls["rising_in_all"], ls["label_flips"]))
    pe = report["per_galaxy_errV"]
    print("measurability cut: %d pass with the registered fixed errV=4 km/s, "
          "%d with per-galaxy errV" % (pe["registered_fixed_4kms"],
                                       pe["per_galaxy_errV"]))
    print("\nwritten to %s" % (RESULTS / "robustness.json"))


if __name__ == "__main__":
    main()
