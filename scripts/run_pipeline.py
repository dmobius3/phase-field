#!/usr/bin/env python3
"""Run the frozen coherence pipeline over the 27-cell sensitivity grid.

Phase 1 code: requires SPARC data in data/ (download with
``scripts/fetch_data.py --confirm-data-contact``). Writes results/, which
is gitignored. The registered claim is that every prediction's pass/fail
verdict is stable across all 27 grid cells, so the output records all 27.
"""

import json
import sys
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

from coherence import loaders  # noqa: E402
from coherence.nulls import (beats_nulls, odr_slope, ols_fit,  # noqa: E402
                             roc_auc, spearman)
from coherence.radii import (flat_onset_radius, rflat_measurable,  # noqa: E402
                             transition_radius)
from coherence.ratios import median_within, scatter_chi2  # noqa: E402
from coherence.scales import (coherence_length_kpc, g_baryonic,  # noqa: E402
                              g_observed)
from coherence.threshold import trigger_ratio  # noqa: E402

DATA = ROOT / "data"
RESULTS = ROOT / "results"

# Registered grid and acceptance numbers (registration/PREREGISTRATION.md).
RATIO_GRID = (1.1, 1.2, 1.3)
TOL_GRID = (0.03, 0.05, 0.07)
PERSISTENCE_GRID = (0.70, 0.80, 0.90)
PRIMARY = {"ratio": 1.2, "tol": 0.05, "persistence": 0.80}

SIGMA_PRED = 0.443
SLOPE_INTERVAL = (0.7, 1.3)
PRECISION_WINDOW = (0.75, 1.25)
COARSE_GATE = (0.5, 2.0)
AUC_MIN = 0.7
CLOSURE_MAX_BELOW = 0.05  # max fraction of flat-curve galaxies with T/T_c < 1
ERRV_KMS = 4.0  # registered fixed representative per-point error (km/s).
# The R_flat measurability cut uses this fixed value, not per-galaxy errV,
# so the sample selection stays pre-data (see PREREGISTRATION.md sample cuts).
MIN_FIT = 10    # minimum galaxies for a correlation fit

# SPARC master-table column names. VERIFY against the .mrt header in
# Phase 1; the machine-readable table fixes the exact spellings.
COL = {"name": "Galaxy", "vflat": "Vflat", "quality": "Q", "inc": "Inc",
       "rdisk": "Rdisk", "reff": "Reff", "mhi": "MHI", "lum": "L[3.6]"}


def load_sample():
    """Load SPARC into per-galaxy records."""
    master = loaders.load_master(DATA / "SPARC_Lelli2016c.mrt")
    by_name = {str(r[COL["name"]]).strip(): r for _, r in master.iterrows()}
    records = []
    for path in loaders.rotation_curve_files(DATA):
        name = path.name.replace("_rotmod.dat", "")
        row = by_name.get(name)
        if row is None:
            continue
        rc = loaders.load_rotation_curve(path)
        rad = rc["rad"].to_numpy(dtype=float)
        records.append({
            "name": name, "rad": rad,
            "vobs": rc["vobs"].to_numpy(dtype=float),
            "vgas": rc["vgas"].to_numpy(dtype=float),
            "vdisk": rc["vdisk"].to_numpy(dtype=float),
            "vbul": rc["vbul"].to_numpy(dtype=float),
            "v_c": float(row[COL["vflat"]]),
            "quality": int(row[COL["quality"]]),
            "inclination": float(row[COL["inc"]]),
            "R_d": float(row[COL["rdisk"]]),
            "R_eff": float(row[COL["reff"]]),
            "R_last": float(rad[-1]),
            # M_b: gas (1.33 M_HI) plus a single-Upsilon stellar mass
            # (0.5 x total 3.6um light). A luminosity-based baryonic-mass
            # PROXY used only for the M_b null-model scalings, not the
            # g_bar mass model, which decomposes disk and bulge with
            # Y_disk = 0.5 and Y_bul = 0.7 separately.
            "M_b": 1.33 * float(row[COL["mhi"]]) + 0.5 * float(row[COL["lum"]]),
        })
    return records


def quality_filtered(records):
    """Quality flag 1 or 2, inclination >= 30 degrees, positive v_c."""
    return [g for g in records if g["quality"] in (1, 2)
            and g["inclination"] >= 30.0 and g["v_c"] > 0]


def measure(galaxy, ratio, tol, persistence):
    """Per-galaxy radii and ratios at one grid cell."""
    rad, vobs = galaxy["rad"], galaxy["vobs"]
    L_f = float(coherence_length_kpc(galaxy["v_c"]))
    g_obs = g_observed(vobs, rad)
    g_bar = g_baryonic(galaxy["vgas"], galaxy["vdisk"], galaxy["vbul"], rad)
    r_t = transition_radius(rad, g_obs, g_bar, ratio, persistence)
    r_flat = flat_onset_radius(rad, vobs, galaxy["v_c"], tol, persistence)
    return {"L_f": L_f, "r_t": r_t, "r_flat": r_flat,
            "eta_t": None if r_t is None else r_t / L_f,
            "eta_flat": None if r_flat is None else r_flat / L_f}


def kinematic_labels(sample):
    """flat (True) if a finite R_flat exists at the primary parameters.

    Frozen at the primary cell and independent of the trigger index, so
    prediction 4 is not circular (PREREGISTRATION.md).
    """
    return {g["name"]: measure(g, PRIMARY["ratio"], PRIMARY["tol"],
                               PRIMARY["persistence"])["r_flat"] is not None
            for g in sample}


def _line_prediction(x, y):
    """Slope in the registered interval, intercept ~0 at 2 sigma,
    significant Spearman correlation. OLS is the registered gate; the
    ODR slope is recorded alongside as a robustness check."""
    if len(x) < MIN_FIT:
        return {"passed": False, "slope": None, "odr_slope": None}
    x, y = np.asarray(x, float), np.asarray(y, float)
    fit = ols_fit(x, y)
    _, p = spearman(x, y)
    slope_ok = SLOPE_INTERVAL[0] <= fit["slope"] <= SLOPE_INTERVAL[1]
    intercept_ok = abs(fit["intercept"]) <= 2.0 * fit["intercept_se"]
    return {"passed": bool(slope_ok and intercept_ok and p < 0.01),
            "slope": fit["slope"], "intercept": fit["intercept"],
            "spearman_p": float(p), "odr_slope": odr_slope(x, y)}


def _ratio_prediction(eta_t, eta_flat):
    """Location within the precision window and scatter consistent with
    sigma_pred, for both eta_t and eta_flat. The coarse gate [0.5, 2.0]
    is recorded but does not affect the pass/fail verdict; the precision
    window carries the falsification weight."""
    if len(eta_t) < MIN_FIT or len(eta_flat) < MIN_FIT:
        return {"passed": False}
    loc_t = median_within(eta_t, PRECISION_WINDOW)[0]
    loc_f = median_within(eta_flat, PRECISION_WINDOW)[0]
    coarse_t = median_within(eta_t, COARSE_GATE)[0]
    coarse_f = median_within(eta_flat, COARSE_GATE)[0]
    p_t = scatter_chi2(eta_t, SIGMA_PRED)[2]
    p_f = scatter_chi2(eta_flat, SIGMA_PRED)[2]
    return {"passed": bool(loc_t and loc_f and p_t >= 0.05 and p_f >= 0.05),
            "location_eta_t": bool(loc_t), "location_eta_flat": bool(loc_f),
            "coarse_eta_t": bool(coarse_t), "coarse_eta_flat": bool(coarse_f),
            "scatter_p_eta_t": float(p_t), "scatter_p_eta_flat": float(p_f)}


def evaluate_cell(sample, labels, triggers, ratio, tol, persistence):
    """Evaluate the four predictions at one grid cell."""
    rows = {g["name"]: measure(g, ratio, tol, persistence) for g in sample}
    flat = [g for g in sample if labels[g["name"]]]

    p1_x = [rows[g["name"]]["L_f"] for g in flat
            if rows[g["name"]]["r_t"] is not None]
    p1_y = [rows[g["name"]]["r_t"] for g in flat
            if rows[g["name"]]["r_t"] is not None]
    pred1 = _line_prediction(p1_x, p1_y)

    measurable = [g for g in flat if rows[g["name"]]["r_flat"] is not None
                  and rflat_measurable(ERRV_KMS, g["v_c"])]
    p2_x = [rows[g["name"]]["L_f"] for g in measurable]
    p2_y = [rows[g["name"]]["r_flat"] for g in measurable]
    line2 = _line_prediction(p2_x, p2_y)
    if len(measurable) >= MIN_FIT:
        catalog = {k: np.array([g[k] for g in measurable])
                   for k in ("v_c", "M_b", "R_d", "R_eff", "R_last")}
        nulls_ok = beats_nulls(np.array(p2_y), np.array(p2_x),
                               catalog)["_overall"]
    else:
        nulls_ok = False
    pred2 = bool(line2["passed"] and nulls_ok)

    eta_t = [rows[g["name"]]["eta_t"] for g in flat
             if rows[g["name"]]["eta_t"] is not None]
    eta_flat = [rows[g["name"]]["eta_flat"] for g in measurable]
    pred3 = _ratio_prediction(eta_t, eta_flat)

    scores = [triggers[g["name"]] for g in sample]
    labs = [labels[g["name"]] for g in sample]
    auc = roc_auc(scores, labs)
    pred4 = bool(auc >= AUC_MIN)

    return {"ratio": ratio, "tol": tol, "persistence": persistence,
            "prediction_1": pred1["passed"], "prediction_2": pred2,
            "prediction_3": pred3["passed"], "prediction_4": pred4,
            "n_flat": len(flat), "n_rflat_measurable": len(measurable),
            "slope_r_t": pred1["slope"], "slope_R_flat": line2["slope"],
            "odr_slope_r_t": pred1["odr_slope"],
            "odr_slope_R_flat": line2["odr_slope"],
            "coarse_eta_t": pred3.get("coarse_eta_t"),
            "coarse_eta_flat": pred3.get("coarse_eta_flat"),
            "auc": auc}


def write_primary_detail(sample, labels, triggers):
    """Per-galaxy detail at the primary cell, for the figures."""
    rows = []
    for g in sample:
        m = measure(g, PRIMARY["ratio"], PRIMARY["tol"], PRIMARY["persistence"])
        rows.append({
            "name": g["name"], "v_c": g["v_c"], "L_f": m["L_f"],
            "r_t": m["r_t"], "R_flat": m["r_flat"], "eta_t": m["eta_t"],
            "eta_flat": m["eta_flat"], "trigger_ratio": triggers[g["name"]],
            "flat": labels[g["name"]],
            "rflat_measurable": rflat_measurable(ERRV_KMS, g["v_c"]),
            "R_d": g["R_d"], "R_eff": g["R_eff"], "R_last": g["R_last"],
            "M_b": g["M_b"]})
    (RESULTS / "primary_cell.json").write_text(json.dumps(rows, indent=2))


def closure_and_contingency(sample, labels, triggers):
    """Closure-identity test and the prediction-4 contingency table.

    Neither depends on the sensitivity-grid cell: the trigger index and
    the kinematic labels are fixed. The closure test counts flat-curve
    galaxies below T/T_c = 1 (registered: no more than 5%). The 2x2
    contingency table at the T/T_c = 1 boundary and its Fisher exact
    test are reported alongside the prediction-4 AUC.
    """
    names = [g["name"] for g in sample]
    trig = np.array([triggers[n] for n in names], dtype=float)
    flat = np.array([labels[n] for n in names], dtype=bool)

    flat_trig = trig[flat]
    frac_below = (float(np.mean(flat_trig < 1.0)) if flat_trig.size
                  else float("nan"))
    closure = {"n_flat_curve": int(flat.sum()),
               "fraction_below_unity": frac_below,
               "max_allowed": CLOSURE_MAX_BELOW,
               "passed": bool(frac_below <= CLOSURE_MAX_BELOW)}

    predicted_flat = trig >= 1.0
    table = [[int(np.sum(~predicted_flat & ~flat)),
              int(np.sum(predicted_flat & ~flat))],
             [int(np.sum(~predicted_flat & flat)),
              int(np.sum(predicted_flat & flat))]]
    _, fisher_p = fisher_exact(table)
    contingency = {"table": table, "rows": ["rising", "flat"],
                   "columns": ["predicted_rising", "predicted_flat"],
                   "fisher_exact_p": float(fisher_p),
                   "auc": float(roc_auc(trig, flat))}
    return closure, contingency


def main():
    if not (DATA / "SPARC_Lelli2016c.mrt").exists():
        print("No SPARC data in data/. run_pipeline.py is Phase 1 code.")
        print("Run: python scripts/fetch_data.py --confirm-data-contact")
        sys.exit(1)

    RESULTS.mkdir(exist_ok=True)
    sample = quality_filtered(load_sample())
    print(f"quality-filtered sample: {len(sample)} galaxies")
    labels = kinematic_labels(sample)
    triggers = {g["name"]: trigger_ratio(
        g["rad"], g["vobs"], coherence_length_kpc(g["v_c"]), g["v_c"])
        for g in sample}
    closure, contingency = closure_and_contingency(sample, labels, triggers)

    grid = [evaluate_cell(sample, labels, triggers, r, t, p)
            for r, t, p in product(RATIO_GRID, TOL_GRID, PERSISTENCE_GRID)]
    keys = ("prediction_1", "prediction_2", "prediction_3", "prediction_4")
    stable = {k: len({c[k] for c in grid}) == 1 for k in keys}

    summary = {"n_cells": len(grid), "primary": PRIMARY,
               "verdict_stable": stable, "all_stable": all(stable.values()),
               "closure_identity": closure,
               "prediction_4_contingency": contingency,
               "cells": grid}
    (RESULTS / "sensitivity_grid.json").write_text(json.dumps(summary, indent=2))
    write_primary_detail(sample, labels, triggers)

    primary = next(c for c in grid if c["ratio"] == PRIMARY["ratio"]
                   and c["tol"] == PRIMARY["tol"]
                   and c["persistence"] == PRIMARY["persistence"])
    print("primary cell: "
          + ", ".join(f"P{i + 1}={primary[keys[i]]}" for i in range(4)))
    print(f"closure identity: {closure['fraction_below_unity']:.3f} of "
          f"flat-curve galaxies below T/T_c=1 (pass={closure['passed']})")
    print(f"prediction 4: AUC={contingency['auc']:.3f}, "
          f"Fisher exact p={contingency['fisher_exact_p']:.3g}")
    print(f"all 27-cell verdicts stable: {summary['all_stable']}")
    print(f"results written to {RESULTS}")


if __name__ == "__main__":
    main()
