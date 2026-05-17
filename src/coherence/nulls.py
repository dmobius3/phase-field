"""Null-model comparison and partial-correlation controls.

L_f must beat a family of generic galactic scalings (PREREGISTRATION.md,
prediction 2). ``v_c**2`` is the decisive null: since L_f is proportional
to v_c^2, a correlation with L_f no tighter than the correlation with
v_c^2 carries no information about a0.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

NULL_MODELS = ("R_d", "R_eff", "R_last", "v_c", "v_c**1.5", "v_c**2",
               "M_b**(1/3)", "M_b**(1/2)")


def null_predictor(name, catalog):
    """Predictor array for a named null model, from a catalog mapping
    with keys R_d, R_eff, R_last, v_c, M_b."""
    v_c = np.asarray(catalog["v_c"], dtype=float)
    m_b = np.asarray(catalog["M_b"], dtype=float)
    table = {
        "R_d": np.asarray(catalog["R_d"], dtype=float),
        "R_eff": np.asarray(catalog["R_eff"], dtype=float),
        "R_last": np.asarray(catalog["R_last"], dtype=float),
        "v_c": v_c,
        "v_c**1.5": v_c**1.5,
        "v_c**2": v_c**2,
        "M_b**(1/3)": np.cbrt(m_b),
        "M_b**(1/2)": np.sqrt(m_b),
    }
    if name not in table:
        raise KeyError(f"unknown null model: {name!r}")
    return table[name]


def ols_fit(x, y):
    """Ordinary-least-squares fit y = slope * x + intercept."""
    res = stats.linregress(np.asarray(x, dtype=float),
                           np.asarray(y, dtype=float))
    return {"slope": float(res.slope), "intercept": float(res.intercept),
            "slope_se": float(res.stderr),
            "intercept_se": float(res.intercept_stderr),
            "r": float(res.rvalue), "p": float(res.pvalue)}


def spearman(x, y):
    """Spearman rank correlation, returns ``(rho, p)``."""
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def partial_spearman(x, y, z):
    """Partial correlation of x and y controlling for z, on ranks.

    Ranks of x and y are linearly residualised on the rank of z and the
    residuals correlated. Returns ``(rho, p)``.
    """
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    rho, p = stats.pearsonr(ex, ey)
    return float(rho), float(p)


def beats_nulls(radius, L_f, catalog, partial_p_max: float = 0.01):
    """Prediction-2 null-model criterion.

    For every null model: the partial correlation of ``radius`` with
    ``L_f`` controlling for that null must stay significant, and the
    Spearman rho of ``radius`` with ``L_f`` must exceed the rho with the
    null. Returns a per-null report plus an overall pass flag.
    """
    radius = np.asarray(radius, dtype=float)
    L_f = np.asarray(L_f, dtype=float)
    rho_lf, _ = spearman(radius, L_f)
    report = {"_rho_Lf": rho_lf}
    overall = True
    for name in NULL_MODELS:
        z = null_predictor(name, catalog)
        rho_null, _ = spearman(radius, z)
        prho, pp = partial_spearman(radius, L_f, z)
        passed = (pp < partial_p_max) and (rho_lf > rho_null)
        report[name] = {"rho_null": rho_null, "partial_rho": prho,
                        "partial_p": pp, "passed": bool(passed)}
        overall = overall and passed
    report["_overall"] = bool(overall)
    return report
