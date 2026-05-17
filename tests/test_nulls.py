import numpy as np

from coherence.nulls import beats_nulls, ols_fit, partial_spearman, spearman
from coherence.scales import A0


def test_ols_recovers_known_line():
    x = np.linspace(1.0, 10.0, 50)
    fit = ols_fit(x, 2.0 * x + 1.0)
    assert np.isclose(fit["slope"], 2.0)
    assert np.isclose(fit["intercept"], 1.0)


def test_partial_spearman_removes_spurious_correlation():
    # x and y are both driven by z; controlling for z removes the link.
    rng = np.random.default_rng(3)
    z = rng.normal(0.0, 1.0, 300)
    x = z + rng.normal(0.0, 0.1, 300)
    y = z + rng.normal(0.0, 0.1, 300)
    rho_raw, _ = spearman(x, y)
    rho_partial, _ = partial_spearman(x, y, z)
    assert rho_raw > 0.9
    assert abs(rho_partial) < 0.3


def test_vc_squared_null_is_decisive():
    # radius == L_f exactly: L_f is proportional to v_c^2, so the v_c**2
    # null ties with L_f and must NOT be counted as beaten.
    v_c = np.linspace(50.0, 250.0, 40)
    L_f = v_c**2 / A0
    catalog = {
        "v_c": v_c,
        "M_b": (v_c / 100.0) ** 4,
        "R_d": np.linspace(1.0, 8.0, 40),
        "R_eff": np.linspace(2.0, 14.0, 40),
        "R_last": np.linspace(5.0, 30.0, 40),
    }
    report = beats_nulls(L_f.copy(), L_f, catalog)
    assert "_overall" in report and "_rho_Lf" in report
    assert report["v_c**2"]["passed"] is False
    assert report["_overall"] is False
