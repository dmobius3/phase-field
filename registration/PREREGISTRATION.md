# Pre-Registration

**Status:** DRAFT. Finalized and frozen at git tag `v1.0-preregistration`.

This document fixes the analysis before any SPARC data is downloaded. The
pipeline code, the thresholds below, the sample cuts, and the acceptance
scatter bound `sigma_pred` are all frozen at the tag. Nothing in the
pipeline is fitted to data.

This is a locked-pipeline blind analysis of pre-existing data. SPARC has
been public since 2016, so the tag defends against one objection only:
that algorithm thresholds were tuned to the data. It does not claim the
data was unseen. The genuine pre-registered test against unseen data is
Euclid DR1.

## Hypothesis

The Mode Identity Theory coherence scale `L_f = v_c^2/a_0`, with
`a_0 = 1.2e-10 m/s^2`, is a real galactic coherence radius: the
gravitational transition radius and the flat-onset radius track `L_f`
specifically, more tightly than they track generic galactic size scalings.

## Sample cuts

Applied uniformly, from the SPARC catalog:

- Include: quality flag 1 or 2.
- Include: inclination >= 30 degrees.
- Exclude: rising or falling rotation curves at the last measured point
  (for the radius-correlation tests; rising-curve galaxies are retained
  as a separate sub-population for the threshold-prediction test).

## Operational definitions

- **Transition radius `r_t`**: smallest radius such that
  `g_obs(r)/g_bar(r) >= 1.2` for all measured points at `r >= r_t`.
  `r_t` is NOT defined by `g_obs = a_0` (that would make `r_t = L_f`
  an arithmetic identity for flat curves).
- **Flat-onset radius `R_flat`**: smallest radius such that
  `|v(r) - v_c|/v_c <= 0.05` for all measured points at `r >= R_flat`.
- **Trigger index** `T/T_c`: computed for every galaxy from its rotation
  curve, with no flat/rising pre-classification.

## Registered predictions

1. **Transition radius.** `r_t` correlates linearly with `L_f`, slope
   approximately 1.
2. **Flat-onset radius.** `R_flat` correlates linearly with `L_f`, slope
   approximately 1, intercept consistent with zero, AND tracks `L_f` more
   tightly than it tracks any null model below (partial correlation).
3. **Primary outcome: dimensionless ratios.** `eta_t = r_t/L_f` and
   `eta_flat = R_flat/L_f` cluster around order unity with scatter
   consistent with `sigma_pred`.
4. **Trigger index as predictor.** `T/T_c`, computed without kinematic
   pre-classification, predicts which galaxies have flat versus rising
   rotation curves.

## Null models

`r_t` and `R_flat` are fit against each of the following with the same
procedure as `L_f`. `L_f` must beat them, `v_c^2` in particular, since
`L_f` is proportional to `v_c^2` and only the constant `1/a_0` carries
physical content.

- Size proxies: `R_d` (disk scale length), `R_eff` (effective radius),
  `R_last` (last measured rotation-curve point).
- Velocity scalings: `v_c`, `v_c^1.5`, `v_c^2`.
- Mass scalings: `M_b^(1/3)`, `M_b^(1/2)` (baryonic mass).

## Acceptance criteria

- **Ratios (primary).** `eta_t` and `eta_flat` are consistent with
  clustering at unity within `sigma_pred`.
  `sigma_pred = TBD` (computed in Phase 0 from a representative SPARC
  error budget; see `registration/sigma_pred.json`, frozen at the tag).
- **Threshold.** Among quality-filtered flat-curve galaxies, no more than
  5% may fall below `T/T_c = 1`. Exceeding 5% falsifies the closure
  identity.

## Sensitivity grid

The two analyst-chosen thresholds are swept; the registered conclusions
must be stable across the grid:

- `r_t` divergence ratio: 1.1, 1.2 (primary), 1.3.
- `R_flat` flatness tolerance: 3%, 5% (primary), 7%.

## Falsification

Component-level: a failure isolates whether `L_f`, the closure identity,
the value of `a_0`, or the single-`L_f` picture is the broken piece. The
failure-mode table is in the companion working note, section V.
