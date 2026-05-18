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
- Exclude from the radius-correlation tests (predictions 1-3): galaxies
  with no finite `R_flat`, i.e. those not labeled `flat` by the
  kinematic-label definition below. Rising-curve galaxies are retained
  as a separate sub-population for the threshold-prediction test
  (prediction 4).
- R_flat analysis only: include a galaxy in the R_flat / eta_flat
  analysis only if `errV/v_c < 0.03`, where `errV = 4 km/s` is a fixed
  representative per-point velocity error (Lelli et al. 2016). With the
  fixed `errV` this is a pre-data cut, equivalent to `v_c > 133 km/s`,
  that selects the registered sub-sample without touching the data.
  Per-galaxy `errV` is deliberately not used for the registered cut:
  that would make sample selection data-dependent and unwind the blind
  analysis. The per-galaxy-`errV` version is run only as a post-hoc
  robustness check, labeled as such. R_flat is biased and unreliable
  when per-point noise approaches the 5% flatness tolerance, as the
  sigma_pred Monte Carlo showed; SPARC errV/v_c is largest in the
  low-v_c tail, so this cut trims low-mass dwarfs from the eta_flat
  analysis only. It does not affect r_t or the threshold-prediction
  test, which use the full quality-filtered sample. The trimmed
  eta_flat sub-sample size is reported.

## Operational definitions

- **Baryonic mass model**: the baryonic rotation velocity is
  `v_bar^2 = |V_gas|V_gas + Y_disk|V_disk|V_disk + Y_bul|V_bul|V_bul`,
  combining the SPARC rotmod components (tabulated at M/L = 1) with
  stellar mass-to-light ratios `Y_disk = 0.5` and `Y_bul = 0.7`
  Msun/Lsun at 3.6 um (McGaugh et al. 2016). The signed-square form
  preserves the sign convention of the catalogue components. The
  baryonic acceleration is `g_bar = v_bar^2 / r`, which feeds the `r_t`
  definition below. `Y_disk` and `Y_bul` are registered values, not
  fitted.
- **Transition radius `r_t`**: smallest radius such that
  `g_obs(r)/g_bar(r) >= 1.2` holds at r_t and at no less than the
  persistence fraction `f` of measured points at `r >= r_t`. `r_t` is
  NOT defined by `g_obs = a_0` (that would make `r_t = L_f` an
  arithmetic identity for flat curves).
- **Flat-onset radius `R_flat`**: smallest radius such that
  `|v(r) - v_c|/v_c <= 0.05` holds at R_flat and at no less than the
  persistence fraction `f` of measured points at `r >= R_flat`.
- **Persistence rule**: the onset is the smallest radius where the
  condition holds at that radius and at a fraction `f` of the points
  beyond it, not at every point. A strict "all points beyond" rule is a
  one-sided ratchet: a single point displaced by measurement error trips
  it and pushes the onset outward, never inward, biasing r_t and R_flat
  high. The persistence rule removes the ratchet while keeping the
  conceptual definition of a sustained transition. Primary `f = 0.80`,
  swept in the sensitivity grid; the same `f` applies to r_t and R_flat.
- **Trigger index** `T/T_c`: computed for every galaxy from its rotation
  curve, with no flat/rising pre-classification.
- **Kinematic label (flat vs. rising)**: a galaxy is labeled `flat` if a
  finite `R_flat` exists within the measured radial range, evaluated at
  the primary-cell parameters (5% flatness tolerance, persistence
  `f = 0.80`), and its SPARC quality flag is 1 or 2; otherwise `rising`.
  This label is the ground truth for prediction 4. It is derived from
  the rotation curve and the SPARC quality flag only, never from
  `T/T_c`, and is assigned before any ROC curve is computed. Using the
  trigger index to set the label would make prediction 4 circular.

## Registered predictions

1. **Transition radius.** `r_t` correlates linearly with `L_f`.
2. **Flat-onset radius.** `R_flat` correlates linearly with `L_f`, and
   tracks `L_f` more tightly than it tracks any null model below.
3. **Primary outcome: dimensionless ratios.** `eta_t = r_t/L_f` and
   `eta_flat = R_flat/L_f` cluster near unity, with scatter explained by
   measurement error.
4. **Trigger index as predictor.** `T/T_c`, computed without kinematic
   pre-classification, predicts which galaxies have flat versus rising
   rotation curves.

Each prediction has a quantitative pass/fail criterion in the next
section. The words "approximately," "near," and "stable" are defined
numerically there and nowhere else.

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

Every interval and test statistic below is frozen at the tag. The numeric
values are provisional in this draft and may be adjusted before the tag;
once frozen, nothing changes after data contact.

**Prediction 1 (transition radius).** Ordinary-least-squares fit of `r_t`
on `L_f`: registered slope interval `[0.7, 1.3]`, registered intercept
consistent with zero at 2 sigma of the fit. Spearman correlation
significant at `p < 0.01`. Both `r_t` and `L_f` carry measurement error,
so OLS is subject to regression dilution that biases the slope toward
zero (typical SPARC uncertainties pull a true unit slope to roughly
0.85). OLS is the registered method, and the `[0.7, 1.3]` interval is set
wide enough to absorb this bias. Orthogonal-distance regression is
reported as a robustness check, not as the registered gate.

**Prediction 2 (flat-onset radius).** OLS fit of `R_flat` on `L_f`:
registered slope interval `[0.7, 1.3]`, intercept consistent with zero at
2 sigma. AND, for the null-model comparison: the partial correlation of
`R_flat` with `L_f`, controlling for each null model in turn, remains
significant at `p < 0.01`; and the Spearman rho of `R_flat` with `L_f`
exceeds that of every null model. The OLS regression-dilution note in
prediction 1 applies here as well.

**Prediction 3 (dimensionless ratios, primary outcome).** Two independent
criteria, both required to pass:

- *Location.* Primary precision window: the sample median of `eta_t`,
  and of `eta_flat`, lies within `[0.75, 1.25]`, and the 95% bootstrap
  confidence interval of each median also lies within `[0.75, 1.25]`.
  This tests that the velocity-to-length constant is `1/a_0`
  specifically. A coarse sanity gate, median within `[0.5, 2.0]`, is
  also recorded; passing the coarse gate but failing the precision gate
  means the ratio is order unity but not `1/a_0`. The precision gate
  carries the falsification weight.

  The window is asymmetric, wider below 1.0, because the measurement
  bias is one-sided: the persistence-rule onset is pulled inward for
  gradual (slow) turnovers and is barely shifted for sharp ones. The
  broadened-population Monte Carlo characterized this pre-data, finding
  an `eta_flat` bias of -16.6% for gradual turnovers, -3.2% for bumpy
  curves, and +2.4% for sharp turnovers. The 0.75 floor absorbs the
  worst-case single-shape bias with margin; the 1.25 ceiling absorbs
  only the small positive sharp-turnover bias.
- *Scatter.* A one-sided chi-squared test of the sample variance of each
  `eta` against `sigma_pred^2`, with `sigma_pred = 0.443`. The criterion
  fails if the observed scatter exceeds `sigma_pred` at `p < 0.05`. This
  tests that measurement error explains the spread. `sigma_pred` is the
  pooled `eta_flat` scatter from the broadened-population Monte Carlo,
  equal-weighted across the sharp, gradual, and bumpy shapes (per-shape
  scatter 0.20, 0.30, 0.67; all recorded in
  `registration/sigma_pred.json`). The equal-weight pooling is the
  documented pre-data choice; a bumpy-heavy real sample marginally
  exceeding 0.443 is a shape-mix discussion, not a theory failure.

**Prediction 4 (trigger index as predictor).** The continuous index
`T/T_c` is scored against the binary kinematic label (flat vs. rising;
defined operationally above and frozen independently of `T/T_c`) by
ROC AUC. Registered criterion: `AUC >= 0.7`. Reported alongside, but not
as the pass/fail metric: the 2x2 contingency table at the `T/T_c = 1`
decision boundary with a Fisher exact test.

**Threshold (closure identity).** Among quality-filtered flat-curve
galaxies, no more than 5% may fall below `T/T_c = 1`. Exceeding 5%
falsifies the closure identity.

`sigma_pred = 0.443`, computed in Phase 0 from the broadened-population
representative SPARC error budget (four rotation-curve morphologies over
a turnover-scale grid) and frozen in `registration/sigma_pred.json`. The
value is rounded to three significant figures: the Monte Carlo produced
0.44256, and its roughly 10% uncertainty from the finite synthetic
population and the equal-weight shape assumption does not warrant more
precision.

## Sensitivity grid

The three analyst-chosen parameters are swept over a 3x3x3 grid:

- `r_t` divergence ratio: 1.1, 1.2 (primary), 1.3.
- `R_flat` flatness tolerance: 3%, 5% (primary), 7%.
- persistence fraction `f`: 0.70, 0.80 (primary), 0.90.

Registered stability criterion: the pass/fail verdict of every prediction
above is unchanged across all twenty-seven cells of the grid. If any
prediction flips its verdict between cells, the result is reported as
threshold-sensitive and is not claimed as a confirmation, regardless of
what the primary cell (1.2, 5%, 0.80) shows.

Kinematic labels (flat vs. rising) are frozen at the primary-cell
parameters and held constant across all twenty-seven cells. `R_flat`
existence depends on the flatness tolerance and the persistence fraction;
letting the labels move with the grid would shift prediction 4's ground
truth, so the ROC AUC could change because the labels moved rather than
because the trigger index did. With the labels fixed, the grid tests the
trigger index, predictions 1-3, and the closure-identity count against a
stable target.

## Falsification

Component-level: a failure isolates whether `L_f`, the closure identity,
the value of `a_0`, or the single-`L_f` picture is the broken piece. The
failure-mode table is in the companion working note, section V.
