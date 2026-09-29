# FSF Package Test Plan

## Test strategy

Level 1 uses small hand-calculated fixtures to protect formulas, exact
boundaries, exact tie semantics, grouping, validation, output schema, and
invariants. Level 2 compares safe numeric quantities against committed
historical references without importing historical 0.60 class semantics.

Safe Level 2 quantities are effects, states, state counts, probabilities,
Signal Identity, SSI, Stability Deviation, and SSI-only ECDF inputs. Historical
regions, classes, architectures, class-specific enrichment, and benchmark
predictions are legacy references, not current goldens. Large datasets are not
copied into this scaffold.

## Grouping contracts

Implemented tests establish that:

- The same `feature_id + perturbation_id` in distinct conditions is accepted
  with `group_cols = "condition"`.
- A duplicate `condition + feature_id + perturbation_id` within one condition
  errors.
- The same feature in two conditions produces two independent identities.
- Each identity uses its own actual perturbation count as denominator.
- Multiple grouping columns are returned before `feature_id` in user order.
- Missing columns, missing or blank argument names, duplicated names, reserved
  identifiers, missing grouping values, and blank character grouping values
  error. `character(0)` may normalize to `NULL`.

## Numerical contracts

Each probability is finite and in `[0,1]`. Probability-sum error at or below
`1e-8` is accepted; error above `1e-8` is rejected. This is validation-only:
there is no normalization, clipping, rounding, boundary snapping, fuzzy tie, or
classification tolerance.

Classification uses raw validated probabilities. Tests cover effect values at
and immediately outside `+/-tau`, and SSI values `1/3`, `0.49`, `0.50`, just
above `0.50`, just below and exactly `0.75`, just below and exactly `0.90`, and
`1`.

Tie tests establish that `(0.50, 0.50, 0)` is tied Low Stability, while
`(0.50, 0.30, 0.20)` reports dominant `up` but remains Low Stability.
`(0.5000000001, 0.4999999999, 0)` is not tied and is Transitional Up. No
rounding occurs before dominant-state, SSI, region, or class evaluation.
Mathematically, two equal maxima cannot exceed `0.50` because their sum would
exceed one.

## Probability and count contracts

Tests assert that counts sum to `n_perturbations`, probabilities sum to one
under the validation rule, SSI is the maximum raw probability, Stability
Deviation is `1 - SSI`, and denominators may differ by feature and grouping
stratum.

## Single-perturbation contracts

Implemented fixtures cover one Up, Constant, and Down perturbation. Each
produces
`n_perturbations = 1`, a one-hot identity, SSI 1, and the corresponding Highly
Stable directional class.

**With one perturbation, SSI = 1 is structurally determined and does not
constitute evidence of cross-perturbation stability.**

These are mathematically valid classifications, but they do not demonstrate
cross-perturbation reproducibility. No class, formula, or output column is
added. Architecture tests permit such identities while requiring documentation
that an exclusively single-perturbation architecture describes directional
composition rather than demonstrated reproducibility.

## Current-lock regression fixtures

Fixtures include an exact-half tie, exact-half unique plurality, UV-like `4/7`,
exact `0.75`, exact `0.90`, pure state, and near-uniform `1/3` identity.
Directional variants cover Up, Down, and Constant where applicable.

## Validation contracts

Implemented tests require errors for `NA`, `NaN`, `Inf`, `-Inf`, character
effects,
blank or missing IDs, duplicate composite keys, invalid `tau`, empty input,
missing required columns, invalid probability ranges or sums, inconsistent
counts, invalid grouping specifications and values, and invalid architecture
condition specifications and duplicate condition-feature rows. Warnings are
insufficient for invalid scientific inputs.

## Architecture contracts

For each value of the supplied `condition_col`, output preserves that column
name and includes all ten signal classes in fixed order, including zero counts.
Counts equal the number of classified features and proportions equal one within
tolerance. Conditions are never pooled and duplicate condition-feature rows
error.

## Output contracts

Tests fix column names and ordering, including grouping columns before
`feature_id`; dominant-state values; ordered region and class levels; integer
counts; identifier preservation; and deterministic ordering. No evidence-level
or replication flag is part of this contract.

## Current implementation behavior

The independent hand-calculated computational contracts are executable for the
implemented mathematical core. The Level 2 safe historical golden-reference
tests are implemented with committed, classification-excluding fixtures.

## Release-hardening contracts

Implemented package tests additionally protect behavior that was previously
implicit:

- supplied count metadata requires a positive `n_perturbations`, integer-like
  nonnegative counts, correct component sums, and count-derived probabilities
  agreeing with supplied probabilities within `1e-8`;
- probability-only `fsf_classify()` input remains valid and classification is
  row-wise on caller-aggregated identities;
- when `state` and `effect` coexist, validated `state` takes precedence in
  `fsf_signal_identity()`;
- extra-column preservation and discard behavior follows each documented
  public output schema;
- identity and condition outputs preserve first-appearance order and are
  deterministic for fixed input order;
- all five exported functions have executable, external-file-free examples.

Isolated build, install, fresh-session API smoke, and `R CMD check` remain
release-process gates. Their presence in this plan does not claim that they
have been run for the current working tree.
