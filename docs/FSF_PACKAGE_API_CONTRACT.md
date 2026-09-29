# FSF Package API Contract

## Authority and scope

This contract is subordinate to `docs/FSF_V1_CURRENT_SCIENTIFIC_LOCK.md`.
FSF v1 performs signal stratification, Signal Identity characterization,
SSI-based response-consistency characterization, signal classification, and
condition-level architecture calculation. It does not perform effect
estimation, differential-expression analysis, feature selection, feature
ranking, enrichment, plotting, prediction, or species-specific analysis.

## Package identity

- Package name: `FSF`.
- Current package version: `0.1.0`.
- License: MIT, as recorded by `DESCRIPTION` and the package `LICENSE` file.
- Computational dependency target: base R only, with no Imports.
- Testing dependency: `testthat (>= 3.0.0)` in Suggests.

The package version matches the `fsf-r-v0.1.0` release-candidate tag and makes
no CRAN-release claim.

## Canonical input

The package-core input is a long-form data frame with these required semantic
columns:

| Column | Contract |
|---|---|
| `feature_id` | Nonmissing, nonblank feature identifier |
| `perturbation_id` | Nonmissing, nonblank perturbation identifier |
| `effect` | Numeric, nonmissing, and finite signed effect |

Arbitrary character effects are not coerced to numeric. Missing rows are not
dropped. `tau` must be one numeric, finite scalar strictly greater than zero;
its default is `0.5`. A zero-row input is invalid.

Features and grouping strata need not share a common number of perturbations.
Each Signal Identity uses its own observed and validated `n_perturbations` as
the probability denominator.

## Generic grouping contract

`group_cols = NULL` calculates one Signal Identity per `feature_id` over the
supplied input. `character(0)` may be normalized to `NULL`.

Otherwise, `group_cols` is a character vector naming one or more input columns.
Signal Identity is calculated independently within `group_cols + feature_id`.
For example, `group_cols = "condition"` treats the same feature under UV and
DES as two independent Signal Identities. This exposes the condition-wise
behavior of the historical real-data pipeline generically; it is not a new FSF
scientific rule. No stress-condition name is hard-coded.

The order of `group_cols` supplied by the user is retained in output. The
argument must not contain missing or blank names, duplicate names,
`feature_id`, or `perturbation_id`; every requested column must exist. Grouping
values must be nonmissing, and character grouping values must be nonblank.

The unique observation key is:

- `feature_id + perturbation_id` when `group_cols = NULL`;
- `group_cols + feature_id + perturbation_id` when grouping is supplied.

Duplicate keys error and are never aggregated or double-weighted. The same
`feature_id` and `perturbation_id` may legitimately occur in different groups.

Identity rows are returned in first-appearance order: grouping values and
features retain the order in which their unique combinations first occur in
the input. Fixed input order therefore produces fixed output order; FSF does
not impose a lexical feature or condition sort.

## Public API

### `fsf_assign_states(data, tau = 0.5)`

Validates canonical signed-effect input and returns one state per input row.
State values are `up`, `down`, and `constant`. Up is `effect > tau`, Down is
`effect < -tau`, and Constant is `-tau <= effect <= tau`; exact boundaries are
Constant.

The function preserves all additional input columns and replaces an existing
`state` column with the state calculated from `effect`.

### `fsf_signal_identity(data, tau = 0.5, group_cols = NULL)`

Returns one row per `group_cols + feature_id` identity, or per `feature_id`
when ungrouped, with `n_perturbations`, state counts, and
`(p_up, p_down, p_const)`. It may accept validated states or canonical signed
effects; the implementation must distinguish these forms explicitly and apply
the same validation rules.

When both `state` and `effect` are present, the validated `state` column is the
operative input and `effect` is ignored. Only grouping columns, `feature_id`,
and computed count/probability columns are propagated; other row-level columns
are discarded during aggregation.

### `fsf_classify(identity)`

Validates Signal Identity and derives `dominant_state`, `ssi`,
`stability_region`, `signal_class`, and `stability_deviation`.

This function operates row-wise on an already aggregated identity table. The
caller is responsible for supplying one row per intended identity because the
function cannot infer which arbitrary columns constitute grouping keys. It
retains additional input columns and appends or replaces the derived
classification columns.

Probability-only identity input remains valid. If any count metadata is
supplied, all of `n_perturbations`, `n_up`, `n_down`, and `n_const` must be
present. Counts must be finite nonnegative integers, `n_perturbations` must be
strictly positive, component counts must sum to it, and count-derived
probabilities must agree with supplied probabilities within `1e-8`.

### `fsf_analyze(data, tau = 0.5, group_cols = NULL)`

The generic high-level composition of state assignment, grouped Signal
Identity, and current-lock classification. It does not estimate effects.

Because its output is aggregated, non-grouping row-level columns are not
propagated. Identity rows use the same first-appearance ordering as
`fsf_signal_identity()`.

### `fsf_architecture(data, condition_col = "condition")`

Consumes already classified feature-level results and calculates signal-class
counts and proportions independently within the column named by
`condition_col`.

The output contains only the condition column, `signal_class`, `class_count`,
and `class_proportion`. Conditions retain first-appearance order; signal
classes always use the fixed ten-class order.

## Probability and numerical contract

Each probability must be numeric, finite, and in `[0, 1]`. A probability
triple is accepted only when:

`abs(p_up + p_down + p_const - 1) <= 1e-8`

The `1e-8` tolerance is validation-only. Probabilities are not renormalized,
clipped, rounded, or snapped to a boundary. The tolerance is not used to create
ties or change a classification.

Classification uses the raw, unrounded validated values. SSI is exactly
`max(p_up, p_down, p_const)`, followed by these exact comparisons:

- `ssi <= 0.50`: Low Stability.
- `0.50 < ssi < 0.75`: Transitional.
- `0.75 <= ssi < 0.90`: Stable.
- `ssi >= 0.90`: Highly Stable.

No classification tolerance applies. Tie detection uses equality of validated
probability values, with no fuzzy or epsilon rule. Thus `(0.50, 0.50, 0)` is a
tie and Low Stability; `(0.50, 0.30, 0.20)` has dominant state `up` but remains
Low Stability; and `(0.5000000001, 0.4999999999, 0)` has unique dominant state
`up` and class `Transitional Up`.

## Single-perturbation interpretation

For `n_perturbations = 1`, state assignment and Signal Identity remain valid,
the probabilities are necessarily one-hot, SSI is necessarily 1, and the
mathematical class is Highly Stable in the corresponding direction.

**With one perturbation, SSI = 1 is structurally determined and does not
constitute evidence of cross-perturbation stability.**

This interpretation safeguard does not alter SSI, suppress classification,
create a class, or convert output to missing. `n_perturbations` remains visible
so evidence depth is explicit. No additional evidence or confidence output
column is introduced in this contract.

## Internal helpers

A later implementation may introduce unexported helpers limited to canonical
input validation, grouping validation, `tau` validation, probability
validation, composite-key detection, exact unique-maximum/tie detection,
stability-region and signal-class mapping, deterministic factors, and
floating-point invariant checks. These remain internal. Historical 0.60 modes,
ranking, selection, plotting, enrichment, and species-specific helpers are
outside package core.

## Feature-level output

Without grouping, stable output columns are:

1. `feature_id`
2. `n_perturbations`
3. `n_up`
4. `n_down`
5. `n_const`
6. `p_up`
7. `p_down`
8. `p_const`
9. `dominant_state`
10. `ssi`
11. `stability_region`
12. `signal_class`
13. `stability_deviation`

With grouping, the user-supplied `group_cols` appear first in supplied order,
then the columns above. `dominant_state` values are `up`, `down`, `constant`,
and `tied`.

Ordered stability-region levels are `Low Stability`, `Transitional`, `Stable`,
and `Highly Stable`. Ordered signal-class levels are:

1. `Low Stability`
2. `Transitional Up`
3. `Transitional Constant`
4. `Transitional Down`
5. `Stable Up`
6. `Stable Constant`
7. `Stable Down`
8. `Highly Stable Up`
9. `Highly Stable Constant`
10. `Highly Stable Down`

The factors provide deterministic presentation order, not feature ranking.
Counts remain integer. Probabilities, SSI, and Stability Deviation remain
unrounded; display rounding is outside scientific computation.

## Architecture contract

`condition_col` must be one nonblank character scalar naming an existing
column whose values are all nonmissing. Character condition values must also be
nonblank. Input must contain `feature_id` and `signal_class`, and every
`condition_col + feature_id` combination must be unique.

The output preserves the supplied condition-column name and contains:

1. `<condition_col>`
2. `signal_class`
3. `class_count`
4. `class_proportion`

Every condition includes all ten current FSF classes in deterministic order,
including zero counts. Within each condition, counts sum to its classified
feature count and proportions sum to one within floating-point tolerance.

Architectures containing `n_perturbations = 1` are permitted. **An architecture
derived exclusively from single-perturbation identities describes directional
class composition, not demonstrated reproducibility across repeated
perturbations.** This affects interpretation, not mathematical assignment.

## Error contract

Invalid scientific inputs produce deterministic errors, not warnings. Required
categories include missing required columns; nonnumeric, missing, or nonfinite
effects; missing or blank identifiers; duplicate mathematical keys; invalid
`tau`; zero rows; invalid probabilities or counts; invalid `group_cols` names
or grouping values; and invalid `condition_col`, condition values, or duplicate
condition-feature classifications. No invalid input is silently repaired.

## Historical behavior

The intentional historical 0.60 definition is provenance, not a selectable
package mode. It must not be exposed as an FSF v1 option. Historical
classifications are not current-lock golden results.

## Installed documentation

The five exported functions have roxygen-generated installed help pages with
executable, dataset-independent examples. Documentation is generated from the
public API comments in `R/`; the generated `man/` files and `NAMESPACE` are
release artifacts and must remain synchronized with those comments.
