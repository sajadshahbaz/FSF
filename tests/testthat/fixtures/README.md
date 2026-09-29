# Safe golden fixtures

These compact fixtures contain only historically committed numerical invariants.

`safe-golden-effects.tsv` comes from `data/processed/effect_estimates.tsv`, with
state references joined exactly by feature and perturbation from
`legacy-bundle:results/signal_states/feature_signal_states.tsv`. Seven complete identities are
selected deterministically: lexicographically first pure Up, pure Down, and pure
Constant; lexicographically first minimum-SSI identity; the next distinct
minimum-SSI tied identity; the lexicographically first identity with `0.50 < SSI < 0.75` and a unique dominant probability;
and the identity nearest SSI 0.90 within `[0.90, 1)`. Five additional rows are
the closest committed effects to `+0.5`, `-0.5`, `+1.5`, `-1.5`, and zero, with
feature and perturbation identifiers breaking ties. No exact `+/-0.5` effect
exists in the source. Historical `const` is renamed to current `constant`; this
is a semantics-preserving spelling transformation. `tau = 0.5`.

`safe-golden-identities.tsv` uses the seven corresponding rows from
`legacy-bundle:results/fsf_metrics/fsf_ssi_stability_deviation.tsv`. Probabilities, SSI,
Stability Deviation, and `n_perturbations` are copied from that committed table.
State counts are deterministic aggregations of the independently committed
row-level state table; they were not stored as columns in the metric table. All
synthetic identities have 100 perturbations, so denominator variation is not
available from this chosen independently state-backed source.

Deliberately excluded are all classification, region, dominant-state,
architecture, enrichment, and manuscript-interpretation fields.
