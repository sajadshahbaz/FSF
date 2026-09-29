# FSF v1 Current Scientific Lock

## 1. Framework Scope

FSF v1 is:

- Feature Signal Stratification.
- A framework for characterizing perturbation-derived Signal Identity.
- A framework for quantifying response consistency using SSI.
- A framework for assigning interpretable signal classes.
- A framework for constructing condition-level signal architectures.

FSF v1 is not:

- Feature selection.
- Feature ranking.
- A replacement for differential-expression analysis.
- Causal inference.
- Predictive classification.

Signal-aware feature selection is explicitly future work and outside FSF v1.

## 2. Mathematical Input

For feature \(i\) and perturbation \(j\),

\[
E(i,j)
\]

is a signed perturbation-derived effect estimate.

Examples may include:

- Log fold change.
- Regression coefficient.
- Standardized signed effect.
- An analogous signed quantitative response.

FSF is estimator-independent once the signed effect estimates have been constructed.

## 3. Directional State Assignment

Let \(\tau\) be a positive finite scalar. The default is:

\[
\tau = 0.5
\]

State assignment is:

**Up:**

\[
E(i,j) > \tau
\]

**Down:**

\[
E(i,j) < -\tau
\]

**Constant:**

\[
-\tau \le E(i,j) \le \tau
\]

Values exactly equal to \(+\tau\) or \(-\tau\) are Constant.

## 4. State Probabilities

For feature \(i\) over \(m\) perturbations:

\[
P_{up}(i) = \frac{N_{up}(i)}{m}
\]

\[
P_{down}(i) = \frac{N_{down}(i)}{m}
\]

\[
P_{const}(i) = \frac{N_{const}(i)}{m}
\]

with:

\[
P_{up} + P_{down} + P_{const} = 1
\]

## 5. Signal Identity

Signal Identity of feature \(i\) is:

\[
\left(P_{up}(i), P_{down}(i), P_{const}(i)\right)
\]

Signal Identity is the fundamental FSF representation.

## 6. Dominant State

The dominant state is the unique argmax of:

- \(P_{up}\)
- \(P_{down}\)
- \(P_{const}\)

when a unique maximum exists.

If the maximum is tied, the dominant state is recorded as tied/undefined.

No artificial precedence such as Up > Down > Constant is permitted.

Dominant-state reporting and directional signal-class assignment are separate concepts.

A feature may have a unique plurality at SSI = 0.50 but still lack strict majority support.

## 7. Signal Stratification Index

\[
SSI(i) = \max\left(P_{up}(i), P_{down}(i), P_{const}(i)\right)
\]

Range:

\[
\frac{1}{3} \le SSI \le 1
\]

SSI quantifies concentration of Signal Identity in its largest component.

## 8. Current Stability Regions

The current authoritative rule is:

**Low Stability:**

\[
SSI \le 0.50
\]

**Transitional:**

\[
0.50 < SSI < 0.75
\]

**Stable:**

\[
0.75 \le SSI < 0.90
\]

**Highly Stable:**

\[
SSI \ge 0.90
\]

Scientific interpretation:

- **Low Stability:** No signal state possesses strict majority support.
- **Transitional:** A signal state possesses strict majority support, but alternative states remain substantial.
- **Stable:** A dominant state has strong support.
- **Highly Stable:** Signal Identity is strongly concentrated near a simplex vertex.

Low Stability denotes low concentration of Signal Identity. It must not be
interpreted as synonymous with biological instability or with the separate
Instability scientific framework or metric.

The theoretical SSI minimum remains \(1/3\), attained when:

\[
P_{up} = P_{down} = P_{const} = \frac{1}{3}
\]

The exact SSI = 0.50 boundary belongs to Low Stability. This includes both:

- Tied maxima at 0.50.
- Unique 0.50 pluralities.

Neither possesses strict majority support.

## 9. Ten Current FSF Signal Classes

The only current FSF v1 signal classes are:

1. Highly Stable Up
2. Highly Stable Constant
3. Highly Stable Down
4. Stable Up
5. Stable Constant
6. Stable Down
7. Transitional Up
8. Transitional Constant
9. Transitional Down
10. Low Stability

Low Stability is deliberately non-directional.

No Low Stability Up, Low Stability Down, Low Stability Constant, Low Stability Tied, Mixed, or Transitional Mixed class exists in current FSF v1.

All records with \(SSI \le 0.50\) belong to the single non-directional signal
class Low Stability. `dominant_state` may independently record `up`, `down`,
`constant`, or `tied`, but it does not create a directional Low-Stability
signal class. Directional FSF signal classes are authorized only when
\(SSI > 0.50\).

The current v1 identifier policy is:

- Display/category label: `Low Stability`.
- Case-preserving gene-set identifier token: `Low_Stability`.
- Lowercase `low_stability` is not a current v1 gene-set identifier token.

This preserves the capitalization convention of the current gene-set
sanitizer. The schema column names remain `stability_region` and
`signal_class`.

The term Instability remains valid for the separate Instability scientific
framework or metric, for historical FSF provenance where that terminology was
actually used, and for any future SIAFS Instability objective or module. Those
uses must not be renamed merely because the current FSF categorical region is
named Low Stability.

## 10. Condition-Level Signal Architecture

For condition \(c\) and signal class \(k\):

\[
A(c,k) = \frac{N(c,k)}{\sum_k N(c,k)}
\]

The architecture vector is the set of class proportions within a condition.

Condition-level architecture is derived from current FSF class membership.

## 11. Stability Deviation

Historical FSF workflows report:

\[
\text{Stability Deviation} = 1 - SSI
\]

This quantity may remain available as an auxiliary/descriptive output for compatibility and interpretation.

It is not:

- An independent classifier.
- A separate stability definition.
- A criterion for signal-class assignment.

SSI remains authoritative for stability-region assignment.

## 12. Input Validation Policy for Package Core

Current package-core implementation must reject:

- Non-numeric effects.
- `NA`.
- `NaN`.
- `Inf`.
- `-Inf`.
- \(\tau \le 0\).
- Non-finite \(\tau\).
- Duplicate feature-perturbation observations.

No missing or non-numeric value may silently become Constant.

A package implementation must not silently double-weight duplicate feature-perturbation pairs.

## 13. Historical 0.60 Definition

Historical FSF used:

**Instable:**

\[
SSI < 0.60
\]

**Weakly Stable / Transitional:**

\[
0.60 \le SSI < 0.75
\]

**Stable:**

\[
0.75 \le SSI < 0.90
\]

**Highly Stable:**

\[
SSI \ge 0.90
\]

This was an intentional historical scientific design. It was documented, implemented, benchmarked, and used in historical manuscript outputs. It must not be described as a coding typo.

For current computation it is superseded by the stability rules in Section 8.

Historical outputs must be retained as legacy provenance and must not be silently overwritten or presented as current-lock results.

## 14. Migration History

Repository audit established that later manuscript-facing code adopted modern terminology and a 0.50 conceptual boundary without fully recomputing historical 0.60-derived memberships.

In particular, presentation code mapped:

- `weakly_stable_transitional` to Transitional.
- `instable` to Instability.

without recomputing the full underlying class membership.

This is classified as:

**INCOMPLETE IMPLEMENTATION MIGRATION**

It does not establish that the historical definition was erroneous.

## 15. Reason for the <= 0.50 Clarification

The current manuscript used:

- Instability: \(SSI < 0.50\).
- Transitional: \(0.50 \le SSI < 0.75\).

while describing Instability as absence of majority support and Transitional as emergence of a dominant directional signal.

Package-conversion auditing identified an exact-boundary ambiguity.

At SSI = 0.50:

- No state possesses strict majority support.
- Exact two-way ties may occur.
- Unique pluralities may also occur.
- Arbitrary tie precedence would violate the intended probabilistic interpretation.

Therefore the current rule is clarified to:

- Low Stability: \(SSI \le 0.50\).
- Transitional: \(SSI > 0.50\) and \(SSI < 0.75\).

This is a formal current-definition clarification discovered during package validation.

Do not claim that this exact `<=` rule was historically used.

## 16. Empirical Consequence of the Current Rule

Repository audit established:

**Real dataset:** 91,122 condition-feature classifications.

**Exactly SSI = 0.50:** 10,065 records, of which:

- 8,588 are two-way top ties.
- 1,477 have a unique plurality but no strict majority.

The strict interval \(0.50 < SSI < 0.60\) contains 2,451 records.

All 2,451 are UV features with:

\[
SSI = \frac{4}{7} = 0.5714285714285714
\]

and current destinations:

- Transitional Up: 432.
- Transitional Constant: 1,624.
- Transitional Down: 395.

No affected record is directionally unresolved.

Only UV condition-level class architecture changes relative to the historical 0.60 implementation.

## 17. Current Condition-Level Stability Architecture

Under the current rule:

| Condition | Low Stability | Transitional | Stable | Highly Stable |
|---|---:|---:|---:|---:|
| DES | 5,837 (38.434%) | 0 | 0 | 9,350 (61.566%) |
| GAM | 1,085 (7.144%) | 2,106 (13.867%) | 6,127 (40.344%) | 5,869 (38.645%) |
| HT | 0 | 0 | 0 | 15,187 (100%) |
| LT | 1,211 (7.974%) | 0 | 0 | 13,976 (92.026%) |
| OSM | 0 | 0 | 0 | 15,187 (100%) |
| UV | 2,948 (19.411%) | 6,697 (44.096%) | 3,151 (20.748%) | 2,391 (15.744%) |

Interpretation:

- DES remains a mixed Highly-Stable/Low-Stability architecture.
- LT remains Highly-Stable dominated.
- HT remains entirely Highly Stable.
- OSM remains entirely Highly Stable.
- GAM remains multi-layer.
- UV remains multi-layer but is Transitional-dominated.

Do not retain the historical manuscript claim that UV has the highest Low-Stability or historically named categorical Instability burden.

Do not retain the manuscript claim that Low Stability is the largest UV region/class.

This terminology change preserves exactly:

- Feature and condition membership.
- SSI and Stability Deviation.
- \(P_{up}\), \(P_{down}\), and \(P_{const}\).
- `dominant_state` and signal-state assignments.
- All threshold boundaries.
- Class and architecture counts.
- Gene-set member vectors and annotation joins.
- GO and KEGG statistics.

The Low-Stability counts remain DES 5,837; GAM 1,085; HT 0; LT 1,211;
OSM 0; and UV 2,948, for a total of 11,081. UV remains numerically:
Low Stability 2,948; Transitional 6,697; Stable 3,151; and Highly Stable
2,391. UV is therefore Transitional-dominated under the current FSF
classification.

## 18. Synthetic Benchmark Migration

Do not redesign synthetic generating distributions.

Historical truth labels include:

- `weakly_stable`.
- `instable`.

These labels are historical provenance. Current evaluation must use the current FSF classifier.

Audit identified 71 synthetic features with \(0.50 < SSI < 0.60\). All have:

- Historical truth scenario: `weakly_stable`.
- Unique dominant state: Up.
- Historical predicted class: `instable`.
- Current predicted class: Transitional Up.

One synthetic feature has exactly SSI = 0.50 and is currently Low Stability.

Current benchmark migration therefore requires classifier/evaluator recomputation, not redesign of the generating distributions.

## 19. Reference Output Policy

Safe numerical references include, where applicable:

- Signed effect estimates.
- State assignments.
- State counts.
- \(P_{up}\).
- \(P_{down}\).
- \(P_{const}\).
- Signal Identity.
- SSI.
- Stability Deviation = \(1 - SSI\).
- SSI-only ECDF results.

Legacy references include:

- Historical `stability_level`.
- Historical dominant/class outputs where tied behavior is involved.
- Historical 0.60-derived signal classes.
- Historical Figure 3 architecture.
- Historical class-specific enrichment.
- Historical benchmark class predictions.

Legacy references must not be used as current class goldens.

Current-lock class goldens must be generated only after implementation of this lock and independent validation.

## 20. Manuscript/Result Recomputation Policy

Class-dependent products must be regenerated from the current classifier.

This includes at minimum:

- Figure 2.
- Figure 3.
- Figure 5.
- Figure 6.
- Figure S1.
- Figure S3.
- Figure S4.
- Class-dependent main tables.
- Complete class annotation atlas.
- Representative Transitional/Low-Stability gene outputs.
- GO enrichment.
- KEGG enrichment.
- Synthetic benchmark summaries.
- Noise-gradient results.
- Tau-sensitivity results.

Figure 4 SSI ECDF is numerically unaffected because SSI does not change.

Conceptual Figures 1 and 7 require terminology/definition review but not feature-level recomputation.

## 21. Supplementary Table Identities

For the current manuscript:

- Supplementary Table S2: Complete FSF gene annotation atlas.
- Supplementary Table S3: Representative FSF genes.

These identities are authoritative for current manuscript reconstruction.

## 22. Supplementary Figure S2

The existing Supplementary Figure S2 is not accepted as a current result.

Audit determined that its plotting script can select tau itself as both x and y because the script chooses the first matching threshold column and then the first numeric column.

It therefore requires a controlled rebuild with an explicitly selected, scientifically meaningful metric after corrected tau-sensitivity results are produced.

No automatic "first numeric column" metric selection is permitted in the replacement.

## 23. Governance

This lock governs:

- FSF current-package implementation.
- Current-package tests.
- Regeneration of current manuscript results.
- Current figures and tables.
- Publication-associated reproducibility claims.

It does not rewrite historical artifacts.

Any future change to:

- State rules.
- Tau semantics.
- Signal Identity.
- SSI definition.
- Stability boundaries.
- Ten-class structure.
- Condition-architecture definition.

requires an explicit new scientific-lock revision before implementation.

## 24. Authority Order

For current FSF package implementation, authority is:

1. This current scientific lock.
2. Current mathematical definitions that do not conflict with this lock.
3. Safe numerical reference outputs.
4. Historical implementation only for provenance and compatibility study.

Historical 0.60 class outputs are not current computational authority.
