# FSF v0.1.0 Final Reproducibility Report

Generated: 2026-09-29T09:18:52+00:00

## Release conclusion

The FSF v0.1.0 clean publication tree satisfies the release gates established
for the manuscript and reusable R package.

**Release status: READY FOR NEW PUBLIC GITHUB REPOSITORY.**

The scientific implementation is frozen at this point. No further scientific
or computational changes are required before public release.

## R package

Package: FSF  
Version: 0.1.0

Public API:

- `fsf_analyze`
- `fsf_architecture`
- `fsf_assign_states`
- `fsf_classify`
- `fsf_signal_identity`

The package source was successfully built, installed into an isolated library,
loaded with `library(FSF)`, exercised through its exported API, and passed its
package tests, examples, documentation checks, and release checks.

## Scientific contract

The final scientific audit confirmed:

- default signal-state threshold: tau = 0.5
- Up: effect > tau
- Down: effect < -tau
- Constant: -tau <= effect <= tau
- SSI = max(p_up, p_down, p_const)
- theoretical SSI range: [1/3, 1]
- Low Stability: SSI <= 0.50
- Transitional: 0.50 < SSI < 0.75
- Stable: 0.75 <= SSI < 0.90
- Highly Stable: SSI >= 0.90
- ten current FSF signal classes
- Stability Deviation = 1 - SSI

The exact SSI=0.50 boundary remains Low Stability.

## Real-data integrity

The governed real-data authority contains 91,122 condition-feature records,
corresponding to 15,187 features across six conditions.

Observed SSI range:

- minimum: 0.357142857142857
- maximum: 1.0
- out-of-range records: 0

The independently recalculated condition architecture matched the frozen
manuscript authority.

HT and OSM each contain a single contrast per feature. Their SSI=1 values are
therefore structural consequences of single-contrast analysis and must not be
interpreted as evidence of cross-perturbation reproducibility.

## Synthetic validation

Synthetic generator scenario labels are provenance/truth labels only.
Predicted states, probabilities, SSI values, stability regions, and signal
classes are computed by the installed FSF package.

Validated components include:

- baseline synthetic benchmark
- noise-gradient validation
- tau sensitivity at 0.25, 0.50, 0.75, and 1.00
- exact boundary behavior
- tie handling
- single-perturbation behavior
- unrounded probability/SSI behavior

No substantive circular-validation mechanism was identified.

## Representative features and enrichment

Representative-feature selection is deterministic and based on FSF-derived
stability, dominant direction, SSI, annotation coverage where applicable,
Stability Deviation, and feature identifier ordering.

GO enrichment uses:

- unique feature-term mappings
- the governed annotation universe
- upper-tail hypergeometric testing
- BH multiple-testing adjustment
- adjusted P <= 0.05

KEGG KO and PATHWAY analyses use separate backgrounds and separate
multiple-testing families.

Frozen expected counts used in workflow code are validation guards; they do
not generate the corresponding scientific results.

## Frozen noncomputational figures

Figure 1 and Figure 7 are intentionally frozen publication assets and are not
script-generated.

SHA-256:

- Figure1.png: `9e3796b99bc0fd1b36d8c5029ee286ebfe76c9bf018d6d763dab596fd83be0f2`
- Figure1.pdf: `cc295cd90aeb44efc2fe6a3a3143ed79040c3ffe410147b152a246b30b5d1f6c`
- Figure7.png: `6ab8be9958a8073303661814fa7def0de9e1a5d9ece745adb23dc1bf377771a0`
- Figure7.pdf: `c5d7af446b5f13e66d7437f0d0197f2400d33c9c3ba28ca7fc502cb2fc184468`

Their identities match the accepted publication authorities.

## Governed biological-review authorities

Two final biological interpretation artifacts are explicitly classified as
frozen governed author-review authorities rather than computationally
regenerated outputs:

- deferred_profile_term_evidence.tsv  
  SHA-256: `e03c94ef535209fa76a587694c536381eb73a69472d693ef61b9c2ae48b613be`

- final_biological_interpretation_summary.md  
  SHA-256: `eeeeaa1180bac85e2de5417b43248c37c84280199e6d44a830204dc3979ed95e`

This boundary is deliberate because no truthful computational producer exists
for those final author-review decisions.

## Manuscript tables

The final Stage12 manuscript table package contains exactly 15 expected files.

Text outputs reproduce the frozen authorities exactly where byte identity is
appropriate. Accepted exceptions are limited to previously documented
serialization behavior such as gzip wrapper metadata and machine-precision
floating-point textual serialization where scientific identity was separately
demonstrated.

The Figure 2, Figure 4, and Figure S4 source-table producer now uses canonical
numeric serialization that reproduces the frozen byte identities directly.

## Portability

The release tree has been checked for:

- machine-specific runtime paths
- dependencies on the historical source repository
- active archive dependencies
- escaping symlinks
- Python cache files
- temporary files
- obsolete runtime paths

All release-sanitation gates passed.

Temporary working directories are allocated portably through `${TMPDIR:-/tmp}`
rather than a fixed machine-specific temporary path.

Public reproduction runner SHA-256:

`63cac139b59f26fecf9db36031705cab7313f9c8f183257119bede8b94ed3947`

## Reproduction evidence

The release decision combines:

1. a complete disposable 0-to-100 reproduction from the clean publication
   boundary;
2. successful package build/install/load and public-API validation;
3. complete workflow-test execution;
4. frozen-reference comparisons;
5. targeted CLEAN integration tests following scientifically neutral
   portability corrections;
6. deterministic canonical serialization verification on the real data;
7. final release-sanitation checks; and
8. a final independent scientific-integrity audit.

A second redundant full scientific recomputation was intentionally not
performed after the final non-scientific portability integrations because each
affected component was validated directly and the complete computational
workflow had already passed.

## Historical provenance note

Historical development commit/tag fields retained inside frozen manuscript
provenance are preserved as historical provenance. They should not be
interpreted as the Git commit identifier of the new clean publication
repository.

## Final status

R package release ready: **PASS**  
Full manuscript reproduction: **PASS**  
Scientific integrity audit: **PASS**  
Release sanitation: **PASS**  
Ready for new GitHub repository: **YES**
