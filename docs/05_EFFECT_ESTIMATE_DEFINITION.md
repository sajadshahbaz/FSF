## 4. Effect Estimate Definition

For each feature *i* and perturbation *p*, FSF requires a directional effect estimate:

[
E(i,p)
]

where:

* *i* denotes the feature.
* *p* denotes the perturbation instance.

### Definition

E(i,p) is the directional effect estimate obtained from the comparison of the target condition against its corresponding reference condition after applying perturbation *p*.

The effect estimate must be calculated using the same analytical procedure across all perturbations within a given analysis.

### Properties

[
E(i,p) > 0
]

indicates activation of the feature.

[
E(i,p) < 0
]

indicates suppression of the feature.

[
E(i,p) \approx 0
]

indicates no substantial directional effect.

### Default Implementation

The default implementation of FSF uses:

[
E(i,p)=\log FC(i,p)
]

where logFC represents the log fold change between the target and reference groups after perturbation.

### Framework Independence

FSF is intentionally agnostic to the specific effect estimator.

Alternative effect estimates may be used provided that:

1. The estimate is directional.
2. The estimate is computed consistently across all perturbations.
3. The same estimator is used throughout the analysis.

Examples include:

* log fold change (logFC)
* standardized effect size
* mean difference
* moderated effect estimates

### State Assignment

Given threshold:

[
\tau > 0
]

the feature state is assigned as:

UP:

[
E(i,p) > \tau
]

DOWN:

[
E(i,p) < -\tau
]

CONST:

[
-\tau \le E(i,p) \le \tau
]

### Default Threshold

The default threshold for Version 1.0 is:

[
\tau = 0.5
]

although user-defined thresholds may be specified.

### Design Principle

FSF does not generate effect estimates.

FSF consumes externally computed effect estimates and evaluates the reproducibility of their directional behavior across perturbations.

The methodological contribution of FSF lies in signal stratification through state frequencies, probability estimation, and stability-based classification rather than in the computation of the effect estimate itself.

