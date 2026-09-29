# Claims and Non-Claims

## Core Claims

FSF is designed to:

1. Stratify candidate signals according to reproducible behavior under perturbation.

2. Characterize stability patterns in heterogeneous transcriptomic datasets.

3. Distinguish reproducible signal classes from inconsistent signal behavior.

4. Provide interpretable signal categories for downstream biological analysis.

5. Complement existing statistical and machine-learning approaches.

## Secondary Claims

FSF may assist in:

* Identifying robust condition markers.
* Characterizing heterogeneous biological responses.
* Comparing stability behavior across datasets and conditions.
* Supporting downstream biological interpretation.

## Explicit Non-Claims

FSF does not:

* Replace DESeq2.
* Replace limma.
* Replace edgeR.
* Replace feature-selection methods.
* Replace machine-learning classifiers.
* Replace statistical hypothesis testing.

FSF is not:

* A differential expression algorithm.
* A feature-selection algorithm.
* A supervised predictive classifier.
* A clustering algorithm.
* A network-inference algorithm.

FSF does assign descriptive stability regions and signal classes from observed
state probabilities. Those assignments are signal stratification, not trained
prediction: they do not by themselves predict future observations or infer
causal, mechanistic, fitness, adaptive, or biological-importance labels.

## Publication Position

FSF should be presented as a signal stratification framework operating on candidate features obtained from external analytical methods.

The primary contribution is signal characterization rather than signal discovery.

