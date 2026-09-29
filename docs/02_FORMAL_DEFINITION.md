# FSF Formal Definition (Version 1.0)

## 1. Authority and purpose

This active formal definition is subordinate to
`docs/FSF_V1_CURRENT_SCIENTIFIC_LOCK.md`, the primary scientific authority for
FSF v1.

The Feature Stability Framework (FSF) is a signal stratification framework
that characterizes candidate biological signals according to their observed
behavior across repeated dataset perturbations. FSF evaluates signal behavior,
not statistical significance, and operates downstream of effect estimation or
candidate-feature selection.

## 2. Inputs

FSF requires perturbation-derived signed effects for identified features. The
same analytical procedure must produce the effect estimate for every
perturbation. The effect may be a log fold change or another consistently
defined signed quantity; effect estimation itself is outside the package core.

## 3. Perturbation definition

A perturbation is a controlled modification of sample composition while
preserving the biological condition being studied. Examples include bootstrap
resampling, random sample removal, leave-one-sample-out,
leave-one-subgroup-out, and leave-one-timepoint-out. Artificial noise
injection is not part of the core framework.

## 4. Signal-state assignment

For feature \(i\) and perturbation \(p\), let \(E_{i,p}\) be the signed effect
and let \(\tau > 0\) be the state threshold. The default is \(\tau = 0.5\), and
the threshold is user-configurable.

- Up: \(E_{i,p} > \tau\)
- Down: \(E_{i,p} < -\tau\)
- Constant: \(-\tau \le E_{i,p} \le \tau\)

Exact \(+\tau\) and \(-\tau\) values are therefore Constant.

## 5. State probabilities

For feature \(i\), with \(N_i\) observed perturbations:

\[
P_{up,i} = \frac{N_{up,i}}{N_i}, \qquad
P_{down,i} = \frac{N_{down,i}}{N_i}, \qquad
P_{const,i} = \frac{N_{const,i}}{N_i}.
\]

The state probabilities satisfy:

\[
P_{up,i} + P_{down,i} + P_{const,i} = 1.
\]

## 6. Signal Stratification Index

For each feature:

\[
SSI_i = \max(P_{up,i}, P_{down,i}, P_{const,i}).
\]

Because three nonnegative state probabilities sum to one, the theoretical
range is:

\[
\frac{1}{3} \le SSI_i \le 1.
\]

SSI describes the concentration of observed perturbation states. It is not a
supervised-prediction confidence score and is not an instability score.

The dominant state is the state attaining SSI. Equal maxima produce a tied
dominant state. No rounding occurs before dominant-state, region, or class
assignment.

## 7. Stability-region assignment

The four current stability regions are assigned from raw, unrounded SSI alone:

- \(SSI \le 0.50\): Low Stability
- \(0.50 < SSI < 0.75\): Transitional
- \(0.75 \le SSI < 0.90\): Stable
- \(SSI \ge 0.90\): Highly Stable

No tolerance changes a region boundary.

## 8. Signal-class assignment

Low Stability is one non-directional signal class. For SSI strictly greater
than 0.50, the stability region combines with the unique dominant state (Up,
Constant, or Down) to produce a directional class.

The ten signal classes are:

1. Low Stability
2. Transitional Up
3. Transitional Constant
4. Transitional Down
5. Stable Up
6. Stable Constant
7. Stable Down
8. Highly Stable Up
9. Highly Stable Constant
10. Highly Stable Down

Thus, `stability_region` has four categories and `signal_class` has ten
categories. A dominant state may still be recorded for a Low-Stability
feature, but it does not create a directional Low-Stability class.

## 9. Stability Deviation

Stability Deviation is defined as:

\[
Stability\ Deviation_i = 1 - SSI_i.
\]

It is an auxiliary descriptive transform only. It does not determine
`stability_region`, `signal_class`, or an independent stability definition.

## 10. Framework outputs

For every feature identity, the current FSF package returns state counts and
probabilities, the dominant state, SSI, stability region, signal class, and
Stability Deviation. The package may also summarize signal-class architecture
within user-defined conditions.

## 11. Scope and interpretation

FSF assigns descriptive signal strata and classes from observed state
probabilities. It is not a supervised predictive classifier and does not, from
an FSF class alone, infer future outcomes, statistical significance, causal or
mechanistic relationships, fitness, adaptation, or biological importance.

FSF does not replace differential-expression analysis, feature selection,
machine-learning prediction, or statistical hypothesis testing. Biological
interpretation requires external evidence and study context.
