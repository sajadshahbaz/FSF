# Feature Signal Stratification (FSF)

**FSF** is a lightweight R package for quantifying how consistently a feature exhibits an **Up**, **Down**, or **Constant** response across perturbations.

FSF converts signed effect estimates into discrete signal states, calculates state probabilities, derives the **Signal Stratification Index (SSI)**, assigns stability and signal classes, and summarizes condition-level signal architecture.

The package is dataset-independent and can be applied to any analysis in which signed feature-level effects are available across one or more perturbations.

---

## Installation

### Install the release from GitHub

The current release is **FSF 0.1.0**.

```r
install.packages("remotes")

remotes::install_github(
  "sajadshahbaz/FSF",
  ref = "fsf-r-v0.1.0"
)
```

Load the package:

```r
library(FSF)
```

### Install from source

```sh
git clone https://github.com/sajadshahbaz/FSF.git
cd FSF
git checkout fsf-r-v0.1.0

R CMD build .
R CMD INSTALL FSF_0.1.0.tar.gz
```

FSF is not currently distributed through CRAN.

---

## Quick start

FSF requires signed effect estimates for features observed across perturbations.

A minimal input data frame contains:

| Column | Description |
|---|---|
| `feature_id` | Nonmissing, nonblank feature identifier |
| `perturbation_id` | Nonmissing, nonblank perturbation identifier |
| `effect` | Finite numeric signed effect estimate |

For example:

```r
library(FSF)

effects <- data.frame(
  feature_id = c(
    "gene1", "gene1", "gene1",
    "gene2", "gene2", "gene2"
  ),
  perturbation_id = c(
    "p1", "p2", "p3",
    "p1", "p2", "p3"
  ),
  effect = c(
    1.2, 0.9, 0.1,
    -1.1, -0.8, 0.0
  )
)

result <- fsf_analyze(
  effects,
  tau = 0.5
)

result
```

`fsf_analyze()` returns one classified row per feature identity.

---

## Analysis with conditions or other strata

If features should be analyzed independently within conditions, tissues, treatments, or other strata, specify the corresponding existing column or columns using `group_cols`.

```r
effects <- data.frame(
  condition = rep(c("control", "treated"), each = 6),
  feature_id = rep(rep(c("f1", "f2"), each = 3), 2),
  perturbation_id = rep(paste0("p", 1:3), 4),
  effect = c(
    1, 1, 0,
    0, 0, 0,
    1, 1, 1,
    -1, -1, 0
  )
)

classified <- fsf_analyze(
  effects,
  tau = 0.5,
  group_cols = "condition"
)

classified
```

Condition-level signal architecture can then be summarized with:

```r
architecture <- fsf_architecture(
  classified,
  condition_col = "condition"
)

architecture
```

`group_cols` is a function argument naming columns already present in the input data. It is not an additional required input column.

---

# How FSF works

## 1. Signal-state assignment

FSF assigns each signed effect to one of three states using the threshold `tau`.

The default is:

```text
tau = 0.5
```

The state rules are:

```text
effect > tau            -> Up
effect < -tau           -> Down
-tau <= effect <= tau   -> Constant
```

Therefore values exactly equal to `+tau` or `-tau` are classified as **Constant**.

With the default threshold:

```text
effect > 0.5            -> Up
effect < -0.5           -> Down
-0.5 <= effect <= 0.5  -> Constant
```

---

## 2. Signal Identity

For each feature identity, FSF counts the number of observed states:

```text
n_up
n_down
n_const
```

where:

```text
n_perturbations = n_up + n_down + n_const
```

The corresponding state probabilities are:

```text
p_up    = n_up / n_perturbations
p_down  = n_down / n_perturbations
p_const = n_const / n_perturbations
```

and:

```text
p_up + p_down + p_const = 1
```

---

## 3. Signal Stratification Index

The **Signal Stratification Index (SSI)** is:

```text
SSI = max(p_up, p_down, p_const)
```

SSI quantifies how strongly the observed perturbation responses concentrate in one signal state.

Because the three state probabilities sum to one:

```text
1/3 <= SSI <= 1
```

The state attaining SSI is reported as the `dominant_state`.

Possible values are:

```text
up
down
constant
tied
```

If two or more states share the maximum probability:

```text
dominant_state = "tied"
```

No rounding is applied before SSI classification.

---

## 4. Stability regions

FSF uses four SSI stability regions:

| SSI | Stability region |
|---|---|
| `SSI <= 0.50` | Low Stability |
| `0.50 < SSI < 0.75` | Transitional |
| `0.75 <= SSI < 0.90` | Stable |
| `SSI >= 0.90` | Highly Stable |

The exact boundary:

```text
SSI = 0.50
```

belongs to **Low Stability**.

Directional signal classes are assigned only when:

```text
SSI > 0.50
```

---

## 5. FSF signal classes

FSF defines ten signal classes.

### Low Stability

```text
Low Stability
```

Low Stability is non-directional.

### Transitional

```text
Transitional Up
Transitional Constant
Transitional Down
```

### Stable

```text
Stable Up
Stable Constant
Stable Down
```

### Highly Stable

```text
Highly Stable Up
Highly Stable Constant
Highly Stable Down
```

A Low-Stability feature can still have a recorded dominant state, but FSF does not create directional Low-Stability subclasses such as:

```text
Low Stability Up
Low Stability Down
Low Stability Constant
```

---

## Stability Deviation

FSF also reports:

```text
Stability Deviation = 1 - SSI
```

This is mathematically complementary to SSI and is provided as an auxiliary descriptive quantity.

It is not an independent stability endpoint.

---

# Output from `fsf_analyze()`

A standard analysis returns fields including:

```text
feature_id
n_perturbations
n_up
n_down
n_const
p_up
p_down
p_const
dominant_state
ssi
stability_region
signal_class
stability_deviation
```

When `group_cols` is supplied, the grouping columns are retained in the result.

Example:

```r
result <- fsf_analyze(
  effects,
  tau = 0.5,
  group_cols = "condition"
)

head(result)
```

---

# Public API

FSF exports five functions.

## `fsf_assign_states()`

Assign Up, Down, and Constant states from signed effects.

```r
states <- fsf_assign_states(
  effects,
  tau = 0.5
)
```

---

## `fsf_signal_identity()`

Calculate state counts and state probabilities for each feature identity.

```r
identity <- fsf_signal_identity(
  effects,
  tau = 0.5,
  group_cols = NULL
)
```

---

## `fsf_classify()`

Classify an already aggregated Signal Identity table.

```r
classified <- fsf_classify(identity)
```

---

## `fsf_analyze()`

Calculate Signal Identity and perform FSF classification in one step.

```r
classified <- fsf_analyze(
  effects,
  tau = 0.5,
  group_cols = NULL
)
```

---

## `fsf_architecture()`

Calculate condition-level composition of the ten FSF signal classes.

```r
architecture <- fsf_architecture(
  classified,
  condition_col = "condition"
)
```

The architecture output includes zero-count classes so that conditions share the same ten-class structure.

Full function documentation is available from R:

```r
?fsf_assign_states
?fsf_signal_identity
?fsf_classify
?fsf_analyze
?fsf_architecture
```

---

# Single-perturbation interpretation

FSF supports feature identities represented by a single perturbation.

With only one perturbation, exactly one of the three state probabilities must equal 1. Therefore:

```text
SSI = 1
```

is mathematically inevitable.

An SSI of 1 based on a single perturbation should therefore **not** be interpreted as evidence of reproducibility or consistency across multiple perturbations.

Users should consider the number of contributing perturbations when interpreting SSI.

---

# Input requirements

The minimum required columns are:

```text
feature_id
perturbation_id
effect
```

Requirements:

- `feature_id` must be nonmissing and nonblank.
- `perturbation_id` must be nonmissing and nonblank.
- `effect` must contain finite numeric signed values.
- duplicate `feature_id + perturbation_id` combinations are not permitted within the same analysis stratum.
- independent strata should be represented using one or more grouping columns and supplied through `group_cols`.

For example:

```r
effects <- data.frame(
  tissue = rep(c("brain", "gut"), each = 6),
  feature_id = rep(c("gene1", "gene2"), each = 3, times = 2),
  perturbation_id = rep(c("p1", "p2", "p3"), 4),
  effect = c(
    1.1, 0.8, 0.2,
    -0.9, -1.2, 0.1,
    0.0, 0.2, 0.1,
    1.0, 1.1, 0.9
  )
)

result <- fsf_analyze(
  effects,
  tau = 0.5,
  group_cols = "tissue"
)
```

---

# What FSF does not do

The FSF package assumes that signed effect estimates have already been obtained from an appropriate upstream analysis.

FSF itself does not:

- estimate effect sizes;
- perform differential-expression testing;
- normalize raw count data;
- select or rank biological features;
- perform functional enrichment;
- perform biological interpretation;
- infer causal biological mechanisms.

This separation allows FSF to operate independently of the method used to generate the original signed effects.

---

# Reproducibility

This repository contains the computational workflows and governed resources associated with development and validation of the FSF framework.

The release-level reproduction workflow can be run from the repository root with:

```sh
./run_fsf_publication_reproduction.sh
```

The runner validates the required release inputs, builds and installs the FSF package in an isolated environment, and executes the associated computational validation workflows.

Detailed provenance and reproducibility records are available under:

```text
docs/
reproducibility/
```

These resources are not required for ordinary use of the installed FSF package.

The current FSF scientific and computational authority is documented in `docs/CURRENT_FSF_AUTHORITY.md`, with current release results under `results/current_fsf_v1/`. Historical development material is documented in `docs/FSF_LEGACY_PROVENANCE.md` and is retained for provenance only; it is not current scientific or runtime authority.

---

# Package versus repository

The installed FSF package contains the reusable, dataset-independent FSF computational core.

The GitHub repository additionally contains development, testing, validation, documentation, and reproducibility resources.

Users who only want to apply FSF to their own data need only install the R package and use the public API described above.

---

# Version

Current software release:

```text
FSF 0.1.0
```

Release tag:

```text
fsf-r-v0.1.0
```

Repository:

```text
https://github.com/sajadshahbaz/FSF
```

---

# Citation

A formal publication citation or DOI will be added when authoritative publication metadata becomes available.

For the installed R package, citation information can be obtained with:

```r
citation("FSF")
```

---

# Author

**Sajad Shahbazi**

Department of Animal Physiology and Development<br>
Faculty of Biology<br>
Adam Mickiewicz University<br>
Poznań, Poland

---

# License

See the `LICENSE` file for the software license governing this repository.
