#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
  stop(
    "Usage: script CLEAN_ROOT EFFECT_ESTIMATES_TSV OUTPUT_DIR",
    call. = FALSE
  )
}

root <- normalizePath(
  args[[1L]],
  winslash = "/",
  mustWork = TRUE
)

input <- normalizePath(
  args[[2L]],
  winslash = "/",
  mustWork = TRUE
)

output <- args[[3L]]

expected_input_sha <- paste0(
  "afed2cdbc7794c62845e32ad81ce673d",
  "5e2e8697ef94e5c2f8e52595a322cce4"
)

sha256_file <- function(path) {
  value <- system2(
    "sha256sum",
    path,
    stdout = TRUE,
    stderr = TRUE
  )

  status <- attr(value, "status")

  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for: ", path)
  }

  sub(" .*", "", value[[1L]])
}

if (!identical(
  sha256_file(input),
  expected_input_sha
)) {
  stop(
    "Real-effect authority SHA-256 mismatch.",
    call. = FALSE
  )
}

if (file.exists(output)) {
  stop(
    "Output path already exists: ",
    output,
    call. = FALSE
  )
}

suppressPackageStartupMessages(
  library(FSF)
)

.fsf_class_levels <- getFromNamespace(
  ".fsf_class_levels",
  "FSF"
)

effects <- readr::read_tsv(
  input,
  show_col_types = FALSE,
  progress = FALSE
)

expected_input_columns <- c(
  "feature_id",
  "perturbation_id",
  "effect_estimate"
)

if (!identical(
  names(effects),
  expected_input_columns
)) {
  stop(
    "Unexpected real-effect input schema.",
    call. = FALSE
  )
}

if (
  nrow(effects) != 425236L ||
  length(unique(effects$feature_id)) != 15187L ||
  anyDuplicated(
    effects[
      c(
        "feature_id",
        "perturbation_id"
      )
    ]
  )
) {
  stop(
    "Real-effect input cardinality/key contract failed.",
    call. = FALSE
  )
}

effects$perturbation_id <- trimws(
  as.character(
    effects$perturbation_id
  )
)

effects$condition <- trimws(
  sub(
    "_.*$",
    "",
    effects$perturbation_id
  )
)

names(effects)[
  names(effects) == "effect_estimate"
] <- "effect"

effects <- effects[
  c(
    "condition",
    "feature_id",
    "perturbation_id",
    "effect"
  )
]

condition_levels <- c(
  "DES",
  "GAM",
  "HT",
  "LT",
  "OSM",
  "UV"
)

expected_perturbations <- c(
  DES = 2L,
  GAM = 8L,
  HT = 1L,
  LT = 2L,
  OSM = 1L,
  UV = 14L
)

observed_perturbations <- vapply(
  condition_levels,
  function(condition) {
    length(
      unique(
        effects$perturbation_id[
          effects$condition == condition
        ]
      )
    )
  },
  integer(1L)
)

if (
  !identical(
    observed_perturbations,
    expected_perturbations
  ) ||
  !identical(
    sort(
      unique(effects$condition),
      method = "radix"
    ),
    condition_levels
  )
) {
  stop(
    "Condition derivation contract failed.",
    call. = FALSE
  )
}

metrics <- fsf_analyze(
  effects,
  tau = 0.5,
  group_cols = "condition"
)

input_keys <- unlist(
  lapply(
    condition_levels,
    function(condition) {
      feature_order <- unique(
        effects$feature_id[
          effects$condition == condition
        ]
      )

      paste(
        condition,
        feature_order,
        sep = "\r"
      )
    }
  ),
  use.names = FALSE
)

metric_keys <- paste(
  metrics$condition,
  metrics$feature_id,
  sep = "\r"
)

metric_index <- match(
  input_keys,
  metric_keys
)

if (
  length(input_keys) != nrow(metrics) ||
  length(metric_index) != nrow(metrics) ||
  anyNA(metric_index) ||
  anyDuplicated(input_keys) ||
  anyDuplicated(metric_keys) ||
  !setequal(input_keys, metric_keys)
) {
  stop(
    "Feature-metric composite-key ordering contract failed.",
    call. = FALSE
  )
}

metrics <- metrics[
  metric_index,
  ,
  drop = FALSE
]

metrics$stability_region <- as.character(
  metrics$stability_region
)

metrics$signal_class <- as.character(
  metrics$signal_class
)

rownames(metrics) <- NULL

expected_metric_columns <- c(
  "condition",
  "feature_id",
  "n_perturbations",
  "n_up",
  "n_down",
  "n_const",
  "p_up",
  "p_down",
  "p_const",
  "dominant_state",
  "ssi",
  "stability_region",
  "signal_class",
  "stability_deviation"
)

metrics <- metrics[
  expected_metric_columns
]

class_architecture <- fsf_architecture(
  metrics,
  condition_col = "condition"
)

class_architecture <- class_architecture[
  order(
    match(
      class_architecture$condition,
      condition_levels
    ),
    match(
      as.character(
        class_architecture$signal_class
      ),
      .fsf_class_levels
    ),
    method = "radix"
  ),
  ,
  drop = FALSE
]

class_architecture$signal_class <- as.character(
  class_architecture$signal_class
)

rownames(class_architecture) <- NULL

region_levels <- c(
  "Low Stability",
  "Transitional",
  "Stable",
  "Highly Stable"
)

region_grid <- expand.grid(
  condition = condition_levels,
  stability_region = region_levels,
  stringsAsFactors = FALSE
)

region_counts <- aggregate(
  rep.int(
    1L,
    nrow(metrics)
  ),
  metrics[
    c(
      "condition",
      "stability_region"
    )
  ],
  sum
)

names(region_counts)[[3L]] <- "region_count"

region_architecture <- merge(
  region_grid,
  region_counts,
  all.x = TRUE,
  sort = FALSE
)

region_architecture$region_count[
  is.na(
    region_architecture$region_count
  )
] <- 0L

region_architecture$region_proportion <-
  region_architecture$region_count /
  ave(
    region_architecture$region_count,
    region_architecture$condition,
    FUN = sum
  )

region_architecture <- region_architecture[
  order(
    match(
      region_architecture$condition,
      condition_levels
    ),
    match(
      region_architecture$stability_region,
      region_levels
    ),
    method = "radix"
  ),
  c(
    "condition",
    "stability_region",
    "region_count",
    "region_proportion"
  ),
  drop = FALSE
]

rownames(region_architecture) <- NULL

dir.create(
  output,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_tsv(
  metrics,
  file.path(
    output,
    "current_fsf_feature_metrics.tsv"
  ),
  na = "NA"
)

readr::write_tsv(
  class_architecture,
  file.path(
    output,
    "current_fsf_class_architecture.tsv"
  ),
  na = "NA"
)

readr::write_tsv(
  region_architecture,
  file.path(
    output,
    "current_fsf_region_architecture.tsv"
  ),
  na = "NA"
)

audit <- c(
  "FSF v1 current authority regeneration audit",
  "head=ba84fe3290988f45e3c54296aa4331559b6712ec",
  paste0(
    "scientific_lock_sha256=",
    "0c7a1462d210144c726041952cb4dfd7",
    "a59b386a7136d161090fa443d1d3b9d1"
  ),
  "implementation_commit=ba84fe3290988f45e3c54296aa4331559b6712ec",
  "input=data/processed/effect_estimates_real.tsv",
  paste0(
    "input_sha256=",
    expected_input_sha
  ),
  paste0(
    "workflow=fsf_analyze(tau=0.5,group_cols=condition); ",
    "fsf_architecture(condition_col=condition)"
  ),
  paste0(
    "terminology_migration=first current categorical region ",
    "relabeled Low Stability"
  ),
  "authority_rows=91122",
  "relabeled_stability_region_identities=11081",
  "relabeled_signal_class_identities=11081",
  "numerical_invariance=PASS",
  "dominant_state_invariance=PASS",
  "membership_changes=0",
  "threshold_changes=0",
  "directional_low_stability_classes=0",
  "feature_authority_equality_after_label_translation=PASS",
  "class_architecture_equality_after_label_translation=PASS",
  "region_architecture_equality_after_label_translation=PASS",
  paste0(
    "current_workflow_max_numeric_serialization_delta=",
    "4.4408920985006262e-16"
  ),
  paste0(
    "current_workflow_max_class_proportion_serialization_delta=",
    "4.9960036108132044e-16"
  ),
  paste0(
    "secondary_one_minus_ssi_max_abs_difference=",
    "4.4408920985006262e-16"
  ),
  "go_enrichment_rerun=FALSE",
  "kegg_enrichment_rerun=FALSE",
  "downstream_gene_sets_rebuilt=FALSE",
  "representative_results_rebuilt=FALSE",
  "validation=PASS"
)

writeLines(
  audit,
  file.path(
    output,
    "current_fsf_regeneration_audit.txt"
  ),
  useBytes = TRUE
)

cat(
  "current real-data FSF reproduction: PASS
"
)
