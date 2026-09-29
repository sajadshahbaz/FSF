#!/usr/bin/env Rscript

generator <- file.path(
  "scripts",
  "current",
  "build_current_biological_review_package.R"
)
source(generator)

root <- file.path("results", "current_fsf_v1", "manuscript", "review")
stopifnot(audit_current_biological_review_package(root, review_state = "final"))

required <- c(
  "current_go_evidence.tsv",
  "current_kegg_evidence.tsv",
  "go_theme_review.tsv",
  "kegg_theme_review.tsv",
  "combined_profile_review.tsv",
  "uv_specific_review.tsv",
  "condition_architecture_review.tsv",
  "representative_feature_review.tsv",
  "review_summary.md"
)

stopifnot(length(required) == 9L)
stopifnot(all(file.exists(file.path(root, required))))
historical_inventory <- file.path(
  root,
  "historical_interpretation_inventory.tsv"
)
stopifnot(file.exists(historical_inventory))
stopifnot(file.info(historical_inventory)$size > 0)
historical_inventory_sha256 <- sub(
  " .*",
  "",
  system2("sha256sum", historical_inventory, stdout = TRUE)[1]
)
stopifnot(identical(
  historical_inventory_sha256,
  "95a4c7f4144b191eabb2854031b8364910b272c8d4b95ddf8af853fa84016cea"
))

go_themes <- readr::read_tsv(
  file.path(root, "go_theme_review.tsv"),
  show_col_types = FALSE
)
kegg_themes <- readr::read_tsv(
  file.path(root, "kegg_theme_review.tsv"),
  show_col_types = FALSE
)
go_evidence <- readr::read_tsv(
  file.path(root, "current_go_evidence.tsv"),
  show_col_types = FALSE
)
kegg_evidence <- readr::read_tsv(
  file.path(root, "current_kegg_evidence.tsv"),
  show_col_types = FALSE
)
profile <- readr::read_tsv(
  file.path(root, "combined_profile_review.tsv"),
  show_col_types = FALSE
)
uv <- readr::read_tsv(
  file.path(root, "uv_specific_review.tsv"),
  show_col_types = FALSE
)
architecture <- readr::read_tsv(
  file.path(root, "condition_architecture_review.tsv"),
  show_col_types = FALSE
)
reps <- readr::read_tsv(
  file.path(root, "representative_feature_review.tsv"),
  show_col_types = FALSE
)

theme_schema <- c(
  "theme_rule_id",
  "rule_definition_source",
  "theme_label",
  "theme_target"
)

stopifnot(all(theme_schema %in% names(go_themes)))
stopifnot(all(theme_schema %in% names(kegg_themes)))
stopifnot(!anyDuplicated(go_themes$theme_rule_id))
stopifnot(!anyDuplicated(kegg_themes$theme_rule_id))

stopifnot(nrow(profile) == 40L)
stopifnot("gene_set_id" %in% names(profile))
stopifnot(!anyDuplicated(profile$gene_set_id))
stopifnot(length(unique(go_evidence$gene_set_id)) == 40L)
stopifnot(length(unique(kegg_evidence$gene_set_id)) == 40L)
stopifnot(nrow(uv) == 13L)
stopifnot(nrow(architecture) == 6L)
stopifnot("condition" %in% names(architecture))
stopifnot(!anyDuplicated(architecture$condition))

current_tables <- list(
  go_themes,
  kegg_themes,
  profile,
  architecture,
  reps
)
stopifnot(all(vapply(
  current_tables,
  function(x) !any(startsWith(names(x), "historical_")),
  logical(1L)
)))

forbidden_generated_schema <- c(
  "historical_architecture_label",
  "historical_profile_status",
  "historical_selection_reason"
)
stopifnot(all(vapply(
  current_tables,
  function(x) !any(forbidden_generated_schema %in% names(x)),
  logical(1L)
)))

obsolete_profile_fields <- c(
  "GO_evidence_change",
  "KEGG_evidence_change",
  "theme_evidence_change"
)
stopifnot(!any(obsolete_profile_fields %in% names(profile)))

for (x in list(profile, reps)) {
  stopifnot(all(c("author_decision", "author_notes") %in% names(x)))
  stopifnot(all(is.na(x$author_decision)))
  stopifnot(all(is.na(x$author_notes)))
}
stopifnot(all(c("author_decision", "author_notes") %in% names(architecture)))
stopifnot(nrow(architecture) == 6L)
stopifnot(setequal(
  architecture$condition,
  c("DES", "GAM", "HT", "LT", "OSM", "UV")
))
stopifnot(!anyDuplicated(architecture$condition))
stopifnot(all(
  !is.na(architecture$author_decision) &
  architecture$author_decision == "RETAIN_DESCRIPTIVE_ONLY"
))
stopifnot(all(
  !is.na(architecture$author_notes) &
  nzchar(trimws(architecture$author_notes))
))

generator_text <- paste(readLines(generator, warn = FALSE), collapse = "\n")

forbidden_generator_text <- c(
  "Table_SignalClass_Biology_CURATED.tsv",
  "historical_interpretation_inventory.tsv",
  "legacy-bundle:scripts/pipeline/",
  "historical_architecture_label",
  "historical_profile_status",
  "historical_selection_reason"
)
stopifnot(!any(vapply(
  forbidden_generator_text,
  function(x) grepl(x, generator_text, fixed = TRUE),
  logical(1L)
)))

required_generator_text <- c(
  ".current_theme_rules <- function()",
  "CURRENT_MATCHING_EVIDENCE",
  "NO_CURRENT_EVIDENCE"
)
stopifnot(all(vapply(
  required_generator_text,
  function(x) grepl(x, generator_text, fixed = TRUE),
  logical(1L)
)))

cat("current biological review package tests: PASS\n")
