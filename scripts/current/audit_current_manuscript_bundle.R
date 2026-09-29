#!/usr/bin/env Rscript

.bundle_script_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(.bundle_script_root, "scripts/current/build_current_external_artifact_manifests.R"))

.bundle_sha256 <- function(path) {
  out <- system2(Sys.which("sha256sum"), path, stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1L]])
}

.bundle_required_files <- c(
  "tables/current_region_architecture.tsv",
  "tables/current_signal_class_architecture.tsv",
  "tables/current_enrichment_gene_set_summary.tsv",
  "enrichment/go/go_enrichment_all.tsv",
  "enrichment/go/go_enrichment_significant.tsv",
  "enrichment/go/go_gene_set_status.tsv",
  "enrichment/go/go_mapping_summary.tsv",
  "enrichment/go/go_significant_top50.tsv",
  "enrichment/kegg/kegg_enrichment_all.tsv",
  "enrichment/kegg/kegg_enrichment_significant.tsv",
  "enrichment/kegg/kegg_gene_set_status.tsv",
  "enrichment/kegg/kegg_mapping_summary.tsv",
  "enrichment/kegg/kegg_significant_top50.tsv",
  "benchmarks/baseline_feature_metrics.tsv",
  "benchmarks/noise_feature_metrics.tsv",
  "benchmarks/tau_sensitivity_tidy.tsv",
  "source_data/figure2_simplex_source.tsv",
  "source_data/figure3_architecture_source.tsv",
  "source_data/figure4_ssi_ecdf_source.tsv",
  "source_data/figure5_go_source.tsv",
  "source_data/figure5_kegg_source.tsv",
  "source_data/figure6_baseline_source.tsv",
  "source_data/figure6_noise_source.tsv",
  "source_data/figureS2_tau_source.tsv",
  "source_data/figureS3_architecture_source.tsv",
  "source_data/figureS4_ssi_distribution_source.tsv",
  "audit/provenance.tsv",
  "audit/output_checksums.tsv"
)

audit_current_manuscript_bundle <- function(root) {
  external_manifest_path <- file.path(root, "audit", "external_artifacts.tsv")
  external <- NULL
  if (file.exists(external_manifest_path)) {
    external <- readr::read_tsv(external_manifest_path, show_col_types = FALSE)
    required_external <- c(
      "relative_repository_path", "external_artifact_path", "bytes", "sha256",
      "artifact_role", "generator_script", "package_tag", "package_commit",
      "scientific_lock_sha256"
    )
    if (!identical(names(external), required_external) ||
        anyDuplicated(external$relative_repository_path)) {
      stop("External artifact manifest is invalid.", call. = FALSE)
    }
    external_root <- .fsf_external_root()
    external$resolved_path <- vapply(
      external$external_artifact_path, .fsf_resolve_external_locator,
      character(1L), root = external_root
    )
    observed_bytes <- file.info(external$resolved_path)$size
    observed_hashes <- vapply(external$resolved_path, .bundle_sha256, character(1L))
    if (!identical(as.numeric(observed_bytes), as.numeric(external$bytes)) ||
        !identical(unname(observed_hashes), external$sha256)) {
      stop("External artifact size or SHA-256 verification failed.", call. = FALSE)
    }
  }
  resolve <- function(relative) {
    local <- file.path(root, relative)
    if (file.exists(local)) return(local)
    if (!is.null(external)) {
      index <- match(relative, external$relative_repository_path)
      if (!is.na(index)) return(external$resolved_path[[index]])
    }
    NA_character_
  }
  paths <- vapply(.bundle_required_files, resolve, character(1L))
  if (anyNA(paths)) {
    stop(paste("Missing bundle files:", paste(.bundle_required_files[is.na(paths)], collapse = ", ")))
  }
  region <- readr::read_tsv(file.path(root, "tables/current_region_architecture.tsv"),
                            show_col_types = FALSE)
  expected <- data.frame(
    condition = rep(c("DES", "GAM", "HT", "LT", "OSM", "UV"), each = 4L),
    stability_region = rep(c("Low Stability", "Transitional", "Stable", "Highly Stable"), 6L),
    region_count = c(5837,0,0,9350, 1085,2106,6127,5869, 0,0,0,15187,
                     1211,0,0,13976, 0,0,0,15187, 2948,6697,3151,2391),
    stringsAsFactors = FALSE
  )
  key <- paste(region$condition, region$stability_region)
  index <- match(paste(expected$condition, expected$stability_region), key)
  if (anyNA(index) || !identical(as.integer(region$region_count[index]), as.integer(expected$region_count)) ||
      any(aggregate(region_count ~ condition, region, sum)$region_count != 15187L)) {
    stop("Region architecture authority mismatch.", call. = FALSE)
  }
  go <- readr::read_tsv(resolve("enrichment/go/go_enrichment_all.tsv"), show_col_types = FALSE)
  kegg <- readr::read_tsv(resolve("enrichment/kegg/kegg_enrichment_all.tsv"), show_col_types = FALSE)
  if (nrow(go) != 637452L || sum(go$significant) != 11119L ||
      nrow(kegg) != 126493L || sum(kegg$significant[kegg$enrichment_type == "KO"]) != 1120L ||
      sum(kegg$significant[kegg$enrichment_type == "PATHWAY"]) != 734L) {
    stop("Serialized enrichment aggregate mismatch.", call. = FALSE)
  }
  tau <- readr::read_tsv(file.path(root, "benchmarks/tau_sensitivity_tidy.tsv"), show_col_types = FALSE)
  if (!identical(names(tau), c("tau", "truth_scenario", "metric", "category", "value")) ||
      !identical(sort(unique(tau$tau)), c(.25, .50, .75, 1.00))) {
    stop("Tau sensitivity tidy contract failed.", call. = FALSE)
  }
  baseline <- readr::read_tsv(file.path(root, "benchmarks/baseline_feature_metrics.tsv"), show_col_types = FALSE)
  migrated <- baseline[baseline$ssi > .5 & baseline$ssi < .6, ]
  half <- baseline[baseline$feature_id == "IN_168", ]
  if (nrow(migrated) != 71L || any(migrated$signal_class != "Transitional Up") ||
      nrow(half) != 1L || half$stability_region != "Low Stability") {
    stop("Baseline synthetic regression failed.", call. = FALSE)
  }
  current_fields <- c("stability_region", "signal_class", "dominant_recovered_signal_class")
  tsv <- list.files(root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  stale <- character()
  for (path in tsv) {
    x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
    fields <- intersect(names(x), current_fields)
    if (length(fields) && any(grepl("Instability|instable|weakly_stable|weakly stable",
                                    unlist(x[fields]), ignore.case = TRUE))) stale <- c(stale, path)
  }
  if (length(stale)) stop("Obsolete terminology in current predicted fields.", call. = FALSE)
  provenance <- readr::read_tsv(file.path(root, "audit/provenance.tsv"), show_col_types = FALSE)
  required_keys <- c("package_tag", "package_commit", "manuscript_branch",
                     "scientific_lock_sha256", "feature_authority_sha256",
                     "annotation_authority_sha256", "annotation_master_sha256", "annotation_helper_sha256",
                     "gene_set_builder_sha256", "go_engine_sha256", "kegg_engine_sha256",
                     "representative_workflow_sha256", "synthetic_workflow_sha256",
                     "generated_utc", "r_version", "scripts")
  if (!all(required_keys %in% provenance$key) || anyDuplicated(provenance$key)) {
    stop("Provenance manifest is incomplete.", call. = FALSE)
  }
  invisible(list(files = length(tsv), go_rows = nrow(go), kegg_rows = nrow(kegg)))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: script BUNDLE_ROOT")
  audit_current_manuscript_bundle(args[[1L]])
  cat("current manuscript bundle audit: PASS\n")
}
