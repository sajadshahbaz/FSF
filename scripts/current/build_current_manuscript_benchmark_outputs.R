#!/usr/bin/env Rscript

# Serialize manuscript benchmark data produced by the frozen current synthetic
# workflow. Historical scenario names are retained only in truth_scenario.

source(file.path("scripts", "current", "run_current_synthetic_validation.R"))

.normalize_baseline_effects <- function(path) {
  x <- readr::read_tsv(path, show_col_types = FALSE)
  names(x)[names(x) == "gene"] <- "feature_id"
  names(x)[names(x) == "perturbation"] <- "perturbation_id"
  names(x)[names(x) == "truth_class"] <- "truth_scenario"
  x
}

.benchmark_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path, na = "NA")
}

.region_counts <- function(x, setting = character()) {
  keys <- c(setting, "truth_scenario", "stability_region")
  out <- aggregate(rep.int(1L, nrow(x)), x[keys], sum)
  names(out)[ncol(out)] <- "n_features"
  denom <- interaction(out[c(setting, "truth_scenario")], drop = TRUE, lex.order = TRUE)
  out$proportion <- out$n_features / ave(out$n_features, denom, FUN = sum)
  out[do.call(order, c(out[keys], list(method = "radix"))), , drop = FALSE]
}

.tau_tidy <- function(tau_result) {
  s <- tau_result$summary
  metric_names <- c("mean_ssi", "median_ssi", "mean_stability_deviation",
                    "median_stability_deviation")
  numeric_rows <- do.call(rbind, lapply(metric_names, function(metric) {
    data.frame(tau = s$tau, truth_scenario = s$truth_scenario,
               metric = metric, category = NA_character_, value = s[[metric]],
               stringsAsFactors = FALSE)
  }))
  c <- tau_result$class_counts
  class_rows <- data.frame(
    tau = c$tau, truth_scenario = c$truth_scenario,
    metric = "signal_class_proportion", category = c$signal_class,
    value = c$proportion, stringsAsFactors = FALSE
  )
  out <- rbind(numeric_rows, class_rows)
  out <- out[order(out$tau, out$truth_scenario, out$metric, out$category,
                   method = "radix", na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

build_current_manuscript_benchmark_outputs <- function(effect_path, output_root) {
  if (!requireNamespace("readr", quietly = TRUE)) stop("readr is required.")
  effects <- .normalize_baseline_effects(effect_path)
  baseline <- run_current_baseline_synthetic(effects, validate_real = TRUE)
  noise <- run_current_noise_gradient()
  tau <- run_current_tau_sensitivity(effects)
  if (nrow(effects) != 500000L || nrow(baseline$feature_metrics) != 5000L ||
      nrow(noise$feature_metrics) != 17500L || nrow(tau$feature_metrics) != 20000L) {
    stop("Synthetic dimensions failed the frozen contract.", call. = FALSE)
  }
  migrated <- baseline$feature_metrics[
    baseline$feature_metrics$ssi > .50 & baseline$feature_metrics$ssi < .60, ]
  half <- baseline$feature_metrics[baseline$feature_metrics$feature_id == "IN_168", ]
  if (nrow(migrated) != 71L || any(migrated$signal_class != "Transitional Up") ||
      nrow(half) != 1L || half$ssi != .5 || half$dominant_state != "constant" ||
      half$signal_class != "Low Stability") {
    stop("Synthetic migration regression failed.", call. = FALSE)
  }
  out <- file.path(output_root, "benchmarks")
  .benchmark_write(baseline$feature_metrics, file.path(out, "baseline_feature_metrics.tsv"))
  .benchmark_write(baseline$summary, file.path(out, "baseline_summary.tsv"))
  .benchmark_write(baseline$class_counts, file.path(out, "baseline_truth_by_class.tsv"))
  .benchmark_write(.region_counts(baseline$feature_metrics),
                   file.path(out, "baseline_truth_by_region.tsv"))
  .benchmark_write(baseline$dominant_recovery,
                   file.path(out, "baseline_dominant_recovery.tsv"))
  .benchmark_write(noise$feature_metrics, file.path(out, "noise_feature_metrics.tsv"))
  .benchmark_write(noise$summary, file.path(out, "noise_summary.tsv"))
  .benchmark_write(noise$class_counts, file.path(out, "noise_truth_by_class.tsv"))
  .benchmark_write(.region_counts(noise$feature_metrics, "noise_sd"),
                   file.path(out, "noise_truth_by_region.tsv"))
  .benchmark_write(noise$dominant_recovery,
                   file.path(out, "noise_dominant_recovery.tsv"))
  .benchmark_write(tau$summary, file.path(out, "tau_summary.tsv"))
  .benchmark_write(tau$class_counts, file.path(out, "tau_truth_by_class.tsv"))
  .benchmark_write(.region_counts(tau$feature_metrics, "tau"),
                   file.path(out, "tau_truth_by_region.tsv"))
  .benchmark_write(tau$dominant_recovery, file.path(out, "tau_dominant_recovery.tsv"))
  .benchmark_write(.tau_tidy(tau), file.path(out, "tau_sensitivity_tidy.tsv"))
  invisible(list(baseline = baseline, noise = noise, tau = tau))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: script EFFECTS_TSV OUTPUT_ROOT")
  build_current_manuscript_benchmark_outputs(args[[1L]], args[[2L]])
}
