# Current FSF v1 synthetic validation workflow.
#
# Historical truth labels are generator-provenance scenario identifiers only.
# All predicted fields are produced by the current package implementation.
# The workflow returns in-memory objects and intentionally has no writers or
# figure logic. The historical Supplementary Figure S2 issue is manuscript-side.

suppressPackageStartupMessages(library(FSF))

.fsf_synthetic_noise_grid <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00)
.fsf_synthetic_tau_grid <- c(0.25, 0.50, 0.75, 1.00)
.fsf_synthetic_scenarios <- c(
  "stable_up", "stable_down", "stable_constant", "weakly_stable", "instable"
)

.fsf_synthetic_truth_scenario <- function(feature_id) {
  feature_id <- as.character(feature_id)
  scenario <- ifelse(
    grepl("^(stable_up_)?SU_", feature_id), "stable_up",
    ifelse(
      grepl("^(stable_down_)?SD_", feature_id), "stable_down",
      ifelse(
        grepl("^(stable_constant_)?SC_", feature_id), "stable_constant",
        ifelse(
          grepl("^(weakly_stable_)?WT_", feature_id), "weakly_stable",
          ifelse(grepl("^(instable_)?IN_", feature_id), "instable", NA_character_)
        )
      )
    )
  )
  if (anyNA(scenario)) {
    stop("Synthetic feature IDs contain an unknown truth scenario.", call. = FALSE)
  }
  scenario
}

.fsf_prepare_synthetic_effects <- function(data) {
  if (!is.data.frame(data) ||
      !all(c("feature_id", "perturbation_id") %in% names(data))) {
    stop("Synthetic effects require feature_id and perturbation_id.", call. = FALSE)
  }
  effect_name <- if ("effect" %in% names(data)) {
    "effect"
  } else if ("effect_estimate" %in% names(data)) {
    "effect_estimate"
  } else {
    stop("Synthetic effects require effect or effect_estimate.", call. = FALSE)
  }
  out <- data
  out$feature_id <- as.character(out$feature_id)
  out$perturbation_id <- as.character(out$perturbation_id)
  out$effect <- as.numeric(out[[effect_name]])
  if (anyNA(out$effect) || any(!is.finite(out$effect))) {
    stop("Synthetic effects must be finite and non-missing.", call. = FALSE)
  }
  if (!"truth_scenario" %in% names(out)) {
    out$truth_scenario <- .fsf_synthetic_truth_scenario(out$feature_id)
  }
  out$truth_scenario <- as.character(out$truth_scenario)
  if (anyNA(out$truth_scenario) ||
      any(!out$truth_scenario %in% .fsf_synthetic_scenarios)) {
    stop("truth_scenario contains an unknown generator scenario.", call. = FALSE)
  }
  out
}

.fsf_synthetic_summaries <- function(feature_metrics, setting_cols) {
  group_key <- interaction(
    feature_metrics[c(setting_cols, "truth_scenario")],
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(seq_len(nrow(feature_metrics)), group_key)
  summary_rows <- lapply(groups, function(index) {
    x <- feature_metrics[index, , drop = FALSE]
    row <- x[1L, setting_cols, drop = FALSE]
    row$truth_scenario <- x$truth_scenario[[1L]]
    row$n_features <- nrow(x)
    row$mean_ssi <- mean(x$ssi)
    row$median_ssi <- stats::median(x$ssi)
    row$mean_stability_deviation <- mean(x$stability_deviation)
    row$median_stability_deviation <- stats::median(x$stability_deviation)
    row
  })
  summary <- do.call(rbind, summary_rows)

  count_key <- interaction(
    feature_metrics[c(setting_cols, "truth_scenario", "signal_class")],
    drop = TRUE, lex.order = TRUE
  )
  count_groups <- split(seq_len(nrow(feature_metrics)), count_key)
  count_rows <- lapply(count_groups, function(index) {
    x <- feature_metrics[index, , drop = FALSE]
    row <- x[1L, setting_cols, drop = FALSE]
    row$truth_scenario <- x$truth_scenario[[1L]]
    row$signal_class <- as.character(x$signal_class[[1L]])
    row$n_features <- nrow(x)
    row
  })
  class_counts <- do.call(rbind, count_rows)
  denominator_key <- interaction(
    class_counts[c(setting_cols, "truth_scenario")],
    drop = TRUE, lex.order = TRUE
  )
  totals <- ave(class_counts$n_features, denominator_key, FUN = sum)
  class_counts$proportion <- class_counts$n_features / totals

  order_cols <- c(setting_cols, "truth_scenario")
  summary <- summary[do.call(order, c(summary[order_cols], list(method = "radix"))), , drop = FALSE]
  class_counts <- class_counts[do.call(order, c(
    class_counts[order_cols], list(-class_counts$n_features),
    list(class_counts$signal_class), list(method = "radix")
  )), , drop = FALSE]
  rownames(summary) <- NULL
  rownames(class_counts) <- NULL

  dominant_key <- interaction(
    class_counts[c(setting_cols, "truth_scenario")],
    drop = TRUE, lex.order = TRUE
  )
  dominant <- class_counts[!duplicated(dominant_key), , drop = FALSE]
  names(dominant)[names(dominant) == "signal_class"] <- "dominant_recovered_signal_class"
  names(dominant)[names(dominant) == "proportion"] <- "recovery_proportion"
  list(summary = summary, class_counts = class_counts, dominant_recovery = dominant)
}

.fsf_run_synthetic_effects <- function(data, tau, setting_cols = character()) {
  prepared <- .fsf_prepare_synthetic_effects(data)
  truth_cols <- unique(c(setting_cols, "feature_id", "truth_scenario"))
  truth <- unique(prepared[truth_cols])
  analysis <- fsf_analyze(
    prepared[c(setting_cols, "feature_id", "perturbation_id", "effect")],
    tau = tau,
    group_cols = setting_cols
  )
  by <- c(setting_cols, "feature_id")
  key <- do.call(paste, c(analysis[by], sep = "\r"))
  truth_key <- do.call(paste, c(truth[by], sep = "\r"))
  truth_index <- match(key, truth_key)
  if (anyNA(truth_index)) stop("Synthetic truth join lost feature identities.", call. = FALSE)
  analysis$truth_scenario <- truth$truth_scenario[truth_index]
  analysis <- analysis[do.call(order, c(analysis[by], list(method = "radix"))), , drop = FALSE]
  rownames(analysis) <- NULL
  analysis
}

run_current_baseline_synthetic <- function(effect_data, validate_real = FALSE) {
  feature_metrics <- .fsf_run_synthetic_effects(effect_data, tau = 0.5)
  summaries <- .fsf_synthetic_summaries(feature_metrics, character())
  if (validate_real) {
    migrated <- feature_metrics[feature_metrics$ssi > 0.50 & feature_metrics$ssi < 0.60, ]
    expected_distribution <- c(`0.53` = 2L, `0.54` = 5L, `0.55` = 5L,
                               `0.56` = 6L, `0.57` = 11L, `0.58` = 19L,
                               `0.59` = 23L)
    observed_distribution <- table(sprintf("%.2f", migrated$ssi))
    exact_half <- feature_metrics[feature_metrics$feature_id == "IN_168", ]
    valid <- nrow(feature_metrics) == 5000L && nrow(migrated) == 71L &&
      all(migrated$dominant_state == "up") &&
      all(as.character(migrated$stability_region) == "Transitional") &&
      all(as.character(migrated$signal_class) == "Transitional Up") &&
      identical(as.integer(observed_distribution[names(expected_distribution)]),
                as.integer(expected_distribution)) &&
      nrow(exact_half) == 1L && exact_half$ssi == 0.50 &&
      exact_half$p_up == 0.26 && exact_half$p_down == 0.24 &&
      exact_half$p_const == 0.50 && exact_half$dominant_state == "constant" &&
      as.character(exact_half$stability_region) == "Low Stability" &&
      as.character(exact_half$signal_class) == "Low Stability"
    if (!valid) stop("Real baseline synthetic regression failed.", call. = FALSE)
  }
  c(list(mode = "BASELINE_SYNTHETIC", feature_metrics = feature_metrics), summaries)
}

.fsf_generate_noise_effects <- function(
    noise_grid = .fsf_synthetic_noise_grid, seed = 12345L,
    n_per_scenario = 500L, n_perturbations = 100L) {
  set.seed(seed)
  output <- vector("list", length(noise_grid))
  for (i in seq_along(noise_grid)) {
    noise <- noise_grid[[i]]
    grid <- expand.grid(
      feature_index = seq_len(n_per_scenario),
      perturbation = seq_len(n_perturbations)
    )
    make <- function(prefix, scenario, mean) {
      data.frame(
        noise_sd = noise, truth_scenario = scenario,
        feature_id = paste0(scenario, "_", prefix, "_", grid$feature_index),
        perturbation_id = paste0("P", grid$perturbation),
        effect = stats::rnorm(nrow(grid), mean = mean, sd = noise),
        stringsAsFactors = FALSE
      )
    }
    stable_up <- make("SU", "stable_up", 1.5)
    stable_down <- make("SD", "stable_down", -1.5)
    stable_constant <- make("SC", "stable_constant", 0)
    weak <- make("WT", "weakly_stable", 0.8)
    hidden <- sample(c("up", "down", "const"), nrow(grid), replace = TRUE)
    up <- stats::rnorm(nrow(grid), 1.5, noise)
    down <- stats::rnorm(nrow(grid), -1.5, noise)
    constant <- stats::rnorm(nrow(grid), 0, noise)
    unstable_effect <- ifelse(hidden == "up", up, ifelse(hidden == "down", down, constant))
    unstable <- data.frame(
      noise_sd = noise, truth_scenario = "instable",
      feature_id = paste0("instable_IN_", grid$feature_index),
      perturbation_id = paste0("P", grid$perturbation), effect = unstable_effect,
      stringsAsFactors = FALSE
    )
    output[[i]] <- rbind(stable_up, stable_down, stable_constant, weak, unstable)
  }
  do.call(rbind, output)
}

run_current_noise_gradient <- function(
    noise_grid = .fsf_synthetic_noise_grid, seed = 12345L,
    n_per_scenario = 500L, n_perturbations = 100L) {
  if (!identical(as.numeric(noise_grid), .fsf_synthetic_noise_grid)) {
    stop("Current noise grid must match the locked historical grid.", call. = FALSE)
  }
  effects <- .fsf_generate_noise_effects(
    noise_grid, seed, n_per_scenario, n_perturbations
  )
  feature_metrics <- .fsf_run_synthetic_effects(
    effects, tau = 0.5, setting_cols = "noise_sd"
  )
  summaries <- .fsf_synthetic_summaries(feature_metrics, "noise_sd")
  c(list(
    mode = "NOISE_GRADIENT", config = list(
      noise_grid = noise_grid, seed = seed, tau = 0.5,
      n_per_scenario = n_per_scenario, n_perturbations = n_perturbations
    ), feature_metrics = feature_metrics
  ), summaries)
}

run_current_tau_sensitivity <- function(
    effect_data, tau_grid = .fsf_synthetic_tau_grid) {
  if (!identical(as.numeric(tau_grid), .fsf_synthetic_tau_grid)) {
    stop("Current tau grid must match the locked historical grid.", call. = FALSE)
  }
  feature_metrics <- do.call(rbind, lapply(tau_grid, function(tau) {
    out <- .fsf_run_synthetic_effects(effect_data, tau = tau)
    out$tau <- tau
    out
  }))
  feature_metrics <- feature_metrics[order(
    feature_metrics$tau, feature_metrics$feature_id, method = "radix"
  ), , drop = FALSE]
  rownames(feature_metrics) <- NULL
  summaries <- .fsf_synthetic_summaries(feature_metrics, "tau")
  c(list(
    mode = "TAU_SENSITIVITY",
    config = list(tau_grid = tau_grid),
    feature_metrics = feature_metrics
  ), summaries)
}

run_current_synthetic_validation <- function(mode, effect_data = NULL) {
  mode <- match.arg(mode, c(
    "BASELINE_SYNTHETIC", "NOISE_GRADIENT", "TAU_SENSITIVITY"
  ))
  if (mode == "BASELINE_SYNTHETIC") {
    return(run_current_baseline_synthetic(effect_data))
  }
  if (mode == "NOISE_GRADIENT") {
    return(run_current_noise_gradient())
  }
  run_current_tau_sensitivity(effect_data)
}
