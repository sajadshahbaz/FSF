source(file.path("scripts", "current", "run_current_synthetic_validation.R"))

fixture <- data.frame(
  feature_id = rep(c("SU_1", "SD_1", "SC_1", "WT_1", "IN_1"), each = 4L),
  perturbation_id = rep(paste0("P", 1:4), 5L),
  effect_estimate = c(
    rep(1.5, 4), rep(-1.5, 4), rep(0, 4),
    c(1, 1, 1, 0), c(1, -1, 0, 0)
  ),
  stringsAsFactors = FALSE
)

baseline_a <- run_current_baseline_synthetic(fixture)
baseline_b <- run_current_baseline_synthetic(fixture[sample(seq_len(nrow(fixture))), ])
stopifnot(nrow(baseline_a$feature_metrics) == 5L)
stopifnot(identical(baseline_a$feature_metrics, baseline_b$feature_metrics))
stopifnot(identical(sort(unique(baseline_a$feature_metrics$truth_scenario)),
                    sort(.fsf_synthetic_scenarios)))
stopifnot(!any(as.character(baseline_a$feature_metrics$stability_region) %in%
                  c("instable", "weakly_stable", "weakly_stable_transitional")))
stopifnot(!any(as.character(baseline_a$feature_metrics$signal_class) %in%
                  .fsf_synthetic_scenarios))

boundary <- fsf_classify(data.frame(
  feature_id = paste0("b", 1:7),
  p_up = c(1/3, .49, .50, .50 + 1e-12, .75, .90, 1),
  p_down = c(1/3, .30, .30, .30 - 1e-12, .15, .05, 0),
  p_const = c(1/3, .21, .20, .20, .10, .05, 0)
))
stopifnot(identical(
  as.character(boundary$stability_region),
  c("Low Stability", "Low Stability", "Low Stability", "Transitional",
    "Stable", "Highly Stable", "Highly Stable")
))

low <- fsf_classify(data.frame(
  feature_id = c("up", "down", "constant", "tied"),
  p_up = c(.50, .30, .30, .50), p_down = c(.30, .50, .20, .50),
  p_const = c(.20, .20, .50, 0)
))
stopifnot(identical(as.character(low$signal_class), rep("Low Stability", 4L)))

stopifnot(identical(.fsf_synthetic_noise_grid,
                    c(.05, .10, .20, .35, .50, .75, 1.00)))
stopifnot(identical(.fsf_synthetic_tau_grid, c(.25, .50, .75, 1.00)))

noise_a <- run_current_noise_gradient(n_per_scenario = 2L, n_perturbations = 4L)
noise_b <- run_current_noise_gradient(n_per_scenario = 2L, n_perturbations = 4L)
stopifnot(identical(noise_a, noise_b))
stopifnot(identical(sort(unique(noise_a$feature_metrics$noise_sd)),
                    .fsf_synthetic_noise_grid))

tau <- run_current_tau_sensitivity(fixture)
stopifnot(identical(sort(unique(tau$feature_metrics$tau)), .fsf_synthetic_tau_grid))
stopifnot(identical(
  tau$feature_metrics,
  tau$feature_metrics[order(tau$feature_metrics$tau,
                            tau$feature_metrics$feature_id,
                            method = "radix"), , drop = FALSE]
))

workflow_text <- readLines(file.path(
  "scripts", "current", "run_current_synthetic_validation.R"
))
stopifnot(any(grepl("fsf_analyze\\(", workflow_text)))
stopifnot(!any(grepl("write_tsv|write\\.table|saveRDS|ggsave|sink\\(|file\\.create",
                    workflow_text)))
boundary_060_lines <- workflow_text[grepl(
  "ssi.*0\\.60|0\\.60.*ssi", workflow_text, ignore.case = TRUE
)]
stopifnot(length(boundary_060_lines) == 1L)
stopifnot(grepl(
  "feature_metrics\\$ssi > 0\\.50.*feature_metrics\\$ssi < 0\\.60",
  boundary_060_lines
))
stopifnot(!any(grepl("stability_level|fsf_signal_class", workflow_text)))
stopifnot(!any(grepl(
  paste(
    "(predicted_)?stability_region.*(instable|weakly_stable)",
    "(predicted_)?signal_class.*(instable|weakly_stable)",
    "weakly_stable_transitional", sep = "|"
  ),
  workflow_text
)))
old_truth_lines <- workflow_text[grepl("instable|weakly_stable", workflow_text)]
stopifnot(all(grepl(
  "synthetic_scenarios|stable_up.*stable_down.*stable_constant|grepl\\(|make\\(|truth_scenario|feature_id = paste0",
  old_truth_lines
)))

real_effect_path <- file.path("data", "processed", "effect_estimates.tsv")
stopifnot(file.exists(real_effect_path))
real_effects <- read.delim(
  real_effect_path, stringsAsFactors = FALSE, check.names = FALSE
)
real_baseline <- run_current_baseline_synthetic(real_effects, validate_real = TRUE)
stopifnot(nrow(real_baseline) == 5000L)

cat("current synthetic validation workflow tests: PASS\n")
