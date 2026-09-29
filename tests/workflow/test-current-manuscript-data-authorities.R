#!/usr/bin/env Rscript

source(file.path("scripts", "current", "build_current_manuscript_enrichment_outputs.R"))
source(file.path("scripts", "current", "build_current_manuscript_benchmark_outputs.R"))
source(file.path("scripts", "current", "build_current_manuscript_tables.R"))

# Ranking is deterministic and uses adjusted p, raw p, decreasing enrichment,
# then term identifier. Fewer than N significant rows are retained in full.
x <- data.frame(gene_set_id = rep("x", 4), significant = c(TRUE, TRUE, TRUE, FALSE),
                adjusted_p_value = c(.01, .01, .01, .001), p_value = c(.02, .02, .01, .001),
                enrichment_ratio = c(2, 3, 1, 99), go_term = c("GO:3", "GO:2", "GO:1", "GO:0"))
r <- .manuscript_rank(x, "go_term", 2L)
stopifnot(identical(r$go_term, c("GO:1", "GO:2")))

metrics <- data.frame(condition = rep("X", 4),
                      stability_region = c("Highly Stable", "Low Stability", "Stable", "Highly Stable"))
a <- .ms_architecture(metrics, c("Low Stability", "Transitional", "Stable", "Highly Stable"))
stopifnot(identical(a$stability_region, c("Low Stability", "Transitional", "Stable", "Highly Stable")))
stopifnot(identical(a$region_count, c(1L, 0L, 1L, 2L)))

tau <- list(
  summary = data.frame(tau = .25, truth_scenario = "stable_up", n_features = 1,
                       mean_ssi = .9, median_ssi = .9,
                       mean_stability_deviation = .1, median_stability_deviation = .1),
  class_counts = data.frame(tau = .25, truth_scenario = "stable_up",
                            signal_class = "Highly Stable Up", n_features = 1, proportion = 1)
)
tidy <- .tau_tidy(tau)
stopifnot(identical(names(tidy), c("tau", "truth_scenario", "metric", "category", "value")))
stopifnot(!any(grepl("Instability|weakly_stable_transitional", tidy$metric)))

files <- c("scripts/current/build_current_manuscript_enrichment_outputs.R",
           "scripts/current/build_current_manuscript_benchmark_outputs.R",
           "scripts/current/build_current_manuscript_tables.R",
           "scripts/current/audit_current_manuscript_bundle.R")
text <- unlist(lapply(files, readLines, warn = FALSE))
stopifnot(!any(grepl("SSI *(>=|>|<) *0\\.60", text)))
stopifnot(!any(grepl("results/manuscript|results/real_data", text)))
cat("current manuscript data-authority workflow tests: PASS\n")
