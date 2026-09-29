#!/usr/bin/env Rscript

# Transactional orchestrator for the deterministic current manuscript data
# bundle. It creates no figures and no biological interpretation products.

source(file.path("scripts", "current", "build_current_manuscript_enrichment_outputs.R"))
source(file.path("scripts", "current", "build_current_manuscript_benchmark_outputs.R"))
source(file.path("scripts", "current", "audit_current_manuscript_bundle.R"))

.ms_hash <- function(path) {
  out <- system2(Sys.which("sha256sum"), path, stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1L]])
}
.ms_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path, na = "NA")
}

.ms_write_canonical_numeric <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

.ms_architecture <- function(metrics, regions) {
  grid <- expand.grid(condition = sort(unique(metrics$condition), method = "radix"),
                      stability_region = regions, stringsAsFactors = FALSE)
  counts <- aggregate(rep.int(1L, nrow(metrics)), metrics[c("condition", "stability_region")], sum)
  names(counts)[3L] <- "region_count"
  out <- merge(grid, counts, all.x = TRUE, sort = FALSE)
  out$region_count[is.na(out$region_count)] <- 0L
  out$region_proportion <- out$region_count / ave(out$region_count, out$condition, FUN = sum)
  out <- out[order(out$condition, match(out$stability_region, regions), method = "radix"), ]
  rownames(out) <- NULL
  out
}

.ms_provenance <- function(root, annotation_path, effect_path, generated_utc) {
  inputs <- c(
    scientific_lock = "docs/FSF_V1_CURRENT_SCIENTIFIC_LOCK.md",
    annotation_authority = "docs/FSF_V1_ANNOTATION_AUTHORITY.md",
    feature_authority = "results/current_fsf_v1/current_fsf_feature_metrics.tsv",
    annotation_master = annotation_path,
    annotation_helper = "scripts/lib/current_fsf_annotation_authority.R",
    gene_set_builder = "scripts/current/build_current_enrichment_gene_sets.R",
    go_engine = "scripts/current/run_current_go_enrichment.R",
    kegg_engine = "scripts/current/run_current_kegg_enrichment.R",
    representative_workflow = "scripts/current/build_current_representative_features.R",
    synthetic_workflow = "scripts/current/run_current_synthetic_validation.R",
    synthetic_effect_input = effect_path
  )
  values <- vapply(inputs, .ms_hash, character(1L))
  keys <- paste0(names(values), "_sha256")
  fixed <- c(package_tag = "fsf-r-v0.1.0",
             package_commit = "ddb347a19e9877032a49aec2ebfeacc551b0b315",
             manuscript_branch = "fsf-manuscript-revision")
  scripts <- paste(c("build_current_manuscript_enrichment_outputs.R",
                     "build_current_manuscript_benchmark_outputs.R",
                     "build_current_manuscript_tables.R",
                     "audit_current_manuscript_bundle.R"), collapse = ";")
  data.frame(key = c(names(fixed), keys, "generated_utc", "r_version", "scripts",
                     "authority_scope"),
             value = c(fixed, values, generated_utc, R.version.string, scripts,
                       "frozen computational authorities plus manuscript presentation data"),
             stringsAsFactors = FALSE)
}

build_current_manuscript_bundle <- function(annotation_master_path, effect_path,
                                            output_root = file.path("results", "current_fsf_v1", "manuscript")) {
  if (!requireNamespace("readr", quietly = TRUE)) stop("readr is required.")
  parent <- dirname(output_root)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  transaction <- tempfile(".manuscript-txn-", tmpdir = parent)
  dir.create(transaction)
  on.exit(if (dir.exists(transaction)) unlink(transaction, recursive = TRUE), add = TRUE)
  for (d in c("tables", "enrichment/go", "enrichment/kegg", "benchmarks", "source_data", "audit")) {
    dir.create(file.path(transaction, d), recursive = TRUE, showWarnings = FALSE)
  }

  enrichment <- build_current_manuscript_enrichment_outputs(annotation_master_path, transaction)
  benchmark <- build_current_manuscript_benchmark_outputs(effect_path, transaction)
  metrics <- readr::read_tsv("results/current_fsf_v1/current_fsf_feature_metrics.tsv", show_col_types = FALSE)
  regions <- c("Low Stability", "Transitional", "Stable", "Highly Stable")
  region <- .ms_architecture(metrics, regions)
  signal <- aggregate(rep.int(1L, nrow(metrics)),
                      metrics[c("condition", "stability_region", "dominant_state", "signal_class")], sum)
  names(signal)[5L] <- "class_count"
  signal$class_proportion <- signal$class_count / ave(signal$class_count, signal$condition, FUN = sum)
  signal <- signal[order(signal$condition, match(signal$stability_region, regions),
                         signal$dominant_state, signal$signal_class, method = "radix"), ]
  manifest <- enrichment$go$results[0, 0, drop = FALSE]
  gene_sets <- build_current_enrichment_gene_sets(annotation_master_path)$current_gene_set_manifest
  gene_summary <- gene_sets[setdiff(names(gene_sets), c("all_features", "annotated_features"))]
  .ms_write(region, file.path(transaction, "tables/current_region_architecture.tsv"))
  .ms_write(signal, file.path(transaction, "tables/current_signal_class_architecture.tsv"))
  .ms_write(gene_summary, file.path(transaction, "tables/current_enrichment_gene_set_summary.tsv"))
  .ms_write(enrichment$go$results[enrichment$go$results$significant, ],
            file.path(transaction, "tables/current_go_statistical_summary.tsv"))
  .ms_write(enrichment$kegg$results[enrichment$kegg$results$significant, ],
            file.path(transaction, "tables/current_kegg_statistical_summary.tsv"))
  .ms_write(benchmark$baseline$summary,
            file.path(transaction, "tables/current_synthetic_baseline_summary.tsv"))
  .ms_write(readr::read_tsv(file.path(transaction, "benchmarks/tau_sensitivity_tidy.tsv"), show_col_types = FALSE),
            file.path(transaction, "tables/current_tau_sensitivity_summary.tsv"))

  simplex <- metrics[c("condition", "feature_id", "p_up", "p_down", "p_const", "ssi",
                       "dominant_state", "stability_region", "signal_class")]
  ssi <- metrics[c("condition", "feature_id", "ssi", "stability_deviation", "stability_region")]
  .ms_write_canonical_numeric(simplex, file.path(transaction, "source_data/figure2_simplex_source.tsv"))
  .ms_write(region, file.path(transaction, "source_data/figure3_architecture_source.tsv"))
  .ms_write_canonical_numeric(ssi, file.path(transaction, "source_data/figure4_ssi_ecdf_source.tsv"))
  .ms_write(enrichment$go$results[enrichment$go$results$significant, ],
            file.path(transaction, "source_data/figure5_go_source.tsv"))
  .ms_write(enrichment$kegg$results[enrichment$kegg$results$significant, ],
            file.path(transaction, "source_data/figure5_kegg_source.tsv"))
  .ms_write(benchmark$baseline$class_counts,
            file.path(transaction, "source_data/figure6_baseline_source.tsv"))
  .ms_write(benchmark$noise$dominant_recovery,
            file.path(transaction, "source_data/figure6_noise_source.tsv"))
  file.copy(file.path(transaction, "benchmarks/tau_sensitivity_tidy.tsv"),
            file.path(transaction, "source_data/figureS2_tau_source.tsv"))
  .ms_write(region, file.path(transaction, "source_data/figureS3_architecture_source.tsv"))
  .ms_write_canonical_numeric(ssi, file.path(transaction, "source_data/figureS4_ssi_distribution_source.tsv"))

  provenance <- .ms_provenance(transaction, annotation_master_path, effect_path,
                               format(Sys.time(), tz = "UTC", usetz = TRUE))
  .ms_write(provenance, file.path(transaction, "audit/provenance.tsv"))
  output_files <- list.files(transaction, recursive = TRUE, full.names = TRUE)
  output_files <- output_files[file.info(output_files)$isdir == FALSE]
  output_files <- output_files[!grepl("audit/(output_checksums|provenance)\\.tsv$", output_files)]
  rel <- substring(output_files, nchar(transaction) + 2L)
  checksums <- data.frame(path = rel, sha256 = vapply(output_files, .ms_hash, character(1L)),
                          stringsAsFactors = FALSE)
  checksums <- checksums[order(checksums$path, method = "radix"), ]
  .ms_write(checksums, file.path(transaction, "audit/output_checksums.tsv"))
  audit_current_manuscript_bundle(transaction)

  if (dir.exists(output_root)) stop("Current manuscript output root already exists; refusing partial replacement.")
  if (!file.rename(transaction, output_root)) stop("Atomic manuscript bundle publication failed.")
  invisible(output_root)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L || length(args) > 3L) {
    stop("Usage: script ANNOTATION_MASTER EFFECTS_TSV [OUTPUT_ROOT]")
  }
  output <- if (length(args) == 3L) args[[3L]] else file.path("results", "current_fsf_v1", "manuscript")
  build_current_manuscript_bundle(args[[1L]], args[[2L]], output)
}
