#!/usr/bin/env Rscript

# Serialize deterministic manuscript-facing enrichment products from the
# frozen current in-memory engines. No biological themes are assigned here.

source(file.path("scripts", "current", "run_current_go_enrichment.R"))
source(file.path("scripts", "current", "run_current_kegg_enrichment.R"))

.manuscript_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path, na = "NA")
}

.manuscript_rank <- function(x, term, n = 50L) {
  x <- x[x$significant, , drop = FALSE]
  if (!nrow(x)) return(x)
  group <- interaction(x[c("gene_set_id")], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(seq_len(nrow(x)), group), function(i) {
    y <- x[i, , drop = FALSE]
    y <- y[order(y$adjusted_p_value, y$p_value, -y$enrichment_ratio,
                 y[[term]], method = "radix"), , drop = FALSE]
    y[seq_len(min(n, nrow(y))), , drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$gene_set_id, out$adjusted_p_value, out$p_value,
            -out$enrichment_ratio, out[[term]], method = "radix"), , drop = FALSE]
}

.validate_enrichment_authority <- function(go, kegg) {
  go_observed <- c(
    background = length(go$background),
    mapped = length(unique(go$combined_mapping$feature_id)),
    terms = length(unique(go$combined_mapping$go_term)),
    combined_pairs = nrow(go$combined_mapping),
    source_pairs = nrow(go$source_specific_mapping),
    rows = nrow(go$results), significant = sum(go$results$significant),
    significant_sets = length(unique(go$results$gene_set_id[go$results$significant]))
  )
  go_expected <- c(13467L, 9321L, 20502L, 1131086L, 1159053L,
                   637452L, 11119L, 47L)
  if (!identical(as.integer(go_observed), go_expected)) {
    stop("Frozen GO aggregate contract failed.", call. = FALSE)
  }
  kr <- kegg$results
  ko <- kr[kr$enrichment_type == "KO", , drop = FALSE]
  pathway <- kr[kr$enrichment_type == "PATHWAY", , drop = FALSE]
  kegg_observed <- c(
    ko_universe = length(kegg$backgrounds$KO),
    pathway_universe = length(kegg$backgrounds$PATHWAY),
    ko_terms = length(unique(kegg$KO_mapping$kegg_term)),
    pathway_terms = length(unique(kegg$PATHWAY_mapping$kegg_term)),
    ko_pairs = nrow(kegg$KO_mapping), pathway_pairs = nrow(kegg$PATHWAY_mapping),
    ko_rows = nrow(ko), pathway_rows = nrow(pathway), combined_rows = nrow(kr),
    significant_ko = sum(ko$significant),
    significant_pathway = sum(pathway$significant)
  )
  kegg_expected <- c(7178L, 4680L, 5755L, 786L, 10146L, 44072L,
                     89789L, 36704L, 126493L, 1120L, 734L)
  if (!identical(as.integer(kegg_observed), kegg_expected)) {
    stop("Frozen KEGG aggregate contract failed.", call. = FALSE)
  }
  list(go = go_observed, kegg = kegg_observed)
}

build_current_manuscript_enrichment_outputs <- function(
    annotation_master_path, output_root, top_n = 50L) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr is required.", call. = FALSE)
  }
  go <- run_current_go_enrichment(annotation_master_path)
  kegg <- run_current_kegg_enrichment(annotation_master_path)
  aggregates <- .validate_enrichment_authority(go, kegg)

  go_dir <- file.path(output_root, "enrichment", "go")
  kegg_dir <- file.path(output_root, "enrichment", "kegg")
  .manuscript_write_tsv(go$results, file.path(go_dir, "go_enrichment_all.tsv"))
  .manuscript_write_tsv(go$results[go$results$significant, , drop = FALSE],
                        file.path(go_dir, "go_enrichment_significant.tsv"))
  .manuscript_write_tsv(go$gene_set_status, file.path(go_dir, "go_gene_set_status.tsv"))
  .manuscript_write_tsv(.manuscript_rank(go$results, "go_term", top_n),
                        file.path(go_dir, paste0("go_significant_top", top_n, ".tsv")))
  go_summary <- data.frame(
    metric = names(aggregates$go), value = as.numeric(aggregates$go),
    stringsAsFactors = FALSE
  )
  .manuscript_write_tsv(go_summary, file.path(go_dir, "go_mapping_summary.tsv"))

  .manuscript_write_tsv(kegg$results, file.path(kegg_dir, "kegg_enrichment_all.tsv"))
  .manuscript_write_tsv(kegg$results[kegg$results$significant, , drop = FALSE],
                        file.path(kegg_dir, "kegg_enrichment_significant.tsv"))
  .manuscript_write_tsv(kegg$gene_set_status,
                        file.path(kegg_dir, "kegg_gene_set_status.tsv"))
  .manuscript_write_tsv(.manuscript_rank(kegg$results, "kegg_term", top_n),
                        file.path(kegg_dir, paste0("kegg_significant_top", top_n, ".tsv")))
  kegg_summary <- data.frame(
    metric = names(aggregates$kegg), value = as.numeric(aggregates$kegg),
    stringsAsFactors = FALSE
  )
  .manuscript_write_tsv(kegg_summary, file.path(kegg_dir, "kegg_mapping_summary.tsv"))
  invisible(list(go = go, kegg = kegg, aggregates = aggregates))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: script ANNOTATION_MASTER OUTPUT_ROOT")
  build_current_manuscript_enrichment_outputs(args[[1L]], args[[2L]])
}
