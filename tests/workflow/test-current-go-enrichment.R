source(file.path("scripts", "current", "run_current_go_enrichment.R"))

annotation <- data.frame(
  feature_id = paste0("f", 1:10),
  eggnog_go = c(
    "prefix GO:0000001 suffix,malformed", "GO:0000001", "GO:0000001",
    "GO:0000002", "GO:0000002", "GO:0000002", "GO:0000003",
    "GO:0000003", "-", NA
  ),
  interpro_go = c(
    "GO:0000001;bad", "GO:0000002", "GO:0000002", NA, "GO:0000002",
    "bad;GO:0000003 trailing", NA, "GO:0000003", "-", ""
  ),
  stringsAsFactors = FALSE
)

make_manifest_row <- function(id, features) {
  row <- data.frame(
    gene_set_id = id, gene_set_family = "main_class", condition = "TEST",
    stability_region = "Stable", dominant_state = "up",
    signal_class = "Stable Up", all_feature_count = length(features),
    annotated_feature_count = length(features), stringsAsFactors = FALSE
  )
  row$all_features <- I(list(sort(unique(features), method = "radix")))
  row$annotated_features <- I(list(sort(unique(features), method = "radix")))
  row
}

manifest <- rbind(
  make_manifest_row("set_four", c("f1", "f2", "f3", "f4")),
  make_manifest_row("set_five_a", c("f1", "f2", "f3", "f4", "f5", "outside")),
  make_manifest_row("set_five_b", c("f1", "f2", "f6", "f7", "f8"))
)
gene_sets <- list(
  current_gene_set_manifest = manifest,
  universes = list(annotation_universe = sort(paste0("f", 1:10), method = "radix"))
)

out <- .run_current_go_enrichment(gene_sets, annotation, validate_real = FALSE)
stopifnot(nrow(out$source_specific_mapping) == 14L)
stopifnot(nrow(out$combined_mapping) == 11L)
stopifnot(!anyDuplicated(out$source_specific_mapping[c("feature_id", "go_term", "source")]))
stopifnot(!anyDuplicated(out$combined_mapping[c("feature_id", "go_term")]))
stopifnot(!any(grepl("malformed|bad", out$combined_mapping$go_term)))
stopifnot(all(grepl("^GO:[0-9]{7}$", out$combined_mapping$go_term)))

status_four <- out$gene_set_status[out$gene_set_status$gene_set_id == "set_four", ]
stopifnot(status_four$query_size == 4L)
stopifnot(status_four$status == "SKIPPED_QUERY_SIZE_LT_5")
stopifnot(!any(out$results$gene_set_id == "set_four"))

a <- out$results[out$results$gene_set_id == "set_five_a", , drop = FALSE]
b <- out$results[out$results$gene_set_id == "set_five_b", , drop = FALSE]
stopifnot(unique(a$query_size) == 5L, unique(a$background_size) == 10L)
term1 <- a[a$go_term == "GO:0000001", , drop = FALSE]
expected_p <- stats::phyper(3L - 1L, 3L, 10L - 3L, 5L, lower.tail = FALSE)
stopifnot(identical(term1$p_value, expected_p))
stopifnot(identical(term1$query_fraction, 3 / 5))
stopifnot(identical(term1$background_fraction, 3 / 10))
stopifnot(identical(term1$enrichment_ratio, 2))
stopifnot(identical(term1$contributing_feature_ids, "f1;f2;f3"))
stopifnot(identical(a$adjusted_p_value, stats::p.adjust(a$p_value, method = "BH")))
stopifnot(identical(b$adjusted_p_value, stats::p.adjust(b$p_value, method = "BH")))

# BH is scoped independently: duplicating another set cannot change set A.
single <- .run_current_go_enrichment(
  list(current_gene_set_manifest = manifest[2, , drop = FALSE], universes = gene_sets$universes),
  annotation,
  validate_real = FALSE
)
stopifnot(identical(a$adjusted_p_value, single$results$adjusted_p_value))
stopifnot(identical(a$significant, a$adjusted_p_value <= 0.05))

ordered <- with(out$results, order(
  gene_set_id, adjusted_p_value, p_value, -enrichment_ratio, go_term,
  method = "radix"
))
stopifnot(identical(ordered, seq_len(nrow(out$results))))
shuffled <- annotation[rev(seq_len(nrow(annotation))), , drop = FALSE]
out_shuffled <- .run_current_go_enrichment(gene_sets, shuffled, validate_real = FALSE)
stopifnot(identical(out$source_specific_mapping, out_shuffled$source_specific_mapping))
stopifnot(identical(out$combined_mapping, out_shuffled$combined_mapping))
stopifnot(identical(out$results, out_shuffled$results))

text <- readLines(file.path("scripts", "current", "run_current_go_enrichment.R"))
stopifnot(sum(grepl("build_current_enrichment_gene_sets\\(annotation_master_path\\)", text)) == 1L)
stopifnot(any(grepl("universes\\$annotation_universe", text)))
stopifnot(any(grepl("stats::phyper", text)))
stopifnot(any(grepl("lower.tail = FALSE", text, fixed = TRUE)))
stopifnot(any(grepl('method = "BH"', text, fixed = TRUE)))
stopifnot(!any(grepl(
  "FSF_condition_class_annotation_master|enrichment_gene_set_manifest|results/real_data",
  text
)))
stopifnot(!any(grepl("stability_level|fsf_signal_class|weakly_stable|instable", text)))
stopifnot(!any(grepl("write_tsv|write\\.table|saveRDS|sink\\(|file\\.create", text)))
stopifnot(!any(grepl("eggnog_kegg|kegg_term|run_current_kegg", text, ignore.case = TRUE)))
stopifnot(!any(grepl("0\\.60|0\\.6[^0-9]", text)))

cat("current GO enrichment workflow tests: PASS\n")
