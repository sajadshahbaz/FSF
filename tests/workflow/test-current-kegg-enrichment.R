source(file.path("scripts", "current", "run_current_kegg_enrichment.R"))

annotation <- data.frame(
  feature_id = paste0("f", 1:10),
  eggnog_kegg_ko = c(
    "ko:K00001, K00001", "ko:K00001;K00002", "K00001", "K00002",
    "K00002", "K00002", "-", NA, "", "K00003"
  ),
  eggnog_kegg_pathway = c(
    "map00010, map00010", "map00010; map00020", "map00010", "map00020",
    "map00020", "map00020", "-", NA, "", "map00030"
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
  universes = list(
    ko_universe = sort(paste0("f", 1:10), method = "radix"),
    pathway_universe = sort(paste0("f", 1:8), method = "radix")
  )
)

out <- .run_current_kegg_enrichment(gene_sets, annotation, validate_real = FALSE)
stopifnot(!anyDuplicated(out$KO_mapping[c("feature_id", "kegg_term")]))
stopifnot(!anyDuplicated(out$PATHWAY_mapping[c("feature_id", "kegg_term")]))
stopifnot(!any(grepl("^ko:", out$KO_mapping$kegg_term)))
stopifnot(all(nzchar(out$KO_mapping$kegg_term)))
stopifnot(all(nzchar(out$PATHWAY_mapping$kegg_term)))
stopifnot("K00001" %in% out$KO_mapping$kegg_term)
stopifnot("map00010" %in% out$PATHWAY_mapping$kegg_term)

four <- out$gene_set_status[out$gene_set_status$gene_set_id == "set_four", ]
stopifnot(nrow(four) == 2L)
stopifnot(all(four$query_size == 4L))
stopifnot(all(four$status == "SKIPPED_QUERY_SIZE_LT_5"))
stopifnot(!any(out$results$gene_set_id == "set_four"))

ko_a <- out$results[
  out$results$gene_set_id == "set_five_a" & out$results$enrichment_type == "KO",
  , drop = FALSE
]
path_a <- out$results[
  out$results$gene_set_id == "set_five_a" & out$results$enrichment_type == "PATHWAY",
  , drop = FALSE
]
stopifnot(unique(ko_a$query_size) == 5L, unique(ko_a$background_size) == 10L)
stopifnot(unique(path_a$query_size) == 5L, unique(path_a$background_size) == 8L)
term <- ko_a[ko_a$kegg_term == "K00001", , drop = FALSE]
expected_p <- stats::phyper(3L - 1L, 3L, 10L - 3L, 5L, lower.tail = FALSE)
stopifnot(identical(term$p_value, expected_p))
stopifnot(identical(term$query_fraction, 3 / 5))
stopifnot(identical(term$background_fraction, 3 / 10))
stopifnot(identical(term$enrichment_ratio, 2))
stopifnot(identical(term$contributing_feature_ids, "f1;f2;f3"))

for (id in unique(out$results$gene_set_id)) {
  for (type in unique(out$results$enrichment_type)) {
    rows <- out$results$gene_set_id == id & out$results$enrichment_type == type
    if (any(rows)) stopifnot(identical(
      out$results$adjusted_p_value[rows],
      stats::p.adjust(out$results$p_value[rows], method = "BH")
    ))
  }
}

# Each gene set and each enrichment type forms an independent BH family.
single <- .run_current_kegg_enrichment(
  list(current_gene_set_manifest = manifest[2, , drop = FALSE], universes = gene_sets$universes),
  annotation,
  validate_real = FALSE
)
single_ko <- single$results[single$results$enrichment_type == "KO", ]
single_path <- single$results[single$results$enrichment_type == "PATHWAY", ]
stopifnot(identical(ko_a$adjusted_p_value, single_ko$adjusted_p_value))
stopifnot(identical(path_a$adjusted_p_value, single_path$adjusted_p_value))
stopifnot(identical(ko_a$significant, ko_a$adjusted_p_value <= 0.05))

ordered <- with(out$results, order(
  gene_set_id, enrichment_type, adjusted_p_value, p_value,
  -enrichment_ratio, kegg_term, method = "radix"
))
stopifnot(identical(ordered, seq_len(nrow(out$results))))
shuffled <- annotation[rev(seq_len(nrow(annotation))), , drop = FALSE]
out_shuffled <- .run_current_kegg_enrichment(gene_sets, shuffled, validate_real = FALSE)
stopifnot(identical(out$KO_mapping, out_shuffled$KO_mapping))
stopifnot(identical(out$PATHWAY_mapping, out_shuffled$PATHWAY_mapping))
stopifnot(identical(out$results, out_shuffled$results))

text <- readLines(file.path("scripts", "current", "run_current_kegg_enrichment.R"))
stopifnot(sum(grepl("build_current_enrichment_gene_sets\\(annotation_master_path\\)", text)) == 1L)
stopifnot(any(grepl("universes\\$ko_universe", text)))
stopifnot(any(grepl("universes\\$pathway_universe", text)))
stopifnot(any(grepl("stats::phyper", text)))
stopifnot(any(grepl("lower.tail = FALSE", text, fixed = TRUE)))
stopifnot(any(grepl('method = "BH"', text, fixed = TRUE)))
stopifnot(!any(grepl(
  "FSF_condition_class_annotation_master|enrichment_gene_set_manifest|results/real_data",
  text
)))
stopifnot(!any(grepl("stability_level|fsf_signal_class|weakly_stable|instable", text)))
stopifnot(!any(grepl("run_current_go|eggnog_go|interpro_go", text)))
stopifnot(!any(grepl("https?://|KEGGREST|enrichKEGG", text)))
stopifnot(!any(grepl("write_tsv|write\\.table|saveRDS|sink\\(|file\\.create", text)))
stopifnot(!any(grepl("0\\.60|0\\.6[^0-9]", text)))

cat("current KEGG enrichment workflow tests: PASS\n")
