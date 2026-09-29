source(file.path("scripts", "current", "build_current_enrichment_gene_sets.R"))

make_enrichment_joined <- function() {
  current <- data.frame(
    condition = c("A", "A", "A", "A", "B", "B", "B", "B"),
    feature_id = c("z", "a", "c", "d", "z", "a", "c", "d"),
    n_perturbations = 4L,
    n_up = c(4L, 0L, 0L, 2L, 4L, 0L, 0L, 2L),
    n_down = c(0L, 4L, 0L, 2L, 0L, 4L, 0L, 2L),
    n_const = c(0L, 0L, 4L, 0L, 0L, 0L, 4L, 0L),
    p_up = c(1, 0, 0, .5, 1, 0, 0, .5),
    p_down = c(0, 1, 0, .5, 0, 1, 0, .5),
    p_const = c(0, 0, 1, 0, 0, 0, 1, 0),
    dominant_state = rep(c("up", "down", "constant", "tied"), 2L),
    ssi = rep(c(1, .8, .6, .5), 2L),
    stability_region = rep(c("Highly Stable", "Stable", "Transitional", "Low Stability"), 2L),
    signal_class = rep(c("Highly Stable Up", "Stable Down", "Transitional Constant", "Low Stability"), 2L),
    stability_deviation = rep(c(0, .2, .4, .5), 2L),
    stringsAsFactors = FALSE
  )
  annotation <- setNames(vector("list", length(.fsf_annotation_columns) - 1L), .fsf_annotation_columns[-1L])
  for (name in names(annotation)) annotation[[name]] <- rep(NA_character_, nrow(current))
  annotation$locus_tag <- paste0("locus_", current$feature_id)
  annotation$protein_id <- paste0("protein_", current$feature_id)
  annotation$protein_length <- rep(c(10L, 20L, 30L, 40L), 2L)
  annotation$translation_status <- "fixture"
  annotation$eggnog_annotated <- rep(c(TRUE, TRUE, FALSE, FALSE), 2L)
  annotation$eggnog_go <- rep(c("GO:0000001,GO:0000002", "-", NA, "bad"), 2L)
  annotation$eggnog_kegg_ko <- rep(c("ko:K00001", "-", NA, "K00002"), 2L)
  annotation$eggnog_kegg_pathway <- rep(c("map00010", "", NA, "map00020"), 2L)
  annotation$interpro_annotated <- rep(c(FALSE, TRUE, FALSE, FALSE), 2L)
  annotation$interpro_go <- rep(c(NA, "GO:0000003", "-", NA), 2L)
  annotation$n_interpro_rows <- 0L
  annotation$pfam_annotated <- FALSE
  annotation$n_pfam_hits <- 0L
  annotation$annotation_source_count <- rep(c(1L, 2L, 0L, 0L), 2L)
  annotation$annotation_status <- "fixture"
  cbind(current, as.data.frame(annotation, stringsAsFactors = FALSE, check.names = FALSE))
}

x <- make_enrichment_joined()
out <- .build_current_enrichment_gene_sets(x, validate_real = FALSE)

stopifnot(identical(names(out), c(
  "main_class_sets", "collapsed_stability_sets", "directional_stable_sets",
  "current_gene_set_manifest", "universes", "uv_migration"
)))
stopifnot(nrow(out$main_class_sets) == 8L)
stopifnot(nrow(out$collapsed_stability_sets) == 8L)
stopifnot(nrow(out$directional_stable_sets) == 4L)
stopifnot(all(is.na(out$collapsed_stability_sets$dominant_state)))
stopifnot(all(is.na(out$collapsed_stability_sets$signal_class)))
stopifnot(all(out$directional_stable_sets$dominant_state %in% c("up", "down")))
stopifnot(!any(out$directional_stable_sets$dominant_state == "constant"))

for (family in out[c("main_class_sets", "collapsed_stability_sets", "directional_stable_sets")]) {
  stopifnot(all(vapply(family$all_features, function(v) {
    identical(v, sort(unique(v), method = "radix"))
  }, logical(1L))))
  stopifnot(all(vapply(family$annotated_features, function(v) {
    identical(v, sort(unique(v), method = "radix"))
  }, logical(1L))))
}
stopifnot(!anyDuplicated(out$current_gene_set_manifest$gene_set_id))
stopifnot(identical(
  out$current_gene_set_manifest$gene_set_id,
  sort(out$current_gene_set_manifest$gene_set_id, method = "radix")
) == FALSE)
ordered <- with(out$current_gene_set_manifest, order(
  gene_set_family, condition, stability_region, dominant_state,
  signal_class, gene_set_id, method = "radix", na.last = TRUE
))
stopifnot(identical(ordered, seq_len(nrow(out$current_gene_set_manifest))))

stopifnot(identical(out$universes$all_feature_universe, c("a", "c", "d", "z")))
stopifnot(identical(out$universes$annotation_universe, c("a", "z")))
stopifnot(identical(out$universes$go_mapped_features, c("a", "z")))
stopifnot(identical(out$universes$ko_universe, c("d", "z")))
stopifnot(identical(out$universes$pathway_universe, c("d", "z")))

main_a_up <- out$main_class_sets[
  out$main_class_sets$condition == "A" &
    out$main_class_sets$dominant_state == "up",
  ,
  drop = FALSE
]
stopifnot(identical(main_a_up$all_features[[1L]], "z"))
stopifnot(identical(main_a_up$annotated_features[[1L]], "z"))

collision_fixture <- x[1:2, , drop = FALSE]
collision_fixture$condition <- c("A B", "A-B")
collision_fixture$stability_region <- "Stable"
collision_fixture$dominant_state <- "up"
collision_fixture$signal_class <- "Stable Up"
collision_error <- tryCatch({
  .build_current_enrichment_gene_sets(collision_fixture, validate_real = FALSE)
  FALSE
}, error = function(e) grepl("collision", conditionMessage(e), fixed = TRUE))
stopifnot(collision_error)

workflow_text <- readLines(file.path("scripts", "current", "build_current_enrichment_gene_sets.R"))
stopifnot(sum(grepl("load_current_fsf_annotation\\(annotation_master_path\\)", workflow_text)) == 1L)
stopifnot(!any(grepl(
  "results/real_data/(condition_fsf|class_annotation|enrichment_gene_sets)|FSF_condition_class_annotation_master",
  workflow_text
)))
stopifnot(!any(grepl("write_tsv|write\\.table|saveRDS|sink\\(|file\\.create", workflow_text)))
stopifnot(!any(grepl("weakly_stable|instable|stability_level|fsf_signal_class", workflow_text)))
stopifnot(sum(grepl("joined\\$ssi < 0\\.60", workflow_text)) == 1L)

cat("current enrichment gene-set workflow tests: PASS\n")
