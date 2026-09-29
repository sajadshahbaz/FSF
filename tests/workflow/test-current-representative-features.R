source(file.path("scripts", "current", "build_current_representative_features.R"))

annotation_names <- .fsf_annotation_columns[-1L]
make_joined <- function() {
  rows <- data.frame(
    condition = rep("TEST", 8L),
    feature_id = c("z_tie", "a_tie", "constant", "transitional", "down", "stable", "unannotated", "low_stability_case"),
    n_perturbations = 4,
    n_up = c(4, 4, 0, 3, 0, 3, 4, 2),
    n_down = c(0, 0, 0, 0, 4, 1, 0, 2),
    n_const = c(0, 0, 4, 1, 0, 0, 0, 0),
    p_up = c(1, 1, 0, .6, 0, .75, 1, .5),
    p_down = c(0, 0, 0, .2, 1, .25, 0, .5),
    p_const = c(0, 0, 1, .2, 0, 0, 0, 0),
    dominant_state = c("up", "up", "constant", "up", "down", "up", "up", "tied"),
    ssi = c(1, 1, 1, .6, 1, .75, 1, .5),
    stability_region = c("Highly Stable", "Highly Stable", "Highly Stable", "Transitional", "Highly Stable", "Stable", "Highly Stable", "Low Stability"),
    signal_class = c("Highly Stable Up", "Highly Stable Up", "Highly Stable Constant", "Transitional Up", "Highly Stable Down", "Stable Up", "Highly Stable Up", "Low Stability"),
    stability_deviation = c(0, 0, 0, .4, 0, .25, 0, .5),
    stringsAsFactors = FALSE
  )
  annotation <- setNames(vector("list", length(annotation_names)), annotation_names)
  for (name in annotation_names) annotation[[name]] <- rep(NA_character_, nrow(rows))
  annotation$locus_tag <- paste0("locus_", rows$feature_id)
  annotation$protein_id <- paste0("protein_", rows$feature_id)
  annotation$protein_length <- seq_len(nrow(rows)) + 100L
  annotation$translation_status <- "translated_clean"
  annotation$eggnog_annotated <- TRUE
  annotation$eggnog_evalue <- seq_len(nrow(rows)) * 1e-20
  annotation$eggnog_description <- 'quoted "annotation"; unchanged'
  annotation$interpro_annotated <- FALSE
  annotation$interpro_signature_descriptions <- 'literal ""quoted"" text'
  annotation$n_interpro_rows <- 0L
  annotation$pfam_annotated <- FALSE
  annotation$n_pfam_hits <- 0L
  annotation$annotation_source_count <- c(2L, 2L, 1L, 3L, 2L, 1L, 0L, 1L)
  annotation$annotation_status <- "fixture"
  cbind(rows, as.data.frame(annotation, stringsAsFactors = FALSE, check.names = FALSE))
}

x <- make_joined()
stopifnot(identical(
  names(dimnames(.fsf_representative_expected_directional)),
  c("condition", "dominant_state")
))
out <- .build_current_representative_features(x, validate_real = FALSE)
stopifnot(identical(names(out), c(
  "current_fsf_stable_signal_catalog",
  "current_fsf_stable_directional_signals",
  "current_fsf_stable_constant_signals",
  "current_fsf_stable_signal_summary",
  "current_fsf_top100_stable_directional_signals",
  "current_fsf_stable_signal_annotated_catalog",
  "current_fsf_stable_signal_annotation_coverage",
  "current_fsf_top100_annotated_stable_directional_signals",
  "current_fsf_top50_annotated_stable_directional_signals_compact"
)))
stable <- out$current_fsf_stable_signal_catalog
stopifnot(setequal(stable$feature_id, c("z_tie", "a_tie", "constant", "down", "stable", "unannotated")))
stopifnot(!any(stable$stability_region %in% c("Transitional", "Low Stability")))
stopifnot("constant" %in% stable$dominant_state)
stopifnot(all(out$current_fsf_stable_directional_signals$dominant_state %in% c("up", "down")))
top_directional_up <- out$current_fsf_top100_stable_directional_signals
top_directional_up <- top_directional_up[top_directional_up$dominant_state == "up", , drop = FALSE]
stopifnot(identical(head(top_directional_up$feature_id, 3L), c("a_tie", "unannotated", "z_tie")))
annotated <- out$current_fsf_stable_signal_annotated_catalog
stopifnot(identical(annotated$eggnog_description, rep('quoted "annotation"; unchanged', nrow(annotated))))
stopifnot(is.numeric(annotated$protein_length), is.numeric(annotated$eggnog_evalue))
stopifnot(!anyDuplicated(annotated[c("condition", "feature_id")]))
stopifnot(!"dominant_signal_identity" %in% names(annotated))
stopifnot(!"fsf_signal_class" %in% names(annotated))
stopifnot(!"stability_level" %in% names(annotated))

# Higher annotation_source_count precedes otherwise stronger SSI; equal rank
# resolves by feature_id.
top_annotated <- out$current_fsf_top100_annotated_stable_directional_signals
top_annotated_up <- top_annotated[top_annotated$dominant_state == "up", , drop = FALSE]
stopifnot(identical(top_annotated_up$feature_id[1L], "a_tie"))

# Public workflow source authority is structurally restricted to the loader.
workflow_text <- readLines(file.path("scripts", "current", "build_current_representative_features.R"))
stopifnot(sum(grepl("load_current_fsf_annotation\\(annotation_master_path\\)", workflow_text)) == 1L)
stopifnot(!any(grepl("results/real_data/(condition_fsf|stable_signal_catalog)|FSF_condition_class_annotation_master", workflow_text)))
stopifnot(!any(grepl("write_tsv|writeLines|write_lines|saveRDS|save\\(", workflow_text)))

cat("current representative-feature workflow tests: PASS\n")
