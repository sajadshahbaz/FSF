# Current FSF v1 real-data annotation authority join layer.
#
# This workflow helper is deliberately outside the FSF package core. It joins
# frozen real-data FSF classifications to the separately frozen biological
# annotation authority. It never calculates or modifies an FSF field.

.fsf_current_columns <- c(
  "condition", "feature_id", "n_perturbations", "n_up", "n_down",
  "n_const", "p_up", "p_down", "p_const", "dominant_state", "ssi",
  "stability_region", "signal_class", "stability_deviation"
)

.fsf_annotation_columns <- c(
  "feature_id", "locus_tag", "protein_id", "protein_length",
  "translation_status", "eggnog_annotated", "eggnog_seed_ortholog",
  "eggnog_evalue", "eggnog_score", "eggnog_ogs", "eggnog_max_annot_lvl",
  "eggnog_cog_category", "eggnog_description", "eggnog_preferred_name",
  "eggnog_go", "eggnog_ec", "eggnog_kegg_ko", "eggnog_kegg_pathway",
  "eggnog_kegg_module", "eggnog_pfam", "interpro_annotated",
  "interpro_analyses", "interpro_signature_accessions",
  "interpro_signature_descriptions", "interpro_ids",
  "interpro_descriptions", "interpro_go", "interpro_pathways",
  "n_interpro_rows", "pfam_annotated", "pfam_domains",
  "pfam_accessions", "pfam_descriptions", "n_pfam_hits",
  "annotation_source_count", "annotation_status"
)

.fsf_current_sha256 <- c(
  frozen = "e7d21744681122f041ad38623bb4cb4807ab3d254132fffd6269ce0ef0b1ab89",
  regenerated = "d4384daa37a2532cbb601bab87556c2a5e8ad3d87a73690f6ddc4f04cce95a57"
)
.fsf_annotation_sha256 <-
  "25627cfcaf544dd2793d9fd2462772d9cb76f4728bbcab8c9d29599520f5d584"

.fsf_stop <- function(message) {
  stop(message, call. = FALSE)
}

.fsf_sha256 <- function(path) {
  executable <- Sys.which("sha256sum")
  if (!nzchar(executable)) {
    .fsf_stop("Cannot validate authority: sha256sum is unavailable.")
  }
  output <- system2(executable, shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    .fsf_stop(paste0("Cannot calculate SHA-256 for: ", path))
  }
  sub("[[:space:]].*$", "", output[[1L]])
}

.fsf_validate_sha256 <- function(path, expected, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || !file.exists(path)) {
    .fsf_stop(paste0(label, " file does not exist: ", path))
  }
  if (!isTRUE(file.info(path)$isdir == FALSE)) {
    .fsf_stop(paste0(label, " path is not a regular file: ", path))
  }
  observed <- .fsf_sha256(path)
  if (!observed %in% unname(expected)) {
    .fsf_stop(paste0(
      label, " SHA-256 mismatch: expected one of ",
      paste(unname(expected), collapse = ", "),
      ", observed ", observed, "."
    ))
  }
  invisible(observed)
}

.fsf_validate_ids <- function(x, label) {
  if (anyNA(x) || any(trimws(as.character(x)) == "")) {
    .fsf_stop(paste0(label, " contains missing or blank feature_id values."))
  }
}

.fsf_validate_current <- function(data, expected_rows = 91122L) {
  if (!is.data.frame(data)) {
    .fsf_stop("Current FSF authority must parse as a data.frame.")
  }
  if (!identical(names(data), .fsf_current_columns)) {
    .fsf_stop("Current FSF authority does not have the exact required 14-column schema.")
  }
  if (nrow(data) != expected_rows) {
    .fsf_stop(paste0("Current FSF authority must contain ", expected_rows, " rows."))
  }
  .fsf_validate_ids(data$feature_id, "Current FSF authority")
  if (anyNA(data$condition) || any(trimws(as.character(data$condition)) == "")) {
    .fsf_stop("Current FSF authority contains missing or blank condition values.")
  }
  if (any(duplicated(data[c("condition", "feature_id")]))) {
    .fsf_stop("Current FSF authority contains duplicate condition + feature_id keys.")
  }
  invisible(data)
}

.fsf_validate_annotation <- function(data, expected_rows = 15187L) {
  if (!is.data.frame(data)) {
    .fsf_stop("Annotation authority must parse as a data.frame.")
  }
  if (!identical(names(data), .fsf_annotation_columns)) {
    .fsf_stop("Annotation authority does not have the exact required 36-column schema.")
  }
  if (nrow(data) != expected_rows) {
    .fsf_stop(paste0("Annotation authority must contain ", expected_rows, " rows."))
  }
  .fsf_validate_ids(data$feature_id, "Annotation authority")
  if (any(duplicated(data$feature_id))) {
    .fsf_stop("Annotation authority contains duplicate feature_id keys.")
  }
  invisible(data)
}

.fsf_validate_join_rows <- function(joined, current) {
  if (nrow(joined) > nrow(current)) {
    .fsf_stop("Annotation join caused row multiplication.")
  }
  if (nrow(joined) < nrow(current)) {
    .fsf_stop("Annotation join caused row loss.")
  }
  invisible(joined)
}

.fsf_join_current_annotation <- function(
    current,
    annotation,
    expected_current_rows = 91122L,
    expected_annotation_rows = 15187L) {
  .fsf_validate_current(current, expected_current_rows)
  .fsf_validate_annotation(annotation, expected_annotation_rows)

  annotation_index <- match(current$feature_id, annotation$feature_id)
  if (anyNA(annotation_index)) {
    missing_features <- unique(current$feature_id[is.na(annotation_index)])
    .fsf_stop(paste0(
      "Annotation authority has unmatched current feature_id values: ",
      paste(head(missing_features, 5L), collapse = ", "), "."
    ))
  }

  joined <- cbind(
    current,
    annotation[annotation_index, .fsf_annotation_columns[-1L], drop = FALSE]
  )
  rownames(joined) <- NULL

  .fsf_validate_join_rows(joined, current)
  if (ncol(joined) != 49L) {
    .fsf_stop("Annotation join did not produce the required 49 columns.")
  }
  if (any(duplicated(joined[c("condition", "feature_id")]))) {
    .fsf_stop("Annotation join produced duplicate condition + feature_id keys.")
  }
  fsf_preserved <- vapply(
    .fsf_current_columns,
    function(name) identical(joined[[name]], current[[name]]),
    logical(1L)
  )
  if (!all(fsf_preserved)) {
    .fsf_stop("Annotation join altered authoritative FSF fields.")
  }
  expected_annotation <- annotation[
    annotation_index,
    .fsf_annotation_columns[-1L],
    drop = FALSE
  ]
  rownames(expected_annotation) <- NULL
  annotation_preserved <- vapply(
    .fsf_annotation_columns[-1L],
    function(name) identical(joined[[name]], expected_annotation[[name]]),
    logical(1L)
  )
  if (!all(annotation_preserved)) {
    .fsf_stop("Annotation join altered authoritative annotation fields.")
  }

  joined
}

#' Load the current FSF annotation view
#'
#' Validates the frozen current FSF table and an explicitly supplied,
#' hash-identified annotation master, then returns their validated in-memory
#' left join by `feature_id`. Nothing is written. The helper does not compute
#' or modify any FSF field.
#'
#' @param annotation_master_path Explicit path to `FSF_v1_annotation_master.tsv`.
#' @param current_fsf_path Path to the frozen current FSF feature authority.
#' @return A 91,122-row data frame with the 14 FSF columns followed by the 35
#'   non-key annotation columns.
load_current_fsf_annotation <- function(
    annotation_master_path,
    current_fsf_path = file.path(
      "results", "current_fsf_v1", "current_fsf_feature_metrics.tsv"
    )) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    .fsf_stop("The workflow helper requires the installed readr package.")
  }

  .fsf_validate_sha256(
    current_fsf_path,
    .fsf_current_sha256,
    "Current FSF authority"
  )
  .fsf_validate_sha256(
    annotation_master_path,
    .fsf_annotation_sha256,
    "Annotation authority"
  )

  current <- readr::read_tsv(current_fsf_path, show_col_types = FALSE)
  annotation <- readr::read_tsv(annotation_master_path, show_col_types = FALSE)

  .fsf_join_current_annotation(current, annotation)
}
