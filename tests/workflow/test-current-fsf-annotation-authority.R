source(file.path("scripts", "lib", "current_fsf_annotation_authority.R"))

current_authority_path <- file.path(
  "results", "current_fsf_v1", "current_fsf_feature_metrics.tsv"
)
stopifnot(.fsf_current_sha256 %in% c(
  "e7d21744681122f041ad38623bb4cb4807ab3d254132fffd6269ce0ef0b1ab89",
  "d4384daa37a2532cbb601bab87556c2a5e8ad3d87a73690f6ddc4f04cce95a57"
))
observed_current_sha256 <- .fsf_validate_sha256(
  current_authority_path,
  .fsf_current_sha256,
  "Current FSF authority"
)
stopifnot(observed_current_sha256 %in% .fsf_current_sha256)

expect_error <- function(code, pattern) {
  error <- tryCatch({
    force(code)
    NULL
  }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

make_current <- function() {
  data.frame(
    condition = c("UV", "DES"),
    feature_id = c("f1", "f1"),
    n_perturbations = c(7L, 2L),
    n_up = c(4L, 0L),
    n_down = c(1L, 0L),
    n_const = c(2L, 2L),
    p_up = c(4 / 7, 0),
    p_down = c(1 / 7, 0),
    p_const = c(2 / 7, 1),
    dominant_state = c("up", "constant"),
    ssi = c(4 / 7, 1),
    stability_region = c("Transitional", "Highly Stable"),
    signal_class = c("Transitional Up", "Highly Stable Constant"),
    stability_deviation = c(3 / 7, 0),
    stringsAsFactors = FALSE
  )
}

make_annotation <- function(ids = "f1") {
  values <- vector("list", length(.fsf_annotation_columns))
  names(values) <- .fsf_annotation_columns
  values$feature_id <- ids
  for (name in setdiff(.fsf_annotation_columns, "feature_id")) {
    values[[name]] <- rep(NA_character_, length(ids))
  }
  values$locus_tag <- paste0("locus_", ids)
  values$protein_id <- paste0("protein_", ids)
  values$protein_length <- as.integer(seq_along(ids) + 100L)
  values$translation_status <- "translated_clean"
  values$eggnog_annotated <- TRUE
  values$eggnog_evalue <- seq_along(ids) * 1e-20
  values$eggnog_description <- 'quoted "annotation"; unchanged'
  values$interpro_annotated <- FALSE
  values$interpro_signature_descriptions <- 'literal ""zincins"" text'
  values$n_interpro_rows <- 0L
  values$pfam_annotated <- FALSE
  values$n_pfam_hits <- 0L
  values$annotation_source_count <- 1L
  values$annotation_status <- "single_source_annotation"
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}

current <- make_current()
annotation <- make_annotation()
joined <- .fsf_join_current_annotation(
  current,
  annotation,
  expected_current_rows = 2L,
  expected_annotation_rows = 1L
)
stopifnot(nrow(joined) == 2L, ncol(joined) == 49L)
stopifnot(identical(joined[.fsf_current_columns], current))
stopifnot(identical(
  joined$interpro_signature_descriptions,
  rep(annotation$interpro_signature_descriptions, 2L)
))
stopifnot(identical(joined$eggnog_evalue, rep(annotation$eggnog_evalue, 2L)))

duplicate_annotation <- rbind(annotation, annotation)
expect_error(
  .fsf_join_current_annotation(current, duplicate_annotation, 2L, 2L),
  "duplicate feature_id"
)

unmatched_annotation <- make_annotation("other")
expect_error(
  .fsf_join_current_annotation(current, unmatched_annotation, 2L, 1L),
  "unmatched current feature_id"
)

wrong_schema <- annotation[-2L]
expect_error(
  .fsf_join_current_annotation(current, wrong_schema, 2L, 1L),
  "exact required 36-column schema"
)

missing_current_column <- current[-3L]
expect_error(
  .fsf_join_current_annotation(missing_current_column, annotation, 2L, 1L),
  "exact required 14-column schema"
)

duplicate_current <- rbind(current, current[1L, , drop = FALSE])
expect_error(
  .fsf_join_current_annotation(duplicate_current, annotation, 3L, 1L),
  "duplicate condition"
)

blank_annotation <- annotation
blank_annotation$feature_id <- ""
expect_error(
  .fsf_join_current_annotation(current, blank_annotation, 2L, 1L),
  "missing or blank feature_id"
)

missing_file <- file.path(tempdir(), "does-not-exist.tsv")
expect_error(
  .fsf_validate_sha256(missing_file, "not-a-hash", "Annotation authority"),
  "does not exist"
)

wrong_hash_file <- tempfile(fileext = ".tsv")
writeLines("fixture", wrong_hash_file)
expect_error(
  .fsf_validate_sha256(wrong_hash_file, paste(rep("0", 64L), collapse = ""),
                       "Annotation authority"),
  "SHA-256 mismatch"
)
unlink(wrong_hash_file)

governed_old_authority_file <- file.path(
  Sys.getenv("FSF_EXTERNAL_ARTIFACT_ROOT"),
  "current_fsf_v1", "manuscript", "historical_authority",
  "current_fsf_feature_metrics.tsv"
)
stopifnot(file.exists(governed_old_authority_file))
stopifnot(identical(
  .fsf_sha256(governed_old_authority_file),
  "f0874a99572884064d52ae0a77af7f4fe95863c1e8776be505fcd358827b095e"
))
old_authority_file <- tempfile(fileext = ".tsv")
stopifnot(file.copy(governed_old_authority_file, old_authority_file))
stopifnot(identical(
  .fsf_sha256(old_authority_file),
  "f0874a99572884064d52ae0a77af7f4fe95863c1e8776be505fcd358827b095e"
))
expect_error(
  .fsf_validate_sha256(
    old_authority_file,
    .fsf_current_sha256,
    "Current FSF authority"
  ),
  "SHA-256 mismatch"
)
unlink(old_authority_file)

expect_error(
  .fsf_validate_join_rows(rbind(joined, joined[1L, , drop = FALSE]), current),
  "row multiplication"
)
expect_error(
  .fsf_validate_join_rows(joined[-1L, , drop = FALSE], current),
  "row loss"
)

cat("current FSF annotation authority helper tests: PASS\n")
