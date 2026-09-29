#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
review <- file.path(root, "results", "current_fsf_v1", "manuscript", "review")
final_packet_path <- file.path(review, "author_decision_packet_FINAL.tsv")
evidence_path <- file.path(review, "deferred_profile_term_evidence.tsv")
summary_path <- file.path(review, "final_biological_interpretation_summary.md")

stopifnot(file.exists(final_packet_path), file.exists(evidence_path), file.exists(summary_path))

read_tsv_text <- function(path) {
  read.delim(
    path,
    check.names = FALSE,
    na.strings = NULL,
    quote = "",
    colClasses = "character",
    stringsAsFactors = FALSE
  )
}

locked <- read_tsv_text(file.path(review, "author_decision_packet_LOCKED.tsv"))
final <- read_tsv_text(final_packet_path)
evidence <- read_tsv_text(evidence_path)

allowed_ids <- c(
  "PROFILE_002", "PROFILE_003", "PROFILE_008",
  "PROFILE_009", "PROFILE_013", "PROFILE_026"
)
allowed_decisions <- c(
  "RETAIN", "RETAIN_WITH_REVISED_WORDING", "PARTIALLY_RETAIN", "DROP"
)

stopifnot(nrow(final) == 107L, anyDuplicated(final$decision_id) == 0L)
stopifnot(identical(final$decision_id, locked$decision_id))
stopifnot(setequal(final$decision_id, locked$decision_id))
stopifnot(sum(final$author_decision == "REVIEW_FURTHER") == 0L)
stopifnot(all(final$author_decision %in% allowed_decisions))

changed <- final$decision_id[
  final$author_decision != locked$author_decision |
    final$author_notes != locked$author_notes
]
stopifnot(identical(changed, allowed_ids))
unchanged <- !final$decision_id %in% allowed_ids
stopifnot(identical(final[unchanged, ], locked[unchanged, ]))
stopifnot(all(final$author_decision[match(allowed_ids, final$decision_id)] == "DROP"))

expected_columns <- c(
  "decision_id", "condition", "stability_region", "dominant_state",
  "source", "term_id", "term_description", "adjusted_p",
  "enrichment_measure", "historical_interpretation", "current_evidence_note"
)
stopifnot(identical(names(evidence), expected_columns))
stopifnot(nrow(evidence) == 507L)
stopifnot(setequal(unique(evidence$decision_id), allowed_ids))
stopifnot(all(evidence$source %in% c("GO", "KEGG")))
stopifnot(all(evidence$term_description == "unavailable in frozen authority"))

expected_go <- c(2L, 0L, 180L, 2L, 0L, 0L)
expected_kegg <- c(74L, 63L, 2L, 78L, 8L, 98L)
observed_go <- table(factor(evidence$decision_id[evidence$source == "GO"], levels = allowed_ids))
observed_kegg <- table(factor(evidence$decision_id[evidence$source == "KEGG"], levels = allowed_ids))
stopifnot(identical(as.integer(observed_go), expected_go))
stopifnot(identical(as.integer(observed_kegg), expected_kegg))
stopifnot(all(table(factor(evidence$decision_id, levels = allowed_ids)) > 0L))

# DROP applies to the unsupported historical interpretation only. The frozen
# quantitative GO/KEGG evidence remains present for every resolved profile.
stopifnot(all(final$author_decision[match(allowed_ids, final$decision_id)] == "DROP"))
stopifnot(all(nzchar(evidence$term_id)))
stopifnot(all(nzchar(evidence$adjusted_p)))
stopifnot(all(nzchar(evidence$enrichment_measure)))

run_and_require_pass <- function(script, marker, args = character()) {
  output <- system2("Rscript", c(script, args), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  stopifnot(is.null(status) || status == 0L)
  stopifnot(any(grepl(marker, output, fixed = TRUE)))
}

run_and_require_pass(
  file.path(root, "scripts", "current", "audit_current_manuscript_bundle.R"),
  "current manuscript bundle audit: PASS",
  file.path(root, "results", "current_fsf_v1", "manuscript")
)
run_and_require_pass(
  file.path(root, "tests", "workflow", "test-current-biological-review-package.R"),
  "current biological review package tests: PASS"
)
run_and_require_pass(
  file.path(root, "tests", "workflow", "test-current-author-decision-lock.R"),
  "current author biological decision lock tests: PASS"
)

cat("current final biological decision lock tests: PASS\n")
