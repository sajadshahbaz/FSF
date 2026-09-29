#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
review_dir <- file.path(root, "results", "current_fsf_v1", "manuscript", "review")
packet_path <- file.path(review_dir, "author_decision_packet.tsv")
locked_path <- file.path(review_dir, "author_decision_packet_LOCKED.tsv")
summary_path <- file.path(review_dir, "author_decision_lock_summary.md")

stopifnot(file.exists(packet_path), file.exists(locked_path), file.exists(summary_path))

read_packet <- function(path) {
  read.delim(
    path,
    check.names = FALSE,
    na.strings = NULL,
    quote = "",
    colClasses = "character",
    stringsAsFactors = FALSE
  )
}

packet <- read_packet(packet_path)
locked <- read_packet(locked_path)
stopifnot(identical(packet, locked))
stopifnot(nrow(locked) == 107L, anyDuplicated(locked$decision_id) == 0L)

allowed <- c(
  "RETAIN",
  "RETAIN_WITH_REVISED_WORDING",
  "PARTIALLY_RETAIN",
  "DROP",
  "REVIEW_FURTHER"
)
stopifnot(all(nzchar(locked$author_decision)))
stopifnot(all(locked$author_decision %in% allowed))
stopifnot(all(nzchar(locked$author_notes)))

expected_counts <- c(
  RETAIN = 2L,
  RETAIN_WITH_REVISED_WORDING = 88L,
  PARTIALLY_RETAIN = 8L,
  DROP = 3L,
  REVIEW_FURTHER = 6L
)
observed_counts <- table(factor(locked$author_decision, levels = names(expected_counts)))
stopifnot(identical(as.integer(observed_counts), unname(expected_counts)))

decision_for <- function(id) locked$author_decision[match(id, locked$decision_id)]
stopifnot(all(decision_for(sprintf("REP_%03d", 1:72)) == "RETAIN_WITH_REVISED_WORDING"))
stopifnot(decision_for("GO_THEME_001") == "DROP")
stopifnot(decision_for("KEGG_THEME_001") == "DROP")
stopifnot(decision_for("KEGG_THEME_002") == "DROP")
stopifnot(decision_for("ARCH_001") == "PARTIALLY_RETAIN")
stopifnot(decision_for("ARCH_004") == "RETAIN")
stopifnot(decision_for("ARCH_006") == "RETAIN")

unresolved <- locked$decision_id[locked$author_decision == "REVIEW_FURTHER"]
stopifnot(identical(
  unresolved,
  c("PROFILE_002", "PROFILE_003", "PROFILE_008", "PROFILE_009", "PROFILE_013", "PROFILE_026")
))

expected_ids <- c(
  sprintf("PROFILE_%03d", 1:26),
  "GO_THEME_001",
  sprintf("KEGG_THEME_%03d", 1:2),
  sprintf("ARCH_%03d", 1:6),
  sprintf("REP_%03d", 1:72)
)
stopifnot(setequal(locked$decision_id, expected_ids))

summary <- readLines(summary_path, warn = FALSE)
required_rules <- c(
  "FSF method is frozen.",
  "No additional threshold redesign.",
  "No new dataset or analytical branch.",
  "Absence of significant enrichment may be reported descriptively.",
  "Absence of enrichment is not evidence of absence of biology.",
  "Broad/manual catch-all themes are not manuscript biological claims.",
  "Partial support means only the currently supported component is retained.",
  "Causal/adaptive/fitness interpretations are not inferred from architecture or enrichment alone.",
  "UV is Transitional-dominant under the current authority.",
  "Representative-feature discovery is closed; only rationale revision is allowed."
)
stopifnot(all(vapply(required_rules, function(x) any(grepl(x, summary, fixed = TRUE)), logical(1))))

review_files <- list.files(review_dir, recursive = TRUE, full.names = TRUE)
stopifnot(!any(grepl("\\.(pdf|png|jpe?g|svg|tiff?)$", review_files, ignore.case = TRUE)))
stopifnot(!any(grepl("manuscript|caption", basename(review_files), ignore.case = TRUE)))

bundle_audit <- system2(
  "Rscript",
  c(file.path(root, "scripts", "current", "audit_current_manuscript_bundle.R"),
    file.path(root, "results", "current_fsf_v1", "manuscript")),
  stdout = TRUE,
  stderr = TRUE
)
stopifnot(is.null(attr(bundle_audit, "status")) || attr(bundle_audit, "status") == 0L)
stopifnot(any(grepl("current manuscript bundle audit: PASS", bundle_audit, fixed = TRUE)))

cat("current author biological decision lock tests: PASS\n")
