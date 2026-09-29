#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
current_manuscript_root <- file.path(
  root, "results", "current_fsf_v1", "manuscript"
)

stopifnot(dir.exists(current_manuscript_root))
stopifnot(!dir.exists(file.path(root, "results", "manuscript")))

all_result_dirs <- list.dirs(
  file.path(root, "results"), recursive = TRUE, full.names = TRUE
)
manuscript_dirs <- all_result_dirs[basename(all_result_dirs) == "manuscript"]
archive_root <- normalizePath(
  file.path(root, "results", "archive"), winslash = "/", mustWork = FALSE
)
normalized_manuscript_dirs <- normalizePath(
  manuscript_dirs, winslash = "/", mustWork = TRUE
)
non_archive_manuscript_dirs <- normalized_manuscript_dirs[
  !startsWith(normalized_manuscript_dirs, paste0(archive_root, "/"))
]
stopifnot(identical(
  non_archive_manuscript_dirs,
  normalizePath(current_manuscript_root, winslash = "/", mustWork = TRUE)
))

current_scripts <- list.files(
  file.path(root, "scripts", "current"), recursive = TRUE,
  full.names = TRUE, include.dirs = FALSE
)
current_script_text <- unlist(
  lapply(current_scripts, readLines, warn = FALSE), use.names = FALSE
)

# One historical provenance locator is intentionally retained as metadata
# in the frozen Stage12 table-package index. It is descriptive provenance,
# not an active runtime dependency.
legacy_s1_provenance <- paste0(
  "Git object d13ab123b021672d20485f14b1381214175ac0e9:",
  "results/manuscript/tables/supplementary/locked_submission_tables/",
  "Supplementary_Table_S1_metadata.csv"
)

stopifnot(
  sum(grepl(legacy_s1_provenance, current_script_text, fixed = TRUE)) == 1L
)

current_script_boundary_text <- gsub(
  legacy_s1_provenance,
  "DECLARED_FROZEN_S1_HISTORICAL_PROVENANCE",
  current_script_text,
  fixed = TRUE
)

stopifnot(!any(grepl("results/manuscript", current_script_boundary_text, fixed = TRUE)))
stopifnot(!any(grepl("results/archive/", current_script_boundary_text, fixed = TRUE)))
stopifnot(!any(grepl("scripts/archive/", current_script_boundary_text, fixed = TRUE)))
stopifnot(!any(grepl("scripts/manuscript/", current_script_boundary_text, fixed = TRUE)))
stopifnot(!any(grepl("historical-archive:", current_script_boundary_text, fixed = TRUE)))

readme <- paste(readLines(file.path(root, "README.md"), warn = FALSE), collapse = "\n")
stopifnot(!grepl("results/manuscript/", readme, fixed = TRUE))
stopifnot(grepl("docs/CURRENT_FSF_AUTHORITY.md", readme, fixed = TRUE))
stopifnot(grepl("results/current_fsf_v1/", readme, fixed = TRUE))
stopifnot(grepl("docs/FSF_LEGACY_PROVENANCE.md", readme, fixed = TRUE))
stopifnot(grepl("provenance only", readme, fixed = TRUE))

authority_path <- file.path(root, "docs", "CURRENT_FSF_AUTHORITY.md")
stopifnot(file.exists(authority_path))
authority <- paste(readLines(authority_path, warn = FALSE), collapse = "\n")
required_authority_text <- c(
  "docs/FSF_V1_CURRENT_SCIENTIFIC_LOCK.md",
  "results/current_fsf_v1/current_fsf_feature_metrics.tsv",
  "results/current_fsf_v1/current_fsf_class_architecture.tsv",
  "results/current_fsf_v1/current_fsf_region_architecture.tsv",
  "results/current_fsf_v1/manuscript/", "Genome Research",
  "docs/FSF_LEGACY_PROVENANCE.md", "Figure 1 and Figure 7",
  "scripts/current/", "FSF_EXTERNAL_ARTIFACT_ROOT"
)
stopifnot(all(vapply(
  required_authority_text, grepl, logical(1L), x = authority, fixed = TRUE
)))

forbidden_active_path <- function(x) {
  x <- gsub("\\\\", "/", as.character(x))
  grepl("(^|/)results/(archive|manuscript)(/|$)", x)
}

source(file.path(root, "scripts/current/build_current_external_artifact_manifests.R"))
external <- read.delim(
  file.path(current_manuscript_root, "audit", "external_artifacts.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
external_active <- c(
  file.path(current_manuscript_root, external$relative_repository_path),
  external$external_artifact_path, external$generator_script
)
stopifnot(!any(forbidden_active_path(external_active)))
stopifnot(all(vapply(external$external_artifact_path,
  .fsf_validate_external_locator, logical(1L))))

figures <- read.delim(
  file.path(current_manuscript_root, "figures", "figure_manifest.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
figure_active <- unlist(figures[c(
  "source_authority", "generation_script", "output_png", "output_pdf"
)], use.names = FALSE)
stopifnot(!any(forbidden_active_path(figure_active)))

figure_root <- file.path(current_manuscript_root, "figures")
main_figure_dir <- file.path(figure_root, "main")
supplementary_figure_dir <- file.path(figure_root, "supplementary")
main_figure_ids <- paste0("Figure", 1:7)
supplementary_figure_ids <- paste0("FigureS", 1:4)
main_outputs <- unlist(lapply(main_figure_ids, function(id) {
  file.path(main_figure_dir, paste0(id, c(".png", ".pdf")))
}))
supplementary_outputs <- unlist(lapply(supplementary_figure_ids, function(id) {
  file.path(supplementary_figure_dir, paste0(id, c(".png", ".pdf")))
}))
flat_outputs <- list.files(
  figure_root, pattern = "^Figure(S[1-4]|[1-7])\\.(png|pdf)$", full.names = TRUE
)
stopifnot(length(flat_outputs) == 0L)
stopifnot(all(file.exists(main_outputs)), all(file.exists(supplementary_outputs)))
stopifnot(setequal(
  list.files(main_figure_dir, pattern = "^Figure[1-7]\\.(png|pdf)$"),
  basename(main_outputs)
))
stopifnot(setequal(
  list.files(supplementary_figure_dir, pattern = "^FigureS[1-4]\\.(png|pdf)$"),
  basename(supplementary_outputs)
))
main_identifiers <- sub("\\.(png|pdf)$", "", basename(main_outputs))
supplementary_identifiers <- sub("\\.(png|pdf)$", "", basename(supplementary_outputs))
stopifnot(length(intersect(unique(main_identifiers), unique(supplementary_identifiers))) == 0L)
stopifnot(all(table(main_identifiers) == 2L), all(table(supplementary_identifiers) == 2L))

semantic <- read.delim(
  file.path(current_manuscript_root, "semantic_authority", "semantic_source_manifest.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
semantic_active <- unlist(semantic[c("source_location", "local_artifact")], use.names = FALSE)
stopifnot(!any(forbidden_active_path(semantic_active)))
stopifnot(all(vapply(semantic$local_artifact,
  .fsf_validate_external_locator, logical(1L))))

legacy_provenance_path <- file.path(root, "docs", "FSF_LEGACY_PROVENANCE.md")
historical_manifest_path <- file.path(root, "docs", "provenance", "historical-archive-files.sha256")
legacy_manifest_path <- file.path(root, "docs", "provenance", "legacy-bundle-files.sha256")
stopifnot(file.exists(legacy_provenance_path), file.exists(historical_manifest_path),
  file.exists(legacy_manifest_path))
legacy_provenance <- paste(readLines(legacy_provenance_path, warn = FALSE), collapse = "\n")
stopifnot(grepl("101 files", legacy_provenance, fixed = TRUE),
  grepl("57 files", legacy_provenance, fixed = TRUE),
  grepl("UNASSIGNED_PENDING_PUBLIC_DEPOSIT", legacy_provenance, fixed = TRUE),
  grepl("docs/FSF_V1_CURRENT_SCIENTIFIC_LOCK.md", legacy_provenance, fixed = TRUE))
stopifnot(length(readLines(historical_manifest_path, warn = FALSE)) == 101L,
  length(readLines(legacy_manifest_path, warn = FALSE)) == 57L)

allowed_regions <- c("Low Stability", "Transitional", "Stable", "Highly Stable")
allowed_classes <- c(
  "Low Stability",
  "Transitional Up", "Transitional Down", "Transitional Constant",
  "Stable Up", "Stable Down", "Stable Constant",
  "Highly Stable Up", "Highly Stable Down", "Highly Stable Constant"
)
class_fields <- c("signal_class", "dominant_recovered_signal_class")

tsv_paths <- list.files(
  current_manuscript_root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE
)
for (path in tsv_paths) {
  data <- read.delim(path, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
  if ("stability_region" %in% names(data)) {
    values <- unique(as.character(data$stability_region))
    values <- values[!is.na(values) & nzchar(values)]
    stopifnot(all(values %in% allowed_regions))
  }
  for (field in intersect(names(data), class_fields)) {
    values <- unique(as.character(data[[field]]))
    values <- values[!is.na(values) & nzchar(values)]
    stopifnot(all(values %in% allowed_classes))
    stopifnot(!any(grepl("^Low Stability .+", values)))
  }
}

if (nzchar(Sys.which("pdftotext"))) {
  pdfs <- c(main_outputs[grepl("[.]pdf$", main_outputs)], supplementary_outputs[grepl("[.]pdf$", supplementary_outputs)])
  pdf_tmp_root <- tempdir()
  for (pdf in pdfs) {
    txt <- tempfile(tmpdir = pdf_tmp_root, fileext = ".txt")
    status <- system2("pdftotext", c(pdf, txt))
    stopifnot(status == 0L, file.exists(txt))
    rendered <- readLines(txt, warn = FALSE)
    unlink(txt)
    stopifnot(!any(grepl("Instability", rendered, fixed = TRUE)))
  }
}

cat("current FSF authority boundary tests: PASS\n")
