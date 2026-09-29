#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source_dir <- file.path(root, "results", "current_fsf_v1", "manuscript", "source_data")
figure_root <- file.path(root, "results", "current_fsf_v1", "manuscript", "figures")
main_figure_dir <- file.path(figure_root, "main")
supplementary_figure_dir <- file.path(figure_root, "supplementary")
builder <- file.path(root, "scripts", "current", "figures", "build_current_manuscript_figures.R")

main_figure_ids <- paste0("Figure", 1:7)
supplementary_figure_ids <- paste0("FigureS", 1:4)
figure_ids <- c(main_figure_ids, supplementary_figure_ids)
source_files <- c(
  "figure2_simplex_source.tsv", "figure3_architecture_source.tsv",
  "figure4_ssi_ecdf_source.tsv", "figure5_go_source.tsv",
  "figure5_kegg_source.tsv", "figure6_baseline_source.tsv",
  "figure6_noise_source.tsv", "figureS2_tau_source.tsv",
  "figureS3_architecture_source.tsv", "figureS4_ssi_distribution_source.tsv"
)
stopifnot(all(file.exists(file.path(source_dir, source_files))))
stopifnot(file.exists(builder))

main_outputs <- unlist(lapply(main_figure_ids, function(id) {
  file.path(main_figure_dir, paste0(id, c(".png", ".pdf")))
}))
supplementary_outputs <- unlist(lapply(supplementary_figure_ids, function(id) {
  file.path(supplementary_figure_dir, paste0(id, c(".png", ".pdf")))
}))
outputs <- c(main_outputs, supplementary_outputs)
stopifnot(all(file.exists(outputs)))
stopifnot(all(file.info(outputs)$size > 1000L))
protected_pngs <- c(
  file.path(main_figure_dir, "Figure1.png"),
  file.path(main_figure_dir, "Figure7.png")
)
protected_sources <- c(
  file.path(figure_root, "new", "Figure_1.png"),
  file.path(figure_root, "new", "Figure_7.png")
)
stopifnot(identical(
  unname(tools::md5sum(protected_pngs)),
  unname(tools::md5sum(protected_sources))
))
stopifnot(setequal(list.files(main_figure_dir, pattern = "^Figure[1-7]\\.(png|pdf)$"), basename(main_outputs)))
stopifnot(setequal(list.files(supplementary_figure_dir, pattern = "^FigureS[1-4]\\.(png|pdf)$"), basename(supplementary_outputs)))
stopifnot(length(list.files(figure_root, pattern = "^Figure(S[1-4]|[1-7])\\.(png|pdf)$")) == 0L)

manifest_path <- file.path(figure_root, "figure_manifest.tsv")
audit_path <- file.path(figure_root, "figure_audit.md")
stopifnot(file.exists(manifest_path), file.exists(audit_path))
manifest <- read.delim(manifest_path, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
expected_manifest_columns <- c(
  "figure_id", "panel", "source_authority", "generation_script",
  "output_png", "output_pdf", "numerical_status",
  "interpretation_status", "notes"
)
stopifnot(identical(names(manifest), expected_manifest_columns))
stopifnot(identical(manifest$figure_id, figure_ids))
stopifnot(anyDuplicated(manifest$figure_id) == 0L)

architecture <- read.delim(
  file.path(source_dir, "figure3_architecture_source.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
expected <- data.frame(
  condition = rep(c("DES", "GAM", "HT", "LT", "OSM", "UV"), each = 4),
  stability_region = rep(c("Low Stability", "Transitional", "Stable", "Highly Stable"), 6),
  region_count = c(
    5837, 0, 0, 9350,
    1085, 2106, 6127, 5869,
    0, 0, 0, 15187,
    1211, 0, 0, 13976,
    0, 0, 0, 15187,
    2948, 6697, 3151, 2391
  ),
  stringsAsFactors = FALSE
)
observed <- architecture[, c("condition", "stability_region", "region_count")]
observed <- observed[order(match(observed$condition, unique(expected$condition)),
                           match(observed$stability_region, unique(expected$stability_region))), ]
rownames(observed) <- NULL
stopifnot(identical(observed$condition, expected$condition))
stopifnot(identical(observed$stability_region, expected$stability_region))
stopifnot(identical(as.integer(observed$region_count), as.integer(expected$region_count)))
stopifnot(all(rowsum(architecture$region_count, architecture$condition)[, 1] == 15187L))
uv <- architecture[architecture$condition == "UV", ]
stopifnot(uv$stability_region[which.max(uv$region_count)] == "Transitional")

tau <- read.delim(
  file.path(source_dir, "figureS2_tau_source.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
stopifnot(identical(names(tau), c("tau", "truth_scenario", "metric", "category", "value")))
builder_text <- readLines(builder, warn = FALSE)
builder_code <- paste(builder_text, collapse = "\n")
stopifnot(grepl(
  "figureS2\\s*<-\\s*ggplot\\s*\\(\\s*s2\\s*,\\s*aes\\s*\\(\\s*tau\\s*,\\s*value\\s*,",
  builder_code, perl = TRUE
))
stopifnot(all(vapply(
  c("color", "linetype", "shape", "group"),
  function(aesthetic) grepl(
    paste0(aesthetic, "\\s*=\\s*generator_display"),
    builder_code, perl = TRUE
  ),
  logical(1)
)))
stopifnot(!grepl(
  "ggplot\\s*\\(\\s*s2\\s*,\\s*aes\\s*\\(\\s*tau\\s*,\\s*tau(?:\\s*[,\\)])",
  builder_code, perl = TRUE
))

pdf_text <- character()
for (pdf in outputs[grepl("[.]pdf$", outputs)]) {
  txt <- tempfile(fileext = ".txt")
  status <- system2("pdftotext", c(pdf, txt))
  stopifnot(status == 0L, file.exists(txt))
  pdf_text <- c(pdf_text, readLines(txt, warn = FALSE))
  unlink(txt)
}
forbidden_terms <- c("Instability", "instable", "weakly stable", "weakly_stable")
pdf_text_lower <- tolower(pdf_text)
stopifnot(!any(vapply(tolower(forbidden_terms), function(x) any(grepl(x, pdf_text_lower, fixed = TRUE)), logical(1))))

figure5_text_file <- tempfile(fileext = ".txt")
stopifnot(system2("pdftotext", c(file.path(main_figure_dir, "Figure5.pdf"), figure5_text_file)) == 0L)
figure5_text <- readLines(figure5_text_file, warn = FALSE)
unlink(figure5_text_file)
dropped <- c(
  "PROFILE_002", "PROFILE_003", "PROFILE_008",
  "PROFILE_009", "PROFILE_013", "PROFILE_026",
  "GO_THEME_001", "KEGG_THEME_001", "KEGG_THEME_002",
  "Other KEGG", "manual review"
)
stopifnot(!any(vapply(dropped, function(x) any(grepl(x, figure5_text, fixed = TRUE)), logical(1))))

audit <- readLines(audit_path, warn = FALSE)
stopifnot(any(grepl("Population: all 15,187 unique analyzed features", audit, fixed = TRUE)))
stopifnot(any(grepl("Endpoint: mean Signal Stratification Index", audit, fixed = TRUE)))
stopifnot(any(grepl("Generator labels describe synthetic distributions", audit, fixed = TRUE)))
stopifnot(any(grepl("Stability Deviation equals 1 - SSI", audit, fixed = TRUE)))
stopifnot(any(grepl("UV shown as Transitional-dominant", audit, fixed = TRUE)))

protected <- c(
  "R", "results/current_fsf_v1/manuscript/source_data",
  "results/current_fsf_v1/manuscript/review",
  "scripts/current/audit_current_manuscript_bundle.R",
  "scripts/current/build_current_manuscript_tables.R"
)
source("tests/workflow/clean_release_hash_guard.R")
assert_paths_match_clean_manifest(protected)
status <- 0L
stopifnot(status == 0L)

cat("current manuscript figure tests: PASS\n")
