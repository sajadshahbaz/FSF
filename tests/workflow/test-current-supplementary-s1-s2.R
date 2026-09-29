#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
builder <- file.path(root, "scripts", "current", "figures", "build_current_manuscript_figures.R")
source_path <- file.path(root, "results", "current_fsf_v1", "manuscript", "source_data", "figureS2_tau_source.tsv")
stopifnot(file.exists(builder), file.exists(source_path))
text <- readLines(builder, warn = FALSE)
joined <- paste(text, collapse = "\n")

# The governed full annotation authority, never the representative subset, supplies S1.
stopifnot(!grepl("representative_feature_review.tsv", joined, fixed = TRUE))
stopifnot(grepl("nrow(annotation_view) != 91122L", joined, fixed = TRUE))
stopifnot(grepl("length(unique(annotation_view$feature_id)) != 15187L", joined, fixed = TRUE))
stopifnot(grepl("source_expected <- c(`0` = 1720L, `1` = 2214L, `2` = 1848L, `3` = 9405L)", joined, fixed = TRUE))
stopifnot(grepl("resource_expected <- c(eggNOG = 9811L, InterPro = 12480L, Pfam = 11834L)", joined, fixed = TRUE))
stopifnot(grepl("Panel A categories are mutually exclusive", joined, fixed = TRUE))
stopifnot(grepl("Panel B resource categories overlap", joined, fixed = TRUE))
s1_plot_start <- grep("^  s1a <- ggplot\\(", text)
s1_plot_end <- grep("^  figureS1 <-", text)
stopifnot(length(s1_plot_start) == 1L, length(s1_plot_end) == 1L, s1_plot_start < s1_plot_end)
s1_plot_text <- paste(text[s1_plot_start:s1_plot_end], collapse = "\n")
common_s1_scale <- paste(c(
  "scale_y_continuous(",
  "limits = c(0, 1.00),",
  "breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),",
  "labels = percent_format(accuracy = 1),",
  "expand = expansion(mult = c(0, 0.02))",
  ")"
), collapse = "\n")
normalized_s1_plot_text <- gsub("^[[:space:]]+", "", strsplit(s1_plot_text, "\n")[[1L]])
normalized_s1_plot_text <- paste(normalized_s1_plot_text, collapse = "\n")
scale_hits <- gregexpr(common_s1_scale, normalized_s1_plot_text, fixed = TRUE)[[1L]]
stopifnot(!identical(scale_hits, -1L), length(scale_hits) == 2L)

annotation_master <- Sys.getenv("FSF_ANNOTATION_MASTER", unset = "")
if (!nzchar(annotation_master)) stop("FSF_ANNOTATION_MASTER is required for the focused S1/S2 test.")
source(file.path(root, "scripts", "lib", "current_fsf_annotation_authority.R"))
annotation_view <- load_current_fsf_annotation(
  annotation_master_path = annotation_master,
  current_fsf_path = file.path(root, "results", "current_fsf_v1", "current_fsf_feature_metrics.tsv")
)
fields <- c("annotation_source_count", "annotation_status", "eggnog_annotated",
            "interpro_annotated", "pfam_annotated")
stopifnot(nrow(annotation_view) == 91122L)
stopifnot(length(unique(annotation_view$feature_id)) == 15187L)
stopifnot(all(table(annotation_view$feature_id) == 6L))
for (field in fields) {
  stopifnot(all(vapply(split(annotation_view[[field]], annotation_view$feature_id),
                         function(value) length(unique(value)) == 1L, logical(1))))
}
feature_view <- annotation_view[!duplicated(annotation_view$feature_id), ]
stopifnot(nrow(feature_view) == 15187L, anyDuplicated(feature_view$feature_id) == 0L)
stopifnot(identical(as.integer(table(factor(feature_view$annotation_source_count, levels = 0:3))),
                    c(1720L, 2214L, 1848L, 9405L)))
stopifnot(identical(c(sum(feature_view$eggnog_annotated),
                      sum(feature_view$interpro_annotated),
                      sum(feature_view$pfam_annotated)),
                    c(9811L, 12480L, 11834L)))

# S2 preserves raw generator IDs but displays generator-specific publication labels.
tau <- read.delim(source_path, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
stopifnot(identical(sort(unique(tau$tau)), c(0.25, 0.50, 0.75, 1.00)))
stopifnot(setequal(unique(tau$truth_scenario),
                   c("instable", "stable_constant", "stable_down", "stable_up", "weakly_stable")))
mean_rows <- tau[tau$metric == "mean_ssi", ]
stopifnot(nrow(mean_rows) == 20L)
stopifnot(grepl('s2 <- s2[s2$metric == "mean_ssi", ]', joined, fixed = TRUE))
stopifnot(!grepl('s2$metric %in% c("mean_ssi", "mean_stability_deviation")', joined, fixed = TRUE))
labels <- c("No-majority generator", "Constant-dominant generator",
            "Down-dominant generator", "Up-dominant generator", "Weak-majority generator")
stopifnot(all(vapply(labels, grepl, logical(1), x = joined, fixed = TRUE)))
stopifnot(grepl("Synthetic generator scenario", joined, fixed = TRUE))
stopifnot(grepl("Mean Signal Stratification Index (SSI)", joined, fixed = TRUE))
stopifnot(grepl("linetype = generator_display", joined, fixed = TRUE))
stopifnot(grepl("shape = generator_display", joined, fixed = TRUE))
s2_plot_start <- grep("^  figureS2 <- ggplot\\($", text)
s2_plot_end <- grep("^  save_figure\\(figureS2,", text)
stopifnot(length(s2_plot_start) == 1L, length(s2_plot_end) == 1L, s2_plot_start < s2_plot_end)
s2_plot_text <- paste(text[s2_plot_start:s2_plot_end], collapse = "\n")
stopifnot(grepl(
  "Supplementary Figure S2. Sensitivity of synthetic scenarios to the signal-state threshold",
  s2_plot_text, fixed = TRUE
))
stopifnot(grepl(
  'x\\s*=\\s*expression\\s*\\(\\s*paste\\s*\\(\\s*"Signal-state threshold, "\\s*,\\s*tau\\s*\\)\\s*\\)',
  s2_plot_text, perl = TRUE
))
stopifnot(grepl('y = "Mean Signal Stratification Index (SSI)"', s2_plot_text, fixed = TRUE))
stopifnot(grepl("limits\\s*=\\s*c\\(\\s*1/3\\s*,\\s*1\\.0\\s*\\)", s2_plot_text, perl = TRUE))
stopifnot(grepl("breaks\\s*=\\s*c\\(\\s*1/3\\s*,\\s*0\\.50\\s*,\\s*0\\.75\\s*,\\s*1\\.00\\s*\\)", s2_plot_text, perl = TRUE))
stopifnot(grepl('labels = c("1/3", "0.50", "0.75", "1.00")', s2_plot_text, fixed = TRUE))
stopifnot(grepl("Synthetic generator scenario", s2_plot_text, fixed = TRUE))
stopifnot(grepl("Up- and down-dominant profiles overlap through tau = 0.75.", s2_plot_text, fixed = TRUE))
stopifnot(!grepl("Scenario names describe synthetic data-generation distributions", s2_plot_text, fixed = TRUE))
stopifnot(!grepl("Stability Deviation is", s2_plot_text, fixed = TRUE))

# The targeted mode encloses all unrelated construction and leaves metadata untouched.
figure2_comment <- match(
  "# Figure 2: conceptual geometry plus observed six-condition distributions.",
  text
)
s1_comment <- match(
  "  # Supplementary Figure S1: complete, feature-invariant annotation authority.",
  text
)
main_start <- which(
  seq_along(text) < figure2_comment &
    text == "if (!supplementary_s1_s2_only) {"
)
s1_guard <- which(
  seq_along(text) < s1_comment &
    text == "if (!figure5_only && !figure4_s4_only) {"
)
s1_start <- s1_comment
full_start <- grep("^if \\(build_all\\) \\{$", text)
stopifnot(
  length(figure2_comment) == 1L, !is.na(figure2_comment),
  length(s1_comment) == 1L, !is.na(s1_comment),
  length(main_start) == 1L, main_start + 1L == figure2_comment,
  length(s1_guard) == 1L, s1_guard + 1L == s1_comment,
  length(s1_start) == 1L, s1_start == s1_comment,
  length(full_start) == 1L
)
stopifnot(main_start < s1_start, s1_start < full_start)
main_text <- paste(text[main_start:(s1_start - 1L)], collapse = "\n")
supp_text <- paste(text[s1_start:(full_start - 1L)], collapse = "\n")
full_text <- paste(text[full_start:length(text)], collapse = "\n")
stopifnot(all(vapply(c("figure2_simplex_source.tsv", "figure3_architecture_source.tsv",
                       "figure4_ssi_ecdf_source.tsv", "figure5_go_source.tsv",
                       "figure5_kegg_source.tsv", "Figure5C_DES_biological_themes_barplot_AUDIT.png",
                       "figure6_baseline_source.tsv", "figure6_noise_source.tsv"),
                     grepl, logical(1), x = main_text, fixed = TRUE)))
stopifnot(!any(vapply(c("Figure5C_DES_biological_themes_barplot_AUDIT.png",
                        "figure2_simplex_source.tsv", "figure3_architecture_source.tsv",
                        "figure4_ssi_ecdf_source.tsv", "figure6_baseline_source.tsv",
                        "figure6_noise_source.tsv", "figureS3_architecture_source.tsv",
                        "figureS4_ssi_distribution_source.tsv", "figure_manifest.tsv", "figure_audit.md"),
                      grepl, logical(1), x = supp_text, fixed = TRUE)))
stopifnot(grepl("figureS3_architecture_source.tsv", full_text, fixed = TRUE))
stopifnot(grepl("figureS4_ssi_distribution_source.tsv", full_text, fixed = TRUE))
metadata_writer_definition <- which(
  trimws(text) ==
    "write_current_figure_metadata <- function(root, figure_root) {"
)
metadata_writer_calls <- which(
  trimws(text) == "write_current_figure_metadata(root, figure_root)"
)
metadata_only_start <- which(trimws(text) == "if (metadata_only) {")
metadata_only_quit <- grep(
  'quit(save = "no", status = 0L)', text, fixed = TRUE
)
annotation_load_start <- which(trimws(text) == "annotation_view <- NULL")
save_figure_definition <- grep(
  "save_figure <- function(", text, fixed = TRUE
)
stopifnot(
  length(metadata_writer_definition) == 1L,
  length(metadata_writer_calls) == 2L,
  length(metadata_only_start) == 1L,
  length(metadata_only_quit) == 1L,
  length(annotation_load_start) == 1L,
  length(save_figure_definition) == 1L
)
metadata_writer_text <- paste(
  text[metadata_writer_definition:(metadata_only_start - 1L)],
  collapse = "\n"
)
stopifnot(grepl(
  'manifest_path <- file.path(figure_root, "figure_manifest.tsv")',
  metadata_writer_text, fixed = TRUE
))
stopifnot(grepl(
  'audit_path <- file.path(figure_root, "figure_audit.md")',
  metadata_writer_text, fixed = TRUE
))
stopifnot(any(
  metadata_writer_calls > metadata_only_start &
    metadata_writer_calls < metadata_only_quit
))
stopifnot(grepl(
  "write_current_figure_metadata(root, figure_root)",
  full_text,
  fixed = TRUE
))
stopifnot(metadata_only_quit < annotation_load_start)
stopifnot(metadata_only_quit < save_figure_definition)
stopifnot(grepl("Supplementary Figures S1 and S2 generated: PASS", joined, fixed = TRUE))
stopifnot(grepl('supplementary_figure_dir <- file.path(figure_root, "supplementary")', joined, fixed = TRUE))
stopifnot(grepl('save_figure(figureS1, "FigureS1", 11.2, 5.6, supplementary_figure_dir)', supp_text, fixed = TRUE))
stopifnot(grepl('save_figure(figureS2, "FigureS2", 10.8, 6.2, supplementary_figure_dir)', supp_text, fixed = TRUE))
stopifnot(!grepl('file.path(figure_root, paste0(id, ".png"))', joined, fixed = TRUE))
stopifnot(!grepl('file.path(figure_root, paste0(id, ".pdf"))', joined, fixed = TRUE))

cat("current supplementary S1/S2 focused tests: PASS\n")
