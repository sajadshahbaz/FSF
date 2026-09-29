#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
packet_path <- file.path(
  root, "results", "current_fsf_v1", "manuscript",
  "biological_theme_decisions", "author_biological_theme_decision_packet.tsv"
)
output_dir <- file.path(
  root, "results", "current_fsf_v1", "manuscript", "figures",
  "biological_theme_prototypes"
)


blank <- function(x) all(is.na(x) | !nzchar(x))

parse_term_count <- function(x, source_name) {
  matched <- grepl("^[0-9]+ terms(?: / [0-9]+ redundancy groups)?$", x)
  if (!all(matched)) stop("Could not parse frozen ", source_name, " support count")
  as.integer(sub(" terms.*$", "", x))
}

strongest_valid_p <- function(go_p, kegg_p) {
  values <- c(go_p, kegg_p)
  values <- values[is.finite(values) & values > 0 & values <= 1]
  if (!length(values)) stop("Approved theme has no valid supporting adjusted P-value")
  min(values)
}

format_p_plotmath <- function(p) {
  exponent <- floor(log10(p))
  coefficient <- p / (10^exponent)
  sprintf("%.2g%%*%%10^{%d}", coefficient, exponent)
}

approved_mappings <- list(
  DES = data.frame(
    FSF_region = c(rep("Low Stability", 3), rep("Highly Stable", 8)),
    candidate_theme = c(
      "Cell adhesion and extracellular organization",
      "Development and morphogenesis",
      "Neural, projection and behavior processes",
      "Translation, ribosome and protein targeting",
      "RNA processing and expression",
      "Vesicle trafficking and organelle organization",
      "Cell cycle and division",
      "Chromosome and nuclear organization",
      "DNA repair and genome maintenance",
      "Energy and core metabolism",
      "Protein turnover and proteostasis"
    ),
    display_theme = c(
      "Cell adhesion & extracellular matrix",
      "Developmental & morphogenetic processes",
      "Neuronal projection & guidance",
      "Translation & ribosome",
      "RNA processing",
      "ER targeting & vesicle trafficking",
      "Cell cycle & division",
      "Chromosome & nuclear organization",
      "DNA repair & genome maintenance",
      "Oxidative phosphorylation & energy metabolism",
      "Protein turnover & proteostasis"
    )
  )
)

excluded_des_themes <- c(
  "Broad regulatory and signaling processes",
  "Cytoskeleton and cell motility",
  "Immune and defense processes",
  "Metabolism and redox processes",
  "Membrane transport and ion homeostasis",
  "Endocrine and physiological regulation",
  "Stress response and signaling"
)

validate_frozen_state <- function() {
  manifest_path <- file.path(root,"results","current_fsf_v1","manuscript","audit","current_release_protected_files.sha256.tsv")
  if (!file.exists(manifest_path)) stop("Current-release protected-file manifest is missing: ",manifest_path)
  frozen <- read.delim(manifest_path,sep="	",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  if (!identical(names(frozen),c("sha256","relative_path"))) stop("Unexpected protected-file manifest schema")
  if (!nrow(frozen)) stop("Protected-file manifest is empty")
  if (anyDuplicated(frozen$relative_path)) stop("Protected-file manifest contains duplicate paths")
  sha256_file <- function(path) { if (!file.exists(path)) stop("Protected file is missing: ",path); value <- system2("sha256sum",path,stdout=TRUE,stderr=TRUE); status <- attr(value,"status"); if (!is.null(status) && status != 0L) stop("sha256sum failed for protected file: ",path); sub(" .*","",value[1]) }
  actual <- vapply(frozen$relative_path,function(rel) sha256_file(file.path(root,rel)),character(1)); mismatch <- which(unname(actual) != unname(frozen$sha256))
  if (length(mismatch)) stop("Current-release protected-file hash changed: ",paste(frozen$relative_path[mismatch],collapse="; "))
  protected_hashes <- c("results/current_fsf_v1/manuscript/figures/new/Figure_1.png"="9e3796b99bc0fd1b36d8c5029ee286ebfe76c9bf018d6d763dab596fd83be0f2","results/current_fsf_v1/manuscript/figures/main/Figure1.png"="9e3796b99bc0fd1b36d8c5029ee286ebfe76c9bf018d6d763dab596fd83be0f2","results/current_fsf_v1/manuscript/figures/new/Figure_7.png"="6ab8be9958a8073303661814fa7def0de9e1a5d9ece745adb23dc1bf377771a0","results/current_fsf_v1/manuscript/figures/main/Figure7.png"="6ab8be9958a8073303661814fa7def0de9e1a5d9ece745adb23dc1bf377771a0")
  actual_figures <- vapply(names(protected_hashes),function(rel) sha256_file(file.path(root,rel)),character(1))
  if (!identical(unname(actual_figures),unname(protected_hashes))) stop("Locked AI-only Figure 1 or Figure 7 hash changed")
  message("Current-release protected-file and figure integrity: PASS")
}

build_plot_data <- function(condition, mapping, packet) {
  if (anyDuplicated(mapping[c("FSF_region", "candidate_theme")])) {
    stop("Approved mapping contains duplicate region/theme keys")
  }
  observed <- packet[packet$condition == condition, ]
  keys <- paste(mapping$FSF_region, mapping$candidate_theme, sep = "\r")
  observed_keys <- paste(observed$FSF_region, observed$candidate_theme, sep = "\r")
  missing <- mapping$candidate_theme[!keys %in% observed_keys]
  if (length(missing)) {
    stop("Exact candidate_theme match failed: ", paste(missing, collapse = "; "))
  }
  selected <- observed[match(keys, observed_keys), ]
  selected$display_theme <- mapping$display_theme
  selected$GO_significant_term_count <- parse_term_count(selected$GO_support, "GO")
  selected$KEGG_significant_term_count <- parse_term_count(selected$KEGG_support, "KEGG")
  selected$total_supporting_term_count <-
    selected$GO_significant_term_count + selected$KEGG_significant_term_count
  selected$strongest_adjusted_p <- mapply(
    strongest_valid_p,
    selected$strongest_GO_adjusted_p,
    selected$strongest_KEGG_adjusted_p
  )
  selected$neg_log10_strongest_adjusted_p <- -log10(selected$strongest_adjusted_p)

  region_order <- c("Low Stability", "Highly Stable")
  selected$FSF_region <- factor(selected$FSF_region, region_order)
  selected <- selected[order(
    selected$FSF_region,
    -selected$total_supporting_term_count,
    -selected$neg_log10_strongest_adjusted_p,
    selected$display_theme
  ), ]
  selected$plot_order <- ave(
    seq_len(nrow(selected)), selected$FSF_region,
    FUN = function(z) seq_along(z)
  )
  selected$FSF_region <- as.character(selected$FSF_region)

  keep <- c(
    "condition", "FSF_region", "candidate_theme", "display_theme",
    "GO_significant_term_count", "KEGG_significant_term_count",
    "total_supporting_term_count", "strongest_GO_adjusted_p",
    "strongest_KEGG_adjusted_p", "strongest_adjusted_p",
    "neg_log10_strongest_adjusted_p", "plot_order"
  )
  selected[keep]
}

validate_plot_data <- function(x) {
  expected_counts <- c("Low Stability" = 3L, "Highly Stable" = 8L)
  actual_counts <- table(factor(x$FSF_region, names(expected_counts)))
  if (!identical(as.integer(actual_counts), unname(expected_counts)) || nrow(x) != 11L) {
    stop("DES plotting data must contain exactly 3 Low Stability and 8 Highly Stable rows")
  }
  if (any(x$candidate_theme %in% excluded_des_themes)) {
    stop("Excluded DES theme appears in plotting data")
  }
  if (any(!is.finite(x$total_supporting_term_count)) ||
      any(!is.finite(x$strongest_adjusted_p)) ||
      any(!is.finite(x$neg_log10_strongest_adjusted_p))) {
    stop("Non-finite plot value detected")
  }
}

make_panel <- function(x, audit = FALSE) {
  x$FSF_region <- factor(x$FSF_region, c("Low Stability", "Highly Stable"))
  x$plot_key <- paste(x$FSF_region, x$display_theme, sep = "\r")
  x$plot_key <- factor(x$plot_key, levels = rev(unique(x$plot_key)))
  label_lookup <- setNames(x$display_theme, as.character(x$plot_key))
  if (audit) {
    x$bar_label <- sprintf(
      "paste('GO=%d | KEGG=%d | total=%d | adj.P=',%s)",
      x$GO_significant_term_count, x$KEGG_significant_term_count,
      x$total_supporting_term_count,
      vapply(x$strongest_adjusted_p, format_p_plotmath, character(1))
    )
  } else {
    x$bar_label <- sprintf(
      "paste('n=%d | adj.P=',%s)", x$total_supporting_term_count,
      vapply(x$strongest_adjusted_p, format_p_plotmath, character(1))
    )
  }
  label_offset <- max(x$total_supporting_term_count) * 0.018
  subtitle <- if (audit) {
    "AUDIT - author-reviewed themes with source-specific supporting-term counts"
  } else {
    "Author-reviewed biological themes supported by functional enrichment"
  }

  ggplot(x, aes(x = total_supporting_term_count, y = plot_key,
                fill = neg_log10_strongest_adjusted_p)) +
    geom_col(width = 0.68) +
    geom_text(
      aes(x = total_supporting_term_count + label_offset, label = bar_label),
      parse = TRUE, hjust = 0, size = if (audit) 3.5 else 3.4,
      color = "#222222", family = "sans"
    ) +
    facet_grid(FSF_region ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = label_lookup) +
    scale_x_continuous(
      expand = expansion(mult = c(0, if (audit) 0.82 else 0.60))
    ) +
    scale_fill_gradient(
      low = "#D8E8F1", high = "#174A6E",
      guide = guide_colourbar(display = "rectangles"),
      name = expression(atop("Evidence strength", -log[10]("strongest adjusted P")))
    ) +
    labs(
      title = "DES: biological differentiation across FSF regions",
      subtitle = subtitle,
      x = "Significant supporting GO + KEGG terms", y = NULL,
      caption = paste(
        "adj.P is the strongest adjusted P-value among significant GO or KEGG terms",
        "supporting each author-reviewed theme; it is not a theme-level P-value."
      )
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 10.5, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, color = "#444444"),
      plot.caption = element_text(size = 9.5, color = "#444444", hjust = 0),
      axis.title.x = element_text(face = "bold", size = 12.5, margin = margin(t = 8)),
      axis.text.y = element_text(size = 10.5, color = "#222222"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#E8E8E8", linewidth = 0.3),
      strip.background = element_rect(fill = "#F1F3F5", color = NA),
      strip.text.y = element_text(face = "bold", size = 11.5, angle = 0),
      strip.placement = "outside",
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.height = grid::unit(0.45, "in"),
      panel.spacing.y = grid::unit(0.55, "lines"),
      plot.margin = margin(12, 30, 12, 12)
    )
}

save_panel <- function(plot, stem, width, height) {
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  ggsave(
    pdf_path, plot, width = width, height = height, units = "in",
    device = grDevices::pdf, bg = "white", useDingbats = FALSE
  )
  png_prefix <- sub("[.]png$", "", png_path)
  status <- system2(
    "pdftocairo", c("-png", "-singlefile", "-r", "600", pdf_path, png_prefix)
  )
  if (status != 0L || !file.exists(png_path)) {
    stop("600 dpi PDF-to-PNG conversion failed for ", stem)
  }
  c(png_path, pdf_path)
}

validate_frozen_state()
if (!file.exists(packet_path)) stop("Missing frozen packet: ", packet_path)
packet <- read.delim(
  packet_path, check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
if (!blank(packet$author_decision) || !blank(packet$author_notes)) {
  stop("Frozen author_decision or author_notes fields are not blank")
}

plot_data <- build_plot_data("DES", approved_mappings$DES, packet)
validate_plot_data(plot_data)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source_path <- file.path(output_dir, "DES_biological_theme_barplot_source.tsv")
write.table(
  plot_data, source_path, sep = "\t", quote = FALSE,
  row.names = FALSE, na = ""
)

clean_plot <- make_panel(plot_data, audit = FALSE)
audit_plot <- make_panel(plot_data, audit = TRUE)
save_panel(clean_plot, "Figure5C_DES_biological_themes_barplot", 12.5, 8.2)
save_panel(audit_plot, "Figure5C_DES_biological_themes_barplot_AUDIT", 14.0, 8.5)

caption <- paste(
  "Author-reviewed biological themes associated with DES FSF regions.",
  "Bar length represents the number of significant GO and KEGG terms supporting",
  "each biological theme. Bar color represents enrichment evidence strength,",
  "expressed as -log10 of the strongest adjusted P-value among the supporting GO or",
  "KEGG terms. Exact supporting-term counts (n) and corresponding strongest adjusted",
  "P-values are shown beside each bar. The displayed adjusted P-value is the strongest",
  "value among supporting enrichment terms, not a P-value from an independent test of",
  "the biological theme. Biological themes are displayed separately according to FSF",
  "region."
)
writeLines(
  caption,
  file.path(output_dir, "Figure5C_DES_biological_themes_caption.txt"),
  useBytes = TRUE
)

message("DES author-reviewed biological theme panels: PASS")
message("Rows: Low Stability 3; Highly Stable 8; total 11")
