#!/usr/bin/env Rscript

# Build current evidence-only biological interpretation review tables.

.review_root <- file.path("results", "current_fsf_v1", "manuscript")
.review_out <- file.path(.review_root, "review")

.write_review <- function(x, name, root = .review_out) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, file.path(root, name), na = "NA")
}

.current_theme_rules <- function() {
  add <- function(prefix, type, labels, rules, targets, mode = "manual") {
    data.frame(
      theme_rule_id = sprintf("%s_%02d", prefix, seq_along(labels)),
      rule_definition_source =
        "scripts/current/build_current_biological_review_package.R", rule_type = type, theme_label = labels,
      mapping_or_rule = rules, theme_target = targets,
      manual_or_deterministic = mode,
      current_review_status = "REQUIRES_AUTHOR_REVIEW",
      stringsAsFactors = FALSE
    )
  }
  go_labels <- c(
    "Translation / ribosome", "Protein targeting / secretion",
    "Genome maintenance / DNA repair", "Chromosome / nuclear organization",
    "Cell cycle / cell division", "Development / morphogenesis",
    "Neural / projection organization", "Movement / behavior / taxis",
    "Cell adhesion / junction", "Metabolism / oxidation-reduction",
    "Broad regulatory / cellular process", "Unclassified / manual review"
  )
  go_rules <- c(
    "GO:0022626;GO:0022625;GO:0022627;GO:0005840;GO:0003735;GO:0002181;GO:0006412;GO:0043043",
    "GO:0006612;GO:0006613;GO:0006614;GO:0045047;GO:0072599;GO:0070972",
    "GO:0006259;GO:0006260;GO:0006271;GO:0006273;GO:0006281;GO:0006974;GO:0000723;GO:0000724;GO:0000725;GO:0000726",
    "GO:0005694;GO:0051276;GO:0000228;GO:0044427;GO:0044454;GO:0098813;GO:0005634;GO:0005654",
    "GO:0007049;GO:0022402;GO:0000278;GO:0000280;GO:0048285;GO:0140013;GO:0051321;GO:1903046;GO:1903047;GO:0000070;GO:0000075;GO:0000077",
    "GO:0035239;GO:0048729;GO:0060562;GO:0009887;GO:0048598;GO:0060541;GO:0001501;GO:0002009",
    "GO:0007411;GO:0097485;GO:0031175;GO:0048812;GO:0061564",
    "GO:0040011;GO:0006935;GO:0050920;GO:0007610;GO:0008045;GO:0008038",
    "GO:0007155;GO:0007156;GO:0098609;GO:0005912;GO:0098742",
    "GO:0016491;GO:0055114;GO:0044281;GO:0006082;GO:0019752;GO:0043436;GO:0071704;GO:0044238",
    "regex:^GO:00", "fallback:TRUE"
  )
  kegg_labels <- c(
    "Translation / ribosome", "Genome maintenance / DNA repair",
    "Cell cycle regulation", "Protein turnover / proteostasis",
    "Energy / core metabolism", "Stress signaling / signal transduction",
    "Cell adhesion / extracellular interaction", "Cytoskeleton / cellular remodeling",
    "Neural signaling", "Endocrine / physiological regulation",
    "Broad proliferation / regulatory pathway", "Other KEGG / manual review"
  )
  kegg_rules <- c(
    "ko03010|map03010|Ribosome|K029|K028",
    "ko03030|map03030|ko034|map034|K107|K108|K109|DNA|repair|replication|recombination",
    "ko04110|map04110|ko04111|map04111|ko04113|map04113|cell cycle",
    "ko03050|map03050|ko04120|map04120|proteasome|ubiquitin|K056",
    "ko00190|map00190|ko01100|map01100|ko01110|map01110|ko01120|map01120|ko01200|map01200|oxidative|metabolism",
    "ko04010|map04010|ko04013|map04013|ko04022|map04022|ko04151|map04151|ko046|map046|signaling|MAPK",
    "ko045|map045|adhesion|junction|ECM", "ko048|map048|cytoskeleton|actin",
    "ko047|map047|synapse|neuro|neuron", "ko049|map049|endocrine|hormone|insulin",
    "ko052|map052|cancer", "fallback:TRUE"
  )
  rbind(
    add("GO", "GO_THEME_MAPPING", go_labels, go_rules, go_labels),
    add(
      "KEGG",
      "KEGG_THEME_MAPPING",
      kegg_labels,
      kegg_rules,
      kegg_labels
    )
  )
}

.rule_evidence <- function(rules, terms, ids, type) {
  out <- rules[rules$rule_type == type, , drop = FALSE]
  out$current_matching_term_count <- 0L
  out$current_matching_gene_set_count <- 0L
  out$current_matching_terms <- NA_character_
  out$evidence_classification <- "AMBIGUOUS_REQUIRES_AUTHOR_REVIEW"
  for (i in seq_len(nrow(out))) {
    rule <- out$mapping_or_rule[[i]]
    hit <- if (startsWith(rule, "fallback:")) rep(FALSE, length(terms)) else if (type == "GO_THEME_MAPPING" && !startsWith(rule, "regex:")) {
      terms %in% strsplit(rule, ";", fixed = TRUE)[[1L]]
    } else {
      pattern <- sub("^regex:", "", rule)
      grepl(pattern, terms, ignore.case = TRUE)
    }
    matched <- sort(unique(terms[hit]), method = "radix")
    out$current_matching_term_count[[i]] <- length(matched)
    out$current_matching_gene_set_count[[i]] <- length(unique(ids[hit]))
    out$current_matching_terms[[i]] <- if (length(matched)) paste(matched, collapse = ";") else NA_character_
    if (!startsWith(rule, "fallback:") && length(matched)) {
      out$evidence_classification[[i]] <- "CURRENT_MATCHING_EVIDENCE"
    } else if (!startsWith(rule, "fallback:")) {
      out$evidence_classification[[i]] <- "NO_CURRENT_EVIDENCE"
    }
  }
  out$author_decision <- NA_character_
  out$author_notes <- NA_character_
  out
}

.compact_evidence <- function(x, term, n = 10L) {
  if (!nrow(x)) return(x)
  groups <- split(seq_len(nrow(x)), x$gene_set_id)
  out <- do.call(rbind, lapply(groups, function(i) {
    y <- x[i, , drop = FALSE]
    y <- y[order(y$adjusted_p_value, y$p_value, -y$enrichment_ratio,
                 y[[term]], method = "radix"), , drop = FALSE]
    y <- y[seq_len(min(n, nrow(y))), , drop = FALSE]
    y$review_rank <- seq_len(nrow(y))
    y$contributing_feature_count <- lengths(strsplit(y$contributing_feature_ids, ";", fixed = TRUE))
    y$term_description <- NA_character_
    y$review_status <- "REQUIRES_AUTHOR_REVIEW"
    y
  }))
  rownames(out) <- NULL
  out
}

.complete_evidence <- function(x, main) {
  missing <- main[!main$gene_set_id %in% x$gene_set_id, , drop = FALSE]
  if (!nrow(missing)) return(x)
  blank <- x[rep(NA_integer_, nrow(missing)), , drop = FALSE]
  common <- intersect(names(blank), names(missing))
  blank[common] <- missing[common]
  blank$review_status <- "REQUIRES_AUTHOR_REVIEW"
  blank$term_description <- NA_character_
  out <- rbind(x, blank)
  out[order(out$gene_set_id, out$review_rank, method = "radix", na.last = TRUE), , drop = FALSE]
}

.evidence_summary <- function(x, term, prefix) {
  sets <- unique(x$gene_set_id)
  rows <- lapply(sets, function(id) {
    y <- x[x$gene_set_id == id, , drop = FALSE]
    data.frame(gene_set_id = id,
      summary = paste0(prefix, " significant terms=", nrow(y), "; top IDs=",
                       paste(head(y[[term]], 5L), collapse = ";"),
                       "; minimum adjusted p=", format(min(y$adjusted_p_value), scientific = TRUE)),
      stringsAsFactors = FALSE)
  })
  if (!length(rows)) return(data.frame(gene_set_id = character(), summary = character()))
  do.call(rbind, rows)
}

.architecture_review <- function(region, signal) {
  conditions <- c("DES", "GAM", "HT", "LT", "OSM", "UV")
  do.call(rbind, lapply(conditions, function(cond) {
    x <- region[region$condition == cond, ]
    dominant <- x$stability_region[x$region_count == max(x$region_count)]
    minor <- x$stability_region[x$region_count > 0 & x$region_count < max(x$region_count)]
    s <- signal[signal$condition == cond, ]
    dir <- aggregate(class_count ~ dominant_state, s, sum)
    data.frame(condition = cond,
      dominant_regions = paste(dominant, collapse = ";"),
      minor_regions = if (length(minor)) paste(minor, collapse = ";") else "none",
      transitional_count = x$region_count[x$stability_region == "Transitional"],
      directional_composition = paste0(dir$dominant_state, "=", dir$class_count, collapse = ";"),
      directly_supported = paste0("observed region counts: ", paste0(x$stability_region, "=", x$region_count, collapse = ";")),
      not_supported = "causal mechanism, biological activation, adaptation, or fitness effect",
      review_status = "REQUIRES_AUTHOR_REVIEW", author_decision = NA_character_, author_notes = NA_character_,
      stringsAsFactors = FALSE)
  }))
}

.representative_review <- function(go, kegg) {
  path <- file.path("results", "current_fsf_v1", "representative_features",
                    "current_fsf_top50_annotated_stable_directional_signals_compact.tsv")
  reps <- readr::read_tsv(path, show_col_types = FALSE)
  membership_count <- function(ids, evidence) {
    vapply(ids, function(id) sum(vapply(strsplit(evidence$contributing_feature_ids, ";", fixed = TRUE),
                                        function(z) id %in% z, logical(1L))), integer(1L))
  }
  reps$current_go_supporting_term_count <- membership_count(reps$feature_id, go)
  reps$current_kegg_supporting_term_count <- membership_count(reps$feature_id, kegg)
  reps$current_evidence_status <- ifelse(
    reps$current_go_supporting_term_count + reps$current_kegg_supporting_term_count > 0,
    "RETAIN", "RETAIN_WITH_REVISED_RATIONALE")
  reps$review_status <- "REQUIRES_AUTHOR_REVIEW"
  reps$author_decision <- NA_character_
  reps$author_notes <- NA_character_
  reps
}

build_current_biological_review_package <- function(root = .review_out) {
  if (!requireNamespace("readr", quietly = TRUE)) stop("readr is required.")
  manifest <- readr::read_tsv(file.path(.review_root, "tables/current_enrichment_gene_set_summary.tsv"), show_col_types = FALSE)
  main <- manifest[manifest$gene_set_family == "main_class", ]
  go_all <- readr::read_tsv(file.path(.review_root, "enrichment/go/go_enrichment_significant.tsv"), show_col_types = FALSE)
  kegg_all <- readr::read_tsv(file.path(.review_root, "enrichment/kegg/kegg_enrichment_significant.tsv"), show_col_types = FALSE)
  go <- go_all[go_all$gene_set_id %in% main$gene_set_id, ]
  kegg <- kegg_all[kegg_all$gene_set_id %in% main$gene_set_id, ]
  go_review <- .complete_evidence(.compact_evidence(go, "go_term", 10L), main)
  kegg_review <- .complete_evidence(.compact_evidence(kegg, "kegg_term", 10L), main)
  rules <- .current_theme_rules()
  go_rules <- .rule_evidence(rules, go$go_term, go$gene_set_id, "GO_THEME_MAPPING")
  kegg_rules <- .rule_evidence(rules, kegg$kegg_term, kegg$gene_set_id, "KEGG_THEME_MAPPING")

  go_sum <- .evidence_summary(go, "go_term", "GO")
  kg_sum <- .evidence_summary(kegg, "kegg_term", "KEGG")
  names(go_sum)[2L] <- "current_GO_evidence_summary"
  names(kg_sum)[2L] <- "current_KEGG_evidence_summary"
  profile <- main[c("condition", "stability_region", "dominant_state", "signal_class", "gene_set_id")]
  profile$region_or_class <- ifelse(is.na(profile$signal_class), profile$stability_region, profile$signal_class)
  profile <- merge(profile, go_sum, by = "gene_set_id", all.x = TRUE, sort = FALSE)
  profile <- merge(profile, kg_sum, by = "gene_set_id", all.x = TRUE, sort = FALSE)
  profile$current_GO_evidence_summary[is.na(profile$current_GO_evidence_summary)] <- "GO significant terms=0"
  profile$current_KEGG_evidence_summary[is.na(profile$current_KEGG_evidence_summary)] <- "KEGG significant terms=0"
  profile$architecture_context <- paste(profile$condition, profile$stability_region, profile$dominant_state, sep = "/")
  go_n <- table(go$gene_set_id)[profile$gene_set_id]; go_n[is.na(go_n)] <- 0
  kg_n <- table(kegg$gene_set_id)[profile$gene_set_id]; kg_n[is.na(kg_n)] <- 0
  profile$current_GO_significant_count <- as.integer(go_n)
  profile$current_KEGG_significant_count <- as.integer(kg_n)
  profile$review_status <- "REQUIRES_AUTHOR_REVIEW"
  profile$author_decision <- NA_character_
  profile$author_notes <- NA_character_

  region <- readr::read_tsv(file.path(.review_root, "tables/current_region_architecture.tsv"), show_col_types = FALSE)
  signal <- readr::read_tsv(file.path(.review_root, "tables/current_signal_class_architecture.tsv"), show_col_types = FALSE)
  architecture <- .architecture_review(region, signal)
  reps <- .representative_review(go_all, kegg_all)

  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  .write_review(go_review, "current_go_evidence.tsv", root)
  .write_review(kegg_review, "current_kegg_evidence.tsv", root)
  .write_review(go_rules, "go_theme_review.tsv", root)
  .write_review(kegg_rules, "kegg_theme_review.tsv", root)
  .write_review(profile, "combined_profile_review.tsv", root)
  .write_review(profile[profile$condition == "UV", ], "uv_specific_review.tsv", root)
  .write_review(architecture, "condition_architecture_review.tsv", root)
  .write_review(reps, "representative_feature_review.tsv", root)

  summary <- c(
    "# FSF biological interpretation review package",
    "Current FSF evidence views and theme-mapping candidates for author review.",
    paste0("- Current main-class gene sets reviewed: ", nrow(main)),
    paste0("- GO mappings with current matching evidence: ",
      sum(go_rules$current_matching_term_count > 0)),
    paste0("- KEGG mappings with current matching evidence: ",
      sum(kegg_rules$current_matching_term_count > 0)),
    "## Current mapping evidence",
    paste0("- GO: ", paste(go_rules$theme_label[
      go_rules$evidence_classification == "CURRENT_MATCHING_EVIDENCE"
    ], collapse = "; ")),
    paste0("- KEGG: ", paste(kegg_rules$theme_label[
      kegg_rules$evidence_classification == "CURRENT_MATCHING_EVIDENCE"
    ], collapse = "; ")),
    paste0("- KEGG without current matches: ", paste(kegg_rules$theme_label[
      kegg_rules$evidence_classification == "NO_CURRENT_EVIDENCE"
    ], collapse = "; ")),
    "Mapping evidence does not constitute author acceptance of a biological theme.",
    "DES, GAM, and LT contain Low Stability rows in the current region authority.",
    "HT and OSM are entirely Highly Stable in the current region authority.",
    "UV current biological profiles require author review.",
    paste0("Current condition architecture -> current signal class -> ",
      "GO/KEGG statistical evidence -> author-reviewed theme/profile."),
    "Author decisions and notes are intentionally blank."
  )
  writeLines(summary, file.path(root, "review_summary.md"), useBytes = TRUE)
  invisible(list(rules = rules, go = go_review, kegg = kegg_review,
                 profiles = profile, architecture = architecture, representatives = reps))
}

audit_current_biological_review_package <- function(
    root = .review_out,
    review_state = c("template", "final")) {
  review_state <- match.arg(review_state)
  required <- c(
    "current_go_evidence.tsv",
    "current_kegg_evidence.tsv",
    "go_theme_review.tsv",
    "kegg_theme_review.tsv",
    "combined_profile_review.tsv",
    "uv_specific_review.tsv",
    "condition_architecture_review.tsv",
    "representative_feature_review.tsv",
    "review_summary.md"
  )
  if (!all(file.exists(file.path(root, required)))) {
    stop("Review package is incomplete.")
  }
  go_themes <- readr::read_tsv(
    file.path(root, "go_theme_review.tsv"),
    show_col_types = FALSE
  )
  kegg_themes <- readr::read_tsv(
    file.path(root, "kegg_theme_review.tsv"),
    show_col_types = FALSE
  )
  go_evidence <- readr::read_tsv(
    file.path(root, "current_go_evidence.tsv"),
    show_col_types = FALSE
  )
  kegg_evidence <- readr::read_tsv(
    file.path(root, "current_kegg_evidence.tsv"),
    show_col_types = FALSE
  )
  profile <- readr::read_tsv(file.path(root, "combined_profile_review.tsv"), show_col_types = FALSE)
  uv <- readr::read_tsv(file.path(root, "uv_specific_review.tsv"), show_col_types = FALSE)
  architecture <- readr::read_tsv(file.path(root, "condition_architecture_review.tsv"), show_col_types = FALSE)
  reps <- readr::read_tsv(file.path(root, "representative_feature_review.tsv"), show_col_types = FALSE)
  theme_schema <- c(
    "theme_rule_id",
    "rule_definition_source",
    "theme_label",
    "theme_target"
  )

  if (!all(theme_schema %in% names(go_themes)) ||
      !all(theme_schema %in% names(kegg_themes)) ||
      !("gene_set_id" %in% names(profile)) ||
      !("condition" %in% names(architecture))) {
    stop("Review schema coverage failed.")
  }

  if (anyDuplicated(go_themes$theme_rule_id) ||
      anyDuplicated(kegg_themes$theme_rule_id) ||
      nrow(profile) != 40L ||
      anyDuplicated(profile$gene_set_id) ||
      length(unique(go_evidence$gene_set_id)) != 40L ||
      length(unique(kegg_evidence$gene_set_id)) != 40L ||
      nrow(uv) != 13L ||
      nrow(architecture) != 6L ||
      anyDuplicated(architecture$condition)) {
    stop("Review key coverage failed.")
  }

  current_tables <- list(
    go_themes,
    kegg_themes,
    profile,
    architecture,
    reps
  )

  has_historical_fields <- vapply(
    current_tables,
    function(x) any(startsWith(names(x), "historical_")),
    logical(1L)
  )

  evidence_change_fields <- c(
    "GO_evidence_change",
    "KEGG_evidence_change",
    "theme_evidence_change"
  )

  if (any(has_historical_fields) ||
      any(evidence_change_fields %in% names(profile))) {
    stop("Current review tables contain obsolete comparison fields.")
  }

  for (x in list(profile, reps)) {
    if (!all(c("author_decision", "author_notes") %in% names(x)) ||
        any(!is.na(x$author_decision)) || any(!is.na(x$author_notes))) {
      stop("Profile and representative-feature author fields must remain blank.")
    }
  }
  if (!all(c("author_decision", "author_notes") %in% names(architecture))) {
    stop("Architecture author-field schema is incomplete.")
  }
  if (review_state == "template") {
    if (any(!is.na(architecture$author_decision)) ||
        any(!is.na(architecture$author_notes))) {
      stop("Template architecture author fields must remain blank.")
    }
  } else {
    expected_conditions <- c("DES", "GAM", "HT", "LT", "OSM", "UV")
    if (nrow(architecture) != 6L ||
        anyDuplicated(architecture$condition) ||
        !setequal(architecture$condition, expected_conditions) ||
        any(is.na(architecture$author_decision)) ||
        any(architecture$author_decision != "RETAIN_DESCRIPTIVE_ONLY") ||
        any(is.na(architecture$author_notes)) ||
        any(!nzchar(trimws(architecture$author_notes)))) {
      stop("Final architecture author decisions failed validation.")
    }
  }
  predicted <- unlist(profile[c("stability_region", "signal_class")])
  if (any(grepl("Instability|instable|weakly_stable|weakly stable", predicted, ignore.case = TRUE))) {
    stop("Obsolete current terminology detected.")
  }
  script <- readLines("scripts/current/build_current_biological_review_package.R", warn = FALSE)
  if (any(grepl("SSI *(>=|>|<) *0\\.60", script))) stop("Obsolete executable boundary detected.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  build_current_biological_review_package()
  audit_current_biological_review_package(review_state = "template")
  cat("current biological review package: PASS\n")
}
