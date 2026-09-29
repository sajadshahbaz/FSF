# Current FSF v1 GO overrepresentation engine.
#
# The workflow consumes the frozen current gene-set builder, parses the
# authoritative annotation fields, and returns deterministic in-memory
# mappings, diagnostics, and results. It writes nothing.

source(file.path("scripts", "current", "build_current_enrichment_gene_sets.R"))

.fsf_go_result_columns <- c(
  "gene_set_id", "gene_set_family", "condition", "stability_region",
  "dominant_state", "signal_class", "go_term", "query_size",
  "background_size", "query_overlap_count", "background_term_count",
  "query_fraction", "background_fraction", "enrichment_ratio", "p_value",
  "adjusted_p_value", "significant", "contributing_feature_ids",
  "all_feature_count", "annotated_feature_count"
)

.fsf_go_empty_results <- function() {
  out <- data.frame(
    gene_set_id = character(), gene_set_family = character(),
    condition = character(), stability_region = character(),
    dominant_state = character(), signal_class = character(),
    go_term = character(), query_size = integer(), background_size = integer(),
    query_overlap_count = integer(), background_term_count = integer(),
    query_fraction = numeric(), background_fraction = numeric(),
    enrichment_ratio = numeric(), p_value = numeric(),
    adjusted_p_value = numeric(), significant = logical(),
    contributing_feature_ids = character(), all_feature_count = integer(),
    annotated_feature_count = integer(), stringsAsFactors = FALSE
  )
  out[.fsf_go_result_columns]
}

.fsf_go_extract_source <- function(feature_id, values, separator, source) {
  rows <- lapply(seq_along(values), function(i) {
    value <- values[[i]]
    if (is.na(value) || !nzchar(value) || identical(value, "-")) return(NULL)
    tokens <- strsplit(as.character(value), separator, fixed = TRUE)[[1L]]
    starts <- regexpr("GO:[0-9]{7}", tokens)
    valid <- starts > 0L
    if (!any(valid)) return(NULL)
    terms <- substring(tokens[valid], starts[valid], starts[valid] + 9L)
    data.frame(
      feature_id = rep(as.character(feature_id[[i]]), length(terms)),
      go_term = terms,
      source = rep(source, length(terms)),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      feature_id = character(), go_term = character(), source = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  out <- unique(out[c("feature_id", "go_term", "source")])
  out[order(out$feature_id, out$go_term, out$source, method = "radix"), , drop = FALSE]
}

.fsf_go_build_mappings <- function(annotation) {
  required <- c("feature_id", "eggnog_go", "interpro_go")
  if (!is.data.frame(annotation) || !all(required %in% names(annotation))) {
    stop("GO annotation input is missing required fields.", call. = FALSE)
  }
  annotation <- annotation[!duplicated(annotation$feature_id), required, drop = FALSE]
  eggnog <- .fsf_go_extract_source(
    annotation$feature_id, annotation$eggnog_go, ",", "eggNOG"
  )
  interpro <- .fsf_go_extract_source(
    annotation$feature_id, annotation$interpro_go, ";", "InterPro"
  )
  source_specific <- rbind(eggnog, interpro)
  source_specific <- unique(source_specific[c("feature_id", "go_term", "source")])
  source_specific <- source_specific[order(
    source_specific$feature_id, source_specific$go_term,
    source_specific$source, method = "radix"
  ), , drop = FALSE]
  rownames(source_specific) <- NULL
  combined <- unique(source_specific[c("feature_id", "go_term")])
  combined <- combined[order(combined$feature_id, combined$go_term, method = "radix"), , drop = FALSE]
  rownames(combined) <- NULL
  list(source_specific = source_specific, combined = combined)
}

.fsf_go_run_one <- function(manifest_row, background, mapping, background_counts) {
  query <- sort(unique(intersect(
    as.character(manifest_row$annotated_features[[1L]]), background
  )), method = "radix")
  query_size <- length(query)
  if (query_size < 5L) return(.fsf_go_empty_results())

  query_mapping <- mapping[mapping$feature_id %in% query, , drop = FALSE]
  if (!nrow(query_mapping)) return(.fsf_go_empty_results())
  term_features <- split(query_mapping$feature_id, query_mapping$go_term)
  terms <- sort(names(term_features), method = "radix")
  overlaps <- lengths(term_features[terms])
  keep <- overlaps > 0L
  terms <- terms[keep]
  overlaps <- overlaps[keep]
  if (!length(terms)) return(.fsf_go_empty_results())

  background_size <- length(background)
  background_term_count <- as.integer(background_counts[terms])
  if (anyNA(background_term_count) || any(background_term_count <= 0L) ||
      any(background_term_count > background_size)) {
    stop("GO term counts contain invalid hypergeometric inputs.", call. = FALSE)
  }
  p_value <- stats::phyper(
    q = overlaps - 1L,
    m = background_term_count,
    n = background_size - background_term_count,
    k = query_size,
    lower.tail = FALSE
  )
  if (any(!is.finite(p_value)) || any(p_value < 0 | p_value > 1)) {
    stop("GO hypergeometric calculation returned invalid probabilities.", call. = FALSE)
  }
  adjusted <- stats::p.adjust(p_value, method = "BH")
  query_fraction <- overlaps / query_size
  background_fraction <- background_term_count / background_size
  contributing <- vapply(term_features[terms], function(ids) {
    paste(sort(unique(ids), method = "radix"), collapse = ";")
  }, character(1L))

  out <- data.frame(
    gene_set_id = manifest_row$gene_set_id,
    gene_set_family = manifest_row$gene_set_family,
    condition = manifest_row$condition,
    stability_region = manifest_row$stability_region,
    dominant_state = manifest_row$dominant_state,
    signal_class = manifest_row$signal_class,
    go_term = terms,
    query_size = query_size,
    background_size = background_size,
    query_overlap_count = as.integer(overlaps),
    background_term_count = background_term_count,
    query_fraction = query_fraction,
    background_fraction = background_fraction,
    enrichment_ratio = query_fraction / background_fraction,
    p_value = p_value,
    adjusted_p_value = adjusted,
    significant = adjusted <= 0.05,
    contributing_feature_ids = contributing,
    all_feature_count = manifest_row$all_feature_count,
    annotated_feature_count = manifest_row$annotated_feature_count,
    stringsAsFactors = FALSE
  )
  out[order(
    out$adjusted_p_value, out$p_value, -out$enrichment_ratio,
    out$go_term, method = "radix"
  ), .fsf_go_result_columns, drop = FALSE]
}

.run_current_go_enrichment <- function(gene_sets, annotation, validate_real = TRUE) {
  required_objects <- c("current_gene_set_manifest", "universes")
  if (!is.list(gene_sets) || !all(required_objects %in% names(gene_sets))) {
    stop("Current gene-set builder output is invalid.", call. = FALSE)
  }
  manifest <- gene_sets$current_gene_set_manifest
  required_manifest <- c(
    "gene_set_id", "gene_set_family", "condition", "stability_region",
    "dominant_state", "signal_class", "all_feature_count",
    "annotated_feature_count", "annotated_features"
  )
  if (!is.data.frame(manifest) || !all(required_manifest %in% names(manifest)) ||
      anyDuplicated(manifest$gene_set_id)) {
    stop("Current gene-set manifest is invalid.", call. = FALSE)
  }
  background <- gene_sets$universes$annotation_universe
  if (is.null(background) || anyNA(background) || anyDuplicated(background) ||
      !identical(background, sort(background, method = "radix"))) {
    stop("GO annotation_universe is invalid.", call. = FALSE)
  }

  mappings <- .fsf_go_build_mappings(annotation)
  mapping <- mappings$combined[
    mappings$combined$feature_id %in% background, , drop = FALSE
  ]
  background_counts <- table(mapping$go_term)

  results <- vector("list", nrow(manifest))
  status <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    query <- sort(unique(intersect(
      as.character(row$annotated_features[[1L]]), background
    )), method = "radix")
    result <- .fsf_go_run_one(row, background, mapping, background_counts)
    results[[i]] <- result
    status[[i]] <- data.frame(
      gene_set_id = row$gene_set_id,
      query_size = length(query),
      status = if (length(query) < 5L) {
        "SKIPPED_QUERY_SIZE_LT_5"
      } else if (!nrow(result)) {
        "EMPTY_NO_MAPPED_GO_TERM"
      } else {
        "TESTED"
      },
      result_rows = nrow(result),
      stringsAsFactors = FALSE
    )
  }
  combined <- do.call(rbind, results)
  if (is.null(combined)) combined <- .fsf_go_empty_results()
  rownames(combined) <- NULL
  if (nrow(combined)) {
    combined <- combined[order(
      combined$gene_set_id, combined$adjusted_p_value, combined$p_value,
      -combined$enrichment_ratio, combined$go_term, method = "radix"
    ), , drop = FALSE]
    rownames(combined) <- NULL
  }
  status <- do.call(rbind, status)
  status <- status[order(status$gene_set_id, method = "radix"), , drop = FALSE]
  rownames(status) <- NULL

  if (validate_real) {
    observed <- c(
      background = length(background),
      mapped_features = length(unique(mappings$combined$feature_id)),
      unique_terms = length(unique(mappings$combined$go_term)),
      combined_pairs = nrow(mappings$combined),
      source_pairs = nrow(mappings$source_specific),
      manifest_rows = nrow(manifest)
    )
    expected <- c(
      background = 13467L, mapped_features = 9321L,
      unique_terms = 20502L, combined_pairs = 1131086L,
      source_pairs = 1159053L, manifest_rows = 70L
    )
    if (!identical(as.integer(observed), as.integer(expected))) {
      stop("Current GO mapping or background counts failed validation.", call. = FALSE)
    }
    if (!nrow(combined)) {
      stop("Current real-data GO enrichment produced no results.", call. = FALSE)
    }
  }

  list(
    source_specific_mapping = mappings$source_specific,
    combined_mapping = mappings$combined,
    background = background,
    gene_set_status = status,
    results = combined
  )
}

#' Run current FSF GO enrichment in memory
#'
#' Uses the frozen current gene-set builder and hash-validated annotation
#' authority. Results and mappings are returned; nothing is written.
run_current_go_enrichment <- function(annotation_master_path) {
  gene_sets <- build_current_enrichment_gene_sets(annotation_master_path)
  joined <- load_current_fsf_annotation(annotation_master_path)
  annotation <- joined[!duplicated(joined$feature_id),
    c("feature_id", "eggnog_go", "interpro_go"), drop = FALSE]
  .run_current_go_enrichment(gene_sets, annotation, validate_real = TRUE)
}
