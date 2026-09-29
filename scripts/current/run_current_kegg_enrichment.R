# Current FSF v1 KEGG overrepresentation engine.
#
# The workflow consumes frozen current gene sets and authoritative eggNOG KO
# and pathway fields. KO and PATHWAY tests use separate backgrounds and
# multiple-testing families. Everything is returned in memory; nothing is
# written and no external database is queried.

source(file.path("scripts", "current", "build_current_enrichment_gene_sets.R"))

.fsf_kegg_result_columns <- c(
  "gene_set_id", "gene_set_family", "condition", "stability_region",
  "dominant_state", "signal_class", "enrichment_type", "kegg_term",
  "query_size", "background_size", "query_overlap_count",
  "background_term_count", "query_fraction", "background_fraction",
  "enrichment_ratio", "p_value", "adjusted_p_value", "significant",
  "contributing_feature_ids", "all_feature_count", "annotated_feature_count"
)

.fsf_kegg_empty_results <- function() {
  out <- data.frame(
    gene_set_id = character(), gene_set_family = character(),
    condition = character(), stability_region = character(),
    dominant_state = character(), signal_class = character(),
    enrichment_type = character(), kegg_term = character(),
    query_size = integer(), background_size = integer(),
    query_overlap_count = integer(), background_term_count = integer(),
    query_fraction = numeric(), background_fraction = numeric(),
    enrichment_ratio = numeric(), p_value = numeric(),
    adjusted_p_value = numeric(), significant = logical(),
    contributing_feature_ids = character(), all_feature_count = integer(),
    annotated_feature_count = integer(), stringsAsFactors = FALSE
  )
  out[.fsf_kegg_result_columns]
}

.fsf_kegg_extract_mapping <- function(feature_id, values, enrichment_type) {
  rows <- lapply(seq_along(values), function(i) {
    value <- values[[i]]
    if (is.na(value) || !nzchar(trimws(value)) || trimws(value) == "-") {
      return(NULL)
    }
    tokens <- trimws(unlist(strsplit(as.character(value), "[,;]")))
    if (enrichment_type == "KO") tokens <- sub("^ko:", "", tokens)
    tokens <- tokens[nzchar(tokens) & tokens != "-"]
    if (!length(tokens)) return(NULL)
    data.frame(
      feature_id = rep(as.character(feature_id[[i]]), length(tokens)),
      kegg_term = tokens,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(feature_id = character(), kegg_term = character()))
  }
  out <- unique(do.call(rbind, rows)[c("feature_id", "kegg_term")])
  out <- out[order(out$feature_id, out$kegg_term, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.fsf_kegg_build_mappings <- function(annotation) {
  required <- c("feature_id", "eggnog_kegg_ko", "eggnog_kegg_pathway")
  if (!is.data.frame(annotation) || !all(required %in% names(annotation))) {
    stop("KEGG annotation input is missing required fields.", call. = FALSE)
  }
  annotation <- annotation[!duplicated(annotation$feature_id), required, drop = FALSE]
  list(
    KO = .fsf_kegg_extract_mapping(
      annotation$feature_id, annotation$eggnog_kegg_ko, "KO"
    ),
    PATHWAY = .fsf_kegg_extract_mapping(
      annotation$feature_id, annotation$eggnog_kegg_pathway, "PATHWAY"
    )
  )
}

.fsf_kegg_run_one <- function(
    manifest_row, enrichment_type, background, mapping, background_counts) {
  query <- sort(unique(intersect(
    as.character(manifest_row$annotated_features[[1L]]), background
  )), method = "radix")
  query_size <- length(query)
  if (query_size < 5L) return(.fsf_kegg_empty_results())

  query_mapping <- mapping[mapping$feature_id %in% query, , drop = FALSE]
  if (!nrow(query_mapping)) return(.fsf_kegg_empty_results())
  term_features <- split(query_mapping$feature_id, query_mapping$kegg_term)
  terms <- sort(names(term_features), method = "radix")
  overlaps <- lengths(term_features[terms])
  keep <- overlaps > 0L
  terms <- terms[keep]
  overlaps <- overlaps[keep]
  if (!length(terms)) return(.fsf_kegg_empty_results())

  background_size <- length(background)
  background_term_count <- as.integer(background_counts[terms])
  if (anyNA(background_term_count) || any(background_term_count <= 0L) ||
      any(background_term_count > background_size)) {
    stop("KEGG term counts contain invalid hypergeometric inputs.", call. = FALSE)
  }
  p_value <- stats::phyper(
    q = overlaps - 1L,
    m = background_term_count,
    n = background_size - background_term_count,
    k = query_size,
    lower.tail = FALSE
  )
  if (any(!is.finite(p_value)) || any(p_value < 0 | p_value > 1)) {
    stop("KEGG hypergeometric calculation returned invalid probabilities.", call. = FALSE)
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
    enrichment_type = enrichment_type,
    kegg_term = terms,
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
    out$kegg_term, method = "radix"
  ), .fsf_kegg_result_columns, drop = FALSE]
}

.run_current_kegg_enrichment <- function(gene_sets, annotation, validate_real = TRUE) {
  if (!is.list(gene_sets) ||
      !all(c("current_gene_set_manifest", "universes") %in% names(gene_sets))) {
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
  backgrounds <- list(
    KO = gene_sets$universes$ko_universe,
    PATHWAY = gene_sets$universes$pathway_universe
  )
  for (type in names(backgrounds)) {
    background <- backgrounds[[type]]
    if (is.null(background) || anyNA(background) || anyDuplicated(background) ||
        !identical(background, sort(background, method = "radix"))) {
      stop(paste0(type, " universe is invalid."), call. = FALSE)
    }
  }

  mappings <- .fsf_kegg_build_mappings(annotation)
  results <- list()
  status <- list()
  result_index <- 0L
  status_index <- 0L
  for (type in c("KO", "PATHWAY")) {
    background <- backgrounds[[type]]
    mapping <- mappings[[type]][
      mappings[[type]]$feature_id %in% background, , drop = FALSE
    ]
    background_counts <- table(mapping$kegg_term)
    for (i in seq_len(nrow(manifest))) {
      row <- manifest[i, , drop = FALSE]
      query <- sort(unique(intersect(
        as.character(row$annotated_features[[1L]]), background
      )), method = "radix")
      result <- .fsf_kegg_run_one(
        row, type, background, mapping, background_counts
      )
      result_index <- result_index + 1L
      results[[result_index]] <- result
      status_index <- status_index + 1L
      status[[status_index]] <- data.frame(
        gene_set_id = row$gene_set_id,
        enrichment_type = type,
        query_size = length(query),
        status = if (length(query) < 5L) {
          "SKIPPED_QUERY_SIZE_LT_5"
        } else if (!nrow(result)) {
          "EMPTY_NO_MAPPED_KEGG_TERM"
        } else {
          "TESTED"
        },
        result_rows = nrow(result),
        stringsAsFactors = FALSE
      )
    }
  }
  combined <- do.call(rbind, results)
  if (is.null(combined)) combined <- .fsf_kegg_empty_results()
  rownames(combined) <- NULL
  if (nrow(combined)) {
    combined <- combined[order(
      combined$gene_set_id, combined$enrichment_type,
      combined$adjusted_p_value, combined$p_value,
      -combined$enrichment_ratio, combined$kegg_term, method = "radix"
    ), , drop = FALSE]
    rownames(combined) <- NULL
  }
  status <- do.call(rbind, status)
  status <- status[order(
    status$gene_set_id, status$enrichment_type, method = "radix"
  ), , drop = FALSE]
  rownames(status) <- NULL

  if (validate_real) {
    observed <- c(
      manifest = nrow(manifest),
      ko_universe = length(backgrounds$KO),
      pathway_universe = length(backgrounds$PATHWAY),
      ko_terms = length(unique(mappings$KO$kegg_term)),
      pathway_terms = length(unique(mappings$PATHWAY$kegg_term))
    )
    expected <- c(
      manifest = 70L, ko_universe = 7178L, pathway_universe = 4680L,
      ko_terms = 5755L, pathway_terms = 786L
    )
    if (!identical(as.integer(observed), as.integer(expected))) {
      stop("Current KEGG mapping or background counts failed validation.", call. = FALSE)
    }
    if (!nrow(combined)) {
      stop("Current real-data KEGG enrichment produced no results.", call. = FALSE)
    }
  }

  list(
    KO_mapping = mappings$KO,
    PATHWAY_mapping = mappings$PATHWAY,
    backgrounds = backgrounds,
    gene_set_status = status,
    results = combined
  )
}

#' Run current FSF KEGG enrichment in memory
#'
#' Uses frozen current gene sets and hash-validated annotation authority.
#' KO and PATHWAY results are returned together; nothing is written.
run_current_kegg_enrichment <- function(annotation_master_path) {
  gene_sets <- build_current_enrichment_gene_sets(annotation_master_path)
  joined <- load_current_fsf_annotation(annotation_master_path)
  annotation <- joined[!duplicated(joined$feature_id), c(
    "feature_id", "eggnog_kegg_ko", "eggnog_kegg_pathway"
  ), drop = FALSE]
  .run_current_kegg_enrichment(gene_sets, annotation, validate_real = TRUE)
}
