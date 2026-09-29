# Current FSF v1 enrichment gene-set construction workflow.
#
# This workflow builds deterministic membership objects in memory only. All
# FSF and annotation fields enter through the frozen annotation-authority join
# helper. No enrichment statistic is calculated and no file is written.

source(file.path("scripts", "lib", "current_fsf_annotation_authority.R"))

.fsf_enrichment_stop <- function(message) {
  stop(message, call. = FALSE)
}

.fsf_enrichment_sanitize_id <- function(parts) {
  gsub("[^A-Za-z0-9_]+", "_", paste(parts, collapse = "__"))
}

.fsf_enrichment_nonempty <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x))) & trimws(as.character(x)) != "-"
}

.fsf_enrichment_has_go <- function(x) {
  .fsf_enrichment_nonempty(x) & grepl("GO:[0-9]{7}", as.character(x))
}

.fsf_enrichment_make_family <- function(data, family, grouping) {
  keys <- do.call(
    paste,
    c(lapply(data[grouping], as.character), list(sep = "\r"))
  )
  indices <- split(seq_len(nrow(data)), keys, drop = TRUE)
  rows <- lapply(indices, function(i) {
    metadata <- data[i[1L], grouping, drop = FALSE]
    all_features <- sort(unique(as.character(data$feature_id[i])), method = "radix")
    annotated_features <- sort(unique(as.character(
      data$feature_id[i][data$annotation_source_count[i] >= 1]
    )), method = "radix")
    values <- as.character(metadata[1L, grouping, drop = TRUE])
    result <- data.frame(
      gene_set_id = .fsf_enrichment_sanitize_id(values),
      gene_set_family = family,
      condition = as.character(metadata$condition),
      stability_region = as.character(metadata$stability_region),
      dominant_state = if ("dominant_state" %in% grouping) {
        as.character(metadata$dominant_state)
      } else {
        NA_character_
      },
      signal_class = if ("signal_class" %in% grouping) {
        as.character(metadata$signal_class)
      } else {
        NA_character_
      },
      all_feature_count = length(all_features),
      annotated_feature_count = length(annotated_features),
      stringsAsFactors = FALSE
    )
    result$all_features <- I(list(all_features))
    result$annotated_features <- I(list(annotated_features))
    result
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(
    result$condition,
    result$stability_region,
    result$dominant_state,
    result$signal_class,
    result$gene_set_id,
    method = "radix",
    na.last = TRUE
  ), , drop = FALSE]
}

.fsf_enrichment_validate_membership <- function(family) {
  valid <- vapply(family$all_features, function(x) {
    identical(x, sort(unique(x), method = "radix"))
  }, logical(1L))
  valid_annotated <- vapply(family$annotated_features, function(x) {
    identical(x, sort(unique(x), method = "radix"))
  }, logical(1L))
  if (!all(valid) || !all(valid_annotated)) {
    .fsf_enrichment_stop("Gene-set memberships are not sorted and unique.")
  }
  if (any(family$all_feature_count != lengths(family$all_features)) ||
      any(family$annotated_feature_count != lengths(family$annotated_features))) {
    .fsf_enrichment_stop("Gene-set manifest counts do not match memberships.")
  }
  invisible(family)
}

.fsf_enrichment_validate_joined <- function(joined) {
  required <- c(.fsf_current_columns, .fsf_annotation_columns[-1L])
  if (!is.data.frame(joined) || !identical(names(joined), required)) {
    .fsf_enrichment_stop(
      "Joined current authority must have the exact required 49-column schema."
    )
  }
  if (anyNA(joined$condition) || anyNA(joined$feature_id) ||
      any(trimws(as.character(joined$condition)) == "") ||
      any(trimws(as.character(joined$feature_id)) == "")) {
    .fsf_enrichment_stop("Joined authority contains missing or blank keys.")
  }
  if (any(duplicated(joined[c("condition", "feature_id")]))) {
    .fsf_enrichment_stop(
      "Joined authority contains duplicate condition + feature_id keys."
    )
  }
  count <- joined$annotation_source_count
  if (!is.numeric(count) || anyNA(count) || any(!is.finite(count)) ||
      any(count < 0) || any(count != as.integer(count))) {
    .fsf_enrichment_stop("annotation_source_count must be finite nonnegative integers.")
  }
  invisible(joined)
}

.fsf_enrichment_expected_uv <- c(
  constant = 1624L,
  down = 395L,
  up = 432L,
  tied = 0L
)

.build_current_enrichment_gene_sets <- function(joined, validate_real = TRUE) {
  .fsf_enrichment_validate_joined(joined)

  main <- .fsf_enrichment_make_family(
    joined,
    "main_class",
    c("condition", "stability_region", "dominant_state", "signal_class")
  )
  collapsed <- .fsf_enrichment_make_family(
    joined,
    "collapsed_stability",
    c("condition", "stability_region")
  )
  directional_data <- joined[
    joined$stability_region %in% c("Stable", "Highly Stable") &
      joined$dominant_state %in% c("up", "down"),
    ,
    drop = FALSE
  ]
  directional <- .fsf_enrichment_make_family(
    directional_data,
    "directional_stable",
    c("condition", "stability_region", "dominant_state")
  )

  lapply(list(main, collapsed, directional), .fsf_enrichment_validate_membership)
  manifest <- rbind(main, collapsed, directional)
  manifest <- manifest[order(
    manifest$gene_set_family,
    manifest$condition,
    manifest$stability_region,
    manifest$dominant_state,
    manifest$signal_class,
    manifest$gene_set_id,
    method = "radix",
    na.last = TRUE
  ), , drop = FALSE]
  rownames(manifest) <- NULL

  collision <- duplicated(manifest$gene_set_id) |
    duplicated(manifest$gene_set_id, fromLast = TRUE)
  if (any(collision)) {
    .fsf_enrichment_stop(paste0(
      "Current gene_set_id collision detected: ",
      paste(unique(manifest$gene_set_id[collision]), collapse = ", "), "."
    ))
  }

  feature_index <- !duplicated(joined$feature_id)
  feature_view <- joined[feature_index, , drop = FALSE]
  universes <- list(
    all_feature_universe = sort(unique(as.character(feature_view$feature_id)), method = "radix"),
    annotation_universe = sort(unique(as.character(
      feature_view$feature_id[feature_view$annotation_source_count >= 1]
    )), method = "radix"),
    go_mapped_features = sort(unique(as.character(feature_view$feature_id[
      .fsf_enrichment_has_go(feature_view$eggnog_go) |
        .fsf_enrichment_has_go(feature_view$interpro_go)
    ])), method = "radix"),
    ko_universe = sort(unique(as.character(feature_view$feature_id[
      .fsf_enrichment_nonempty(feature_view$eggnog_kegg_ko)
    ])), method = "radix"),
    pathway_universe = sort(unique(as.character(feature_view$feature_id[
      .fsf_enrichment_nonempty(feature_view$eggnog_kegg_pathway)
    ])), method = "radix")
  )

  uv <- joined[
    joined$condition == "UV" & joined$ssi > 0.50 & joined$ssi < 0.60,
    ,
    drop = FALSE
  ]
  uv_levels <- names(.fsf_enrichment_expected_uv)
  uv_breakdown <- table(factor(uv$dominant_state, levels = uv_levels))

  if (validate_real) {
    expected_sizes <- c(
      main = 40L,
      collapsed = 14L,
      directional = 16L,
      all = 15187L,
      annotation = 13467L,
      go = 9321L,
      ko = 7178L,
      pathway = 4680L
    )
    observed_sizes <- c(
      main = nrow(main),
      collapsed = nrow(collapsed),
      directional = nrow(directional),
      all = length(universes$all_feature_universe),
      annotation = length(universes$annotation_universe),
      go = length(universes$go_mapped_features),
      ko = length(universes$ko_universe),
      pathway = length(universes$pathway_universe)
    )
    if (!identical(as.integer(observed_sizes), as.integer(expected_sizes))) {
      .fsf_enrichment_stop("Current gene-set family or universe sizes failed validation.")
    }
    if (nrow(uv) != 2451L ||
        any(uv$stability_region != "Transitional") ||
        !identical(as.integer(uv_breakdown), as.integer(.fsf_enrichment_expected_uv))) {
      .fsf_enrichment_stop("Current UV migration contract failed validation.")
    }
  }

  # Frozen background contracts for later statistical phases:
  # GO uses annotation_universe (13,467 genes), not go_mapped_features.
  # KEGG uses separate ko_universe (7,178) and pathway_universe (4,680).
  list(
    main_class_sets = main,
    collapsed_stability_sets = collapsed,
    directional_stable_sets = directional,
    current_gene_set_manifest = manifest,
    universes = universes,
    uv_migration = list(
      feature_ids = sort(uv$feature_id, method = "radix"),
      dominant_state_counts = uv_breakdown
    )
  )
}

#' Build current enrichment gene-set objects in memory
#'
#' @param annotation_master_path Explicit path to the hash-locked biological
#'   annotation master.
#' @return Current membership families, deterministic manifest, background
#'   universes, and UV migration diagnostics. Nothing is written.
build_current_enrichment_gene_sets <- function(annotation_master_path) {
  joined <- load_current_fsf_annotation(annotation_master_path)
  .build_current_enrichment_gene_sets(joined, validate_real = TRUE)
}
