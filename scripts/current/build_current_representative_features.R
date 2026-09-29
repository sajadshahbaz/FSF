# Current FSF v1 representative/stable-feature workflow.
#
# This workflow returns validated in-memory objects only. It obtains all FSF
# and annotation fields through the frozen annotation-authority join helper.

source(file.path("scripts", "lib", "current_fsf_annotation_authority.R"))

.fsf_representative_expected_candidates <- c(
  DES = 9350L, GAM = 11996L, HT = 15187L,
  LT = 13976L, OSM = 15187L, UV = 5542L
)

.fsf_representative_expected_directional <- matrix(
  c(1233L, 678L, 1755L, 1074L, 2564L, 2390L,
    84L, 42L, 238L, 203L, 355L, 458L),
  ncol = 2L,
  byrow = TRUE,
  dimnames = list(
    condition = c("DES", "GAM", "HT", "LT", "OSM", "UV"),
    dominant_state = c("down", "up")
  )
)

.fsf_rep_order <- function(data, annotated = FALSE, prioritize_region = FALSE) {
  keys <- list(data$condition, data$dominant_state)
  if (prioritize_region) {
    keys <- c(keys, list(ifelse(data$stability_region == "Highly Stable", 0L, 1L)))
  }
  if (annotated) keys <- c(keys, list(-data$annotation_source_count))
  keys <- c(keys, list(-data$ssi, data$stability_deviation, data$feature_id))
  data[do.call(order, c(keys, list(method = "radix"))), , drop = FALSE]
}

.fsf_rep_group_head <- function(data, n) {
  if (!nrow(data)) return(data)
  group <- interaction(data$condition, data$dominant_state, drop = TRUE, lex.order = TRUE)
  position <- ave(seq_len(nrow(data)), group, FUN = seq_along)
  data[position <= n, , drop = FALSE]
}

.fsf_rep_summary <- function(data) {
  keys <- c("condition", "stability_region", "dominant_state", "signal_class")
  out <- aggregate(
    rep.int(1L, nrow(data)),
    data[keys],
    FUN = sum
  )
  names(out)[ncol(out)] <- "n_features"
  totals <- ave(out$n_features, out$condition, FUN = sum)
  out$proportion <- out$n_features / totals
  out[do.call(order, c(out[keys], list(method = "radix"))), , drop = FALSE]
}

.fsf_rep_coverage <- function(data) {
  keys <- c("condition", "stability_region", "dominant_state", "signal_class")
  groups <- interaction(data[keys], drop = TRUE, lex.order = TRUE)
  indices <- split(seq_len(nrow(data)), groups)
  rows <- lapply(indices, function(i) {
    x <- data[i, , drop = FALSE]
    result <- x[1L, keys, drop = FALSE]
    result$n_features <- nrow(x)
    result$eggnog_annotated <- sum(x$eggnog_annotated, na.rm = TRUE)
    result$interpro_annotated <- sum(x$interpro_annotated, na.rm = TRUE)
    result$pfam_annotated <- sum(x$pfam_annotated, na.rm = TRUE)
    result$at_least_one_annotation <- sum(x$annotation_source_count >= 1, na.rm = TRUE)
    result$at_least_two_annotation_layers <- sum(x$annotation_source_count >= 2, na.rm = TRUE)
    result$all_three_annotation_layers <- sum(x$annotation_source_count == 3, na.rm = TRUE)
    result$unannotated <- sum(x$annotation_source_count == 0, na.rm = TRUE)
    result$mean_annotation_source_count <- mean(x$annotation_source_count)
    result$mean_ssi <- mean(x$ssi)
    result$median_ssi <- stats::median(x$ssi)
    result$mean_stability_deviation <- mean(x$stability_deviation)
    result
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$eggnog_fraction <- out$eggnog_annotated / out$n_features
  out$interpro_fraction <- out$interpro_annotated / out$n_features
  out$pfam_fraction <- out$pfam_annotated / out$n_features
  out$annotated_fraction <- out$at_least_one_annotation / out$n_features
  out$unannotated_fraction <- out$unannotated / out$n_features
  out[do.call(order, c(out[keys], list(method = "radix"))), , drop = FALSE]
}

.build_current_representative_features <- function(joined, validate_real = TRUE) {
  required <- c(.fsf_current_columns, .fsf_annotation_columns[-1L])
  if (!is.data.frame(joined) || !identical(names(joined), required)) {
    stop("Joined current authority must have the exact required 49-column schema.", call. = FALSE)
  }
  if (any(duplicated(joined[c("condition", "feature_id")]))) {
    stop("Joined current authority contains duplicate condition + feature_id keys.", call. = FALSE)
  }

  eligible <- joined$stability_region %in% c("Stable", "Highly Stable")
  stable_annotated <- joined[eligible, , drop = FALSE]
  stable_annotated <- .fsf_rep_order(stable_annotated)
  stable <- stable_annotated[.fsf_current_columns]

  directional_annotated <- stable_annotated[
    stable_annotated$dominant_state %in% c("up", "down"), , drop = FALSE
  ]
  directional_annotated <- .fsf_rep_order(directional_annotated)
  directional <- directional_annotated[.fsf_current_columns]
  constant <- stable[stable$dominant_state == "constant", , drop = FALSE]

  top100 <- .fsf_rep_group_head(directional, 100L)

  annotated_directional <- directional_annotated[
    directional_annotated$annotation_source_count >= 1, , drop = FALSE
  ]
  annotated_top100 <- .fsf_rep_group_head(
    .fsf_rep_order(annotated_directional, annotated = TRUE, prioritize_region = TRUE),
    100L
  )
  annotated_top50 <- .fsf_rep_group_head(
    .fsf_rep_order(annotated_directional, annotated = TRUE),
    50L
  )

  output <- list(
    current_fsf_stable_signal_catalog = stable,
    current_fsf_stable_directional_signals = directional,
    current_fsf_stable_constant_signals = constant,
    current_fsf_stable_signal_summary = .fsf_rep_summary(stable),
    current_fsf_top100_stable_directional_signals = top100,
    current_fsf_stable_signal_annotated_catalog = stable_annotated,
    current_fsf_stable_signal_annotation_coverage = .fsf_rep_coverage(stable_annotated),
    current_fsf_top100_annotated_stable_directional_signals = annotated_top100,
    current_fsf_top50_annotated_stable_directional_signals_compact = annotated_top50
  )

  if (validate_real) {
    observed <- table(factor(stable$condition, levels = names(.fsf_representative_expected_candidates)))
    if (nrow(stable) != 71238L || !identical(as.integer(observed), unname(.fsf_representative_expected_candidates))) {
      stop("Current stable candidate counts do not match the frozen validation contract.", call. = FALSE)
    }
    directional_counts <- xtabs(~ condition + dominant_state, directional)
    directional_counts <- directional_counts[
      rownames(.fsf_representative_expected_directional),
      colnames(.fsf_representative_expected_directional),
      drop = FALSE
    ]
    directional_counts_valid <-
      identical(
        dim(directional_counts),
        dim(.fsf_representative_expected_directional)
      ) &&
      identical(
        dimnames(directional_counts),
        dimnames(.fsf_representative_expected_directional)
      ) &&
      all(is.finite(directional_counts)) &&
      all(directional_counts == as.integer(directional_counts)) &&
      identical(
        as.integer(directional_counts),
        as.integer(.fsf_representative_expected_directional)
      )
    if (!directional_counts_valid) {
      stop("Current directional candidate counts do not match the frozen validation contract.", call. = FALSE)
    }
    uv_migrated <- joined$condition == "UV" & joined$ssi > 0.50 & joined$ssi < 0.60
    if (sum(uv_migrated) != 2451L ||
        any(joined$stability_region[uv_migrated] != "Transitional") ||
        any(eligible[uv_migrated])) {
      stop("UV Transitional migration guard failed.", call. = FALSE)
    }
    if (nrow(top100) != 1126L || nrow(annotated_top100) != 1115L ||
        nrow(annotated_top50) != 590L) {
      stop("Current representative top-list totals do not match the frozen contract.", call. = FALSE)
    }
  }

  output
}

#' Build current FSF representative-feature objects in memory
#'
#' @param annotation_master_path Explicit path to the hash-locked biological
#'   annotation master.
#' @return A named list of nine current representative/stable-feature objects.
#'   No files are written.
build_current_representative_features <- function(annotation_master_path) {
  joined <- load_current_fsf_annotation(annotation_master_path)
  .build_current_representative_features(joined, validate_real = TRUE)
}
