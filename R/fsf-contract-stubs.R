# Internal constants -------------------------------------------------------

.fsf_probability_tolerance <- 1e-8

.fsf_stability_levels <- c(
  "Low Stability", "Transitional", "Stable", "Highly Stable"
)

.fsf_class_levels <- c(
  "Low Stability",
  "Transitional Up", "Transitional Constant", "Transitional Down",
  "Stable Up", "Stable Constant", "Stable Down",
  "Highly Stable Up", "Highly Stable Constant", "Highly Stable Down"
)

.fsf_stop <- function(message) {
  stop(message, call. = FALSE)
}

.fsf_validate_data_frame <- function(data, argument = "data") {
  if (!is.data.frame(data)) {
    .fsf_stop(sprintf("%s must be a data.frame", argument))
  }
  if (nrow(data) == 0L) {
    .fsf_stop(sprintf("%s must contain at least one row", argument))
  }
  invisible(data)
}

.fsf_require_columns <- function(data, columns) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) {
    .fsf_stop(sprintf("missing required column: %s", absent[[1L]]))
  }
  invisible(data)
}

.fsf_is_blank_character <- function(x) {
  is.character(x) & trimws(x) == ""
}

.fsf_validate_identifier <- function(x, name) {
  if (anyNA(x)) {
    .fsf_stop(sprintf("%s must not contain missing values", name))
  }
  if (is.character(x) && any(.fsf_is_blank_character(x))) {
    .fsf_stop(sprintf("%s must not contain blank values", name))
  }
  invisible(x)
}

.fsf_validate_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1L) {
    .fsf_stop("tau must be one numeric scalar")
  }
  if (is.na(tau) || !is.finite(tau) || tau <= 0) {
    .fsf_stop("tau must be finite and greater than zero")
  }
  invisible(tau)
}

.fsf_validate_group_cols <- function(data, group_cols) {
  if (is.null(group_cols) || (is.character(group_cols) && !length(group_cols))) {
    return(character())
  }
  if (!is.character(group_cols)) {
    .fsf_stop("group_cols must be NULL or a character vector")
  }
  if (anyNA(group_cols) || any(trimws(group_cols) == "")) {
    .fsf_stop("group_cols must not contain missing or blank names")
  }
  if (anyDuplicated(group_cols)) {
    .fsf_stop("group_cols must not contain duplicated names")
  }
  reserved <- intersect(group_cols, c("feature_id", "perturbation_id"))
  if (length(reserved)) {
    .fsf_stop("group_cols must not include feature_id or perturbation_id")
  }
  absent <- setdiff(group_cols, names(data))
  if (length(absent)) {
    .fsf_stop(sprintf("group_cols column not found: %s", absent[[1L]]))
  }
  for (column in group_cols) {
    value <- data[[column]]
    if (anyNA(value)) {
      .fsf_stop(sprintf("grouping column %s must not contain missing values", column))
    }
    if (is.character(value) && any(.fsf_is_blank_character(value))) {
      .fsf_stop(sprintf("grouping column %s must not contain blank values", column))
    }
  }
  group_cols
}

.fsf_validate_duplicate_key <- function(data, key) {
  if (anyDuplicated(data[key])) {
    .fsf_stop(sprintf("duplicate mathematical observation key: %s",
                      paste(key, collapse = " + ")))
  }
  invisible(data)
}

.fsf_validate_effect_input <- function(data, tau, group_cols = NULL) {
  .fsf_validate_data_frame(data)
  .fsf_require_columns(data, c("feature_id", "perturbation_id", "effect"))
  .fsf_validate_tau(tau)
  groups <- .fsf_validate_group_cols(data, group_cols)
  .fsf_validate_identifier(data$feature_id, "feature_id")
  .fsf_validate_identifier(data$perturbation_id, "perturbation_id")
  if (!is.numeric(data$effect)) {
    .fsf_stop("effect must be numeric")
  }
  if (anyNA(data$effect)) {
    .fsf_stop("effect must not contain missing values")
  }
  if (any(!is.finite(data$effect))) {
    .fsf_stop("effect must be finite")
  }
  .fsf_validate_duplicate_key(
    data, c(groups, "feature_id", "perturbation_id")
  )
  groups
}

.fsf_validate_state_input <- function(data, tau, group_cols = NULL) {
  .fsf_validate_data_frame(data)
  .fsf_require_columns(data, c("feature_id", "perturbation_id", "state"))
  .fsf_validate_tau(tau)
  groups <- .fsf_validate_group_cols(data, group_cols)
  .fsf_validate_identifier(data$feature_id, "feature_id")
  .fsf_validate_identifier(data$perturbation_id, "perturbation_id")
  if (!(is.character(data$state) || is.factor(data$state)) || anyNA(data$state)) {
    .fsf_stop("state must contain only up, down, or constant")
  }
  states <- as.character(data$state)
  if (any(!states %in% c("up", "down", "constant"))) {
    .fsf_stop("state must contain only up, down, or constant")
  }
  .fsf_validate_duplicate_key(
    data, c(groups, "feature_id", "perturbation_id")
  )
  groups
}

.fsf_identity_rows <- function(data, identity_cols) {
  keys <- unique(data[identity_cols])
  rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    selected <- rep(TRUE, nrow(data))
    for (column in identity_cols) {
      selected <- selected & data[[column]] == keys[[column]][[i]]
    }
    rows[[i]] <- which(selected)
  }
  list(keys = keys, rows = rows)
}

.fsf_validate_probabilities <- function(identity) {
  .fsf_validate_data_frame(identity, "identity")
  .fsf_require_columns(identity, c("feature_id", "p_up", "p_down", "p_const"))
  .fsf_validate_identifier(identity$feature_id, "feature_id")
  probability_names <- c("p_up", "p_down", "p_const")
  for (column in probability_names) {
    value <- identity[[column]]
    if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
        any(value < 0 | value > 1)) {
      .fsf_stop("Signal Identity probabilities must be numeric, finite, and in [0,1]")
    }
  }
  probability_sum <- identity$p_up + identity$p_down + identity$p_const
  if (any(abs(probability_sum - 1) > .fsf_probability_tolerance)) {
    .fsf_stop("Signal Identity probabilities must sum to one within 1e-8")
  }

  count_names <- c("n_perturbations", "n_up", "n_down", "n_const")
  present <- count_names %in% names(identity)
  if (any(present) && !all(present)) {
    .fsf_stop("Signal Identity counts must include all count columns")
  }
  if (all(present)) {
    for (column in count_names) {
      value <- identity[[column]]
      if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
          any(value < 0 | value != floor(value))) {
        .fsf_stop("Signal Identity counts must be finite nonnegative integers")
      }
    }
    if (any(identity$n_perturbations <= 0)) {
      .fsf_stop("n_perturbations must be greater than zero")
    }
    if (any(identity$n_up + identity$n_down + identity$n_const !=
            identity$n_perturbations)) {
      .fsf_stop("Signal Identity counts must sum to n_perturbations")
    }
    count_probabilities <- cbind(identity$n_up, identity$n_down,
                                 identity$n_const) / identity$n_perturbations
    supplied_probabilities <- as.matrix(identity[probability_names])
    if (any(abs(count_probabilities - supplied_probabilities) >
            .fsf_probability_tolerance)) {
      .fsf_stop("Signal Identity counts and probabilities must agree within 1e-8")
    }
  }
  invisible(identity)
}

.fsf_state_values <- function(effect, tau) {
  ifelse(
    effect > tau,
    "up",
    ifelse(effect < -tau, "down", "constant")
  )
}

# Public API ---------------------------------------------------------------

#' Assign FSF directional states
#'
#' Assigns Up, Down, or Constant to finite signed effects. Up is
#' `effect > tau`, Down is `effect < -tau`, and Constant is
#' `-tau <= effect <= tau`.
#'
#' @param data A data frame containing `feature_id`, `perturbation_id`, and
#'   `effect` columns.
#' @param tau One positive finite numeric scalar. The default is `0.5`.
#'
#' @return The validated input with a character `state` column containing
#'   `up`, `down`, or `constant`.
#' @details Additional input columns are retained. An existing `state` column
#'   is replaced. Invalid input produces an error rather than being repaired.
#'
#' @section Errors:
#' Errors on missing columns, invalid identifiers, duplicate observation keys,
#' nonnumeric or nonfinite effects, empty input, or invalid `tau`.
#'
#' @examples
#' effects <- data.frame(
#'   feature_id = paste0("f", 1:4),
#'   perturbation_id = "p1",
#'   effect = c(0.5, -0.5, 0.6, -0.6)
#' )
#' fsf_assign_states(effects, tau = 0.5)
#' # Exact +tau and -tau are constant; only values beyond them are directional.
#'
#' @export
fsf_assign_states <- function(data, tau = 0.5) {
  .fsf_validate_effect_input(data, tau)
  output <- data
  output$state <- .fsf_state_values(data$effect, tau)
  output
}

#' Calculate FSF Signal Identity
#'
#' Aggregates feature-perturbation states into state counts and Signal Identity
#' probabilities independently within optional grouping strata.
#'
#' @param data A data frame containing canonical signed effects or validated
#'   `up`, `down`, and `constant` states.
#' @param tau One positive finite numeric scalar used when `data` contains
#'   signed effects. The default is `0.5`.
#' @param group_cols `NULL`, `character(0)`, or a character vector naming
#'   columns that define independent grouping strata.
#'
#' @return One row per Signal Identity. Grouping columns occur first in
#'   user-supplied order, followed by `feature_id`, counts, and probabilities.
#'   Identities follow their first appearance in the input.
#' @details If both `state` and `effect` columns are present, the validated
#'   `state` column takes precedence and `effect` is not used. Columns other
#'   than grouping columns and the required identity inputs are not propagated
#'   to the aggregated output. State probabilities are the corresponding state
#'   counts divided by each identity's own `n_perturbations`.
#'
#' @section Errors:
#' Errors on invalid effect/state input, grouping specifications, identifiers,
#' duplicate observation keys, empty input, or invalid `tau`.
#'
#' @examples
#' effects <- data.frame(
#'   condition = c("A", "A", "B", "B"),
#'   feature_id = "gene1",
#'   perturbation_id = c("p1", "p2", "p1", "p2"),
#'   effect = c(1, 0, -1, -1)
#' )
#' fsf_signal_identity(effects, tau = 0.5, group_cols = "condition")
#'
#' @export
fsf_signal_identity <- function(data, tau = 0.5, group_cols = NULL) {
  .fsf_validate_data_frame(data)
  if ("state" %in% names(data)) {
    groups <- .fsf_validate_state_input(data, tau, group_cols)
    states <- as.character(data$state)
  } else {
    groups <- .fsf_validate_effect_input(data, tau, group_cols)
    states <- .fsf_state_values(data$effect, tau)
  }

  identity_cols <- c(groups, "feature_id")
  partitions <- .fsf_identity_rows(data, identity_cols)
  output <- partitions$keys
  n <- lengths(partitions$rows)
  count_state <- function(state) {
    as.integer(vapply(partitions$rows, function(index) {
      sum(states[index] == state)
    }, integer(1L)))
  }
  output$n_perturbations <- as.integer(n)
  output$n_up <- count_state("up")
  output$n_down <- count_state("down")
  output$n_const <- count_state("constant")
  output$p_up <- output$n_up / output$n_perturbations
  output$p_down <- output$n_down / output$n_perturbations
  output$p_const <- output$n_const / output$n_perturbations

  if (any(output$n_up + output$n_down + output$n_const !=
          output$n_perturbations)) {
    .fsf_stop("Signal Identity counts must sum to n_perturbations")
  }
  probability_sum <- output$p_up + output$p_down + output$p_const
  if (any(abs(probability_sum - 1) > .fsf_probability_tolerance)) {
    .fsf_stop("Signal Identity probabilities must sum to one within 1e-8")
  }
  output
}

#' Classify FSF Signal Identities
#'
#' Derives exact dominant state, SSI, current stability region, current signal
#' class, and auxiliary Stability Deviation from validated probabilities.
#'
#' @param identity A data frame containing `feature_id`, `p_up`, `p_down`, and
#'   `p_const`, with optional grouping and state-count columns.
#'
#' @return The identity columns plus `dominant_state`, `ssi`,
#'   `stability_region`, `signal_class`, and `stability_deviation`.
#' @details `identity` must already contain one aggregated row per intended
#'   identity; this row-wise function does not infer grouping columns or merge
#'   duplicate identities. Additional columns are retained. When count metadata
#'   is supplied, all four count columns are required, `n_perturbations` must
#'   be positive, component counts must sum to it, and count-derived
#'   probabilities must agree with supplied probabilities within `1e-8`.
#'
#' SSI is the exact maximum of `p_up`, `p_down`, and `p_const`. Regions are
#' Low Stability for SSI <= 0.50, Transitional for 0.50 < SSI < 0.75,
#' Stable for 0.75 <= SSI < 0.90, and Highly Stable for SSI >= 0.90.
#' Equal maximum probabilities produce `dominant_state = "tied"`. Low
#' Stability is always non-directional, even when its dominant state is unique.
#'
#' @section Errors:
#' Errors on invalid identifiers, probabilities, partial or inconsistent count
#' metadata, or any classification state outside the current FSF contract.
#'
#' @examples
#' identity <- data.frame(
#'   feature_id = c("half", "majority", "stable", "high"),
#'   p_up = c(0.50, 0.51, 0.75, 0.90),
#'   p_down = c(0.30, 0.29, 0.15, 0.05),
#'   p_const = c(0.20, 0.20, 0.10, 0.05)
#' )
#' fsf_classify(identity)
#'
#' @export
fsf_classify <- function(identity) {
  .fsf_validate_probabilities(identity)
  probabilities <- as.matrix(identity[c("p_up", "p_down", "p_const")])
  ssi <- apply(probabilities, 1L, max)
  maximum_count <- rowSums(probabilities == ssi)
  maximum_index <- max.col(probabilities, ties.method = "first")
  state_names <- c("up", "down", "constant")
  dominant_state <- state_names[maximum_index]
  dominant_state[maximum_count > 1L] <- "tied"

  region <- ifelse(
    ssi <= 0.50, "Low Stability",
    ifelse(ssi < 0.75, "Transitional",
           ifelse(ssi < 0.90, "Stable", "Highly Stable"))
  )

  directional <- ssi > 0.50
  if (any(directional & dominant_state == "tied")) {
    .fsf_stop("validated identity has an impossible tied maximum above 0.50")
  }
  state_title <- c(up = "Up", down = "Down", constant = "Constant")
  signal_class <- rep("Low Stability", nrow(identity))
  signal_class[directional] <- paste(
    region[directional], state_title[dominant_state[directional]]
  )
  if (any(!signal_class %in% .fsf_class_levels)) {
    .fsf_stop("classification produced an invalid FSF signal class")
  }

  output <- identity
  output$dominant_state <- dominant_state
  output$ssi <- as.numeric(ssi)
  output$stability_region <- factor(
    region, levels = .fsf_stability_levels, ordered = TRUE
  )
  output$signal_class <- factor(
    signal_class, levels = .fsf_class_levels, ordered = TRUE
  )
  output$stability_deviation <- 1 - output$ssi
  output
}

#' Run the complete FSF feature-level analysis
#'
#' Composes state assignment, grouped Signal Identity calculation, and
#' current-lock classification. It does not perform effect estimation,
#' feature selection, or feature ranking.
#'
#' @param data A long-form data frame containing `feature_id`,
#'   `perturbation_id`, and finite numeric `effect` columns.
#' @param tau One positive finite numeric scalar. The default is `0.5`.
#' @param group_cols Optional grouping columns defining independent strata.
#'
#' @return One row per Signal Identity following the stable feature-level
#'   output contract, with grouping columns first when supplied.
#'
#' @details With one perturbation, SSI = 1 is structurally determined and does
#'   not constitute evidence of cross-perturbation stability.
#'   Non-grouping extra input columns are not propagated because analysis
#'   returns the aggregated Signal Identity schema. Identity rows retain the
#'   first-appearance order of their grouping values and feature identifiers.
#'
#' @section Errors:
#' Propagates deterministic validation errors from Signal Identity calculation
#' and classification; invalid inputs are never silently repaired.
#'
#' @examples
#' effects <- data.frame(
#'   feature_id = rep(c("f1", "f2"), each = 3),
#'   perturbation_id = rep(paste0("p", 1:3), 2),
#'   effect = c(1, 1, 0, -1, -1, -1)
#' )
#' fsf_analyze(effects, tau = 0.5)
#'
#' @export
fsf_analyze <- function(data, tau = 0.5, group_cols = NULL) {
  identity <- fsf_signal_identity(data, tau, group_cols)
  fsf_classify(identity)
}

#' Calculate condition-level FSF signal architecture
#'
#' Counts all ten current signal classes within each supplied condition and
#' calculates directly comparable class proportions.
#'
#' @param data A classified feature data frame containing a condition column,
#'   `feature_id`, and `signal_class`.
#' @param condition_col One nonblank character scalar naming the condition
#'   column.
#'
#' @return All ten signal classes for every condition in deterministic order,
#'   including zero-count classes. The first column retains the name supplied
#'   by `condition_col`.
#'
#' @details An architecture derived exclusively from single-perturbation
#'   identities describes directional class composition, not demonstrated
#'   reproducibility across repeated perturbations.
#'   Conditions follow their first appearance in `data`; within every condition
#'   all ten signal classes use the fixed FSF class order. Other columns are not
#'   propagated to the four-column architecture output.
#'
#' @section Errors:
#' Errors on an invalid condition column, identifiers, duplicate
#' condition-feature keys, empty input, or non-current signal classes.
#'
#' @examples
#' classified <- data.frame(
#'   condition = c("A", "A", "B"),
#'   feature_id = c("f1", "f2", "f3"),
#'   signal_class = c("Low Stability", "Stable Up", "Highly Stable Down")
#' )
#' fsf_architecture(classified)
#'
#' @export
fsf_architecture <- function(data, condition_col = "condition") {
  .fsf_validate_data_frame(data)
  if (!is.character(condition_col) || length(condition_col) != 1L ||
      is.na(condition_col) || trimws(condition_col) == "") {
    .fsf_stop("condition_col must be one nonblank character scalar")
  }
  if (!condition_col %in% names(data)) {
    .fsf_stop(sprintf("condition_col column not found: %s", condition_col))
  }
  .fsf_require_columns(data, c("feature_id", "signal_class"))
  .fsf_validate_identifier(data[[condition_col]], condition_col)
  .fsf_validate_identifier(data$feature_id, "feature_id")
  .fsf_validate_duplicate_key(data, c(condition_col, "feature_id"))

  classes <- as.character(data$signal_class)
  if (anyNA(classes) || any(!classes %in% .fsf_class_levels)) {
    .fsf_stop("signal_class must contain only current FSF classes")
  }

  conditions <- unique(data[condition_col])
  pieces <- vector("list", nrow(conditions))
  for (i in seq_len(nrow(conditions))) {
    selected <- data[[condition_col]] == conditions[[condition_col]][[i]]
    class_count <- as.integer(tabulate(
      match(classes[selected], .fsf_class_levels),
      nbins = length(.fsf_class_levels)
    ))
    condition_values <- conditions[rep(i, length(.fsf_class_levels)), , drop = FALSE]
    pieces[[i]] <- data.frame(
      condition_values,
      signal_class = factor(
        .fsf_class_levels, levels = .fsf_class_levels, ordered = TRUE
      ),
      class_count = class_count,
      class_proportion = class_count / sum(class_count),
      check.names = FALSE
    )
  }
  output <- do.call(rbind, pieces)
  rownames(output) <- NULL
  output
}
