test_that("invalid scientific inputs fail deterministically", {
  valid <- data.frame(
    feature_id = "f1",
    perturbation_id = "p1",
    effect = 0
  )

  invalid_effects <- list(NA_real_, NaN, Inf, -Inf, "0.5")
  for (value in invalid_effects) {
    bad <- valid
    bad$effect <- value
    expect_error(fsf_assign_states(bad), "effect")
  }

  for (value in list(0, -0.5, Inf, NA_real_, c(0.5, 1))) {
    expect_error(fsf_assign_states(valid, tau = value), "tau")
  }
})

test_that("identifiers, keys, rows, and columns are validated", {
  valid <- data.frame(
    feature_id = "f1",
    perturbation_id = "p1",
    effect = 0
  )

  for (column in names(valid)) {
    expect_error(fsf_assign_states(valid[setdiff(names(valid), column)]), column)
  }

  for (column in c("feature_id", "perturbation_id")) {
    bad_blank <- valid
    bad_blank[[column]] <- ""
    expect_error(fsf_assign_states(bad_blank), column)

    bad_missing <- valid
    bad_missing[[column]] <- NA_character_
    expect_error(fsf_assign_states(bad_missing), column)
  }

  duplicated <- rbind(valid, valid)
  expect_error(fsf_assign_states(duplicated), "duplicate")
  expect_error(fsf_assign_states(valid[FALSE, ]), "row")
})
