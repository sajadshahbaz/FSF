test_that("duplicate keys are scoped by grouping strata", {
  cross_group <- data.frame(
    condition = c("UV", "DES"),
    feature_id = c("f1", "f1"),
    perturbation_id = c("p1", "p1"),
    effect = c(1, -1)
  )
  expect_no_error(fsf_signal_identity(cross_group, group_cols = "condition"))

  within_group <- rbind(cross_group[1, ], cross_group[1, ])
  expect_error(
    fsf_signal_identity(within_group, group_cols = "condition"),
    "duplicate"
  )
})

test_that("features are analyzed independently within groups", {
  data <- data.frame(
    condition = c("UV", "UV", "DES", "DES", "DES"),
    feature_id = "f1",
    perturbation_id = c("u1", "u2", "d1", "d2", "d3"),
    effect = c(1, -1, 0, 0, 1)
  )
  observed <- fsf_signal_identity(data, group_cols = "condition")
  expect_equal(nrow(observed), 2L)
  expect_equal(names(observed)[1:2], c("condition", "feature_id"))
  expect_equal(sort(observed$n_perturbations), c(2L, 3L))
  expect_equal(observed$p_up, observed$n_up / observed$n_perturbations)
})

test_that("multiple grouping columns preserve user order", {
  data <- data.frame(
    batch = "b1", condition = "UV", feature_id = "f1",
    perturbation_id = "p1", effect = 1
  )
  observed <- fsf_analyze(data, group_cols = c("condition", "batch"))
  expect_equal(names(observed)[1:3], c("condition", "batch", "feature_id"))
})

test_that("group_cols and grouping values are validated", {
  valid <- data.frame(
    condition = "UV", feature_id = "f1", perturbation_id = "p1", effect = 1
  )
  invalid_args <- list("absent", NA_character_, "", c("condition", "condition"),
                       "feature_id", "perturbation_id")
  for (value in invalid_args) {
    expect_error(fsf_signal_identity(valid, group_cols = value), "group_cols")
  }
  expect_no_error(fsf_signal_identity(valid, group_cols = character(0)))

  missing_group <- valid
  missing_group$condition <- NA_character_
  expect_error(fsf_signal_identity(missing_group, group_cols = "condition"),
               "condition")
  blank_group <- valid
  blank_group$condition <- ""
  expect_error(fsf_signal_identity(blank_group, group_cols = "condition"),
               "condition")
})
