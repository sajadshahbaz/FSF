test_that("architectures include all ten classes in fixed order", {
  classified <- data.frame(
    feature_id = c("a", "b", "c"),
    condition = c("one", "one", "two"),
    signal_class = c("Low Stability", "Stable Up", "Highly Stable Constant")
  )
  observed <- fsf_architecture(classified)

  class_order <- c(
    "Low Stability",
    "Transitional Up", "Transitional Constant", "Transitional Down",
    "Stable Up", "Stable Constant", "Stable Down",
    "Highly Stable Up", "Highly Stable Constant", "Highly Stable Down"
  )
  expect_equal(nrow(observed), 2L * 10L)
  expect_equal(as.character(unique(observed$signal_class)), class_order)

  count_sums <- tapply(observed$class_count, observed$condition, sum)
  proportion_sums <- tapply(observed$class_proportion, observed$condition, sum)
  expect_equal(as.vector(count_sums), c(2, 1))
  expect_equal(as.vector(proportion_sums), c(1, 1), tolerance = 1e-15)
})


test_that("architecture preserves condition_col and validates identities", {
  classified <- data.frame(
    feature_id = c("a", "b"),
    stratum = c("one", "two"),
    signal_class = c("Low Stability", "Stable Up")
  )
  observed <- fsf_architecture(classified, condition_col = "stratum")
  expect_identical(names(observed)[1], "stratum")

  for (value in list(NA_character_, "", c("stratum", "other"), "absent")) {
    expect_error(fsf_architecture(classified, condition_col = value), "condition")
  }
  missing <- classified
  missing$stratum[1] <- NA_character_
  expect_error(fsf_architecture(missing, condition_col = "stratum"), "stratum")
  duplicated <- rbind(classified[1, ], classified[1, ])
  expect_error(fsf_architecture(duplicated, condition_col = "stratum"), "duplicate")
})

test_that("single-perturbation architectures remain descriptive", {
  classified <- data.frame(
    feature_id = c("a", "b"), condition = "one", n_perturbations = 1L,
    signal_class = c("Highly Stable Up", "Highly Stable Constant")
  )
  expect_no_error(fsf_architecture(classified))
})
