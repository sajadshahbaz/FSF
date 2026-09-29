test_that("current-lock integer-count fixtures classify exactly", {
  fixtures <- data.frame(
    feature_id = c("tie_half", "plurality_half", "uv_four_sevenths",
                   "stable_boundary", "high_boundary", "pure", "uniform"),
    p_up = c(1 / 2, 1 / 2, 4 / 7, 3 / 4, 9 / 10, 1, 1 / 3),
    p_down = c(1 / 2, 3 / 10, 2 / 7, 1 / 8, 1 / 20, 0, 1 / 3),
    p_const = c(0, 1 / 5, 1 / 7, 1 / 8, 1 / 20, 0, 1 / 3)
  )

  observed <- fsf_classify(fixtures)
  expect_equal(
    as.character(observed$signal_class),
    c("Low Stability", "Low Stability", "Transitional Up", "Stable Up",
      "Highly Stable Up", "Highly Stable Up", "Low Stability")
  )
  expect_equal(
    observed$dominant_state,
    c("tied", "up", "up", "up", "up", "up", "tied")
  )
})
