test_that("Signal Identity counts and probabilities are coherent", {
  states <- data.frame(
    feature_id = c(rep("a", 3), rep("b", 4)),
    perturbation_id = c(paste0("p", 1:3), paste0("p", 1:4)),
    state = c("up", "down", "constant", "up", "up", "up", "constant")
  )
  observed <- fsf_signal_identity(states)
  expect_equal(observed$n_up + observed$n_down + observed$n_const,
               observed$n_perturbations)
  expect_equal(observed$p_up + observed$p_down + observed$p_const,
               rep(1, nrow(observed)), tolerance = 1e-15)
  expect_equal(observed$n_perturbations, c(3L, 4L))
})

test_that("SSI and Stability Deviation retain unrounded values", {
  identity <- data.frame(
    feature_id = "x",
    p_up = 4 / 7,
    p_down = 2 / 7,
    p_const = 1 / 7
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$ssi, 4 / 7, tolerance = 1e-15)
  expect_equal(observed$stability_deviation, 3 / 7, tolerance = 1e-15)
})
