test_that("single perturbations retain mathematically valid classes", {
  data <- data.frame(
    feature_id = c("up_one", "constant_one", "down_one"),
    perturbation_id = "p1",
    effect = c(1, 0, -1)
  )
  observed <- fsf_analyze(data, tau = 0.5)
  expect_equal(observed$n_perturbations, rep(1L, 3))
  expect_equal(observed$p_up, c(1, 0, 0))
  expect_equal(observed$p_down, c(0, 0, 1))
  expect_equal(observed$p_const, c(0, 1, 0))
  expect_equal(observed$ssi, rep(1, 3))
  expect_equal(
    as.character(observed$signal_class),
    c("Highly Stable Up", "Highly Stable Constant", "Highly Stable Down")
  )
})
