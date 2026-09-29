test_that("effect-state boundaries follow the scientific lock", {
  tau <- 0.5
  fixture <- data.frame(
    feature_id = paste0("f", 1:4),
    perturbation_id = "p1",
    effect = c(tau, -tau, tau + 1e-12, -tau - 1e-12)
  )
  expect_equal(
    fsf_assign_states(fixture, tau)$state,
    c("constant", "constant", "up", "down")
  )
})

test_that("SSI region boundaries follow the scientific lock", {
  cases <- data.frame(
    p_up = c(1 / 3, 0.49, 0.50, 0.50 + 1e-12, 0.75 - 1e-12,
             0.75, 0.90 - 1e-12, 0.90, 1),
    p_down = c(1 / 3, 0.30, 0.30, 0.30, 0.20,
               0.15, 0.05, 0.05, 0),
    p_const = c(1 / 3, 0.21, 0.20, 0.20 - 1e-12, 0.05 + 1e-12,
                0.10, 0.05 + 1e-12, 0.05, 0)
  )
  cases$feature_id <- paste0("f", seq_len(nrow(cases)))

  observed <- fsf_classify(cases)
  expect_equal(
    as.character(observed$stability_region),
    c("Low Stability", "Low Stability", "Low Stability", "Transitional",
      "Transitional", "Stable", "Stable", "Highly Stable", "Highly Stable")
  )
})

test_that("classification does not use rounded probabilities", {
  identity <- data.frame(
    feature_id = c("below", "above"),
    p_up = c(0.749999999999, 0.750000000001),
    p_down = c(0.20, 0.20),
    p_const = c(0.050000000001, 0.049999999999)
  )
  observed <- fsf_classify(identity)
  expect_equal(
    as.character(observed$stability_region),
    c("Transitional", "Stable")
  )
})
