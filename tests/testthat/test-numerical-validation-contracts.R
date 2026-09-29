test_that("probability sum tolerance is validation-only", {
  accepted <- data.frame(
    feature_id = "accepted",
    p_up = 0.5000000001,
    p_down = 0.4999999999,
    p_const = 5e-9
  )
  observed <- fsf_classify(accepted)
  expect_equal(observed$p_up, accepted$p_up)
  expect_equal(observed$p_down, accepted$p_down)
  expect_equal(observed$p_const, accepted$p_const)

  rejected <- accepted
  rejected$p_const <- 2e-8
  expect_error(fsf_classify(rejected), "probabilities")
})

test_that("individual probabilities must be finite and within zero and one", {
  for (value in c(NA_real_, NaN, Inf, -Inf, -1e-12, 1 + 1e-12)) {
    identity <- data.frame(feature_id = "x", p_up = value,
                           p_down = 0.5, p_const = 0.5)
    expect_error(fsf_classify(identity), "probabilities")
  }
})
