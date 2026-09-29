test_that("safe committed effects reproduce historical states", {
  fixture <- read.delim(
    testthat::test_path("fixtures", "safe-golden-effects.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  observed <- fsf_assign_states(
    fixture[c("feature_id", "perturbation_id", "effect")], tau = 0.5
  )
  expect_identical(observed$state, fixture$historical_state)
})

test_that("safe committed identities reproduce invariant metrics", {
  effects <- read.delim(
    testthat::test_path("fixtures", "safe-golden-effects.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  expected <- read.delim(
    testthat::test_path("fixtures", "safe-golden-identities.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  effects <- effects[effects$include_identity, ]
  observed <- fsf_analyze(
    effects[c("feature_id", "perturbation_id", "effect")], tau = 0.5
  )
  index <- match(expected$feature_id, observed$feature_id)
  expect_false(anyNA(index))
  observed <- observed[index, ]

  for (field in c("n_perturbations", "n_up", "n_down", "n_const")) {
    expect_identical(observed[[field]], as.integer(expected[[field]]))
  }
  for (field in c("p_up", "p_down", "p_const", "ssi",
                  "stability_deviation")) {
    exact <- observed[[field]] == expected[[field]]
    if (!all(exact)) {
      expect_lte(max(abs(observed[[field]] - expected[[field]])), 1e-12)
    }
    expect_true(all(exact |
      abs(observed[[field]] - expected[[field]]) <= 1e-12))
  }
})
