test_that("count metadata is positive and agrees with probabilities", {
  valid <- data.frame(
    feature_id = "f1",
    n_perturbations = 4L,
    n_up = 2L,
    n_down = 1L,
    n_const = 1L,
    p_up = 0.5,
    p_down = 0.25,
    p_const = 0.25
  )
  expect_no_error(fsf_classify(valid))

  zero <- valid
  zero[c("n_perturbations", "n_up", "n_down", "n_const")] <- 0L
  expect_error(fsf_classify(zero), "greater than zero")

  bad_sum <- valid
  bad_sum$n_const <- 0L
  expect_error(fsf_classify(bad_sum), "sum to n_perturbations")

  bad_probability <- valid
  bad_probability$p_up <- 0.4
  bad_probability$p_down <- 0.35
  expect_error(fsf_classify(bad_probability), "counts and probabilities")

  within_tolerance <- valid
  within_tolerance$p_up <- within_tolerance$p_up + 5e-9
  within_tolerance$p_down <- within_tolerance$p_down - 5e-9
  expect_no_error(fsf_classify(within_tolerance))
})

test_that("state takes precedence when state and effect coexist", {
  both <- data.frame(
    feature_id = c("f1", "f1"),
    perturbation_id = c("p1", "p2"),
    effect = c(-10, -10),
    state = c("up", "constant")
  )
  observed <- fsf_signal_identity(both, tau = 0.5)
  expect_equal(observed$n_up, 1L)
  expect_equal(observed$n_down, 0L)
  expect_equal(observed$n_const, 1L)
})

test_that("extra-column behavior follows each public output schema", {
  effects <- data.frame(
    feature_id = c("f1", "f1"),
    perturbation_id = c("p1", "p2"),
    effect = c(1, 0),
    note = c("first", "second")
  )
  states <- fsf_assign_states(effects)
  expect_identical(states$note, effects$note)

  identity <- fsf_signal_identity(effects)
  expect_false("note" %in% names(identity))

  identity$annotation <- "retained"
  classified <- fsf_classify(identity)
  expect_identical(classified$annotation, "retained")

  analyzed <- fsf_analyze(effects)
  expect_false("note" %in% names(analyzed))

  architecture_input <- data.frame(
    condition = "A", feature_id = "f1",
    signal_class = "Stable Up", note = "discarded"
  )
  architecture <- fsf_architecture(architecture_input)
  expect_identical(
    names(architecture),
    c("condition", "signal_class", "class_count", "class_proportion")
  )
})

test_that("identity and architecture preserve first-appearance order", {
  effects <- data.frame(
    condition = c("B", "B", "A", "A"),
    feature_id = c("z", "z", "y", "x"),
    perturbation_id = c("p1", "p2", "p1", "p1"),
    effect = c(1, 0, -1, 0)
  )
  first <- fsf_signal_identity(effects, group_cols = "condition")
  second <- fsf_signal_identity(effects, group_cols = "condition")
  expect_identical(first, second)
  expect_identical(first$condition, c("B", "A", "A"))
  expect_identical(first$feature_id, c("z", "y", "x"))

  classified <- data.frame(
    condition = c("B", "A", "B"),
    feature_id = c("z", "y", "x"),
    signal_class = c("Stable Up", "Low Stability", "Highly Stable Down")
  )
  architecture_first <- fsf_architecture(classified)
  architecture_second <- fsf_architecture(classified)
  expect_identical(architecture_first, architecture_second)
  expect_identical(
    as.character(unique(architecture_first$condition)), c("B", "A")
  )
})

test_that("fsf_classify remains row-wise for already aggregated identities", {
  duplicated_identity <- data.frame(
    feature_id = c("same", "same"),
    p_up = c(0.6, 0.2),
    p_down = c(0.2, 0.6),
    p_const = c(0.2, 0.2)
  )
  observed <- fsf_classify(duplicated_identity)
  expect_equal(nrow(observed), 2L)
  expect_equal(as.character(observed$signal_class),
               c("Transitional Up", "Transitional Down"))
})
