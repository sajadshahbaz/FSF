test_that("an exact 0.50 tie is non-directional Low Stability", {
  identity <- data.frame(
    feature_id = "tie",
    p_up = 0.50,
    p_down = 0.50,
    p_const = 0
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, "tied")
  expect_equal(observed$ssi, 0.50)
  expect_equal(as.character(observed$stability_region), "Low Stability")
  expect_equal(as.character(observed$signal_class), "Low Stability")
})

test_that("an exact 0.50 unique plurality remains Low Stability", {
  identity <- data.frame(
    feature_id = "plurality",
    p_up = 0.50,
    p_down = 0.30,
    p_const = 0.20
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, "up")
  expect_equal(observed$ssi, 0.50)
  expect_equal(as.character(observed$stability_region), "Low Stability")
  expect_equal(as.character(observed$signal_class), "Low Stability")
})

test_that("Low Stability remains non-directional for every dominant-state value", {
  identity <- data.frame(
    feature_id = c("up", "down", "constant", "tied"),
    p_up = c(0.50, 0.30, 0.30, 0.50),
    p_down = c(0.30, 0.50, 0.20, 0.50),
    p_const = c(0.20, 0.20, 0.50, 0)
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, c("up", "down", "constant", "tied"))
  expect_equal(as.character(observed$stability_region), rep("Low Stability", 4L))
  expect_equal(as.character(observed$signal_class), rep("Low Stability", 4L))
})

test_that("strict majority produces a unique Transitional class", {
  identity <- data.frame(
    feature_id = "majority",
    p_up = 0.51,
    p_down = 0.49,
    p_const = 0
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, "up")
  expect_equal(as.character(observed$stability_region), "Transitional")
  expect_equal(as.character(observed$signal_class), "Transitional Up")
})


test_that("raw near-half probabilities do not create a fuzzy tie", {
  identity <- data.frame(
    feature_id = "raw_majority",
    p_up = 0.5000000001,
    p_down = 0.4999999999,
    p_const = 0
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, "up")
  expect_equal(observed$ssi, 0.5000000001)
  expect_equal(as.character(observed$stability_region), "Transitional")
  expect_equal(as.character(observed$signal_class), "Transitional Up")
})
