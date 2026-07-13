test_that("predict_components returns correct structure", {
  rp <- rpf(mpg ~ ., data = mtcars, max_interaction = 2)

  m <- predict_components(rp, mtcars)

  expect_named(m$m, c("cyl", "disp", "hp", "drat", "wt", "qsec", "vs",
                    "am", "gear", "carb", "cyl:disp", "cyl:hp", "cyl:drat", "cyl:wt",
                    "cyl:qsec", "cyl:vs", "am:cyl", "cyl:gear", "carb:cyl", "disp:hp",
                    "disp:drat", "disp:wt", "disp:qsec", "disp:vs", "am:disp", "disp:gear",
                    "carb:disp", "drat:hp", "hp:wt", "hp:qsec", "hp:vs", "am:hp",
                    "gear:hp", "carb:hp", "drat:wt", "drat:qsec", "drat:vs", "am:drat",
                    "drat:gear", "carb:drat", "qsec:wt", "vs:wt", "am:wt", "gear:wt",
                    "carb:wt", "qsec:vs", "am:qsec", "gear:qsec", "carb:qsec", "am:vs",
                    "gear:vs", "carb:vs", "am:gear", "am:carb", "carb:gear"))

  expect_s3_class(m$m, "data.frame")
  expect_equal(nrow(mtcars), nrow(m$m))

})

test_that("predict_components returns requested predictors, in order", {
  rp <- rpf(mpg ~ ., data = mtcars, max_interaction = 2)

  m <- predict_components(rp, mtcars, predictors = c("cyl", "am", "vs"))

  expect_named(m$m, c("cyl", "am", "vs", "am:cyl", "cyl:vs", "am:vs"))

  expect_s3_class(m$m, "data.frame")
  expect_length(m$m, 6)
  expect_equal(nrow(mtcars), nrow(m$m))

})

test_that("predict_components purifies if needed", {
  rp <- rpf(mpg ~ ., data = mtcars[, c(1:4)], max_interaction = 2)

  expect_false(is_purified(rp))

  predict_components(rp, mtcars)

  expect_true(is_purified(rp))
})


test_that("predict_components sums to prediction", {
  rp <- rpf(mpg ~ cyl + am + gear, data = mtcars, max_interaction = 3, purify = TRUE)

  m <- predict_components(rp, mtcars)

  expect_equal(
    predict(rp, mtcars)[[1]],
    rowSums(m$m) + m$intercept
  )

})

test_that("predict_components errors appropriately", {
  rp <- rpf(mpg ~ cyl + am + gear, data = mtcars, max_interaction = 3)

  # Unknown predictors requested
  expect_error(predict_components(rp, mtcars, predictors = c("x1", "x2")))

  # Predictors not in fit on
  expect_error(predict_components(rp, mtcars, predictors = c("cyl", "disp")))

  # new_data does not contain all original predictors
  expect_error(predict_components(rp, mtcars[, c("cyl", "am")]))

  # max_interaction higher than original fit
  expect_error(predict_components(rp, mtcars, max_interaction = 30))

  # max_interaction implausible
  expect_error(predict_components(rp, mtcars, max_interaction = -1))

  # predictors is wrong type
  expect_error(predict_components(rp, mtcars, predictors = list("cyl")))
})

test_that(".predict_single_component is consistent with predictor order", {
  rp <- rpf(mpg ~ cyl + am + gear, data = mtcars, max_interaction = 3, purify = TRUE)

  # Internal data preprocessing only done in predict_components to save time
  processed <- hardhat::forge(mtcars, rp$blueprint)
  new_data <- preprocess_predictors_predict(rp, processed$predictors)

  expect_equal(
    .predict_single_component(rp, new_data, c("cyl", "am")),
    .predict_single_component(rp, new_data, c("am", "cyl"))
  )

  expect_equal(
    Reduce(expect_equal, list(
      .predict_single_component(rp, new_data, c("cyl", "am", "gear")),
      .predict_single_component(rp, new_data, c("gear", "am", "cyl")),
      .predict_single_component(rp, new_data, c("am", "cyl", "gear")),
      .predict_single_component(rp, new_data, c("gear", "cyl", "am")),
      .predict_single_component(rp, new_data, c("am", "gear", "cyl")),
      .predict_single_component(rp, new_data, c("cyl", "gear", "am"))
    )),
    .predict_single_component(rp, new_data, c("cyl", "gear", "am"))
  )

})
