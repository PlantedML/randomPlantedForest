test_that("extract_components returns correct structure", {
  rp <- rpf(mpg ~ ., data = mtcars, max_interaction = 2)

  m <- extract_components(rp, mtcars)

  expect_named(m, c("intercept", "cyl", "disp", "hp", "drat", "wt", "qsec", "vs",
                    "am", "gear", "carb", "cyl:disp", "cyl:hp", "cyl:drat", "cyl:wt",
                    "cyl:qsec", "cyl:vs", "am:cyl", "cyl:gear", "carb:cyl", "disp:hp",
                    "disp:drat", "disp:wt", "disp:qsec", "disp:vs", "am:disp", "disp:gear",
                    "carb:disp", "drat:hp", "hp:wt", "hp:qsec", "hp:vs", "am:hp",
                    "gear:hp", "carb:hp", "drat:wt", "drat:qsec", "drat:vs", "am:drat",
                    "drat:gear", "carb:drat", "qsec:wt", "vs:wt", "am:wt", "gear:wt",
                    "carb:wt", "qsec:vs", "am:qsec", "gear:qsec", "carb:qsec", "am:vs",
                    "gear:vs", "carb:vs", "am:gear", "am:carb", "carb:gear"))

  expect_s3_class(m, "data.frame")
  expect_equal(nrow(mtcars), nrow(m))

})

test_that("extract_components returns requested predictors, in order", {
  rp <- rpf(mpg ~ ., data = mtcars, max_interaction = 2)

  m <- extract_components(rp, mtcars, predictors = c("cyl", "am", "vs"))

  expect_named(m, c("intercept", "cyl", "am", "vs", "am:cyl", "cyl:vs", "am:vs"))

  expect_s3_class(m, "data.frame")
  expect_length(m, 7)
  expect_equal(nrow(mtcars), nrow(m))

})

test_that("extract_components purifies if needed", {
  rp <- rpf(mpg ~ ., data = mtcars[, c(1:4)], max_interaction = 2)

  expect_false(is_purified(rp))

  extract_components(rp, mtcars)

  expect_true(is_purified(rp))
})


test_that("extract_components sums to prediction", {
  rp <- rpf(mpg ~ cyl + am + gear, data = mtcars, max_interaction = 3, purify = TRUE)

  m <- extract_components(rp, mtcars)

  expect_equal(
    predict(rp, mtcars)[[1]],
    rowSums(m)
  )

})

test_that("extract_component is consistent with predictor order", {
  rp <- rpf(mpg ~ cyl + am + gear, data = mtcars, max_interaction = 3, purify = TRUE)

  # Internal data preprocessing only done in extract_components to save time
  processed <- hardhat::forge(mtcars, rp$blueprint)
  new_data <- preprocess_predictors_predict(rp, processed$predictors)

  expect_equal(
    extract_component(rp, new_data, c("cyl", "am")),
    extract_component(rp, new_data, c("am", "cyl"))
  )

  expect_equal(
    extract_component(rp, new_data, c("cyl", "am", "gear")),
    extract_component(rp, new_data, c("gear", "am", "cyl"))
  )

  expect_equal(
    extract_component(rp, new_data, c("am", "cyl", "gear")),
    extract_component(rp, new_data, c("gear", "am", "cyl"))
  )
})
