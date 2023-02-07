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

test_that("extract_components returns requested predictors", {
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

test_that("extract_component is consistent with predictor order", {
  rp <- rpf(mpg ~ ., data = mtcars, max_interaction = 2, purify = TRUE)

  expect_equal(
    extract_component(rp, mtcars, c("cyl", "am")),
    extract_component(rp, mtcars, c("am", "cyl"))
  )
})
