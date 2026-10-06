# Regression -----------------------------------------------------------------
test_that("Basic fit: All numeric", {
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})

test_that("Prediction: All numeric", {
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)

  pred <- predict(rpf_fit, mtcars[, c(2, 6)], type = "numeric")

  expect_s3_class(pred, "tbl_df")
  expect_named(pred, ".pred")
  expect_equal(nrow(pred), nrow(mtcars))
})

test_that("Fit + predict: Categorical features", {
  mtcars_cat <- mtcars
  mtcars_cat$cyl <- factor(mtcars$cyl)

  # Coercible to integer
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars_cat)
  pred <- predict(rpf_fit, mtcars_cat[, c(2, 6)], type = "numeric")

  expect_s3_class(pred, "tbl_df")

  # Not coercible to integer
  mtcars_cat$wt_cat <- ifelse(mtcars$wt > 3.2, "heavy", "light")
  rpf_fit <- rpf(mpg ~ wt_cat, data = mtcars_cat)

  pred <- predict(rpf_fit, mtcars_cat[, c("wt_cat"), drop = FALSE], type = "numeric")

  expect_s3_class(pred, "tbl_df")
})

test_that("Warn for y = 0,1", {
  xdat <- data.frame(
    y01 = sample(c(0L, 1L), 100, replace = TRUE),
    x1 = rnorm(100),
    x2 = rnorm(100)
  )
  bin_fit <- suppressWarnings(rpf(y01 ~ x1 + x2, data = xdat, loss = "L2"))

  expect_warning(
    predict(bin_fit, new_data = xdat, type = "class"),
    regexp = "^Only predict type 'numeric' supported for regression"
  )

  expect_warning(
    predict(bin_fit, new_data = xdat, type = "link"),
    regexp = "^Only predict type 'numeric' supported for regression"
  )
})

test_that("predictions do not depend on nthreads", {
  fit <- rpf(mpg ~ cyl + wt + hp, data = mtcars, max_interaction = 2, nthreads = 2)
  expect_identical(
    predict(fit, mtcars, nthreads = 1L),
    predict(fit, mtcars, nthreads = 3L)
  )
  expect_identical(predict(fit, mtcars), predict(fit, mtcars, nthreads = 1L))
})

test_that("predict() rejects unknown arguments", {
  fit <- rpf(mpg ~ cyl + wt, data = mtcars, ntrees = 5)
  expect_error(predict(fit, mtcars, newdata = mtcars), class = "rlib_error_dots_nonempty")
})
