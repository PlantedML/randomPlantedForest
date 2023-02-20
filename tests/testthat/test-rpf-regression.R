# Regression -----------------------------------------------------------------
test_that("Basic fit: All numeric", {
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})

test_that("Prediction: All numeric", {
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)

  # check default behavior: components = 0
  pred_default <- predict(rpf_fit, mtcars[, c(2, 6)], type = "numeric")
  pred_0 <- predict(rpf_fit, mtcars[, c(2, 6)], type = "numeric", components = 0)

  expect_identical(pred_default, pred_0)
  expect_s3_class(pred_default, "tbl_df")

  # FIXME: Test components
  pred_1 <- predict(rpf_fit, mtcars[, c(2, 6)], type = "numeric", components = c(0, 1))
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

# FIXME: Since L1 loss not supported in regression, this test is somewhat weird
test_that("Setting L1 or L2 loss is saved in rpf object", {
  fit_l1 <- rpf(mpg ~ wt + cyl, data = mtcars, loss = "L1")
  fit_l2 <- rpf(mpg ~ wt + cyl, data = mtcars, loss = "L2")
  fit_default <- rpf(mpg ~ wt + cyl, data = mtcars)

  expect_identical(fit_l1$params$loss, "L1")
  expect_identical(fit_l2$params$loss, "L2")
  expect_identical(fit_default$params$loss, "L2")
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

