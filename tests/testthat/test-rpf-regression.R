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
