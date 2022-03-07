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
