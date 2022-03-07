# Regression -----------------------------------------------------------------
test_that("Regression fit", {

  # XY: data.frame / matrix
  pred_df <-  mtcars[, c(2, 6)]
  rpf_fit1 <- rpf(x = pred_df, y = mtcars$mpg)
  rpf_fit2 <- rpf(x = as.matrix(pred_df), y = mtcars$mpg)

  # Formula interface
  rpf_fit3 <- rpf(mpg ~ wt + cyl, data = mtcars)

  expect_s3_class(rpf_fit1, "rpf")
  expect_s3_class(rpf_fit2, "rpf")
  expect_s3_class(rpf_fit3, "rpf")
  expect_s4_class(rpf_fit1$fit, "Rcpp_RandomPlantedForest")
  expect_s4_class(rpf_fit2$fit, "Rcpp_RandomPlantedForest")
  expect_s4_class(rpf_fit3$fit, "Rcpp_RandomPlantedForest")
})

test_that("Regression prediction", {
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)

  # check default behavior: components = 0
  pred_default <- predict(rpf_fit, mtcars[, c(2, 6)], type = "numeric")
  pred_0 <- predict(rpf_fit, mtcars[, c(2, 6)], type = "numeric", components = 0)

  expect_identical(pred_default, pred_0)
  expect_s3_class(pred_default, "tbl_df")

  # FIXME: Test components
  pred_1 <- predict(rpf_fit, mtcars[, c(2, 6)], type = "numeric", components = c(0, 1))

})


# Classification ----------------------------------------------------------
