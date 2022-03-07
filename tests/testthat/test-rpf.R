test_that("Formula interface", {
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})

test_that("XY data.frame interface", {

  pred_df <-  mtcars[, c(2, 6)]
  rpf_fit <- rpf(x = pred_df, y = mtcars$mpg)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})

test_that("XY matrix interface", {

  pred_mat <-  as.matrix(mtcars[, c(2, 6)])
  rpf_fit <- rpf(x = pred_mat, y = mtcars$mpg)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})
