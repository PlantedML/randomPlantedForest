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

test_that("Setting seed in R works", {
  set.seed(13)
  rpf_fit1 <- rpf(mpg ~ wt + cyl, data = mtcars)
  pred1 <- predict(rpf_fit1, mtcars[1:5, ])
  
  set.seed(13)
  rpf_fit2 <- rpf(mpg ~ wt + cyl, data = mtcars)
  pred2 <- predict(rpf_fit2, mtcars[1:5, ])
  
  # No seed set, so should be different
  rpf_fit3 <- rpf(mpg ~ wt + cyl, data = mtcars)
  pred3 <- predict(rpf_fit3, mtcars[1:5, ])
  
  expect_equal(pred1, pred2)
  expect_failure(expect_equal(pred1, pred3))
})

test_that("Rcpp RNG/R RNG interference", {
  set.seed(1)
  r11 <- runif(1)
  r12 <- runif(1)
  
  set.seed(1)
  r21 <- runif(1)
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)
  r22 <- runif(1)
  
  # If this fails R is broken
  expect_equal(r11, r21)
  # If this fails Rcpp does not properly reset R's RNG
  expect_equal(r12, r22) 
})
