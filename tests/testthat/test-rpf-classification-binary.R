xdat <- data.frame(
  y01 = sample(c(0L, 1L), 100, replace = TRUE),
  y12 = sample(c(1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("pos", "neg"), 100, replace = TRUE)),
  ychar = sample(c("pos", "neg"), 100, replace = TRUE),
  ylogi = sample(c(TRUE, FALSE), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100)
)

# Basic model creation ----------------------------------------------------
test_that("Binary: All numeric", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat)

  expect_s3_class(bin_fit, "rpf")
  expect_s4_class(bin_fit$fit, "Rcpp_ClassificationRPF")
})

test_that("Binary: All numeric, logit loss", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "logit")

  expect_s3_class(bin_fit, "rpf")
  expect_s4_class(bin_fit$fit, "Rcpp_ClassificationRPF")
  expect_identical(bin_fit$loss, "logit")
})

test_that("Binary: All numeric, exponential loss", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "exponential")
  
  expect_s3_class(bin_fit, "rpf")
  expect_s4_class(bin_fit$fit, "Rcpp_ClassificationRPF")
  expect_identical(bin_fit$loss, "exponential")
})


# Binary task detection ---------------------------------------------------

test_that("Binary detection: factor", {
  y_fact <- rpf(yfact ~ x1 + x2, xdat)
  expect_s4_class(y_fact$fit, "Rcpp_ClassificationRPF")
})

# Ambiguous, use regression and warn
test_that("Binary detection: Regression for 0,1", {
  # y in 0, 1: Ambiguous, expect warning, but should classifybrows
  expect_warning(rpf(y01 ~ x1 + x2, xdat), regexp = "^y is.*assuming")
  y_01 <- suppressWarnings(rpf(y01 ~ x1 + x2, xdat))
  expect_s4_class(y_01$fit, "Rcpp_RandomPlantedForest")
})

test_that("Binary detection: Regression for 1,2", {
  # y in 1, 2: See 0,1
  expect_warning(rpf(y12 ~ x1 + x2, xdat), regexp = "^y is.*assuming")
  y_12 <- suppressWarnings(rpf(y12 ~ x1 + x2, xdat))
  expect_s4_class(y_12$fit, "Rcpp_RandomPlantedForest")
})

test_that("Binary detection: Fail for character, logical", {
  # y two-level character: should fail because ambiguous
  # similar problem as factor but w/o levels no order can be assumed
  expect_error(rpf(ychar ~ x1 + x2, xdat), regexp = "^y should be")

  # y logical: should error and note what it expects
  expect_error(rpf(ylogi ~ x1 + x2, xdat), regexp = "^y should be")
})
