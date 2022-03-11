# Binary classification ----------------------------------------------------------
test_that("Binary: All numeric", {
  bin_fit <- rpf(am ~ hp + disp, data = mtcars)

  expect_s3_class(bin_fit, "rpf")
  expect_s4_class(bin_fit$fit, "Rcpp_ClassificationRPF")
})

test_that("Binary: All numeric, logit loss", {
  expect_silent(rpf(am ~ hp + disp, data = mtcars, loss = "logit"))

  bin_fit <- rpf(am ~ hp + disp, data = mtcars, loss = "logit")

  expect_s3_class(bin_fit, "rpf")
  expect_s4_class(bin_fit$fit, "Rcpp_ClassificationRPF")
})


test_that("Binary: Detection works", {
  xdat <- data.frame(
    y01 = sample(c(0, 1), 100, replace = TRUE),
    yfact = factor(sample(c("pos", "neg"), 100, replace = TRUE)),
    ychar = sample(c("pos", "neg"), 100, replace = TRUE),
    ylogi = sample(c(TRUE, FALSE), 100, replace = TRUE),
    x1 = rnorm(100),
    x2 = rnorm(100)
  )

  # y in 0, 1
  y_01 <- rpf(y01 ~ x1 + x2, xdat)
  expect_s4_class(y_01$fit, "Rcpp_ClassificationRPF")

  # FIXME: decide factor order assumption?
  # y two-level factor: should work assuming which level is positive?
  y_fact <- rpf(yfact ~ x1 + x2, xdat)
  expect_s4_class(y_fact$fit, "Rcpp_ClassificationRPF")

  # y two-level character: should fail because ambiguous
  # similar problem as factor but w/o levels no order can be assumed
  expect_error(rpf(ychar ~ x1 + x2, xdat))

  # y logical: should warn but work, user should supply numerical/factor
  y_logi <- expect_warning(rpf(ylogi ~ x1 + x2, xdat))
  expect_s4_class(y_logi$fit, "Rcpp_ClassificationRPF")
})

# Multiclass classification -----------------------------------------------
# Using iris is suboptimal but it's easy and comes with 3-class y
# penguins would require package in Suggests and na.omit() preprocessing
test_that("Multiclass: All numeric", {
  classif_fit <- rpf(Species ~ ., data = iris)

  expect_s3_class(classif_fit, "rpf")
  expect_s4_class(classif_fit$fit, "Rcpp_ClassificationRPF")
})

test_that("Multiclass: All numeric, class predictions", {
  classif_fit <- rpf(Species ~ ., data = iris)

  classif_pred_class <- predict(classif_fit, iris[, -5], type = "class")
})

test_that("Multiclass: All numeric, probability predictions", {
  classif_fit <- rpf(Species ~ ., data = iris)

  classif_pred_prob <- predict(classif_fit, iris[, -5], type = "prob")
})

test_that("Multiclass: Detection works", {
  xdat <- data.frame(
    yint = sample(c(0L, 1L, 2L), 100, replace = TRUE),
    yfact = factor(sample(c("hi", "mid", "lo"), 100, replace = TRUE)),
    ychar = sample(c("hi", "mid", "lo"), 100, replace = TRUE),
    x1 = rnorm(100),
    x2 = rnorm(100)
  )

  # y is integer: should _not_ be treated as classif task
  y_int <- rpf(yint ~ x1 + x2, xdat)
  expect_failure(expect_s4_class(y_int$fit, "Rcpp_ClassificationRPF"))

  # y 3-level factor
  y_fact <- rpf(yfact ~ x1 + x2, xdat)
  expect_s4_class(y_fact$fit, "Rcpp_ClassificationRPF")

  # y 3-level character
  y_char <- rpf(ychar ~ x1 + x2, xdat)
  expect_s4_class(y_char$fit, "Rcpp_ClassificationRPF")
})
