# Binary classification ----------------------------------------------------------
test_that("Binary: All numeric", {
  bin_fit <- rpf(am ~ hp + disp, data = mtcars)

  expect_s3_class(bin_fit, "rpf")
  expect_s4_class(bin_fit$fit, "Rcpp_ClassificationRPF")
})

test_that("Binary: All numeric, logit loss", {
  # Currently prints "test" (3x) in console, not sure how to capture
  # Hopefully capturable with switch to Rcout?
  expect_silent(rpf(am ~ hp + disp, data = mtcars, loss = "logit"))

  bin_fit <- rpf(am ~ hp + disp, data = mtcars, loss = "logit")

  expect_s3_class(bin_fit, "rpf")
  expect_s4_class(bin_fit$fit, "Rcpp_ClassificationRPF")
})

test_that("Binary: All numeric, class predictions", {
  bin_fit <- rpf(am ~ hp + disp, data = mtcars)

  bin_pred_class <- predict(bin_fit, mtcars[c("hp", "disp")], type = "class")
})

test_that("Binary: All numeric, probability predictions", {
  bin_fit <- rpf(am ~ hp + disp, data = mtcars)

  bin_pred_prob <- predict(bin_fit, mtcars[c("hp", "disp")], type = "prob")
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
