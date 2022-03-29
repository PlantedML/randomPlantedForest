xdat <- data.frame(
  yint = sample(c(0L, 1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("hi", "mid", "lo"), 100, replace = TRUE)),
  ychar = sample(c("hi", "mid", "lo"), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100)
)

# Basic model creation ----------------------------------------------------
test_that("Multiclass: All numeric", {
  classif_fit <- rpf(yfact ~ ., data = xdat)
  
  expect_s3_class(classif_fit, "rpf")
  expect_s4_class(classif_fit$fit, "Rcpp_ClassificationRPF")
})

# Classif task detection ---------------------------------------------------
test_that("Multiclass: Detection works", {
  # y 3-level factor
  y_fact <- rpf(yfact ~ x1 + x2, xdat)
  expect_s4_class(y_fact$fit, "Rcpp_ClassificationRPF")
  
  # y is integer: should _not_ be treated as classif task
  y_int <- rpf(yint ~ x1 + x2, xdat)
  expect_failure(expect_s4_class(y_int$fit, "Rcpp_ClassificationRPF"))
  
  # y 3-level character
  expect_error(rpf(ychar ~ x1 + x2, xdat), regexp = "^y should be")
})