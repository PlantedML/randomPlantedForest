# Using iris is suboptimal but it's easy and comes with 3-class y
# penguins would require package in Suggests and na.omit() preprocessing
test_that("Multiclass: All numeric", {
  classif_fit <- rpf(Species ~ ., data = iris)

  expect_s3_class(classif_fit, "rpf")
  expect_s4_class(classif_fit$fit, "Rcpp_ClassificationRPF")
})

test_that("Multiclass: All numeric, class predictions", {
  classif_fit <- rpf(Species ~ ., data = iris)

  classif_pred_class <- predict(classif_fit, iris, type = "class")
  expect_equal(dim(classif_pred_class), c(nrow(iris), 1))
  expect_equal(levels(classif_pred_class$.pred_class), levels(iris$Species))
})

test_that("Multiclass: All numeric, probability predictions", {
  classif_fit <- rpf(Species ~ ., data = iris)

  classif_pred_prob <- predict(classif_fit, iris, type = "prob")
  expect_equal(dim(classif_pred_prob), c(nrow(iris), nlevels(iris$Species)))
  expect_gte(min(classif_pred_prob), 0)
  expect_lte(min(classif_pred_prob), 1)
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
  expect_error(rpf(ychar ~ x1 + x2, xdat))
})
