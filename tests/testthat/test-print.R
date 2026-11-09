test_that("print() summarises a regression fit and its purification state", {
  fit <- rpf_050(mpg ~ cyl + wt, data = mtcars, n_families = 3)
  expect_snapshot(print(fit))

  purify(fit)
  expect_snapshot(print(fit))
})

test_that("print() covers x/y fits, classification losses and deterministic fits", {
  fit_logit <- rpf_050(x = iris[, 1:4], y = iris$Species, n_families = 1, loss = "logit", deterministic = TRUE)
  expect_snapshot(print(fit_logit))

  fit_l2 <- rpf_050(Species ~ ., data = iris, n_families = 2, max_interaction = 1, loss = "L2")
  expect_snapshot(print(fit_l2))
})

test_that("print() returns its input invisibly and matches format()", {
  fit <- rpf_050(mpg ~ cyl + wt, data = mtcars, n_families = 2)

  output <- capture.output(result <- withVisible(print(fit)))

  expect_false(result$visible)
  expect_identical(result$value, fit)
  expect_identical(output, format(fit))
})

test_that("an exported forest prints compactly", {
  fit <- rpf_050(mpg ~ cyl + wt, data = mtcars, n_families = 3, export_forest = TRUE)

  expect_snapshot({
    print(fit$forest)
    str(fit$forest)
  })
})

test_that("print() counts the fit's own columns, not the formula's raw ones", {
  set.seed(4)
  d <- data.frame(x1 = runif(60, 0.1, 1), x2 = runif(60), x3 = runif(60))
  d$y <- log(d$x1) + d$x2 + rnorm(60, sd = 0.1)

  # y ~ . - x3: two columns, all interactions
  fit <- rpf_050(y ~ . - x3, data = d, max_interaction = 2, ntrees = 2, splits = 5)
  out <- paste(format(fit), collapse = "\n")
  expect_match(out, "on 2 predictors, all possible interactions", fixed = TRUE)

  # log(x1) + I(x2^2) + x2 + x3: four columns from three raw ones
  fit <- rpf_050(y ~ log(x1) + I(x2^2) + x2 + x3, data = d, max_interaction = 4, ntrees = 2, splits = 5)
  out <- paste(format(fit), collapse = "\n")
  expect_match(out, "on 4 predictors, all possible interactions", fixed = TRUE)
})
