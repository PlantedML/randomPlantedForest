test_that("print() summarises a regression fit and its purification state", {
  fit <- rpf(mpg ~ cyl + wt, data = mtcars, ntrees = 3)
  expect_snapshot(print(fit))

  purify(fit)
  expect_snapshot(print(fit))
})

test_that("print() covers x/y fits, classification losses and deterministic fits", {
  fit_logit <- rpf(x = iris[, 1:4], y = iris$Species, ntrees = 1, loss = "logit", deterministic = TRUE)
  expect_snapshot(print(fit_logit))

  fit_l2 <- rpf(Species ~ ., data = iris, ntrees = 2)
  expect_snapshot(print(fit_l2))
})

test_that("print() returns its input invisibly and matches format()", {
  fit <- rpf(mpg ~ cyl + wt, data = mtcars, ntrees = 2)

  output <- capture.output(result <- withVisible(print(fit)))

  expect_false(result$visible)
  expect_identical(result$value, fit)
  expect_identical(output, format(fit))
})

test_that("an exported forest prints compactly", {
  fit <- rpf(mpg ~ cyl + wt, data = mtcars, ntrees = 3, export_forest = TRUE)

  expect_snapshot({
    print(fit$forest)
    str(fit$forest)
  })
})
