test_that("bundle round-trip", {
  skip_if_not_installed("bundle")
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5)
  b <- bundle::bundle(fit)
  tmp <- tempfile(fileext = ".rds")
  saveRDS(b, tmp)
  restored <- bundle::unbundle(readRDS(tmp))
  expect_identical(predict(restored, mtcars), predict(fit, mtcars))
})
