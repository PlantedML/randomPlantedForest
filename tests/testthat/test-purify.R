set.seed(124)
xdat <- data.frame(
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = cut(runif(100), 3, labels = 1:3),
  x4 = cut(runif(100), 2, labels = 1:2),
  x5 = factor(sample(letters[1:4], 100, replace = TRUE))
)

xdat$yfact <- factor(rbinom(100, 1, prob = plogis(xdat$x1 + xdat$x2)), labels = c("neg", "pos"))

test_that("purify() and is_purified() work", {
  bin_fit <- rpf(yfact ~ x1 + x2 + x3 + x4 + x5, data = xdat, max_interaction = 3)

  expect_false(is_purified(bin_fit))
  expect_identical(is_purified(bin_fit), bin_fit$fit$is_purified())

  purify(bin_fit)

  expect_true(is_purified(bin_fit))
  expect_identical(is_purified(bin_fit), bin_fit$fit$is_purified())
})

test_that("purification does not alter predictions (null effect)", {
  bin_fit <- rpf(yfact ~ x1 + x2 + x3 + x4 + x5, data = xdat, max_interaction = 3)
  pred_pre <- predict(bin_fit, new_data = xdat, type = "numeric")

  expect_false(is_purified(bin_fit))
  purify(bin_fit)
  expect_true(is_purified(bin_fit))

  pred_post <- predict(bin_fit, new_data = xdat, type = "numeric")

  expect_equal(pred_pre, pred_post, tolerance = 1e-15)
})

test_that("purification does not alter predictions (with effect)", {
  train <-  mtcars[1:20, 1:4]
  test <-  mtcars[21:32, 1:4]

  set.seed(23)
  rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, ntrees = 30)
  pred_pre <- predict(rpfit, test)

  # Not purified
  expect_false(is_purified(rpfit))

  purify(rpfit)

  # Purified
  expect_true(is_purified(rpfit))

  pred_post <- predict(rpfit, test)

  expect_equal(pred_pre, pred_post, tolerance = 1e-15)
})
