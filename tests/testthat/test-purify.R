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
  bin_fit <- rpf_050(yfact ~ x1 + x2 + x3 + x4 + x5, data = xdat, max_interaction = 3)

  expect_false(is_purified(bin_fit))
  expect_identical(is_purified(bin_fit), bin_fit$fit$is_purified())

  purify(bin_fit)

  expect_true(is_purified(bin_fit))
  expect_identical(is_purified(bin_fit), bin_fit$fit$is_purified())
})

test_that("purification does not alter predictions (null effect)", {
  bin_fit <- rpf_050(yfact ~ x1 + x2 + x3 + x4 + x5, data = xdat, max_interaction = 3)
  pred_pre <- predict(bin_fit, new_data = xdat, type = "numeric")

  expect_false(is_purified(bin_fit))
  purify(bin_fit)
  expect_true(is_purified(bin_fit))

  pred_post <- predict(bin_fit, new_data = xdat, type = "numeric")

  expect_equal(pred_pre, pred_post, tolerance = 1e-10)
})

test_that("purification does not alter predictions (with effect)", {
  train <- mtcars[1:20, 1:4]
  test <- mtcars[21:32, 1:4]

  set.seed(23)
  rpfit <- rpf_050(mpg ~ ., data = train, max_interaction = 3, ntrees = 30)
  pred_pre <- predict(rpfit, test)

  # Not purified
  expect_false(is_purified(rpfit))

  purify(rpfit)

  # Purified
  expect_true(is_purified(rpfit))

  pred_post <- predict(rpfit, test)

  expect_equal(pred_pre, pred_post, tolerance = 1e-10)
})

test_that("purify() purifies in place, returns invisibly and rejects unknown arguments", {
  rp <- rpf_050(mpg ~ cyl + wt + hp, data = mtcars, max_interaction = 2, ntrees = 5)
  expect_false(is_purified(rp))
  expect_invisible(purify(rp))
  expect_true(is_purified(rp))
  expect_error(purify(rp, nthread = 2), class = "rlib_error_dots_nonempty")
})

test_that("purify(maxp_interaction = k) keeps predict(), and predict_components() sums to it", {
  set.seed(3)
  n <- 400
  X <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n))
  y <- 4 * X$x1 * X$x2 * X$x3 + X$x4 + rnorm(n, sd = 0.3)
  Xt <- X[1:100, ]
  for (mode in 1:2) {
    set.seed(1)
    fit <- rpf_050(x = X, y = y, max_interaction = 4, ntrees = 10, splits = 60)
    p0 <- predict(fit, Xt)
    purify(fit, maxp_interaction = 2, mode = mode)
    # the whole forest's prediction, not the intercept plus the components of order <= 2
    expect_identical(predict(fit, Xt), p0)
    # the components of order 3 and 4 are zero columns; their sum is the remainder
    pc <- predict_components(fit, Xt)
    expect_false(is.null(pc$remainder))
    expect_equal(rowSums(pc$m) + pc$intercept + pc$remainder, p0$.pred, tolerance = 1e-12)
    expect_gt(max(abs(pc$remainder)), 0.01)
    expect_true(fit$fit$is_purify_capped())
    # purifying all orders again clears the cap
    purify(fit, mode = mode)
    expect_false(fit$fit$is_purify_capped())
    expect_equal(predict(fit, Xt)$.pred, p0$.pred, tolerance = 1e-12)
    expect_null(predict_components(fit, Xt)$remainder)
  }
})

test_that("purify(maxp_interaction = k) keeps predict() for multiclass classification", {
  set.seed(5)
  fit <- rpf_050(Species ~ ., data = iris, max_interaction = 3, ntrees = 10, loss = "L2")
  p0 <- predict(fit, iris, type = "numeric")
  purify(fit, maxp_interaction = 1)
  expect_identical(predict(fit, iris, type = "numeric"), p0)
  pc <- predict_components(fit, iris)
  for (level in pc$target_levels) {
    cols <- endsWith(names(pc$m), paste0("__class:", level))
    expect_equal(
      rowSums(pc$m[, cols, with = FALSE]) + pc$intercept[[level]] + pc$remainder[[level]],
      p0[[paste0(".pred_", level)]],
      tolerance = 1e-12
    )
  }
})

test_that("purify() above an earlier cap purifies again, as a fresh purification would", {
  set.seed(3)
  n <- 300
  X <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  y <- X$x1 * X$x2 + X$x2 * X$x3 + X$x1 + rnorm(n, sd = 0.2)
  for (mode in 1:2) {
    set.seed(1)
    fit <- rpf(x = X, y = y, max_interaction = 3, ntrees = 5, splits = 40)
    set.seed(1)
    fresh <- rpf(x = X, y = y, max_interaction = 3, ntrees = 5, splits = 40)
    purify(fresh, mode = mode)

    purify(fit, maxp_interaction = 1, mode = mode)
    expect_identical(fit$fit$get_purify_cap(), 1L)
    # at or below the cap: returned as it is
    purify(fit, maxp_interaction = 1, mode = mode)
    expect_identical(fit$fit$get_purify_cap(), 1L)
    # above it: purified again, up to the order asked for
    purify(fit, maxp_interaction = 2, mode = mode)
    expect_identical(fit$fit$get_purify_cap(), 2L)
    pc2 <- predict_components(fit, X)
    expect_true(any(pc2$m[["x1:x2"]] != 0))
    expect_true(all(pc2$m[["x1:x2:x3"]] == 0))
    # all orders (the default), as a fresh purification of all orders
    purify(fit, mode = mode)
    expect_identical(fit$fit$get_purify_cap(), 0L)
    expect_false(fit$fit$is_purify_capped())
    pc <- predict_components(fit, X)
    expect_identical(pc, predict_components(fresh, X))
    expect_null(pc$remainder)
    expect_identical(predict(fit, X), predict(fresh, X))
    # a full purification covers every cap
    purify(fit, maxp_interaction = 1, mode = mode)
    expect_identical(fit$fit$get_purify_cap(), 0L)
  }
})

test_that("purify() above an earlier cap purifies again for multiclass classification", {
  set.seed(5)
  fit <- rpf(Species ~ ., data = iris, max_interaction = 2, ntrees = 5, loss = "L2")
  set.seed(5)
  fresh <- rpf(Species ~ ., data = iris, max_interaction = 2, ntrees = 5, loss = "L2")
  purify(fit, maxp_interaction = 1)
  purify(fit)
  purify(fresh)
  expect_false(fit$fit$is_purify_capped())
  expect_identical(predict_components(fit, iris), predict_components(fresh, iris))
})

test_that("purify() checks its arguments on a purified forest too", {
  fit <- rpf(mpg ~ ., data = mtcars, max_interaction = 2, ntrees = 5)
  purify(fit)
  expect_error(purify(fit, mode = 99), "mode")
  expect_error(purify(fit, maxp_interaction = 0), "maxp_interaction")
  expect_error(purify(fit, nthreads = 0), "nthreads")
  expect_error(purify(fit, 2), class = "rlib_error_dots_nonempty")
})
