# rpf(prune = TRUE): the forest pruned by prune_components() with its defaults, on the training rows.
# It must equal prune_components(rpf(..., prune = FALSE), data) bit for bit, from the same random number state.

rp_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n))
  x$y <- x$x1 + 2 * x$x2 * x$x3 * x$x4 + rnorm(n, sd = 0.1)
  x
}

# The pruned model without what cannot be compared with identical(): the C++ forest (compared by its model, its
# purified rectangles and its predictions), the caches of the built store and the display's built indices (external
# pointers) and the timings (the display's, and a recipe's own prep times)
rp_strip <- function(p) {
  p$fit$fit <- NULL
  if (!is.null(p$fit$blueprint$recipe)) {
    p$fit$blueprint$recipe$fit_times <- NULL
  }
  p$cache <- NULL
  p$display$cache <- NULL
  p$display$seconds <- NULL
  for (nm in names(p$display$coarse)) {
    p$display$coarse[[nm]]$seconds <- NULL
  }
  p
}

expect_same_pruned <- function(a, b, newdata) {
  expect_s3_class(a, "rpf_pruned")
  expect_s3_class(b, "rpf_pruned")
  expect_identical(rp_strip(a), rp_strip(b))
  expect_identical(a$fit$fit$get_model(), b$fit$fit$get_model())
  expect_identical(a$fit$fit$get_purified_rects(), b$fit$fit$get_purified_rects())
  # the pruning leaves the forest unpurified
  expect_false(is_purified(a$fit))
  expect_identical(predict(a, newdata), predict(b, newdata))
  expect_identical(predict_components(a, newdata), predict_components(b, newdata))
}

# fit both ways from the same random number state (max_panels as given to rpf()); the caller's stream ends in the same
# state too
rp_both <- function(fit_call, data, y = NULL, max_panels = 3) {
  set.seed(7)
  a <- fit_call(prune = TRUE)
  seed_a <- .Random.seed
  set.seed(7)
  b <- prune_components(fit_call(prune = FALSE), data, y, max_panels = max_panels)
  expect_identical(.Random.seed, seed_a)
  list(a = a, b = b)
}

rp_args <- list(max_interaction = 3, ntrees = 12, split_try = 5, t_try = 0.5, max_panels = 3)

test_that("rpf(prune = TRUE) equals prune_components(rpf(prune = FALSE), data): formula and x/y", {
  d <- rp_data()
  x <- d[, 1:4]
  # one formula: its environment is part of the blueprint
  fml <- y ~ .
  for (splits in list(60, NULL)) {
    args <- c(rp_args, list(splits = splits))
    f <- rp_both(function(prune) do.call(rpf, c(list(fml, data = d, prune = prune), args)), d)
    expect_same_pruned(f$a, f$b, d)
    # the display ran: a component of order 3 is coarsened to at most 3 heat maps
    expect_true(f$a$display$second)
    xy <- rp_both(function(prune) do.call(rpf, c(list(x = x, y = d$y, prune = prune), args)), x, d$y)
    expect_same_pruned(xy$a, xy$b, x)
    xm <- as.matrix(x)
    m <- rp_both(function(prune) do.call(rpf, c(list(x = xm, y = d$y, prune = prune), args)), xm, d$y)
    expect_same_pruned(m$a, m$b, xm)
    # the interfaces fit the same forest, so they prune the same
    expect_identical(predict(f$a, d), predict(xy$a, x))
    expect_identical(predict(f$a, d), predict(m$a, xm))
    if (is.null(splits)) {
      expect_false(is.null(f$a$fit$splits_oob))
      expect_identical(f$a$fit$params$splits, f$a$fit$splits_oob$splits)
    }
  }
})

test_that("rpf(prune = TRUE) equals prune_components(rpf(prune = FALSE), data): recipe", {
  skip_if_not_installed("recipes")
  d <- rp_data(n = 250, seed = 2)
  rec <- recipes::step_normalize(recipes::recipe(y ~ ., data = d), recipes::all_numeric_predictors())
  for (splits in list(60, NULL)) {
    args <- c(rp_args, list(splits = splits))
    r <- rp_both(function(prune) do.call(rpf, c(list(x = rec, data = d, prune = prune), args)), d)
    expect_same_pruned(r$a, r$b, d)
  }
})

test_that("max_panels = Inf: no display; purify is not used with prune = TRUE", {
  d <- rp_data(n = 250, seed = 3)
  fit_call <- function(prune, ...) {
    rpf(y ~ ., data = d, max_interaction = 3, ntrees = 12, splits = 60, max_panels = Inf, prune = prune, ...)
  }
  f <- rp_both(fit_call, d, max_panels = Inf)
  expect_same_pruned(f$a, f$b, d)
  expect_null(f$a$display)
  set.seed(7)
  p <- fit_call(prune = TRUE, purify = TRUE)
  expect_same_pruned(p, f$a, d)
  expect_false(p$fit$params$purify)
  # prune = FALSE: the forest, as before
  set.seed(7)
  expect_s3_class(fit_call(prune = FALSE), "rpf")
})

test_that("prune = TRUE prunes classification forests too, is an error for deterministic ones; arguments are checked", {
  d <- rp_data(n = 100)
  set.seed(1)
  pi <- suppressMessages(rpf(Species ~ ., data = iris, ntrees = 2, prune = TRUE))
  expect_s3_class(pi, "rpf_pruned")
  expect_identical(colnames(predict(pi, iris)), paste0(".pred_", levels(iris$Species)))
  expect_error(
    rpf(y ~ ., data = d, ntrees = 1, deterministic = TRUE, prune = TRUE),
    "not available with `deterministic = TRUE`"
  )
  expect_error(rpf(y ~ ., data = d, ntrees = 2, prune = NA), "prune")
  expect_error(rpf_050(y ~ ., data = d, ntrees = 2, max_panels = 0), "max_panels")
  expect_error(rpf_050(y ~ ., data = d, ntrees = 2, max_panels = 2.5), "max_panels")
  # the errors name rpf()
  err <- tryCatch(rpf(y ~ ., data = d, ntrees = 1, deterministic = TRUE, prune = TRUE), error = identity)
  expect_match(conditionMessage(err), "deterministic")
  expect_identical(deparse(conditionCall(err)[[1]]), "rpf")
})
