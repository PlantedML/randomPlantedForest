# prune_components(): stored components only, the purification, the random number stream and the training rows

fix_sim <- function(n = 300, p = 4, seed = 11) {
  set.seed(seed)
  x <- matrix(runif(n * p, -1, 1), ncol = p, dimnames = list(NULL, paste0("x", seq_len(p))))
  y <- 1 + 2 * x[, 1] + 4 * x[, 2] * x[, 3] + rnorm(n, sd = 0.3)
  list(x = x, y = y)
}

fix_fit <- function(d, max_interaction = 2, ntrees = 20, splits = 40, ...) {
  set.seed(5)
  rpf(d$x, d$y, max_interaction = max_interaction, ntrees = ntrees, splits = splits, ...)
}

# The decomposition over every predictor set up to the forest's order, for comparison
all_sets_decompose <- function(fit, x) {
  p <- ncol(x)
  sets <- unlist(
    lapply(seq_len(min(p, fit$params$max_interaction)), function(k) utils::combn(p, k, simplify = FALSE)),
    recursive = FALSE
  )
  M <- fit$fit$predict_components_matrix(x, sets, 1L)
  v <- apply(M, 2, function(c) mean((c - mean(c))^2))
  list(sets = sets[v > 0], M = M[, v > 0, drop = FALSE])
}

test_that("only the stored components are evaluated, with the results of all subsets", {
  ns <- asNamespace("randomPlantedForest")
  for (k in 2:3) {
    d <- fix_sim(p = 5)
    fit <- fix_fit(d, max_interaction = k)
    pr <- prune_components(fit, d$x, d$y, rule = "D", tau = 0)
    ref <- all_sets_decompose(fit, d$x)
    expect_identical(pr$components$component, ns$prune_set_names(ref$sets, colnames(d$x)))
    dec <- ns$prune_decompose(fit, ns$prune_training_data(fit, d$x, d$y)$predictors, rep(TRUE, nrow(d$x)))
    expect_identical(unname(dec$own), unname(ref$M))
    # only the stored sets were evaluated, ordered as utils::combn() orders all sets
    expect_length(dec$sets, length(fit$fit$purified_components()))
    expect_lte(length(dec$sets), sum(choose(5, seq_len(k))))
    expect_identical(dec$sets[dec$nonzero], ref$sets)
  }
})

test_that("many predictors at full order: only the stored components, not 2^p - 1 sets", {
  skip_on_cran()
  # 2^20 - 1 candidate sets, a few hundred stored components
  d <- fix_sim(n = 200, p = 20)
  set.seed(5)
  fit <- rpf(d$x, d$y, max_interaction = 20, ntrees = 3, splits = 15)
  t0 <- proc.time()[["elapsed"]]
  pr <- prune_components(fit, d$x, d$y, rule = "R", tau = c(0, 0.01), folds = 2, seed = 1)
  expect_lt(proc.time()[["elapsed"]] - t0, 60)
  expect_lt(nrow(pr$components), 1000)
  expect_equal(
    predict(prune_components(fit, d$x, d$y, rule = "D", tau = 0), d$x)$.pred,
    predict(fit, d$x)$.pred,
    tolerance = 1e-10
  )
})

test_that("a forest purified below its order is purified again at all orders", {
  d <- fix_sim()
  ref <- prune_components(fix_fit(d, max_interaction = 3), d$x, d$y, folds = 3, seed = 2)
  fit <- fix_fit(d, max_interaction = 3)
  p0 <- predict(fit, d$x)
  purify(fit, maxp_interaction = 1)
  expect_true(fit$fit$is_purify_capped())
  expect_message(pr <- prune_components(fit, d$x, d$y, folds = 3, seed = 2), "below its order")
  expect_false(fit$fit$is_purify_capped())
  expect_identical(pr$cv, ref$cv)
  expect_identical(pr$beta, ref$beta)
  expect_identical(pr$components, ref$components)
  expect_equal(predict(fit, d$x)$.pred, p0$.pred, tolerance = 1e-12)
  # a given rule too
  f2 <- fix_fit(d, max_interaction = 3)
  purify(f2, maxp_interaction = 2, mode = 1)
  expect_message(p2 <- prune_components(f2, d$x, d$y, rule = "R", tau = 0.01), "below its order")
  expect_identical(p2$components$component, ref$components$component)

  # a purification at the forest's order (or above) is fine, and so is a mode-1 purification (equal up to rounding)
  f3 <- fix_fit(d, max_interaction = 3)
  purify(f3, maxp_interaction = 3)
  expect_silent(p3 <- prune_components(f3, d$x, d$y, folds = 3, seed = 2))
  expect_identical(p3$cv, ref$cv)
  f1 <- fix_fit(d, max_interaction = 3)
  purify(f1, mode = 1)
  p1 <- prune_components(f1, d$x, d$y, folds = 3, seed = 2)
  expect_equal(p1$cv, ref$cv, tolerance = 1e-10)
  expect_equal(p1$beta, ref$beta, tolerance = 1e-10)
})

test_that("the caller's random number stream is restored", {
  d <- fix_sim()
  fit <- fix_fit(d)
  set.seed(9)
  before <- .Random.seed
  prune_components(fit, d$x, d$y, folds = 3, seed = 1)
  expect_identical(.Random.seed, before)
  # without a seed, one is drawn: the stream advances by that draw only
  set.seed(9)
  s <- sample.int(.Machine$integer.max, 1L)
  after_draw <- .Random.seed
  set.seed(9)
  pr <- prune_components(fit, d$x, d$y, folds = 3)
  expect_identical(.Random.seed, after_draw)
  expect_identical(pr$seed, s)
  # a given rule draws nothing
  set.seed(9)
  prune_components(fit, d$x, d$y, rule = "R", tau = 0.01)
  expect_identical(.Random.seed, before)
})

test_that("`data` and `y` must be the training rows, also for a given rule", {
  d <- fix_sim()
  fit <- fix_fit(d)
  perm <- sample(nrow(d$x))
  expect_error(prune_components(fit, d$x[perm, ], d$y[perm], folds = 3, seed = 1), "training rows")
  expect_error(prune_components(fit, d$x[perm, ], d$y[perm], rule = "R", tau = 0.01), "training rows")
  expect_error(prune_components(fit, d$x[-1, ], d$y[-1], rule = "R", tau = 0.01), "training rows")
  expect_error(prune_components(fit, d$x, rev(d$y), rule = "R", tau = 0.01), "training rows")
  # a restored forest without its training data cannot be checked
  restored <- rpf_unmarshal(rpf_marshal(fit))
  expect_error(prune_components(restored, d$x, d$y, rule = "R", tau = 0.01), "no training data")
  pr <- prune_components(rpf_unmarshal(rpf_marshal(fit, include_data = TRUE)), d$x, d$y, rule = "R", tau = 0.01)
  expect_identical(pr$beta, prune_components(fit, d$x, d$y, rule = "R", tau = 0.01)$beta)
})
