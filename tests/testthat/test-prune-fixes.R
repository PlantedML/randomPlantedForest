# prune_components(): stored components only, the purification, the random number stream, the training rows, and
# the out-of-bag choice with few families

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
  pr <- prune_components(fit, d$x, d$y, rule = "R", tau = c(0, 0.01), select_by = "cv", folds = 2, seed = 1)
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
  for (by in c("oob", "cv")) {
    ref <- prune_components(fix_fit(d, max_interaction = 3), d$x, d$y, select_by = by, folds = 3, seed = 2)
    expect_false(is.null(ref[[by]]))
    fit <- fix_fit(d, max_interaction = 3)
    p0 <- predict(fit, d$x)
    purify(fit, maxp_interaction = 1)
    expect_true(fit$fit$is_purify_capped())
    expect_message(pr <- prune_components(fit, d$x, d$y, select_by = by, folds = 3, seed = 2), "below its order")
    expect_false(fit$fit$is_purify_capped())
    expect_identical(pr[[by]], ref[[by]])
    expect_identical(pr$beta, ref$beta)
    expect_identical(pr$components, ref$components)
    expect_equal(predict(fit, d$x)$.pred, p0$.pred, tolerance = 1e-12)

    # a purification at the forest's order (or above) is fine, and so is a mode-1 purification (equal up to
    # rounding; the out-of-bag choice purifies in mode 2 again)
    f3 <- fix_fit(d, max_interaction = 3)
    purify(f3, maxp_interaction = 3)
    expect_silent(p3 <- prune_components(f3, d$x, d$y, select_by = by, folds = 3, seed = 2))
    expect_identical(p3[[by]], ref[[by]])
    f1 <- fix_fit(d, max_interaction = 3)
    purify(f1, mode = 1)
    p1 <- prune_components(f1, d$x, d$y, select_by = by, folds = 3, seed = 2)
    expect_equal(p1[[by]], ref[[by]], tolerance = 1e-10)
    expect_equal(p1$beta, ref$beta, tolerance = 1e-10)
  }
  # a given rule too
  f2 <- fix_fit(d, max_interaction = 3)
  purify(f2, maxp_interaction = 2, mode = 1)
  expect_message(p2 <- prune_components(f2, d$x, d$y, rule = "R", tau = 0.01), "below its order")
  expect_identical(p2$components$component, ref$components$component)
})

test_that("the caller's random number stream is restored", {
  d <- fix_sim()
  fit <- fix_fit(d)
  for (by in c("oob", "cv")) {
    set.seed(9)
    before <- .Random.seed
    prune_components(fit, d$x, d$y, select_by = by, folds = 3, seed = 1)
    expect_identical(.Random.seed, before)
    # without a seed, one is drawn: the stream advances by that draw only
    set.seed(9)
    s <- sample.int(.Machine$integer.max, 1L)
    after_draw <- .Random.seed
    set.seed(9)
    pr <- prune_components(fit, d$x, d$y, select_by = by, folds = 3)
    expect_identical(pr$select_by, by)
    expect_identical(.Random.seed, after_draw)
    expect_identical(pr$seed, s)
  }
  # a given rule draws nothing
  set.seed(9)
  before <- .Random.seed
  prune_components(fit, d$x, d$y, rule = "R", tau = 0.01)
  expect_identical(.Random.seed, before)
})

test_that("`data` and `y` must be the training rows, also for a given rule", {
  d <- fix_sim()
  fit <- fix_fit(d)
  perm <- sample(nrow(d$x))
  for (by in c("oob", "cv")) {
    expect_error(
      prune_components(fit, d$x[perm, ], d$y[perm], select_by = by, folds = 3, seed = 1),
      "training rows"
    )
  }
  expect_error(prune_components(fit, d$x[perm, ], d$y[perm], rule = "R", tau = 0.01), "training rows")
  expect_error(prune_components(fit, d$x[-1, ], d$y[-1], rule = "R", tau = 0.01), "training rows")
  expect_error(prune_components(fit, d$x, rev(d$y), rule = "R", tau = 0.01), "training rows")
  # a restored forest without its training data cannot be checked
  restored <- rpf_unmarshal(rpf_marshal(fit))
  expect_error(prune_components(restored, d$x, d$y, rule = "R", tau = 0.01), "no training data")
  pr <- prune_components(rpf_unmarshal(rpf_marshal(fit, include_data = TRUE)), d$x, d$y, rule = "R", tau = 0.01)
  expect_identical(pr$beta, prune_components(fit, d$x, d$y, rule = "R", tau = 0.01)$beta)
})

test_that("with few families the default scores by CV, and an empty half gets the mean variance", {
  d <- fix_sim()
  fit <- fix_fit(d, ntrees = 8)
  expect_message(pr <- prune_components(fit, d$x, d$y, folds = 3, seed = 1), "scored by CV")
  expect_identical(pr$select_by, "cv")
  expect_false(is.null(pr$cv))
  expect_warning(
    po <- prune_components(fit, d$x, d$y, select_by = "oob", rule = "D", folds = 3, seed = 1),
    "fewer than 10"
  )
  expect_identical(po$select_by, "oob")
  # 10 families: out of bag, no message
  expect_silent(p10 <- prune_components(fix_fit(d, ntrees = 10), d$x, d$y, folds = 3, seed = 1))
  expect_identical(p10$select_by, "oob")

  # D's corrected OOB error by hand: rows with OOB families in one half only get the mean v of the other rows
  o <- fix_fit(d, ntrees = 8)$fit$purify_oob(0L, 1L)
  m0 <- o$count[, 1]
  m1 <- o$count[, 2]
  m <- m0 + m1
  use <- m > 0
  both <- m0 > 0 & m1 > 0
  expect_true(any(use & !both))
  nm <- vapply(o$sets, function(s) paste(colnames(d$x)[s], collapse = ":"), "")
  colnames(o$sum0) <- colnames(o$sum1) <- nm
  share <- stats::setNames(po$components$share, po$components$component)
  for (tau in unique(po$oob$tau)) {
    kept <- names(share)[share >= tau]
    g <- (o$intercept0[, 1] +
      o$intercept1[, 1] +
      rowSums(o$sum0[, kept, drop = FALSE]) +
      rowSums(o$sum1[, kept, drop = FALSE])) /
      m
    g0 <- (o$intercept0[, 1] + rowSums(o$sum0[, kept, drop = FALSE])) / pmax(m0, 1)
    g1 <- (o$intercept1[, 1] + rowSums(o$sum1[, kept, drop = FALSE])) / pmax(m1, 1)
    v <- (g0 - g1)^2 / (1 / pmax(m0, 1) + 1 / pmax(m1, 1))
    v[!both] <- mean(v[both & use])
    e2 <- ((d$y - g)^2 - (1 / m - 1 / 8) * v)[use]
    row <- po$oob$tau == tau
    expect_equal(po$oob$oob_rmse[row], sqrt(mean(e2)), tolerance = 1e-12)
    expect_equal(attr(po$oob, "n_rows"), sum(use))
  }
})
