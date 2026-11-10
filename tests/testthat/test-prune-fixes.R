# prune_components(): stored components only, the purification, the random number stream, the training rows, the
# out-of-bag choice with few families, and a corrected OOB error of 0

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

test_that("candidates whose corrected OOB MSE was set to 0 are not chosen, unless all are", {
  ns <- asNamespace("randomPlantedForest")
  cand <- data.frame(rule = rep(c("D", "R"), each = 3), tau = rep(c(0, 0.01, 0.1), 2))
  # R at 0 and 0.01 left out (at 0, below every other RMSE and within any SE): min D 0; 1se within D 0's SE, D 0.01
  rmse <- c(1, 1.1, 2, 0, 0, 3)
  se <- c(0.3, 0.1, 0.1, 0.5, 5, 0.1)
  ch <- ns$prune_choose(rmse, se, cand, left_out = rmse == 0)
  expect_identical(which(ch$min), 1L)
  expect_identical(which(ch$one_se), 2L)
  # all left out: they tie at 0, the largest tau, then D
  ch <- ns$prune_choose(rep(0, 6), rep(0, 6), cand, left_out = rep(TRUE, 6))
  expect_identical(which(ch$min), 3L)
  expect_identical(which(ch$one_se), 3L)
  # none left out: as without left_out
  expect_identical(ns$prune_choose(rmse + 1, se, cand, left_out = rep(FALSE, 6)), ns$prune_choose(rmse + 1, se, cand))
  # without left_out (the CV), an RMSE of 0 is a tie at 0: at the same tau, D
  ch <- ns$prune_choose(c(0, 1, 1, 0, 2, 2), rep(0.1, 6), cand)
  expect_identical(which(ch$min), 1L)
  expect_identical(which(ch$one_se), 1L)
  # above 0, ties are relative as before
  ch <- ns$prune_choose(c(2, 1, 1 + 1e-13, 3, 3, 3), rep(0, 6), cand)
  expect_identical(which(ch$min), 3L)
  expect_identical(which(ch$one_se), 3L)
})

test_that("the adaptive stop, the pruning and the tuner choose the lowest corrected error by one rule", {
  ns <- asNamespace("randomPlantedForest")
  # a corrected MSE below 0 (the path's, unclamped) or clamped at 0 (the pruning's and the tuner's tables): left out
  mse <- c(0.5, -0.01, 0.3, 0.31)
  expect_identical(ns$oob_select(mse), 3L)
  expect_identical(ns$tune_oob_choose(pmax(mse, 0), simpler = 1:4), 3L)
  rmse <- sqrt(pmax(mse, 0))
  cand <- data.frame(rule = rep("D", 4), tau = c(0.03, 0.02, 0.01, 0))
  expect_identical(which(ns$prune_choose(rmse, rep(0, 4), cand, left_out = rmse == 0)$min), 3L)
  # all at or below 0: all tie
  expect_identical(ns$oob_lowest(c(0, -1, 0)), 1:3)
})

test_that("tiny data: a corrected OOB error of 0 does not break the choice", {
  # 32 rows, 10 families: R at a small tau has about as many coefficients as a fold has rows, its predictions differ
  # wildly between the two halves of the families, and the correction exceeds the uncorrected OOB MSE
  set.seed(1)
  pm <- rpf(mpg ~ ., data = mtcars, ntrees = 10, splits = NULL, prune = TRUE, nthreads = 1)
  expect_s3_class(pm, "rpf_pruned")
  oob <- pm$oob
  zero <- oob$oob_rmse == 0
  expect_equal(sum(oob$min), 1)
  expect_equal(sum(oob$one_se), 1)
  # no choice at 0 while a candidate is above 0
  expect_true(all(zero) || !any(zero[oob$min | oob$one_se]))
  expect_identical(c(pm$rule, pm$tau), unlist(oob[oob$one_se, c("rule", "tau")], use.names = FALSE))
  expect_true(all(is.finite(predict(pm, mtcars)$.pred)))
  pmin <- prune_components(pm$fit, mtcars, select = "min", seed = pm$seed, max_panels = Inf, nthreads = 1)
  expect_identical(pmin$oob, oob)
  expect_identical(c(pmin$rule, pmin$tau), unlist(oob[oob$min, c("rule", "tau")], use.names = FALSE))
  # this platform's fit: R at 0 and 0.01% at 0, left out
  skip_on_cran()
  expect_identical(zero, oob$rule == "R" & oob$tau <= 1e-4)
  expect_true(all(oob$oob_rmse_raw > 0))
  expect_identical(c(pmin$rule, pmin$tau), c("D", "0"))
  expect_identical(c(pm$rule, pm$tau), c("D", "0.01"))
})

test_that("noiseless data: R's corrected OOB errors are 0, and the choice keeps the true component", {
  # y = x1 exactly: R's uncorrected OOB error is almost all the OOB means' extra variance, and the correction,
  # estimated from the two halves of the families, exceeds it
  n <- 24
  i <- seq_len(n)
  x <- data.frame(x1 = (i - 0.5) / n, x2 = ((i * 7) %% n) / n)
  y <- x$x1
  set.seed(1)
  fit <- rpf(x, y, max_interaction = 1, ntrees = 10, splits = 10, nthreads = 1)
  pr <- prune_components(fit, x, y, seed = 1, nthreads = 1)
  oob <- pr$oob
  zero <- oob$oob_rmse == 0
  expect_equal(sum(oob$min), 1)
  expect_equal(sum(oob$one_se), 1)
  expect_true(all(zero) || !any(zero[oob$min | oob$one_se]))
  # this platform's fit: every R candidate at 0, left out; D at 10% with x1 alone
  skip_on_cran()
  expect_identical(zero, oob$rule == "R")
  expect_true(all(oob$oob_rmse_raw > 0))
  expect_identical(oob$min, oob$rule == "D" & oob$tau == 0.1)
  expect_identical(oob$one_se, oob$min)
  expect_identical(c(pr$rule, pr$tau), c("D", "0.1"))
  expect_identical(names(pr$beta), "x1")
})
