# prune_components(): stored components only, the purification, the random number stream, the training rows, the
# out-of-bag choice with few families, and a corrected OOB error of 0

fix_sim <- function(n = 300, p = 4, seed = 11) {
  set.seed(seed)
  x <- matrix(runif(n * p, -1, 1), ncol = p, dimnames = list(NULL, paste0("x", seq_len(p))))
  y <- 1 + 2 * x[, 1] + 4 * x[, 2] * x[, 3] + rnorm(n, sd = 0.3)
  list(x = x, y = y)
}

fix_fit <- function(d, max_interaction = 2, n_families = 20, n_splits = 40, ...) {
  set.seed(5)
  rpf_050(d$x, d$y, max_interaction = max_interaction, n_families = n_families, n_splits = n_splits, ...)
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

test_that("only the components the forest has are evaluated, with the results of all subsets", {
  ns <- asNamespace("randomPlantedForest")
  for (k in 2:3) {
    d <- fix_sim(p = 5)
    fit <- fix_fit(d, max_interaction = k)
    pr <- prune_components(fit, d$x, d$y, rule = "D", tau = 0)
    dec <- ns$prune_columns(fit, colnames(d$x), 1L, oob = FALSE)$dec
    expect_false(is_purified(fit))
    # against every predictor set of the purified forest
    purify(fit)
    ref <- all_sets_decompose(fit, d$x)
    expect_identical(pr$components$component, ns$prune_set_names(ref$sets, colnames(d$x)))
    expect_identical(dec$sets[dec$nonzero], ref$sets)
    expect_lt(max(abs(unname(dec$own) - unname(ref$M))), 1e-13)
    # only the sets the purified forest stores were evaluated, ordered as utils::combn() orders all sets
    expect_identical(dec$sets, ns$prune_stored_sets(fit))
    expect_lte(length(dec$sets), sum(choose(5, seq_len(k))))
  }
})

test_that("many predictors at full order: only the stored components, not 2^p - 1 sets", {
  skip_on_cran()
  # 2^20 - 1 candidate sets, a few hundred stored components
  d <- fix_sim(n = 200, p = 20)
  set.seed(5)
  fit <- rpf_050(d$x, d$y, max_interaction = 20, n_families = 3, n_splits = 15)
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

test_that("the forest is not purified, and an earlier purification of it does not matter", {
  d <- fix_sim()
  model <- c("rule", "tau", "intercept", "beta", "components", "store")
  for (by in c("oob", "cv")) {
    f0 <- fix_fit(d, max_interaction = 3)
    ref <- prune_components(f0, d$x, d$y, select_by = by, folds = 3, seed = 2)
    expect_false(is.null(ref[[by]]))
    expect_false(is_purified(f0))
    # purified below the forest's order, at its order, or in mode 1: the same model, the purification left as it is
    for (how in list(list(maxp_interaction = 1), list(maxp_interaction = 3), list(mode = 1))) {
      fit <- fix_fit(d, max_interaction = 3)
      do.call(purify, c(list(fit), how))
      capped <- fit$fit$is_purify_capped()
      p0 <- predict(fit, d$x)
      expect_silent(pr <- prune_components(fit, d$x, d$y, select_by = by, folds = 3, seed = 2))
      expect_identical(pr[c(by, model)], ref[c(by, model)])
      expect_true(is_purified(fit))
      expect_identical(fit$fit$is_purify_capped(), capped)
      expect_identical(predict(fit, d$x), p0)
    }
  }
  # a given rule too
  f2 <- fix_fit(d, max_interaction = 3)
  purify(f2, maxp_interaction = 2, mode = 1)
  p2 <- prune_components(f2, d$x, d$y, rule = "R", tau = 0.01)
  expect_identical(
    p2[model],
    prune_components(fix_fit(d, max_interaction = 3), d$x, d$y, rule = "R", tau = 0.01)[model]
  )
  expect_true(f2$fit$is_purify_capped())
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
  fit <- fix_fit(d, n_families = 8)
  expect_message(pr <- prune_components(fit, d$x, d$y, folds = 3, seed = 1), "scored by CV")
  expect_identical(pr$select_by, "cv")
  expect_false(is.null(pr$cv))
  expect_warning(
    po <- prune_components(fit, d$x, d$y, select_by = "oob", rule = "D", folds = 3, seed = 1),
    "fewer than 10"
  )
  expect_identical(po$select_by, "oob")
  # 10 families: out of bag, no message
  expect_silent(p10 <- prune_components(fix_fit(d, n_families = 10), d$x, d$y, folds = 3, seed = 1))
  expect_identical(p10$select_by, "oob")

  # D's corrected OOB error by hand: rows with OOB families in one half only get the mean v of the other rows
  o <- fix_fit(d, n_families = 8)$fit$purify_oob(0L, 1L)
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

# The choice in an OOB table, by its definition (prune_choose()): a candidate whose corrected OOB error was clamped
# to 0 is left out unless all are; among the others, min has the lowest RMSE, and 1se the largest tau whose MSE is
# within one SE of min's. Which candidates are clamped depends on the forest, and so on the platform's random draws.
expect_oob_choice <- function(oob) {
  zero <- oob$oob_rmse == 0
  use <- if (all(zero)) rep(TRUE, nrow(oob)) else !zero
  expect_equal(sum(oob$min), 1)
  expect_equal(sum(oob$one_se), 1)
  expect_false(any((oob$min | oob$one_se) & !use))
  best <- which(oob$min)
  one <- which(oob$one_se)
  expect_lte(oob$oob_rmse[best], min(oob$oob_rmse[use]) * (1 + 1e-12))
  within <- use & oob$oob_rmse^2 <= (oob$oob_rmse[best]^2 + oob$se_mse[best]) * (1 + 1e-12)
  expect_true(within[one])
  expect_identical(oob$tau[one], max(oob$tau[within]))
}

test_that("tiny data: a corrected OOB error of 0 does not break the choice", {
  # 32 rows, 10 families: R at a small tau has about as many coefficients as a fold has rows, its predictions differ
  # wildly between the two halves of the families, and the correction can exceed the uncorrected OOB MSE
  # (0.5.0's defaults and the exact draw: the fit the clamped candidates below are named for)
  old <- options(randomPlantedForest.exact_draw = TRUE)
  on.exit(options(old))
  set.seed(1)
  pm <- rpf_050(
    mpg ~ .,
    data = mtcars,
    n_families = 10,
    n_splits = NULL,
    prune = TRUE,
    nthreads = 1,
    keep_forest = TRUE
  )
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
  # the clamped candidates (on macOS arm64 R at 0 and 0.01%; none in the Linux forest) are left out of the choice
  expect_true(all(oob$oob_rmse_raw > 0))
  expect_oob_choice(oob)
})

test_that("noiseless data: corrected OOB errors of 0 are left out, and the choice keeps the true component", {
  # y = x1 exactly: R's uncorrected OOB error is almost all the OOB means' extra variance, and the correction,
  # estimated from the two halves of the families, can exceed it
  n <- 24
  i <- seq_len(n)
  x <- data.frame(x1 = (i - 0.5) / n, x2 = ((i * 7) %% n) / n)
  y <- x$x1
  set.seed(1)
  fit <- rpf_050(x, y, max_interaction = 1, n_families = 10, n_splits = 10, nthreads = 1)
  pr <- prune_components(fit, x, y, seed = 1, nthreads = 1)
  oob <- pr$oob
  zero <- oob$oob_rmse == 0
  expect_equal(sum(oob$min), 1)
  expect_equal(sum(oob$one_se), 1)
  expect_true(all(zero) || !any(zero[oob$min | oob$one_se]))
  # the clamped candidates (on macOS arm64 every R candidate; none in the Linux forest) are left out of the choice
  expect_true(all(oob$oob_rmse_raw > 0))
  expect_oob_choice(oob)
  expect_identical(c(pr$rule, pr$tau), unlist(oob[oob$one_se, c("rule", "tau")], use.names = FALSE))
  expect_identical(names(pr$beta), "x1")
})
