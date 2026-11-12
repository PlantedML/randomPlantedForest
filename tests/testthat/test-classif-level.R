# The classification level (classif_level()): one constant per score column, the maximum-likelihood intercept of the
# log-loss given the forest's out-of-bag scores (the balance property), in every family's null leaf. The option
# randomPlantedForest.level (TRUE by default) switches it; off, a fit is the fit of earlier versions.

level_data <- function(n = 400, seed = 3, K = 2) {
  set.seed(seed)
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  if (K == 2) {
    p <- plogis(-2.5 + 2 * x$x1 + sin(2 * pi * x$x2))
    x$y <- factor(ifelse(runif(n) < p, "yes", "no"), levels = c("no", "yes"))
  } else {
    eta <- cbind(0, 2 * x$x1 - 1, -2.5 + 2 * x$x2)
    cls <- apply(eta, 1, function(e) sample.int(3, 1, prob = exp(e)))
    x$y <- factor(c("a", "b", "c")[cls])
  }
  x
}

level_raw <- function(fit, d) {
  X <- preprocess_predictors_predict(fit, hardhat::forge(d, fit$blueprint)$predictors)
  matrix(fit$fit$predict_matrix(X, 0, 0L), nrow = nrow(d))
}

with_level <- function(on, code) {
  old <- options(randomPlantedForest.level = on)
  on.exit(options(old))
  code
}

test_that("classif_level() solves the balance equations (binary, multiclass logit, SAMME)", {
  set.seed(11)
  n <- 500
  for (case in list(c("logit", 2), c("logit", 4), c("exponential", 2), c("exponential", 3))) {
    K <- as.integer(case[2])
    y <- factor(sample(letters[1:K], n, TRUE, prob = c(0.6, rep(0.4 / (K - 1), K - 1))), levels = letters[1:K])
    crit <- oob_criterion("classification", case[1], y)
    G <- matrix(rnorm(n * crit$V, sd = 0.7), n)
    # no correction (B = m): the mean probability of every class is its frequency
    lev <- classif_level(crit, seq_len(n), G, G, G, rep(1L, n), rep(1L, n), 2)
    expect_length(lev, crit$V)
    expect_lte(attr(lev, "iterations"), 10L)
    P <- crit_logloss(crit, G + rep(lev, each = n), crit$cls)$P
    Pf <- if (crit$link == "sigmoid") {
      cbind(1 - P[, 1], P[, 1])
    } else if (crit$link == "softmax_ref") {
      cbind(1 - rowSums(P), P)
    } else {
      P
    }
    expect_equal(colMeans(Pf), as.numeric(table(y)) / n, tolerance = 1e-10)
    if (crit$link == "softmax") {
      expect_equal(sum(lev), 0, tolerance = 1e-12)
    }
    # the level is the maximum-likelihood intercept: the log-loss is lowest there
    ll <- function(c) sum(crit_logloss(crit, G + rep(c, each = n), crit$cls)$e)
    for (v in seq_len(crit$V)) {
      e <- replace(numeric(crit$V), v, 1e-3)
      expect_lt(ll(lev), ll(lev + e))
      expect_lt(ll(lev), ll(lev - e))
    }
  }
  # a class without an out-of-bag row: no finite level, 0
  y <- factor(c(rep("a", 10), "b"), levels = c("a", "b", "c"))
  crit <- oob_criterion("classification", "logit", y)
  G <- matrix(0, 10, 2)
  expect_identical(as.numeric(classif_level(crit, 1:10, G, G, G, rep(1L, 10), rep(1L, 10), 2)), c(0, 0))
})

test_that("the level moves the forest's score by a constant; off, no level and the same splits", {
  d <- level_data()
  fits <- lapply(c(FALSE, TRUE), function(on) {
    with_level(on, {
      set.seed(4)
      rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 20, n_splits = 30, prune = FALSE)
    })
  })
  expect_null(fits[[1]]$level)
  expect_false(fits[[1]]$params$level)
  expect_true(fits[[2]]$params$level)
  lev <- fits[[2]]$level
  expect_length(lev, 1)
  # a rare class: the forest's mean probability is below the frequency, the level raises it
  expect_gt(lev, 0)
  diff <- level_raw(fits[[2]], d) - level_raw(fits[[1]], d)
  expect_equal(diff, matrix(lev, nrow(d), 1), tolerance = 1e-12)
  # the same leaves but the null leaves
  m0 <- fits[[1]]$fit$get_model()
  m1 <- fits[[2]]$fit$get_model()
  for (b in seq_along(m0)) {
    null <- which(vapply(m0[[b]]$variables, function(v) identical(as.numeric(v), 0), TRUE))
    expect_identical(m0[[b]]$values[-null], m1[[b]]$values[-null])
    expect_equal(unlist(m1[[b]]$values[null]), unlist(m0[[b]]$values[null]) + lev, tolerance = 1e-15)
  }
  # the out-of-bag balance at the forest's last split: the OOB mean probability is the frequency up to the
  # second-order term of the OOB correction, which is large with 20 families (1 / m - 1 / B about 0.08)
  path <- fits[[2]]$fit$oob_path()
  m0 <- path$count[, 1]
  m1 <- path$count[, 2]
  ok <- m0 + m1 > 0
  last <- length(path$checkpoints)
  s0 <- path$sum0[ok, last, drop = FALSE]
  s1 <- path$sum1[ok, last, drop = FALSE]
  g <- (s0 + s1) / (m0 + m1)[ok]
  expect_equal(mean(plogis(g + lev)), mean(d$y[ok] == "yes"), tolerance = 0.05)
  crit <- oob_criterion("classification", "logit", d$y)
  expect_equal(
    as.numeric(classif_level(crit, which(ok), g, s0 / pmax(m0[ok], 1), s1 / pmax(m1[ok], 1), m0[ok], m1[ok], 20)),
    lev
  )
  # without the second-order term (B = m: as many families as the OOB means), the plain balance
  lev0 <- classif_level(crit, which(ok), g, g, g, rep(1L, sum(ok)), rep(1L, sum(ok)), 2)
  expect_equal(mean(plogis(g + as.numeric(lev0))), mean(d$y[ok] == "yes"), tolerance = 1e-10)
  expect_match(paste(format(fits[[2]]), collapse = "\n"), "level: +", fixed = TRUE)
  expect_no_match(paste(format(fits[[1]]), collapse = "\n"), "level:", fixed = TRUE)
})

test_that("multiclass logit and exponential get one level per score column", {
  d <- level_data(K = 3)
  for (loss in c("logit", "exponential")) {
    set.seed(5)
    fit <- rpf(y ~ ., data = d, loss = loss, max_interaction = 2, n_families = 20, n_splits = 25, prune = FALSE)
    expect_length(fit$level, if (loss == "logit") 2L else 3L)
    if (loss == "exponential") {
      expect_equal(sum(fit$level), 0, tolerance = 1e-12)
    }
    with_level(FALSE, {
      set.seed(5)
      fit0 <- rpf(y ~ ., data = d, loss = loss, max_interaction = 2, n_families = 20, n_splits = 25, prune = FALSE)
    })
    expect_equal(
      level_raw(fit, d) - level_raw(fit0, d),
      matrix(fit$level, nrow(d), length(fit$level), byrow = TRUE),
      tolerance = 1e-12
    )
  }
})

test_that("L1 / L2 classification and regression get no level", {
  d <- level_data()
  set.seed(6)
  fit <- rpf(y ~ ., data = d, loss = "L2", n_families = 10, n_splits = 10, prune = FALSE)
  expect_null(fit$level)
  expect_null(fit$params$level)
  d$z <- as.numeric(d$y == "yes") + d$x1
  set.seed(6)
  fit <- rpf(z ~ x1 + x2 + x3, data = d, n_families = 10, n_splits = 10, prune = FALSE)
  expect_null(fit$level)
})

test_that("the adaptive fit scores every checkpoint at its level and returns the forest at the chosen one's", {
  d <- level_data(n = 300)
  set.seed(7)
  fit <- rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 20, prune = FALSE)
  # the chosen checkpoint's level, recomputed from the fit's path at that checkpoint (a fit with the chosen splits)
  set.seed(7)
  fixed <- rpf(
    y ~ .,
    data = d,
    loss = "logit",
    max_interaction = 2,
    n_families = 20,
    prune = FALSE,
    n_splits = fit$params$splits
  )
  expect_equal(fit$level, fixed$level, tolerance = 1e-12)
  expect_equal(level_raw(fit, d), level_raw(fixed, d), tolerance = 1e-12)
  # the path of a fit with fixed splits scores each checkpoint at its own level; the last one is the adaptive path's
  set.seed(7)
  path_fit <- rpf(
    y ~ .,
    data = d,
    loss = "logit",
    max_interaction = 2,
    n_families = 20,
    prune = FALSE,
    n_splits = fit$params$splits,
    oob_path = TRUE
  )
  k <- match(fit$params$splits, fit$oob_path$n_splits)
  expect_equal(tail(path_fit$oob_path$oob_loss, 1), fit$oob_path$oob_loss[k], tolerance = 1e-12)
})

test_that("a path of fewer checkpoints than the splits reports those only; the level is the last split's", {
  d <- level_data(n = 300)
  set.seed(8)
  fit <- rpf(
    y ~ .,
    data = d,
    loss = "logit",
    max_interaction = 2,
    n_families = 20,
    n_splits = 30,
    prune = FALSE,
    oob_path = c(5, 10)
  )
  expect_identical(fit$oob_path$n_splits, c(5L, 10L))
  set.seed(8)
  fit30 <- rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 20, n_splits = 30, prune = FALSE)
  expect_identical(fit$level, fit30$level)
})

test_that("purify = TRUE purifies the forest with its level; the pruned model's intercept holds it", {
  d <- level_data()
  set.seed(9)
  fit <- rpf(
    y ~ .,
    data = d,
    loss = "logit",
    max_interaction = 2,
    n_families = 20,
    n_splits = 30,
    prune = FALSE,
    purify = TRUE
  )
  expect_true(fit$fit$is_purified())
  set.seed(9)
  raw <- rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 20, n_splits = 30, prune = FALSE)
  expect_equal(predict(fit, d, type = "link"), predict(raw, d, type = "link"), tolerance = 1e-12)
  # rule D at tau = 0 keeps every component and the intercept: the forest
  pr <- prune_components(raw, d, rule = "D", tau = 0, select = "min")
  expect_equal(predict(pr, d, type = "prob"), predict(raw, d, type = "prob"), tolerance = 1e-10)
})

test_that("deterministic: the level balances the forest's own predictions", {
  d <- level_data()
  fit <- rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 1, n_splits = 15, deterministic = TRUE)
  expect_length(fit$level, 1)
  expect_equal(mean(predict(fit, d, type = "prob")$.pred_yes), mean(d$y == "yes"), tolerance = 1e-10)
})

test_that("the pruning's cross-validation refits each fold with its level; a fit without params$level gets none", {
  d <- level_data(n = 300)
  set.seed(10)
  fit <- rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 10, n_splits = 20, prune = FALSE)
  refit <- prune_refit(fit, d[c("x1", "x2", "x3")], d$y, seq_len(200), seed = 1, nthreads = 1L)
  # its level: the one a fit on those rows gets
  set.seed(1)
  direct <- rpf(
    y ~ .,
    data = d[seq_len(200), ],
    loss = "logit",
    max_interaction = 2,
    n_families = 10,
    n_splits = 20,
    prune = FALSE
  )
  expect_equal(level_raw(refit, d), level_raw(direct, d), tolerance = 1e-12)
  # a fit saved before the level existed has no params$level: its refits get none
  old <- fit
  old$params$level <- NULL
  expect_false(level_on(old$params, "classification"))
  refit_old <- prune_refit(old, d[c("x1", "x2", "x3")], d$y, seq_len(200), seed = 1, nthreads = 1L)
  with_level(FALSE, {
    set.seed(1)
    direct0 <- rpf(
      y ~ .,
      data = d[seq_len(200), ],
      loss = "logit",
      max_interaction = 2,
      n_families = 10,
      n_splits = 20,
      prune = FALSE
    )
  })
  expect_identical(level_raw(refit_old, d), level_raw(direct0, d))
})

test_that("the option must be TRUE or FALSE", {
  d <- level_data(n = 100)
  with_level("yes", {
    expect_error(
      rpf(y ~ ., data = d, loss = "logit", n_families = 2, n_splits = 2, prune = FALSE),
      "randomPlantedForest.level"
    )
  })
})

test_that("families without a split: the level is computed from their fits, and is small", {
  d <- level_data(n = 100)
  # 100 rows: no threshold keeps 51 rows on each side, so no family splits
  fit_path <- function(level) {
    with_level(level, {
      set.seed(2)
      rpf(
        y ~ .,
        data = d,
        loss = "logit",
        n_families = 10,
        n_splits = 5,
        min_leaf_rows = 51,
        prune = FALSE,
        oob_path = TRUE
      )
    })
  }
  on <- fit_path(TRUE)
  expect_true(all(vapply(on$fit$get_model(), function(f) length(f$variables) == 1, TRUE)))
  path <- on$fit$oob_path()
  m <- path$count[, 1] + path$count[, 2]
  ok <- m > 0
  g <- (path$sum0[ok, 1] + path$sum1[ok, 1]) / m[ok]
  # the OOB path holds the families' fits (test-no-split-family.R), without the level
  nulls <- vapply(on$fit$get_model(), function(f) f$values[[1]][[1]][1], 0) - on$level
  expect_true(all(g >= min(nulls) - 1e-12 & g <= max(nulls) + 1e-12))
  expect_lt(max(abs(g - stats::qlogis(mean(d$y == "yes")))), 0.5)
  # the OOB fits of the families are balanced already, up to the bootstrap: a small level
  expect_lt(abs(on$level), 0.2)
  # off: the same path and the same families, without the level
  off <- fit_path(FALSE)
  expect_identical(off$fit$oob_path(), path)
  nulls_off <- vapply(off$fit$get_model(), function(f) f$values[[1]][[1]][1], 0)
  expect_equal(nulls_off, nulls, tolerance = 1e-12)
})
