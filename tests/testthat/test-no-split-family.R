no_split_data <- function(n = 100) {
  set.seed(42)
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  x$y <- 10 + x$x1 + rnorm(n, sd = 0.1)
  x$cl <- factor(ifelse(runif(n) < 0.7, "a", "b"))
  x$cl3 <- factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.5, 0.3, 0.2)))
  x
}

test_that("a family without a split predicts the mean of its rows", {
  dat <- no_split_data()
  constant <- dat
  constant[c("x1", "x2")] <- 1
  for (structure in c("leaf", "leaf_binned", "tree_then_leaf", "tree_then_all_leaves", "target_tree_then_leaf")) {
    fit <- rpf_050(
      y ~ x1 + x2,
      data = constant,
      n_families = 1,
      n_splits = 5,
      split_structure = structure,
      deterministic = TRUE
    )
    expect_equal(predict(fit, dat[1:3, ])$.pred, rep(mean(dat$y), 3), label = structure)
  }

  # a single row: every bootstrap sample is that row
  set.seed(1)
  fit <- rpf_050(y ~ x1 + x2, data = dat[1, ], n_families = 5, n_splits = 5)
  expect_equal(predict(fit, dat[1:3, ])$.pred, rep(dat$y[1], 3))
})

test_that("a classification family without a split predicts its rows' fit", {
  dat <- no_split_data()
  dat[c("x1", "x2")] <- 1
  shares <- function(y) as.numeric(prop.table(table(y)))
  # exponential, and logit with exact_step: the minimiser of the penalised loss from the scores 0, the class shares
  # with one pseudo-row split evenly over the classes
  smoothed <- function(y) {
    n <- as.numeric(table(y))
    (n + 1 / length(n)) / (sum(n) + 1)
  }
  # logit with exact_step = FALSE: one Newton step from the scores 0, the pseudo-rows' Hessian at its bound
  # (1 / 2)(I - 11' / K1); in the coding with the first class at 0 (the step does not depend on the coding)
  newton <- function(y) {
    n <- as.numeric(table(y))
    K1 <- length(n)
    K <- K1 - 1
    H <- sum(n) * (diag(K) / K1 - matrix(1 / K1^2, K, K)) + 0.5 * (diag(K) - matrix(1 / K1, K, K))
    step <- solve(H, n[-1] - sum(n) / K1)
    exp(c(0, step)) / sum(exp(c(0, step)))
  }
  fit <- function(y, ...) {
    rpf_050(
      stats::as.formula(paste(y, "~ x1 + x2")),
      data = dat,
      n_families = 1,
      n_splits = 5,
      deterministic = TRUE,
      ...
    )
  }
  p1 <- function(f) unlist(predict(f, dat[1, ]), use.names = FALSE)
  for (y in c("cl", "cl3")) {
    for (loss in c("L1", "L2")) {
      expect_equal(p1(fit(y, loss = loss)), shares(dat[[y]]), label = paste(loss, y))
    }
    expect_equal(p1(fit(y, loss = "exponential")), smoothed(dat[[y]]), tolerance = 1e-12, label = y)
    expect_equal(p1(fit(y, loss = "logit")), smoothed(dat[[y]]), tolerance = 1e-12, label = y)
    expect_equal(p1(fit(y, loss = "logit", exact_step = FALSE)), newton(dat[[y]]), tolerance = 1e-12, label = y)
  }

  # no clamp: a rare class keeps its smoothed share (3 of 100 rows: 3.5 / 101)
  dat$rare <- factor(ifelse(seq_len(nrow(dat)) <= 97, "a", "b"))
  for (loss in c("exponential", "logit")) {
    expect_equal(p1(fit("rare", loss = loss)), c(97.5, 3.5) / 101, tolerance = 1e-12, label = loss)
  }
})

test_that("in a forest where only some families split, the others predict their rows' mean", {
  # a bootstrap sample of these three rows has one distinct x in 9 of 27 cases: then no split is possible
  x <- data.frame(x = c(0, 0, 1))
  y <- c(1, 10, 100)
  set.seed(3)
  fit <- rpf_050(x, y, n_families = 40, n_splits = 3, export_forest = TRUE)
  null_value <- vapply(
    fit$forest,
    function(family) {
      null_tree <- which(vapply(family$variables, function(v) identical(as.integer(v), 0L), logical(1)))
      family$values[[null_tree]][[1]]
    },
    numeric(1)
  )
  split <- vapply(fit$forest, function(family) length(family$variables) > 1, logical(1))
  expect_true(any(split) && any(!split))
  # families with a split keep their null leaf at 0; the others hold the mean of three rows of one x value
  expect_true(all(null_value[split] == 0))
  expect_true(all(null_value[!split] %in% c(1, 4, 7, 10, 100)))

  # purification puts the null leaves' values into the intercept, in both modes
  pred <- predict(fit, x)$.pred
  for (mode in 1:2) {
    set.seed(3)
    fit <- rpf_050(x, y, n_families = 40, n_splits = 3)
    purify(fit, mode = mode)
    components <- predict_components(fit, x)
    expect_equal(predict(fit, x)$.pred, pred)
    expect_equal(rowSums(components$m) + components$intercept, pred)
  }
})

test_that("min_leaf_rows above half the rows leaves every family without a split: they predict their rows' mean", {
  dat <- no_split_data()
  # 100 rows: no threshold keeps 51 rows on each side, with either draw of candidates
  for (exact in c(TRUE, FALSE)) {
    old <- options(randomPlantedForest.exact_draw = exact)
    set.seed(1)
    fit <- rpf_050(y ~ x1 + x2, data = dat, n_families = 5, n_splits = 10, min_leaf_rows = 51, export_forest = TRUE)
    options(old)
    expect_true(all(vapply(fit$forest, function(family) length(family$variables) == 1, logical(1))))
    pred <- predict(fit, dat[1:3, ])$.pred
    expect_equal(pred, rep(pred[1], 3))
    # the mean of five bootstrap means
    expect_lt(abs(pred[1] - mean(dat$y)), 0.1)
  }
  # without the bootstrap, the mean itself
  fit <- rpf_050(y ~ x1 + x2, data = dat, n_families = 1, n_splits = 10, min_leaf_rows = 51, deterministic = TRUE)
  expect_equal(predict(fit, dat[1:3, ])$.pred, rep(mean(dat$y), 3))
})

test_that("rpf()'s default call (splits chosen out of bag): a family without a split predicts its rows' fit", {
  # The adaptive fit rebuilds the forest from the families' records at the chosen number of splits; a family without
  # a split there gets its rows' fit, as a fit with that number of splits gives it. Checked on the forest itself
  # (`$fit`, kept with keep_forest = TRUE) and on the pruned model.
  dat <- no_split_data()
  constant <- dat
  constant[c("x1", "x2")] <- 1
  shares <- function(y) as.numeric(prop.table(table(y)))
  check_reg <- function(fit, label) {
    for (model in list(fit, fit$fit)) {
      pred <- predict(model, dat[1:3, ])$.pred
      expect_equal(pred, rep(pred[1], 3), label = label)
      # the mean of the families' bootstrap means
      expect_lt(abs(pred[1] - mean(dat$y)), 0.1, label = label)
    }
  }
  check_cls <- function(fit, y, label) {
    for (model in list(fit, fit$fit)) {
      p <- unlist(predict(model, dat[1, ]), use.names = FALSE)
      expect_lt(max(abs(p - shares(y))), 0.1, label = label)
    }
  }
  set.seed(1)
  check_reg(rpf(y ~ x1 + x2, data = constant, keep_forest = TRUE), "constant predictors")
  for (outcome in c("cl", "cl3")) {
    for (loss in c("exponential", "logit")) {
      set.seed(1)
      fit <- rpf(stats::as.formula(paste(outcome, "~ x1 + x2")), data = constant, loss = loss, keep_forest = TRUE)
      check_cls(fit, constant[[outcome]], paste(outcome, loss))
    }
  }

  # min_leaf_rows above half the rows: no family can split
  set.seed(1)
  check_reg(rpf(y ~ x1 + x2, data = dat, min_leaf_rows = 51, keep_forest = TRUE), "min_leaf_rows")
  for (outcome in c("cl", "cl3")) {
    set.seed(1)
    fit <- rpf(stats::as.formula(paste(outcome, "~ x1 + x2")), data = dat, min_leaf_rows = 51, keep_forest = TRUE)
    check_cls(fit, dat[[outcome]], paste(outcome, "min_leaf_rows"))
  }
})

test_that("the pruned model's intercept includes the families that made no split", {
  # one row differs in the only predictor: the families whose bootstrap sample misses it cannot split
  set.seed(1)
  n <- 60
  d <- data.frame(x1 = c(1, rep(0, n - 1)), y = rnorm(n, 5))
  fit <- rpf(y ~ x1, data = d, n_families = 20, n_splits = 3, prune = FALSE, max_interaction = 1)
  families <- fit$fit$get_model()
  no_split <- vapply(families, function(f) identical(as.integer(unlist(f$variables)), 0L), logical(1))
  expect_true(any(no_split) && !all(no_split))
  # rule D at tau = 0 is the forest
  forest <- prune_components(fit, d, rule = "D", tau = 0, max_panels = Inf)
  expect_equal(predict(forest, d)$.pred, predict(fit, d)$.pred, tolerance = 1e-12)
})
