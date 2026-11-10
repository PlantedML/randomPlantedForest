# Regression split structures: deterministic mode evaluates every sampled candidate (cur_trees_1) and spreads its
# thresholds over each leaf (res_trees, cur_trees_2); cur_trees_1 draws leaves in proportion to their available
# thresholds (the incremental sampling weights after a split agree with a full rebuild).

# Low-level parameter vector: max_interaction, ntrees, splits, split_try, t_try, purify, deterministic,
# nthreads, cv, split_decay_rate, max_candidates, delete_leaves, split_structure
low_params <- function(
  mode,
  max_interaction = 2,
  ntrees = 3,
  splits = 25,
  split_try = 5,
  t_try = 0.4,
  deterministic = 0,
  delete_leaves = 1
) {
  c(max_interaction, ntrees, splits, split_try, t_try, 0, deterministic, 1, 0, 0.1, 50, delete_leaves, mode)
}

tree_dims <- function(fit) {
  unique(unlist(lapply(fit$get_model(), function(fam) lapply(fam$variables, as.numeric))))
}

test_that("deterministic cur_trees_1 evaluates every sampled candidate", {
  # y does not depend on x1, the feature of the first candidate
  set.seed(1)
  n <- 100
  X <- matrix(runif(n * 3), n)
  Y <- matrix(sin(2 * pi * X[, 2]) + X[, 3], ncol = 1)
  fit <- new(RandomPlantedForest, Y, X, low_params(2, ntrees = 1, splits = 20, deterministic = 1))
  expect_true(all(c(2, 3) %in% tree_dims(fit)))
  mse <- mean((fit$predict_matrix(X, 0, 1L) - Y)^2)
  expect_lt(mse, 0.5 * var(Y[, 1]))
})

test_that("cur_trees_1 draws only leaves with an available threshold (ties)", {
  # One feature with 10 distinct values, one draw per round (split_try = 1). The root leaf can always be split,
  # so if leaves are weighted by their available thresholds every round finds a split; with delete_leaves = 0
  # each split adds two leaves to the tree {1}.
  set.seed(2)
  X <- matrix(rep(1:10, each = 30), ncol = 1)
  Y <- matrix(X[, 1] + rnorm(nrow(X)), ncol = 1)
  splits <- 30
  set.seed(3)
  fit <- new(
    RandomPlantedForest,
    Y,
    X,
    low_params(2, max_interaction = 1, ntrees = 5, splits = splits, split_try = 1, t_try = 1, delete_leaves = 0)
  )
  for (fam in fit$get_model()) {
    i <- which(vapply(fam$variables, function(v) identical(as.numeric(v), 1), TRUE))
    expect_equal(length(fam$intervals[[i]]), 2 * splits)
  }
})

# Friedman 1 with 5 features; one deterministic family
det_data <- function() {
  set.seed(8)
  gen <- function(n) data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n), x5 = runif(n))
  f <- function(X) 10 * sin(pi * X$x1 * X$x2) + 20 * (X$x3 - 0.5)^2 + 10 * X$x4 + 5 * X$x5
  X <- gen(600)
  Xt <- gen(2000)
  list(X = X, y = f(X) + rnorm(600), Xt = Xt, ft = f(Xt))
}

test_that("deterministic cur_trees_1, res_trees and cur_trees_2 fit as well as the other structures", {
  d <- det_data()
  for (ss in c("cur_trees_1", "res_trees", "cur_trees_2")) {
    fit <- rpf(
      x = d$X,
      y = d$y,
      max_interaction = 2,
      ntrees = 1,
      splits = 100,
      deterministic = TRUE,
      split_structure = ss
    )
    # before the fix: 75-95% of Var f unexplained (cur_trees_1 split only x1; res_trees and cur_trees_2 tried only
    # the lowest thresholds of each leaf)
    mse <- mean((predict(fit, d$Xt)$.pred - d$ft)^2) / var(d$ft)
    expect_lt(mse, 0.15, label = paste(ss, "test MSE / Var f"))
    # the main-effect split points are not piled up at the bottom of each feature
    m <- fit$fit$get_model()[[1]]
    low <- unlist(lapply(seq_along(m$variables), function(t) {
      v <- as.numeric(m$variables[[t]])
      if (length(v) != 1 || v == 0) {
        return(NULL)
      }
      b <- vapply(m$intervals[[t]], function(iv) iv[1, v], numeric(1))
      b <- b[b > min(d$X[[v]])]
      b < stats::quantile(d$X[[v]], 0.05)
    }))
    expect_gt(length(low), 10)
    expect_lt(mean(low), 0.25, label = paste(ss, "share of split points below the 5% quantile"))
  }
})
