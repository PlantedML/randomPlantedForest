# cur_trees_1: deterministic mode evaluates every sampled candidate, and leaves are drawn in proportion to their
# available thresholds (the incremental sampling weights after a split agree with a full rebuild).

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
