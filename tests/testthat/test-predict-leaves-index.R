# predict_matrix() of an unpurified forest indexes the leaves of each tree of more than 64 leaves when it predicts
# 512 rows or more, and scans every leaf otherwise (src/lib/predict.cpp). Both must give bitwise the same predictions.

# the scan: chunks of at most 7 rows
scan_rows <- function(fit, X, components = 0) {
  chunks <- split(seq_len(nrow(X)), ceiling(seq_len(nrow(X)) / 7))
  rows <- lapply(chunks, function(i) fit$fit$predict_matrix(X[i, , drop = FALSE], components, 1L))
  do.call(rbind, rows)
}

# rows to predict (at least 512): the training rows (values on the split points), random rows, the bounds, values
# beyond them, Inf and NaN
leaf_test_rows <- function(x, n_random = 0) {
  X <- as.matrix(x)
  rng <- apply(X, 2, range)
  random <- apply(rng, 2, function(r) runif(n_random, r[1] - 0.1, r[2] + 0.1))
  X <- rbind(X, random, rng[1, ], rng[2, ], rng[1, ] - 1, rng[2, ] + 1, rep(Inf, ncol(X)), rep(-Inf, ncol(X)))
  X[nrow(X), 2] <- NaN
  X
}

# the number of leaves of the forest's largest tree: the index is used for trees of more than 64
max_tree_leaves <- function(fit) {
  max(unlist(lapply(fit$fit$get_model(), function(family) lengths(family$values))))
}

test_that("the leaf index gives the scan's predictions bitwise (regression)", {
  set.seed(3)
  n <- 600
  x <- data.frame(x1 = round(runif(n), 2), x2 = runif(n), x3 = sample(0:3, n, TRUE), x4 = rnorm(n))
  y <- x$x1 + x$x2 * x$x3 + sin(3 * x$x4) + rnorm(n, sd = 0.2)
  X <- leaf_test_rows(x)
  for (st in c("leaves", "cur_trees_1", "cur_trees_2", "res_trees", "hist")) {
    for (dl in c(TRUE, FALSE)) {
      set.seed(5)
      fit <- rpf_050(
        x,
        y,
        max_interaction = 3,
        ntrees = 4,
        splits = 400,
        split_structure = st,
        delete_leaves = dl
      )
      label <- paste(st, dl)
      expect_gt(max_tree_leaves(fit), 64, label = label)
      scan <- scan_rows(fit, X)
      for (nt in 1:2) {
        expect_identical(fit$fit$predict_matrix(X, 0, nt), scan, label = paste(label, nt))
      }
      for (S in list(1, c(1, 3), c(2, 3), c(1, 2, 4), -1)) {
        Xs <- X[, abs(S), drop = FALSE]
        expect_identical(fit$fit$predict_matrix(Xs, S, 2L), scan_rows(fit, Xs, S), label = paste(label, toString(S)))
      }
    }
  }
})

test_that("the leaf index gives the scan's predictions bitwise (classification)", {
  set.seed(4)
  n <- 300
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  f <- x$x1 + x$x2 * x$x3
  yb <- factor(ifelse(f + rnorm(n, sd = 0.2) > 0.75, "a", "b"))
  y3 <- factor(cut(f + rnorm(n, sd = 0.2), 3, labels = c("a", "b", "c")))
  X <- leaf_test_rows(x, n_random = 300)
  for (loss in c("L2", "exponential")) {
    for (y in list(yb, y3)) {
      for (dl in c(TRUE, FALSE)) {
        set.seed(6)
        fit <- rpf_050(x, y, max_interaction = 3, ntrees = 4, splits = 150, loss = loss, delete_leaves = dl)
        label <- paste(loss, nlevels(y), dl)
        expect_gt(max_tree_leaves(fit), 64, label = label)
        expect_identical(fit$fit$predict_matrix(X, 0, 2L), scan_rows(fit, X), label = label)
        expect_identical(
          fit$fit$predict_matrix(X[, 1:2], c(1, 2), 2L),
          scan_rows(fit, X[, 1:2], c(1, 2)),
          label = paste(label, "x1:x2")
        )
      }
    }
  }
})

test_that("predict() gives the same predictions for all rows at once as in chunks", {
  set.seed(7)
  n <- 600
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  x$y <- x$x1 + x$x2 * x$x3 + rnorm(n, sd = 0.1)
  fit <- rpf_050(y ~ ., data = x, max_interaction = 3, ntrees = 4, splits = 300, delete_leaves = FALSE)
  expect_gt(max_tree_leaves(fit), 64)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / 100))
  expect_identical(
    predict(fit, x, nthreads = 2L)$.pred,
    unlist(lapply(chunks, function(i) predict(fit, x[i, ])$.pred), use.names = FALSE)
  )
})
