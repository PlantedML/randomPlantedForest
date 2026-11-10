# Classification split structures behave as their regression versions (same candidate pool,
# sampling and aging; only the loss differs).

structures <- c(res_trees = 0, cur_trees_2 = 1, cur_trees_1 = 2, leaves = 3, hist = 4)

set.seed(101)
n_cs <- 200
X_cs <- cbind(matrix(runif(n_cs * 3), n_cs), rbinom(n_cs, 1, 0.5), sample(1:4, n_cs, TRUE))
eta_cs <- 2 * (X_cs[, 1] - 0.5) + X_cs[, 4] - 0.5 + 3 * (X_cs[, 2] - 0.5) * (X_cs[, 3] - 0.5)
Y_cs <- matrix(rbinom(n_cs, 1, plogis(eta_cs)), ncol = 1)
dat_cs <- data.frame(y = factor(Y_cs[, 1]), X_cs)

# Low-level parameter vectors: max_interaction, n_families, n_splits, n_try_split_points, share_candidates, purify,
# deterministic, nthreads, cv, split_decay_rate, max_candidates, delete_leaves, split_structure (+ delta, epsilon)
reg_params <- function(mode, deterministic = 0, n_families = 3, n_splits = 25) {
  c(2, n_families, n_splits, 5, 0.4, 0, deterministic, 1, 0, 0.1, 50, 1, mode)
}

# Leaves of a tree D that have both bounds inside the data range in some d in D: such a leaf comes from
# refining a leaf of D within D itself.
own_dim_interior <- function(fit, X) {
  lo <- apply(X, 2, min)
  hi <- apply(X, 2, max)
  cnt <- 0
  for (fam in fit$fit$get_model()) {
    for (t in seq_along(fam$variables)) {
      dims <- fam$variables[[t]]
      dims <- dims[dims > 0]
      for (iv in fam$intervals[[t]]) {
        for (d in dims) {
          cnt <- cnt + (iv[1, d] > lo[d] && iv[2, d] < hi[d])
        }
      }
    }
  }
  cnt
}

test_that("classification with loss L2 reproduces regression on a 0/1 response", {
  for (nm in names(structures)) {
    for (det in c(0, 1)) {
      p <- reg_params(structures[[nm]], det)
      set.seed(5)
      reg <- new(RandomPlantedForest, Y_cs, X_cs, p)
      set.seed(5)
      cls <- new(ClassificationRPF, Y_cs, X_cs, "L2", c(p, 0, 0.1))
      d <- max(abs(cls$predict_matrix(X_cs, 0, 1L) - reg$predict_matrix(X_cs, 0, 1L)))
      expect_lt(d, 1e-10, label = paste(nm, "deterministic =", det))
    }
  }
})

test_that("classification cur_trees_1 differs from cur_trees_2", {
  for (loss in c("logit", "exponential", "L2")) {
    set.seed(3)
    f1 <- rpf_050(
      y ~ .,
      data = dat_cs,
      loss = loss,
      n_families = 5,
      max_interaction = 2,
      split_structure = "tree_then_leaf"
    )
    set.seed(3)
    f2 <- rpf_050(
      y ~ .,
      data = dat_cs,
      loss = loss,
      n_families = 5,
      max_interaction = 2,
      split_structure = "tree_then_all_leaves"
    )
    p1 <- predict(f1, dat_cs, type = "prob")
    p2 <- predict(f2, dat_cs, type = "prob")
    expect_gt(max(abs(as.matrix(p1) - as.matrix(p2))), 1e-6, label = loss)
  }
})

test_that("classification cur_trees_2 and res_trees refine leaves within their own tree", {
  for (s in c("tree_then_all_leaves", "target_tree_then_leaf", "tree_then_leaf", "leaf", "leaf_binned")) {
    set.seed(4)
    f <- rpf_050(
      y ~ .,
      data = dat_cs,
      loss = "logit",
      n_families = 5,
      n_splits = 40,
      max_interaction = 2,
      split_structure = s
    )
    expect_gt(own_dim_interior(f, X_cs), 0, label = s)
  }
})

test_that("split_decay_rate changes a res_trees classification fit", {
  fits <- lapply(c(0, 2), function(rate) {
    set.seed(6)
    f <- rpf_050(
      y ~ .,
      data = dat_cs,
      loss = "logit",
      n_families = 5,
      split_structure = "target_tree_then_leaf",
      split_decay_rate = rate
    )
    as.matrix(predict(f, dat_cs, type = "prob"))
  })
  expect_gt(max(abs(fits[[1]] - fits[[2]])), 1e-6)
})

test_that("all classification losses run under every split structure, binary and multiclass", {
  losses <- c(
    "logit",
    "logit_2",
    "logit_3",
    "logit_4",
    "exponential",
    "exponential_2",
    "exponential_3",
    "L1",
    "L2",
    "median"
  )
  cls3 <- cut(eta_cs + rnorm(n_cs), 3, labels = FALSE)
  Y_multi <- outer(cls3, 1:3, "==") * 1
  for (Y in list(binary = Y_cs, multiclass = Y_multi)) {
    for (loss in losses) {
      Yl <- if (startsWith(loss, "exponential")) 2 * Y - 1 else Y
      for (nm in names(structures)) {
        set.seed(8)
        cls <- new(ClassificationRPF, Yl, X_cs, loss, c(reg_params(structures[[nm]], n_splits = 15), 0, 0.1))
        pr <- cls$predict_matrix(X_cs, 0, 1L)
        expect_equal(dim(pr), c(n_cs, ncol(Y)), label = paste(loss, nm))
        expect_true(all(is.finite(pr)), label = paste(loss, nm, ncol(Y)))
      }
    }
  }
})

test_that("deterministic res_trees and cur_trees_2 learn in classification", {
  # Friedman 1 with 5 features as the log-odds; one deterministic family
  set.seed(8)
  gen <- function(n) data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n), x5 = runif(n))
  eta <- function(X) (10 * sin(pi * X$x1 * X$x2) + 20 * (X$x3 - 0.5)^2 + 10 * X$x4 + 5 * X$x5 - 14) / 3
  X <- gen(600)
  Xt <- gen(2000)
  set.seed(9)
  y <- factor(ifelse(runif(nrow(X)) < plogis(eta(X)), "b", "a"))
  pt <- plogis(eta(Xt))
  for (loss in c("logit", "L2")) {
    for (ss in c("target_tree_then_leaf", "tree_then_all_leaves")) {
      fit <- rpf_050(
        x = X,
        y = y,
        max_interaction = 2,
        n_families = 1,
        n_splits = 100,
        loss = loss,
        deterministic = TRUE,
        split_structure = ss
      )
      p <- predict(fit, Xt, type = "prob")$.pred_b
      # before the fix no better than the constant (about 0.25); fixed: as cur_trees_1
      expect_lt(mean(abs(p - pt)), 0.22, label = paste(loss, ss))
    }
  }
})
