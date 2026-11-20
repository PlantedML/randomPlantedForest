# The exponential loss with K >= 3 classes is SAMME (Zhu, Zou, Rosset and Hastie 2009): one loss, exp(-y.f / K),
# gives the split criterion, the leaf step and one weight per row. The forest stores g = f / (K - 1), whose softmax
# predict() returns, so after growing every row's weight must be exp(-g_c), c the row's class, recomputed from the
# scores the forest predicts at the row, and the family's scores must be those predictions. Binary data: the weight
# exp(-y g / 2), y = -1 / 1.

samme_data <- function(K, n = 150) {
  set.seed(3)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n), x3 = sample(0:3, n, TRUE))
  eta <- cbind(0, 3 * d$x1 - 1.5, 2 * (d$x2 > 0.5) + d$x3 / 2 - 1, 1 - 4 * d$x1 * d$x2)
  cls <- apply(eta, 1, function(e) sample.int(4, 1, prob = exp(e)))
  d$y <- factor(letters[pmin(cls, K)], levels = letters[seq_len(K)])
  d
}

# one family grown with its weights kept (internal hook): the data rows of its sample, their classes, their weights,
# their scores in the family and the scores the forest predicts there
samme_family <- function(d, structure, deterministic) {
  set.seed(1)
  m <- rpf(
    y ~ .,
    data = d,
    loss = "exponential",
    max_interaction = 3,
    n_families = 1,
    n_splits = 40,
    split_structure = structure,
    n_try_split_points = 3,
    share_candidates = 0.3,
    deterministic = deterministic,
    prune = FALSE
  )
  m$fit$set_parameters("keep_weights", 1) # refits, keeping the weights
  kw <- m$fit$kept_weights()
  expect_length(kw, 1)
  rows <- kw[[1]]$rows
  list(
    class = as.integer(d$y[rows]),
    w = kw[[1]]$weights,
    scores = kw[[1]]$scores,
    g = as.matrix(predict(m, d[rows, ], type = "link")),
    p = as.matrix(predict(m, d[rows, ], type = "prob"))
  )
}

test_that("multiclass exponential (SAMME): each row's weight is exp(-g_c) of its predicted scores", {
  structures <- c("leaf", "tree_then_leaf", "tree_then_all_leaves", "target_tree_then_leaf", "leaf_binned")
  for (K in 3:4) {
    d <- samme_data(K)
    for (s in structures) {
      for (det in c(FALSE, TRUE)) {
        lab <- sprintf("K = %d, %s, deterministic = %s", K, s, det)
        f <- samme_family(d, s, det)
        # one weight per row
        expect_identical(dim(f$w), c(nrow(d), 1L), label = lab)
        expect_identical(dim(f$scores), c(nrow(d), K), label = lab)
        expect_equal(f$scores, f$g, tolerance = 1e-12, ignore_attr = TRUE, label = lab)
        g_c <- f$g[cbind(seq_along(f$class), f$class)]
        expect_lt(max(abs(f$w[, 1] / exp(-g_c) - 1)), 1e-12, label = lab)
        # the fit moved the scores, which sum to 0 over the classes (every leaf step does), and predict() gives their
        # softmax
        expect_gt(max(abs(g_c)), 0.5, label = lab)
        expect_lt(max(abs(rowSums(f$g))), 1e-12, label = lab)
        expect_equal(f$p, exp(f$g) / rowSums(exp(f$g)), tolerance = 1e-12, ignore_attr = TRUE, label = lab)
      }
    }
  }
})

test_that("binary exponential: the weight is exp(-y g / 2) of the predicted score", {
  d <- samme_data(2)
  for (det in c(FALSE, TRUE)) {
    f <- samme_family(d, "leaf", det)
    expect_identical(dim(f$w), c(nrow(d), 1L))
    expect_equal(f$scores, f$g, tolerance = 1e-12, ignore_attr = TRUE)
    y <- ifelse(f$class == 2, 1, -1)
    expect_lt(max(abs(f$w[, 1] / exp(-y * f$g[, 1] / 2) - 1)), 1e-12)
    expect_gt(max(abs(f$g)), 0.5)
  }
})

test_that("the exponential's first leaf step is log(n_k + a / K) - mean_j log(n_j + a / K), a = 1 pseudo-row", {
  # the leaves' steps, before the level (test-classif-level.R)
  d <- samme_data(3)
  set.seed(1)
  m <- rpf(
    y ~ .,
    data = d,
    loss = "exponential",
    max_interaction = 1,
    n_families = 1,
    n_splits = 1,
    n_try_split_points = 3,
    share_candidates = 1,
    deterministic = TRUE,
    prune = FALSE,
    level = FALSE
  )
  g <- as.matrix(predict(m, d, type = "link"))
  # the one split leaves two leaves, all weights 1 and scores 0 before it: each leaf's step from its class counts n_k
  # and one pseudo-row split over the K classes at its mean score 0, no clamp
  side <- factor(apply(round(g, 12), 1, paste, collapse = " "))
  expect_length(levels(side), 2)
  for (s in levels(side)) {
    in_leaf <- side == s
    v <- tabulate(as.integer(d$y[in_leaf]), 3) + 1 / 3
    expect_equal(unname(g[which(in_leaf)[1], ]), log(v) - mean(log(v)), tolerance = 1e-14)
  }
  # binary: log((n_1 + 1 / 2) / (n_0 + 1 / 2))
  d2 <- samme_data(2)
  set.seed(1)
  m <- rpf(
    y ~ .,
    data = d2,
    loss = "exponential",
    max_interaction = 1,
    n_families = 1,
    n_splits = 1,
    n_try_split_points = 3,
    share_candidates = 1,
    deterministic = TRUE,
    prune = FALSE,
    level = FALSE
  )
  g <- predict(m, d2, type = "link")[[1]]
  for (s in unique(round(g, 12))) {
    in_leaf <- round(g, 12) == s
    n <- tabulate(as.integer(d2$y[in_leaf]), 2)
    expect_equal(g[which(in_leaf)[1]], log((n[2] + 0.5) / (n[1] + 0.5)), tolerance = 1e-14)
  }
})
