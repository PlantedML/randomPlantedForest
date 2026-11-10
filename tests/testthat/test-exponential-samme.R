# The exponential loss with K >= 3 classes is SAMME (Zhu, Zou, Rosset and Hastie 2009): one loss, exp(-y.f / K),
# gives the split criterion, the leaf step and one weight per row. The forest stores g = f / (K - 1), whose softmax
# predict() returns, so after growing every row's weight must be exp(-g_c), c the row's class, recomputed from the
# scores the forest predicts at the row. Binary data keep master's code: the weight exp(-y g / 2), y = -1 / 1.

samme_data <- function(K, n = 150) {
  set.seed(3)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n), x3 = sample(0:3, n, TRUE))
  eta <- cbind(0, 3 * d$x1 - 1.5, 2 * (d$x2 > 0.5) + d$x3 / 2 - 1, 1 - 4 * d$x1 * d$x2)
  cls <- apply(eta, 1, function(e) sample.int(4, 1, prob = exp(e)))
  d$y <- factor(letters[pmin(cls, K)], levels = letters[seq_len(K)])
  d
}

# one family grown with its weights kept (internal hook): the data rows of its sample, their classes, their weights and
# the scores the forest predicts there
samme_family <- function(d, structure, deterministic) {
  set.seed(1)
  m <- rpf(
    y ~ .,
    data = d,
    loss = "exponential",
    max_interaction = 3,
    ntrees = 1,
    splits = 40,
    split_structure = structure,
    split_try = 3,
    t_try = 0.3,
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
    g = as.matrix(predict(m, d[rows, ], type = "link")),
    p = as.matrix(predict(m, d[rows, ], type = "prob"))
  )
}

test_that("multiclass exponential (SAMME): each row's weight is exp(-g_c) of its predicted scores", {
  structures <- c("leaves", "cur_trees_1", "cur_trees_2", "res_trees", "hist")
  for (K in 3:4) {
    d <- samme_data(K)
    for (s in structures) {
      for (det in c(FALSE, TRUE)) {
        lab <- sprintf("K = %d, %s, deterministic = %s", K, s, det)
        f <- samme_family(d, s, det)
        expect_identical(dim(f$w), c(nrow(d), K), label = lab)
        # one weight per row, kept in each of its K columns
        expect_true(all(f$w == f$w[, 1]), label = lab)
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

test_that("binary exponential: the weight is exp(-y g / 2) of the predicted score (master's code)", {
  d <- samme_data(2)
  for (det in c(FALSE, TRUE)) {
    f <- samme_family(d, "leaves", det)
    expect_identical(dim(f$w), c(nrow(d), 1L))
    y <- ifelse(f$class == 2, 1, -1)
    expect_lt(max(abs(f$w[, 1] / exp(-y * f$g[, 1] / 2) - 1)), 1e-12)
    expect_gt(max(abs(f$g)), 0.5)
  }
})

test_that("the multiclass exponential's first leaf step is SAMME's: log p_k - mean_j log p_j of the clamped shares", {
  d <- samme_data(3)
  set.seed(1)
  m <- rpf(
    y ~ .,
    data = d,
    loss = "exponential",
    max_interaction = 1,
    ntrees = 1,
    splits = 1,
    split_try = 3,
    t_try = 1,
    deterministic = TRUE,
    prune = FALSE
  )
  g <- as.matrix(predict(m, d, type = "link"))
  # the one split leaves two leaves, all weights 1 before it: the step from each leaf's class shares, clamped to
  # [epsilon, 1 - epsilon] (epsilon = 0.1)
  side <- factor(apply(round(g, 12), 1, paste, collapse = " "))
  expect_length(levels(side), 2)
  for (s in levels(side)) {
    in_leaf <- side == s
    p <- pmin(pmax(tabulate(as.integer(d$y[in_leaf]), 3) / sum(in_leaf), 0.1), 0.9)
    expect_equal(unname(g[which(in_leaf)[1], ]), log(p) - mean(log(p)), tolerance = 1e-14)
  }
})
