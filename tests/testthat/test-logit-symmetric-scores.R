# Multiclass logit fits K - 1 scores against the first level (whose score is 0). Every output reports one score per
# level, the first level's 0 joined and centred across the K levels, so that they sum to 0, as multiclass exponential
# reports its scores. The fit, its probabilities and classes, and its pruned model are those of the K - 1 scores.

sym_data <- function(n = 300, seed = 7) {
  set.seed(seed)
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  eta <- cbind(0, 3 * x$x1 - 1.5, 2 * (x$x2 - 0.5) + 2 * x$x3 * x$x1)
  cls <- apply(eta, 1, function(e) sample.int(3, 1, prob = exp(e)))
  x$y <- factor(c("a", "b", "c")[cls])
  x
}

# the fit's K - 1 scores S, reported: 0 joined, centred across the K levels
sym_centred <- function(S) {
  S <- cbind(0, as.matrix(S))
  unname(S - rowMeans(S))
}

# the probabilities of the fit's K - 1 scores, as the forest computes them
sym_prob_ref <- function(S) {
  E <- exp(as.matrix(S))
  unname(cbind(1, E) / (1 + rowSums(E)))
}

sym_raw <- function(fit, d) {
  X <- preprocess_predictors_predict(fit, hardhat::forge(d, fit$blueprint)$predictors)
  fit$fit$predict_matrix(X, 0, 0L)
}

test_that("multiclass logit: the link has one score per level, summing to 0; probabilities and classes as before", {
  d <- sym_data()
  set.seed(1)
  fit <- rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 10, n_splits = 40, prune = FALSE)
  raw <- sym_raw(fit, d)
  expect_identical(ncol(raw), 2L)

  link <- predict(fit, d, type = "link")
  expect_named(link, c(".pred_a", ".pred_b", ".pred_c"))
  expect_identical(predict(fit, d, type = "numeric"), link)
  L <- unname(as.matrix(link))
  expect_lt(max(abs(rowSums(L))), 1e-12)
  expect_equal(L, sym_centred(raw), tolerance = 1e-15)

  prob <- unname(as.matrix(predict(fit, d, type = "prob")))
  expect_equal(prob, sym_prob_ref(raw), tolerance = 1e-12)
  expect_equal(prob, exp(L) / rowSums(exp(L)), tolerance = 1e-12)
  expect_identical(
    predict(fit, d, type = "class")$.pred_class,
    factor(c("a", "b", "c")[max.col(prob, ties.method = "first")], levels = c("a", "b", "c"))
  )

  # zero rows: the reported columns
  expect_identical(dim(predict(fit, d[0, ], type = "link")), c(0L, 3L))
})

test_that("multiclass logit: components and intercepts per level sum to 0 and add up to the link", {
  d <- sym_data()
  set.seed(2)
  fit <- rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 10, n_splits = 40, prune = FALSE)
  purify(fit)
  link <- unname(as.matrix(predict(fit, d, type = "link")))
  pc <- predict_components(fit, d)
  expect_identical(pc$target_levels, c("a", "b", "c"))
  expect_named(pc$intercept, c("a", "b", "c"))
  expect_lt(abs(sum(pc$intercept)), 1e-12)
  M <- as.matrix(pc$m)
  expect_identical(colnames(M)[1:3], paste0("x1__class:", c("a", "b", "c")))
  comps <- unique(sub("__class:.*$", "", colnames(M)))
  for (k in seq_along(comps)) {
    expect_lt(max(abs(rowSums(M[, (k - 1) * 3 + 1:3]))), 1e-12)
  }
  recon <- sapply(1:3, function(v) pc$intercept[v] + rowSums(M[, seq(v, ncol(M), by = 3), drop = FALSE]))
  expect_equal(unname(recon), link, tolerance = 1e-10)

  # the K - 1 score columns' components, reported
  X <- preprocess_predictors_predict(fit, hardhat::forge(d, fit$blueprint)$predictors)
  raw <- fit$fit$predict_components_matrix(X, list(1L), 0L)
  expect_equal(unname(M[, 1:3]), sym_centred(raw), tolerance = 1e-15)

  # left-out components: the remainder per level, so that the three add up to the link
  pc1 <- predict_components(fit, d, max_interaction = 1)
  expect_named(pc1$remainder, c("a", "b", "c"))
  M1 <- as.matrix(pc1$m)
  R1 <- as.matrix(pc1$remainder)
  recon1 <- sapply(1:3, function(v) pc1$intercept[v] + rowSums(M1[, seq(v, ncol(M1), by = 3)]) + R1[, v])
  expect_equal(unname(recon1), link, tolerance = 1e-10)
  expect_lt(max(abs(rowSums(R1))), 1e-12)

  # zero rows
  pc0 <- predict_components(fit, d[0, ], max_interaction = 1)
  expect_identical(ncol(pc0$m), 9L)
  expect_named(pc0$remainder, c("a", "b", "c"))
})

test_that("multiclass logit: the pruned model reports one score per level; its probabilities are its K - 1 scores'", {
  d <- sym_data()
  set.seed(3)
  fit <- rpf(y ~ ., data = d, loss = "logit", max_interaction = 2, n_families = 10, n_splits = 40, prune = FALSE)
  pr <- prune_components(fit, d, seed = 1)
  expect_length(pr$intercept, 2L) # stored: the fit's score columns

  raw <- pruned_scores(pr, prune_kept_matrix(pr, d))
  link <- predict(pr, d, type = "link")
  expect_named(link, c(".pred_a", ".pred_b", ".pred_c"))
  L <- unname(as.matrix(link))
  expect_lt(max(abs(rowSums(L))), 1e-12)
  expect_equal(L, sym_centred(raw), tolerance = 1e-15)
  expect_equal(unname(as.matrix(predict(pr, d, type = "prob"))), sym_prob_ref(raw), tolerance = 1e-12)

  pc <- predict_components(pr, d)
  expect_identical(pc$target_levels, c("a", "b", "c"))
  expect_lt(abs(sum(pc$intercept)), 1e-12)
  M <- as.matrix(pc$m)
  expect_identical(ncol(M), 3L * length(pr$beta))
  expect_identical(colnames(M), paste0(rep(names(pr$beta), each = 3), "__class:", c("a", "b", "c")))
  recon <- sapply(1:3, function(v) pc$intercept[v] + rowSums(M[, seq(v, ncol(M), by = 3), drop = FALSE]))
  expect_equal(unname(recon), L, tolerance = 1e-10)

  # only the main effects: the rest is the remainder, per level
  pc1 <- predict_components(pr, d, max_interaction = 1)
  if (!is.null(pc1$remainder)) {
    M1 <- as.matrix(pc1$m)
    R1 <- as.matrix(pc1$remainder)
    expect_named(pc1$remainder, c("a", "b", "c"))
    recon1 <- sapply(1:3, function(v) {
      pc1$intercept[v] + rowSums(M1[, seq(v, ncol(M1), by = 3), drop = FALSE]) + R1[, v]
    })
    expect_equal(unname(recon1), L, tolerance = 1e-10)
  }

  out <- paste(format(pr), collapse = "\n")
  ic <- format(as.numeric(sym_centred(matrix(pr$intercept, 1))), digits = 6)
  expect_match(out, paste0("Intercept per class: a ", ic[1], ", b ", ic[2], ", c ", ic[3]), fixed = TRUE)
})

test_that("binary logit and multiclass exponential report as before", {
  d <- sym_data()
  d2 <- d
  d2$y <- factor(ifelse(d$y == "a", "a", "b"))
  set.seed(4)
  fit2 <- rpf(y ~ ., data = d2, loss = "logit", max_interaction = 2, n_families = 5, n_splits = 20, prune = FALSE)
  link2 <- predict(fit2, d2, type = "link")
  expect_named(link2, ".pred")
  expect_identical(unname(as.matrix(link2)), unname(sym_raw(fit2, d2)))

  set.seed(5)
  fite <- rpf(y ~ ., data = d, loss = "exponential", max_interaction = 2, n_families = 5, n_splits = 20, prune = FALSE)
  linke <- predict(fite, d, type = "link")
  expect_named(linke, c(".pred_a", ".pred_b", ".pred_c"))
  expect_identical(unname(as.matrix(linke)), unname(sym_raw(fite, d)))
})
