# predict_components() evaluates only the components the purified forest stores, all in one pass
# (predict_components_matrix); the others are zero columns. Every column must equal the one-at-a-time extraction.

op_data <- function(n = 250) {
  set.seed(5201)
  x <- data.frame(
    x1 = runif(n),
    x2 = sample(0:6, n, replace = TRUE),
    x3 = rnorm(n),
    x4 = runif(n),
    x5 = runif(n)
  )
  x$y <- x$x1 + 0.5 * x$x2 * x$x3 - 0.3 * x$x4 + rnorm(n, sd = 0.5)
  x
}

# the components one at a time, as predict_components() extracted them before
single_components <- function(fit, new_data) {
  X <- preprocess_predictors_predict(fit, hardhat::forge(new_data, fit$blueprint)$predictors)
  preds <- names(fit$blueprint$ptypes$predictors)
  sets <- unlist(
    lapply(seq_len(fit$params$max_interaction), function(i) utils::combn(preds, i, simplify = FALSE)),
    recursive = FALSE
  )
  do.call(cbind, lapply(sets, function(s) .predict_single_component(fit, X, s)))
}

test_that("one pass equals one component at a time: mode 2, mode 1, unpurified, with absent components", {
  dat <- op_data()
  nd <- dat[sample(nrow(dat), 120), 1:5]
  nd[1:4, "x3"] <- c(-9, 9, -9, 9)
  # few splits: some pairs and triples are never grown
  set.seed(52)
  fit <- rpf_050(y ~ ., data = dat, max_interaction = 3, n_families = 10, n_splits = 10)
  X <- preprocess_predictors_predict(fit, hardhat::forge(nd, fit$blueprint)$predictors)
  sets <- unlist(lapply(1:3, function(i) utils::combn(5, i, simplify = FALSE)), recursive = FALSE)

  # unpurified: the leaf scan of each set
  m0 <- fit$fit$predict_components_matrix(X, sets, 0L)
  for (k in seq_along(sets)) {
    expect_identical(m0[, k], as.numeric(fit$fit$predict_matrix(X[, sets[[k]], drop = FALSE], sets[[k]], 0L)))
  }
  expect_error(fit$fit$purified_components(), "not purified")

  for (mode in c(2L, 1L)) {
    fit$fit$purify_threads(0L, 0L, mode)
    stored <- vapply(fit$fit$purified_components(), paste, character(1), collapse = ",")
    # the stored components are among the candidate sets (which ones the forest grows depends on the platform)
    expect_true(all(stored %in% vapply(sets, paste, character(1), collapse = ",")))

    pc <- predict_components(fit, nd)
    single <- single_components(fit, nd)
    expect_identical(colnames(pc$m), colnames(single))
    expect_identical(unname(as.matrix(pc$m)), unname(single))
    expect_equal(ncol(pc$m), sum(choose(5, 1:3)))

    # every subset the forest does not store is a zero column
    absent <- !vapply(sets, paste, character(1), collapse = ",") %in% stored
    expect_true(any(absent))
    expect_true(all(as.matrix(pc$m)[, absent] == 0))
    expect_equal(rowSums(pc$m) + pc$intercept, predict(fit, nd)$.pred, tolerance = 1e-10)
  }
})

test_that("one pass equals one component at a time for multiclass, with 1 and 2 threads", {
  dat <- op_data()
  dat$y <- factor(cut(dat$y, 3, labels = c("a", "b", "c")))
  for (loss in c("logit", "L2")) {
    set.seed(53)
    fit <- rpf_050(y ~ ., data = dat, max_interaction = 2, n_families = 8, n_splits = 30, loss = loss)
    purify(fit, nthreads = 1L)
    pc1 <- predict_components(fit, dat[, 1:5])
    single <- single_components(fit, dat[, 1:5])
    expect_identical(colnames(pc1$m), colnames(single))
    expect_identical(unname(as.matrix(pc1$m)), unname(single))
    # one column per level (logit's scores centred across the levels)
    expect_equal(ncol(pc1$m), length(pc1$target_levels) * sum(choose(5, 1:2)))
    expect_identical(pc1$target_levels, c("a", "b", "c"))

    big <- dat[rep(seq_len(nrow(dat)), 20), 1:5] # several chunks of rows
    X <- preprocess_predictors_predict(fit, hardhat::forge(big, fit$blueprint)$predictors)
    sets <- unlist(lapply(1:2, function(i) utils::combn(5, i, simplify = FALSE)), recursive = FALSE)
    expect_identical(
      fit$fit$predict_components_matrix(X, sets, 1L),
      fit$fit$predict_components_matrix(X, sets, 2L)
    )
  }
})

test_that("predict_components_matrix() rejects bad sets", {
  dat <- op_data()
  set.seed(55)
  fit <- rpf_050(y ~ ., data = dat, max_interaction = 2, n_families = 5, n_splits = 20)
  X <- preprocess_predictors_predict(fit, hardhat::forge(dat[1:5, 1:5], fit$blueprint)$predictors)
  expect_error(fit$fit$predict_components_matrix(X, list(c(1L, 6L)), 0L), "distinct features")
  expect_error(fit$fit$predict_components_matrix(X, list(c(2L, 2L)), 0L), "distinct features")
  expect_error(fit$fit$predict_components_matrix(X[, 1:4], list(1L), 0L), "wrong dimension")
  expect_identical(dim(fit$fit$predict_components_matrix(X, list(), 0L)), c(5L, 0L))
})
