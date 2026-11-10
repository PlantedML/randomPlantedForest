# prune_components() without the purified store: the components at the training rows from one pass over the families
# (purify_columns()), and the pruned model's store of its kept components only (store_components()). The forest is not
# purified.

pl_data <- function(n = 400, p = 5, seed = 1) {
  set.seed(seed)
  x <- matrix(runif(n * p, -1, 1), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  x[, 2] <- round(x[, 2], 1) # ties
  y <- x[, 1] + 2 * x[, 2] * x[, 3] + 3 * x[, 1] * x[, 4] * x[, 5] + rnorm(n, sd = 0.3)
  list(x = x, y = y, xt = matrix(runif(200 * p, -1.1, 1.1), 200, p, dimnames = list(NULL, colnames(x))))
}

pl_fit <- function(d, max_interaction = ncol(d$x), ntrees = 12) {
  set.seed(3)
  rpf_050(d$x, d$y, max_interaction = max_interaction, ntrees = ntrees, splits = 80, split_try = 5, t_try = 0.5)
}

test_that("purify_columns() gives purify_oob()'s sums and the purified components, and leaves the forest", {
  d <- pl_data()
  for (mi in c(2L, 5L)) {
    fit <- pl_fit(d, mi)
    cc <- fit$fit$purify_columns(0L, 3L, TRUE)
    expect_false(is_purified(fit))
    # threads do not matter
    expect_identical(fit$fit$purify_columns(0L, 1L, TRUE), cc)
    # without the out-of-bag sums: the same columns
    c0 <- fit$fit$purify_columns(0L, 2L, FALSE)
    expect_named(c0, c("sets", "own", "intercept"))
    expect_identical(c0, cc[names(c0)])
    # the out-of-bag sums of purify_oob() up to rounding
    o <- fit$fit$purify_oob(0L, 2L)
    expect_identical(cc$sets, o$sets)
    expect_identical(cc$count, o$count)
    for (k in c("sum0", "sum1", "intercept0", "intercept1")) {
      expect_lt(max(abs(cc[[k]] - o[[k]])), 1e-12)
    }
    # the components of the purified forest up to rounding, its intercept bitwise
    expect_identical(fit$fit$purified_components(), cc$sets)
    M <- fit$fit$predict_components_matrix(d$x, cc$sets, 1L)
    expect_lt(max(abs(cc$own - M)), 1e-13)
    expect_identical(cc$intercept, as.numeric(fit$fit$predict_matrix(d$x[1, 1, drop = FALSE], -1, 0L)))
    # they add up to the forest's prediction
    expect_lt(max(abs(cc$intercept + rowSums(cc$own) - predict(fit, d$x)$.pred)), 1e-12)
  }
})

test_that("store_components() holds the chosen components exactly as purify() stores them", {
  ns <- asNamespace("randomPlantedForest")
  d <- pl_data()
  fit <- pl_fit(d)
  sets <- lapply(fit$fit$purify_columns(0L, 1L, FALSE)$sets, as.integer)
  sets <- sets[order(lengths(sets))]
  pick <- sets[c(1, 3, 7, length(sets))]
  st <- fit$fit$store_components(pick, 2L)
  expect_false(is_purified(fit))
  expect_length(st$comps, length(pick))
  xp <- ns$rpf_store_build(st, 2L)
  expect_true(ns$rpf_store_valid(xp))
  purify(fit)
  for (X in list(d$x, d$xt)) {
    ref <- fit$fit$predict_components_matrix(X, pick, 1L)
    expect_identical(ns$rpf_store_eval(xp, X, pick, 1L), ref)
    expect_identical(ns$rpf_store_eval(xp, X, pick, 3L), ref)
    expect_identical(ns$rpf_store_eval(xp, X, rev(pick), 2L), ref[, rev(seq_along(pick))])
  }
  expect_identical(dim(ns$rpf_store_eval(xp, d$x[0, , drop = FALSE], pick, 1L)), c(0L, length(pick)))
  expect_error(ns$rpf_store_eval(xp, d$x, sets[2], 1L), "not in the store")
  expect_error(ns$rpf_store_eval(xp, d$x[, 1:3], pick, 1L), "wrong dimension")
  # a corrupt store is refused
  bad <- st
  bad$comps[[2]]$groups[[1]]$box[1, 1] <- -1L
  expect_error(ns$rpf_store_build(bad, 1L), "Corrupt")
  bad <- st
  bad$comps[[1]]$groups[[1]]$w <- cbind(bad$comps[[1]]$groups[[1]]$w, 0)
  expect_error(ns$rpf_store_build(bad, 1L), "Corrupt")
})

test_that("the pruned model stores its kept components only and predicts from them", {
  ns <- asNamespace("randomPlantedForest")
  d <- pl_data(n = 500)
  fit <- pl_fit(d, ntrees = 20)
  # the same forest, purified: the reference for the components
  ref_fit <- pl_fit(d, ntrees = 20)
  purify(ref_fit)
  for (P in c(Inf, 3)) {
    pr <- prune_components(fit, d$x, d$y, max_panels = P, seed = 4)
    expect_false(is_purified(fit))
    coarse <- names(pr$beta) %in% names(pr$display$coarse)
    stored <- vapply(pr$store$comps, function(z) paste(z$S, collapse = ","), "")
    expect_setequal(stored, vapply(pr$sets[!coarse], function(s) paste(sort(s), collapse = ","), ""))
    if (is.finite(P)) {
      expect_true(any(coarse))
    }
    # predict(): the stored components (as the purified forest's), the coarsened ones from their stored form
    M <- ref_fit$fit$predict_components_matrix(d$xt, pr$sets, 1L)
    for (k in which(coarse)) {
      M[, k] <- rpf_coarse_predict(pr$display$coarse[[names(pr$beta)[k]]]$stored, d$xt)
    }
    expect_identical(predict(pr, d$xt)$.pred, pr$intercept + rowSums(sweep(M, 2, pr$beta, "*")))
    expect_identical(unname(as.matrix(predict_components(pr, d$xt)$m)), sweep(M, 2, pr$beta, "*"))
    expect_identical(predict(pr, d$xt, nthreads = 3)$.pred, predict(pr, d$xt)$.pred)
    # the store is built once and kept
    expect_true(ns$rpf_store_valid(pr$cache$.store))
  }
  # the same choice and model as from the purified store: purify_oob()'s sums and the store's components
  pr <- prune_components(fit, d$x, d$y, max_panels = Inf, seed = 4)
  tr <- ns$prune_training_data(fit, d$x, d$y)
  old <- ns$prune_oob_score(
    ref_fit$fit$purify_oob(0L, 1L),
    ns$prune_decompose(ref_fit, tr$predictors, rep(TRUE, nrow(d$x))),
    colnames(d$x),
    d$y,
    pr$oob[c("rule", "tau")],
    pr$fold_id,
    20L
  )
  choice <- c("rule", "tau", "n_kept", "min", "one_se")
  expect_identical(pr$oob[choice], old$table[choice])
  expect_equal(pr$oob, old$table, tolerance = 1e-12)
  expect_identical(pr$components$component, names(old$dec$share))
  expect_equal(pr$components$share, unname(old$dec$share), tolerance = 1e-13)
  rf <- ns$prune_fit_rule(pr$rule, pr$tau, old$dec, d$y, oob = old$oob_fit)
  expect_equal(pr$beta, rf$beta, tolerance = 1e-12)
  expect_equal(pr$intercept, rf$alpha, tolerance = 1e-12)
})

test_that("CV, a given rule and deterministic forests are pruned without purifying the forest", {
  d <- pl_data(n = 300)
  fit <- pl_fit(d, 3L)
  for (args in list(
    list(select_by = "cv", folds = 2, seed = 1),
    list(rule = "D", tau = 0.01),
    list(rule = "R", tau = 0.01)
  )) {
    pr <- do.call(prune_components, c(list(fit, d$x, d$y, max_panels = 2), args))
    expect_false(is_purified(fit))
    expect_false(is.null(pr$store))
    expect_length(pr$store$comps, sum(!names(pr$beta) %in% names(pr$display$coarse)))
  }
  set.seed(1)
  det <- rpf_050(d$x, d$y, max_interaction = 3, ntrees = 1, splits = 30, deterministic = TRUE)
  pr <- prune_components(det, d$x, d$y, folds = 2, seed = 1)
  expect_identical(pr$select_by, "cv")
  expect_false(is_purified(det))
  # D at tau = 0 keeps every component: the forest's prediction
  expect_equal(
    predict(prune_components(det, d$x, d$y, rule = "D", tau = 0, max_panels = Inf), d$xt)$.pred,
    predict(det, d$xt)$.pred,
    tolerance = 1e-12
  )
  # without in-bag counts (a blob of an earlier version): CV, and the columns need none
  blob <- rpf_marshal(fit, include_data = TRUE)
  blob$fit_state$inbag <- NULL
  old <- rpf_unmarshal(blob)
  expect_error(old$fit$purify_columns(0L, 1L, TRUE), "no in-bag counts")
  expect_named(old$fit$purify_columns(0L, 1L, FALSE), c("sets", "own", "intercept"))
  expect_message(pr <- prune_components(old, d$x, d$y, folds = 2, seed = 1), "no in-bag counts")
  expect_identical(pr$select_by, "cv")
  expect_identical(pr$beta, prune_components(fit, d$x, d$y, select_by = "cv", folds = 2, seed = 1)$beta)
})
