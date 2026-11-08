# Out-of-bag (OOB) components and prune_components(select_by = "oob")

oob_sim <- function(n = 400, seed = 11) {
  set.seed(seed)
  x <- matrix(runif(n * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
  # a main effect and a pure pair (no main effects of x2, x3); x4 is noise
  y <- 1 + 2 * x[, 1] + 4 * x[, 2] * x[, 3] + rnorm(n, sd = 0.3)
  list(x = x, y = y)
}

oob_fit <- function(d, max_interaction = 2, ntrees = 20, ...) {
  set.seed(5)
  rpf(d$x, d$y, max_interaction = max_interaction, ntrees = ntrees, splits = 40, ...)
}

# One family's purified component f_S at the rows of X, directly from its leaves (mode 2's formula):
# sum over the leaves of trees T containing S of c * sum_{W subset S} (-1)^(|S| - |W|) 1{x_W in box_W} q(T \ W),
# q(U) = the share of the training rows x_train inside the leaf's box on U.
family_component <- function(fam, S, X, x_train) {
  out <- numeric(nrow(X))
  inside <- function(Z, iv, U) {
    if (length(U) == 0) {
      return(rep(TRUE, nrow(Z)))
    }
    ok <- rep(TRUE, nrow(Z))
    for (u in U) {
      ok <- ok & Z[, u] >= iv[1, u] & Z[, u] < iv[2, u]
    }
    ok
  }
  for (t in seq_along(fam$variables)) {
    T <- fam$variables[[t]]
    if (identical(as.integer(T), 0L) || !all(S %in% T)) {
      next
    }
    for (l in seq_along(fam$values[[t]])) {
      iv <- fam$intervals[[t]][[l]]
      for (wm in 0:(2^length(S) - 1)) {
        W <- S[bitwAnd(wm, 2^(seq_along(S) - 1)) > 0]
        q <- mean(inside(x_train, iv, setdiff(T, W)))
        out <- out + (-1)^(length(S) - length(W)) * fam$values[[t]][[l]][1] * q * inside(X, iv, W)
      }
    }
  }
  out
}

# One family's prediction at the rows of X, from its leaves
family_prediction <- function(fam, X) {
  out <- numeric(nrow(X))
  for (t in seq_along(fam$variables)) {
    for (l in seq_along(fam$values[[t]])) {
      iv <- fam$intervals[[t]][[l]]
      ok <- rep(TRUE, nrow(X))
      for (u in seq_len(ncol(X))) {
        ok <- ok & X[, u] >= iv[1, u] & X[, u] < iv[2, u]
      }
      out <- out + fam$values[[t]][[l]][1] * ok
    }
  }
  out
}

test_that("the in-bag counts are the bootstrap samples", {
  d <- oob_sim(n = 150)
  for (structure in c("leaves", "res_trees", "cur_trees_1", "cur_trees_2", "hist")) {
    fit <- oob_fit(d, ntrees = 6, split_structure = structure)
    ib <- fit$fit$inbag()
    expect_identical(dim(ib), c(150L, 6L))
    expect_true(is.integer(ib))
    expect_true(all(colSums(ib) == 150))
    # about 37% of the rows are out of bag in each family
    expect_true(all(abs(colMeans(ib == 0) - exp(-1)) < 0.15))
  }
  det <- oob_fit(d, ntrees = 1, deterministic = TRUE)
  expect_true(all(det$fit$inbag() == 1))
  # classification records them too
  set.seed(5)
  cfit <- rpf(d$x, factor(d$y > 1), ntrees = 4, splits = 10)
  expect_true(all(colSums(cfit$fit$inbag()) == 150))
})

test_that("OOB rows are really out of bag: their outcome does not change the family", {
  d <- oob_sim(n = 150)
  fit <- oob_fit(d, ntrees = 8)
  ib <- fit$fit$inbag()
  # a row that is out of bag in some families and in bag in others
  i <- which(rowSums(ib == 0) >= 2 & rowSums(ib > 0) >= 2)[1]
  d2 <- d
  d2$y[i] <- d2$y[i] + 100
  fit2 <- oob_fit(d2, ntrees = 8)
  expect_identical(fit2$fit$inbag(), ib)
  m1 <- fit$fit$get_model()
  m2 <- fit2$fit$get_model()
  for (f in seq_along(m1)) {
    if (ib[i, f] == 0) {
      expect_identical(m2[[f]], m1[[f]])
    } else {
      expect_false(identical(m2[[f]], m1[[f]]))
    }
  }
})

test_that("the OOB sums equal a per-family evaluation from the leaves", {
  set.seed(3)
  n <- 120
  x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, paste0("x", 1:3)))
  y <- x[, 1] + 2 * x[, 2] * x[, 3] + rnorm(n, sd = 0.3)
  set.seed(11)
  fit <- rpf(x, y, max_interaction = 3, ntrees = 5, splits = 12)
  ib <- fit$fit$inbag()
  mod <- fit$fit$get_model()
  o <- fit$fit$purify_oob(0L, 1L)
  expect_named(o, c("sets", "sum0", "sum1", "intercept0", "intercept1", "count"))
  sets <- o$sets
  expect_equal(dim(o$sum0), c(n, length(sets)))

  brute0 <- brute1 <- matrix(0, n, length(sets))
  pred0 <- pred1 <- numeric(n)
  for (f in seq_along(mod)) {
    out <- ib[, f] == 0
    for (k in seq_along(sets)) {
      v <- family_component(mod[[f]], sets[[k]], x, x)
      if (f %% 2 == 1) {
        brute0[out, k] <- brute0[out, k] + v[out]
      } else {
        brute1[out, k] <- brute1[out, k] + v[out]
      }
    }
    p <- family_prediction(mod[[f]], x)
    if (f %% 2 == 1) {
      pred0[out] <- pred0[out] + p[out]
    } else {
      pred1[out] <- pred1[out] + p[out]
    }
  }
  expect_equal(o$sum0, brute0, tolerance = 1e-12)
  expect_equal(o$sum1, brute1, tolerance = 1e-12)
  expect_identical(o$count[, 1], as.integer(rowSums(ib[, c(1, 3, 5)] == 0)))
  expect_identical(o$count[, 2], as.integer(rowSums(ib[, c(2, 4)] == 0)))
  # intercept + components of each half = the half's family predictions (the families' own decomposition)
  expect_equal(o$intercept0[, 1] + rowSums(o$sum0), pred0, tolerance = 1e-12)
  expect_equal(o$intercept1[, 1] + rowSums(o$sum1), pred1, tolerance = 1e-12)
  # all families together are the forest's components
  set.seed(11)
  fit_all <- rpf(x, y, max_interaction = 3, ntrees = 5, splits = 12)
  purify(fit_all)
  brute_all <- sapply(seq_along(sets), function(k) {
    Reduce(`+`, lapply(mod, function(fam) family_component(fam, sets[[k]], x, x))) / length(mod)
  })
  expect_equal(fit_all$fit$predict_components_matrix(x, sets, 1L), brute_all, tolerance = 1e-12)
})

test_that("purify_oob() purifies the forest exactly as purify(), on any number of threads", {
  d <- oob_sim()
  a <- oob_fit(d, max_interaction = 3)
  b <- oob_fit(d, max_interaction = 3)
  c2 <- oob_fit(d, max_interaction = 3)
  o1 <- a$fit$purify_oob(0L, 1L)
  purify(b)
  o2 <- c2$fit$purify_oob(0L, 2L)
  expect_true(is_purified(a))
  expect_false(a$fit$is_purify_capped())
  sets <- b$fit$purified_components()
  expect_identical(a$fit$purified_components(), sets)
  # the same store (the per-subset counts are the k-d tree's)
  expect_identical(a$fit$get_purified_rects(), b$fit$get_purified_rects())
  expect_identical(a$fit$predict_components_matrix(d$x, sets, 1L), b$fit$predict_components_matrix(d$x, sets, 1L))
  expect_identical(a$fit$predict_matrix(d$x, 0, 1L), b$fit$predict_matrix(d$x, 0, 1L))
  expect_identical(o2, o1)
  expect_identical(c2$fit$get_purified_rects(), b$fit$get_purified_rects())
  # the stored components are the OOB components' sets
  expect_identical(o1$sets, sets)
  # a purification capped by maxp_interaction is capped here too, and replaced by the next one
  o_cap <- a$fit$purify_oob(1L, 1L)
  expect_true(a$fit$is_purify_capped())
  expect_true(all(lengths(o_cap$sets) == 1))
  a$fit$purify_oob(0L, 1L)
  expect_identical(a$fit$get_purified_rects(), b$fit$get_purified_rects())
})

test_that("OOB pruning keeps the true components and scores the candidates as documented", {
  d <- oob_sim()
  # with 20 families (7 OOB families per row) the OOB components are too noisy to tell R(1%) from R(0) here
  oob_fit <- function(d) {
    set.seed(5)
    rpf(d$x, d$y, max_interaction = 2, ntrees = 100, splits = 40)
  }
  fit <- oob_fit(d)
  fold_id <- rep_len(1:5, nrow(d$x))
  # the default: the largest tau within one SE of the lowest corrected OOB MSE
  pr <- prune_components(fit, d$x, d$y, select_by = "oob", fold_id = fold_id, seed = 3)
  expect_s3_class(pr, "rpf_pruned")
  expect_identical(pr$select, "1se")
  expect_identical(pr$select_by, "oob")
  expect_null(pr$cv)
  expect_setequal(names(pr$beta), c("x1", "x2:x3"))
  expect_equal(nrow(pr$oob), 16)
  expect_equal(attr(pr$oob, "n_rows"), nrow(d$x))
  expect_equal(sum(pr$oob$min), 1)
  expect_equal(sum(pr$oob$one_se), 1)
  chosen <- pr$oob[pr$oob$one_se, ]
  expect_equal(c(chosen$rule, chosen$tau), c(pr$rule, pr$tau))
  best <- pr$oob[pr$oob$min, ]
  expect_equal(best$oob_rmse, min(pr$oob$oob_rmse))
  expect_lte(chosen$oob_rmse^2, best$oob_rmse^2 + best$se_mse)
  expect_gte(pr$tau, best$tau)
  expect_true(all(pr$oob$tau[pr$oob$oob_rmse^2 <= best$oob_rmse^2 + best$se_mse] <= pr$tau))

  # the final model is the chosen rule and tau applied to the forest, as without a choice (R's coefficients fitted on
  # the OOB components in both, below)
  direct <- prune_components(oob_fit(d), d$x, d$y, rule = pr$rule, tau = pr$tau)
  expect_identical(pr$components, direct$components)
  expect_identical(pr$sets, direct$sets)
  expect_identical(pr$beta, direct$beta)
  expect_identical(pr$intercept, direct$intercept)
  expect_identical(predict(pr, d$x), predict(direct, d$x))

  # the OOB RMSEs by hand: D = OOB intercept + kept OOB components; R = coefficients fitted on the OOB components
  # outside the fold, applied to the OOB components of the fold. Raw: the mean squared error; corrected: minus
  # (1/m - 1/B) v per row, v from the means over the even and the odd families
  f2 <- oob_fit(d)
  o <- f2$fit$purify_oob(0L, 1L)
  m0 <- o$count[, 1]
  m1 <- o$count[, 2]
  m <- m0 + m1
  expect_true(all(m0 > 0 & m1 > 0))
  nm <- vapply(o$sets, function(s) paste(colnames(d$x)[s], collapse = ":"), "")
  colnames(o$sum0) <- colnames(o$sum1) <- nm
  share <- stats::setNames(pr$components$share, pr$components$component)
  rmse <- function(g, g0, g1) {
    e2 <- (d$y - g)^2
    v <- (g0 - g1)^2 / (1 / m0 + 1 / m1)
    c(raw = sqrt(mean(e2)), corrected = sqrt(mean(e2 - (1 / m - 1 / 100) * v)))
  }
  for (tau in c(0, 0.01)) {
    kept <- names(share)[share >= tau]
    S0 <- o$sum0[, kept, drop = FALSE]
    S1 <- o$sum1[, kept, drop = FALSE]
    g_d <- rmse(
      (o$intercept0[, 1] + o$intercept1[, 1] + rowSums(S0) + rowSums(S1)) / m,
      (o$intercept0[, 1] + rowSums(S0)) / m0,
      (o$intercept1[, 1] + rowSums(S1)) / m1
    )
    row_d <- pr$oob$rule == "D" & pr$oob$tau == tau
    expect_equal(pr$oob$oob_rmse_raw[row_d], unname(g_d["raw"]), tolerance = 1e-12)
    expect_equal(pr$oob$oob_rmse[row_d], unname(g_d["corrected"]), tolerance = 1e-12)
    F_oob <- (S0 + S1) / m
    g <- g0 <- g1 <- numeric(nrow(d$x))
    for (k in 1:5) {
      te <- fold_id == k
      cf <- stats::lm.fit(cbind(1, F_oob[!te, , drop = FALSE]), d$y[!te])$coefficients
      g[te] <- cbind(1, F_oob[te, , drop = FALSE]) %*% cf
      g0[te] <- cbind(1, S0[te, , drop = FALSE] / m0[te]) %*% cf
      g1[te] <- cbind(1, S1[te, , drop = FALSE] / m1[te]) %*% cf
    }
    g_r <- rmse(g, g0, g1)
    row_r <- pr$oob$rule == "R" & pr$oob$tau == tau
    expect_equal(pr$oob$oob_rmse_raw[row_r], unname(g_r["raw"]), tolerance = 1e-10)
    expect_equal(pr$oob$oob_rmse[row_r], unname(g_r["corrected"]), tolerance = 1e-10)
  }
  expect_true(all(pr$oob$oob_rmse <= pr$oob$oob_rmse_raw))
  # here R is chosen: its coefficients are least squares of y on the kept OOB components of every row
  expect_identical(pr$rule, "R")
  kept <- names(pr$beta)
  F_oob <- (o$sum0[, kept, drop = FALSE] + o$sum1[, kept, drop = FALSE]) / m
  cf <- stats::lm.fit(cbind(1, F_oob), d$y)$coefficients
  expect_equal(unname(c(pr$intercept, pr$beta)), unname(cf), tolerance = 1e-10)
  expect_equal(pr$components$beta[pr$components$kept], unname(pr$beta[pr$components$component[pr$components$kept]]))
  expect_output(print(pr), "Chosen by out-of-bag error \\(1-SE rule; R refitted in 5 folds, seed 3\\)")
  expect_output(print(pr), "fitted out of bag")

  # select = "min": the same table, the lowest corrected OOB RMSE (here more components)
  pr_min <- prune_components(oob_fit(d), d$x, d$y, fold_id = fold_id, seed = 3, select = "min")
  expect_identical(pr_min$oob, pr$oob)
  expect_equal(c(pr_min$rule, pr_min$tau), c(best$rule, best$tau))
  expect_gte(length(pr_min$beta), length(pr$beta))
  expect_output(print(pr_min), "lowest OOB RMSE")
})

test_that("the least squares of the OOB folds equal the final model's", {
  set.seed(2)
  X <- cbind(1, matrix(rnorm(300 * 20), 300, 20))
  y <- rnorm(300)
  ns <- asNamespace("randomPlantedForest")
  expect_equal(ns$prune_lstsq_qr(X, y), ns$prune_lstsq(X, y), tolerance = 1e-12)
  # rank-deficient (a column twice) and wide: the minimum-norm solution
  Xd <- cbind(X, X[, 2])
  expect_identical(ns$prune_lstsq_qr(Xd, y), ns$prune_lstsq(Xd, y))
  expect_identical(ns$prune_lstsq_qr(X[1:10, ], y[1:10]), ns$prune_lstsq(X[1:10, ], y[1:10]))
})

test_that("OOB pruning does not depend on nthreads, is the default, and CV is still available", {
  d <- oob_sim()
  fold_id <- rep_len(1:3, nrow(d$x))
  p1 <- prune_components(oob_fit(d), d$x, d$y, select_by = "oob", fold_id = fold_id, seed = 3, nthreads = 1)
  p2 <- prune_components(oob_fit(d), d$x, d$y, select_by = "oob", fold_id = fold_id, seed = 3, nthreads = 2)
  expect_identical(p2$oob, p1$oob)
  expect_identical(p2$beta, p1$beta)
  expect_identical(predict(p2, d$x), predict(p1, d$x))
  # the default is OOB
  def <- prune_components(oob_fit(d), d$x, d$y, fold_id = fold_id, seed = 3)
  expect_identical(def$select_by, "oob")
  expect_identical(def$oob, p1$oob)
  cv <- prune_components(oob_fit(d), d$x, d$y, select_by = "cv", fold_id = fold_id, seed = 3)
  expect_identical(cv$select_by, "cv")
  expect_null(cv$oob)
  expect_equal(nrow(cv$cv), 16)
  # a given rule and tau: no choice
  given <- prune_components(oob_fit(d), d$x, d$y, rule = "D", tau = 0.01)
  expect_null(given$select_by)
  expect_null(given$oob)
  expect_output(print(given), "Given, not chosen")
})

test_that("OOB pruning gives clear errors; a deterministic forest falls back to CV by default", {
  d <- oob_sim(n = 150)
  fit <- oob_fit(d, ntrees = 5)
  expect_error(prune_components(fit, d$x[150:1, ], d$y[150:1], select_by = "oob"), "training rows")
  det <- oob_fit(d, ntrees = 1, deterministic = TRUE)
  expect_error(prune_components(det, d$x, d$y, select_by = "oob"), "deterministic")
  expect_silent(pr <- prune_components(det, d$x, d$y, folds = 3, seed = 1))
  expect_identical(pr$select_by, "cv")
  expect_false(is.null(pr$cv))
  expect_error(prune_components(fit, d$x, d$y, select_by = "bag"))
})

test_that("a forest restored without in-bag counts falls back to CV by default", {
  d <- oob_sim(n = 200)
  fit <- oob_fit(d)
  blob <- rpf_marshal(fit, include_data = TRUE)
  ref <- prune_components(rpf_unmarshal(blob), d$x, d$y, folds = 3, seed = 1)
  expect_identical(ref$select_by, "oob")
  # a blob with data made before the in-bag counts were stored
  blob$fit_state$inbag <- NULL
  expect_message(pr <- prune_components(rpf_unmarshal(blob), d$x, d$y, folds = 3, seed = 1), "no in-bag counts")
  expect_identical(pr$select_by, "cv")
  expect_error(
    prune_components(rpf_unmarshal(blob), d$x, d$y, select_by = "oob", folds = 3, seed = 1),
    "no in-bag counts"
  )
})

test_that("OOB pruning reports its progress", {
  d <- oob_sim(n = 200)
  expect_message(
    prune_components(oob_fit(d), d$x, d$y, folds = 3, seed = 1, verbose = TRUE),
    "Out of bag: [0-9]+ of 200 rows with out-of-bag families \\([0-9.]+ per row\\), [0-9]+ nonzero components"
  )
})
