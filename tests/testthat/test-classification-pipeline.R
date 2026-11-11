# Classification through the pipeline: the OOB path and n_splits = NULL scored by the loss's criterion, the pruning
# (centred shares, R with one beta per component shared by the classes), the display on V score columns, the pruned
# model's predict types and components, and the tuner.

cp_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n))
  eta <- 2 * (x$x1 - 0.5) + 6 * (x$x2 - 0.5) * (x$x3 - 0.5) + 16 * (x$x2 - 0.5) * (x$x3 - 0.5) * (x$x4 - 0.5)
  yb <- factor(ifelse(runif(n) < stats::plogis(2 * eta), "b", "a"))
  ym <- factor(c("a", "b", "c")[1 + (eta + rnorm(n, sd = 0.3) > -0.3) + (eta + rnorm(n, sd = 0.3) > 0.3)])
  list(x = x, yb = yb, ym = ym, y01 = as.numeric(yb == "b"))
}

cp_fit <- function(d, y, loss, n_families = 16, n_splits = 120, ...) {
  set.seed(3)
  rpf_050(d$x, y, loss = loss, max_interaction = 4, n_families = n_families, n_splits = n_splits, ...)
}

cp_ns <- function() asNamespace("randomPlantedForest")

# the scores of a pruned model from its components (per score column), as predict(type = "numeric") gives them
cp_recon <- function(m, x) {
  pc <- predict_components(m, x)
  M <- as.matrix(pc$m)
  V <- length(pc$intercept)
  sapply(seq_len(V), function(v) {
    cols <- if (V == 1) seq_len(ncol(M)) else seq(v, ncol(M), by = V)
    pc$intercept[v] + rowSums(M[, cols, drop = FALSE])
  })
}

test_that("L2 classification on 0/1 is the regression pipeline bitwise (OOB path, n_splits = NULL, pruning, display, tuner)", {
  skip_on_cran()
  d <- cp_data(500, seed = 11)
  # the exact draw (the default is the fast one), as this comparison was written for it
  old <- options(randomPlantedForest.exact_draw = TRUE)
  on.exit(options(old))
  pipe <- function(y, ...) {
    set.seed(3)
    m <- suppressWarnings(rpf(
      d$x,
      y,
      max_interaction = 4,
      n_families = 16,
      n_try_split_points = 3,
      share_candidates = 0.1,
      n_splits = NULL,
      prune = TRUE,
      max_panels = 2,
      ...
    ))
    list(
      n_splits = m$fit$params$splits,
      path = m$fit$oob_path,
      rule = m$rule,
      tau = m$tau,
      beta = m$beta,
      icpt = m$intercept,
      first = m$display$first[c("rule", "tau", "beta", "intercept", "oob")],
      coarse = lapply(m$display$coarse, function(z) z[c("S", "C", "q", "cuts", "var_kept")]),
      oob = m$oob,
      pred = as.numeric(as.matrix(predict(m, d$x, type = "numeric")))
    )
  }
  a <- pipe(d$yb, loss = "L2")
  b <- pipe(d$y01)
  expect_identical(a, b)
  # the display coarsened (4 of the first pruning's 13 components) and pruned again: what this example fit does, on the
  # platform it was chosen on (other platforms draw other forests)
  if (on_macos_arm64()) {
    expect_gt(length(a$coarse), 0)
    expect_identical(a$rule, "R")
  }
  # also with min_leaf_rows and the fast draw
  options(randomPlantedForest.exact_draw = FALSE)
  expect_identical(pipe(d$yb, loss = "L2", min_leaf_rows = 5), pipe(d$y01, min_leaf_rows = 5))
  options(randomPlantedForest.exact_draw = TRUE)
  # the OOB path of a fixed fit, and the tuner
  d <- cp_data(300, seed = 11)
  path <- function(y, ...) {
    set.seed(4)
    suppressWarnings(rpf_050(d$x, y, max_interaction = 3, n_families = 8, n_splits = 80, oob_path = TRUE, ...))$oob_path
  }
  expect_identical(path(d$yb, loss = "L2"), path(d$y01))
  # the tuner (two ladder points, two leaf sizes)
  tun <- function(y, ...) {
    set.seed(4)
    m <- suppressWarnings(rpf(
      d$x,
      y,
      max_interaction = 3,
      n_families = 12,
      tune = tune_control(search = c(5, 7), min_leaf_rows = c(1, 5), neighbours = FALSE),
      ...
    ))
    cols <- c(
      "n_try_split_points",
      "share_candidates",
      "min_leaf_rows",
      "n_splits",
      "rule",
      "tau",
      "n_kept",
      "oob_rmse_forest",
      "oob_rmse"
    )
    list(tab = m$tune$search[, cols], pred = as.numeric(as.matrix(predict(m, d$x, type = "numeric"))))
  }
  expect_identical(tun(d$yb, loss = "L2"), tun(d$y01))
})

test_that("L2 classification on 0/1 is the regression fit bitwise where candidate splits tie in loss", {
  skip_on_cran()
  # Splits of 0/1 data often tie in loss: in a leaf whose rows all have the same residual every threshold scores the
  # same, and others tie to an ulp. Both modes score a split with the same arithmetic, so they break such a tie alike.
  # Before they did not (classification computed the equal -2 M S + n M^2): with these data the OOB paths of the
  # default structure parted for 189 of the seeds 1:300 on macOS arm64, and for other seeds on other platforms. Each
  # case below parted there before (structure, seed, min_leaf_rows, draw).
  d <- cp_data(300, seed = 11)
  path <- function(y, seed, ...) {
    set.seed(seed)
    suppressWarnings(rpf_050(d$x, y, max_interaction = 3, n_families = 8, n_splits = 80, oob_path = TRUE, ...))$oob_path
  }
  cases <- list(
    list("leaf", 6, 1, TRUE),
    list("leaf", 38, 1, TRUE),
    list("leaf", 4, 5, TRUE),
    list("leaf_binned", 33, 1, TRUE),
    list("tree_then_leaf", 70, 1, TRUE),
    list("tree_then_all_leaves", 44, 1, TRUE),
    list("target_tree_then_leaf", 55, 1, TRUE),
    list("leaf", 6, 1, FALSE),
    list("leaf", 38, 5, FALSE)
  )
  old <- options(randomPlantedForest.exact_draw = TRUE)
  on.exit(options(old))
  for (cs in cases) {
    options(randomPlantedForest.exact_draw = cs[[4]])
    mlr <- if (cs[[3]] > 1) list(min_leaf_rows = cs[[3]]) else list()
    a <- do.call(path, c(list(d$yb, cs[[2]], loss = "L2", split_structure = cs[[1]]), mlr))
    b <- do.call(path, c(list(d$y01, cs[[2]], split_structure = cs[[1]]), mlr))
    expect_identical(a, b, info = paste(unlist(cs), collapse = " "))
  }
})

test_that("components + intercept = predict per class; probabilities valid; class = argmax (all losses)", {
  skip_on_cran()
  d <- cp_data()
  for (y in list(d$yb, d$ym)) {
    for (loss in c("L1", "L2", "logit", "exponential")) {
      fit <- cp_fit(d, y, loss)
      first <- prune_components(fit, d$x, y, max_panels = Inf, seed = 1)
      shown <- prune_components(fit, d$x, y, max_panels = 2, seed = 1)
      given <- prune_components(fit, d$x, y, rule = "R", tau = 0.01, max_panels = 2)
      for (m in list(first, shown, given)) {
        raw <- as.matrix(predict(m, d$x, type = "numeric"))
        expect_equal(unname(cp_recon(m, d$x)), unname(raw), tolerance = 1e-12)
        expect_identical(predict(m, d$x, type = "link"), predict(m, d$x, type = "numeric"))
        P <- as.matrix(predict(m, d$x))
        expect_identical(colnames(P), paste0(".pred_", levels(y)))
        expect_true(all(P >= 0 & P <= 1))
        expect_equal(unname(rowSums(P)), rep(1, nrow(P)), tolerance = 1e-12)
        cl <- predict(m, d$x, type = "class")$.pred_class
        expect_identical(as.character(cl), levels(y)[max.col(P, "first")])
        expect_identical(m$criterion, if (loss %in% c("L1", "L2")) "brier" else "logloss")
      }
      pc <- predict_components(shown, d$x)
      # one column per class and component (V = K for multiclass, logit too)
      V <- length(pc$intercept)
      expect_identical(ncol(pc$m), V * length(shown$beta))
      if (V > 1) {
        expect_identical(pc$target_levels, names(pc$intercept))
        expect_true(all(grepl("__class:", names(pc$m), fixed = TRUE)))
        # the remainder per class: the components left out
        pm <- predict_components(shown, d$x, max_interaction = 1)
        rec <- sapply(seq_len(V), function(v) {
          rowSums(as.matrix(pm$m)[, seq(v, ncol(pm$m), by = V), drop = FALSE]) + pm$intercept[v] + pm$remainder[[v]]
        })
        expect_equal(unname(rec), unname(as.matrix(predict(shown, d$x, type = "numeric"))), tolerance = 1e-12)
      }
      out <- format(shown)
      expect_true(any(grepl("Classification", out, fixed = TRUE)))
    }
  }
})

test_that("L1 / L2 multiclass: R's shared beta keeps the scores summing to 1", {
  d <- cp_data()
  for (loss in c("L1", "L2")) {
    fit <- cp_fit(d, d$ym, loss)
    pr <- prune_components(fit, d$x, d$ym, rule = "R", tau = c(0, 0.01), max_panels = Inf, seed = 1)
    expect_identical(pr$rule, "R")
    raw <- as.matrix(predict(pr, d$x, type = "numeric"))
    expect_equal(unname(rowSums(raw)), rep(1, nrow(raw)), tolerance = 1e-10)
    expect_length(pr$beta, sum(pr$components$kept))
  }
})

test_that("max_panels = Inf gives the first pruning's model; coarsened multiclass components have V forms", {
  skip_on_cran()
  d <- cp_data()
  ns <- cp_ns()
  for (y in list(d$yb, d$ym)) {
    for (loss in c("logit", "exponential")) {
      fit <- cp_fit(d, y, loss)
      pi <- prune_components(fit, d$x, y, max_panels = .Machine$integer.max, seed = 7)
      pr <- prune_components(fit, d$x, y, max_panels = Inf, seed = 7)
      expect_false(pi$display$second)
      expect_identical(pi$beta, pr$beta)
      expect_identical(pi$intercept, pr$intercept)
      expect_identical(predict(pi, d$x), predict(pr, d$x))
    }
  }
  # a coarsened 3-way of a multiclass forest: one panel form per component, V stored forms; its variance kept is the
  # centred variance (crit_centred_var(), as the pruning's shares) of the coarsened over the forest's component
  fit <- cp_fit(d, d$ym, "exponential", n_families = 20, n_splits = 300)
  pd <- prune_components(fit, d$x, d$ym, max_panels = 2, tau = c(0, 1e-4), select = "min", seed = 1)
  expect_gt(length(pd$display$coarse), 0)
  for (z in pd$display$coarse) {
    expect_length(z$stored, 3)
    expect_true(z$cells <= 2)
    expect_true(z$var_kept > 0 && z$var_kept <= 1 + 1e-12)
  }
  expect_equal(unname(cp_recon(pd, d$x)), unname(as.matrix(predict(pd, d$x, type = "numeric"))), tolerance = 1e-12)
  for (loss in c("logit", "exponential")) {
    fit <- cp_fit(d, d$ym, loss, n_families = 12, n_splits = 200)
    tr <- ns$prune_training_data(fit, d$x, d$ym)
    pcol <- ns$prune_columns(fit, names(d$x), 1L, oob = FALSE, crit = tr$crit)$dec
    hi <- which(lengths(pcol$sets[pcol$nonzero]) >= 3)
    skip_if(length(hi) == 0)
    S <- pcol$sets[pcol$nonzero][hi[1]]
    fv <- pcol$own[, ns$dec_cols(pcol, names(pcol$share)[hi[1]]), drop = FALSE]
    sel <- fit$fit$coarsen_select(S, 2L, 1L, fv)[[1]]
    expect_identical(sel$form, "coarsened")
    expect_equal(dim(sel$g), dim(fv))
    expect_equal(
      sel$var_kept,
      ns$crit_centred_var(tr$crit, sel$g) / ns$crit_centred_var(tr$crit, fv),
      tolerance = 1e-12
    )
    oo <- fit$fit$coarsen_oob(list(sel[c("S", "C", "cuts")]), 1L)
    expect_identical(dim(oo$sum0), c(nrow(d$x), ncol(fv)))
  }
})

test_that("threads do not change the classification pipeline", {
  skip_on_cran()
  # 600 rows: the three-class logit fit keeps a component of three predictors, which the display coarsens
  d <- cp_data(600)
  # the exact draw: the three-class fit below has coarsened components with it
  old <- options(randomPlantedForest.exact_draw = TRUE)
  on.exit(options(old))
  # n_splits = NULL (or a fixed fit with its OOB path), pruning and the display
  run <- function(y, loss, nt, n_splits = NULL) {
    set.seed(5)
    m <- rpf(
      d$x,
      y,
      loss = loss,
      max_interaction = 3,
      n_families = 10,
      n_splits = n_splits,
      prune = TRUE,
      max_panels = 2,
      nthreads = nt,
      oob_path = TRUE
    )
    list(
      m$fit$params$splits,
      m$fit$oob_path,
      m$rule,
      m$tau,
      m$beta,
      m$intercept,
      m$oob,
      lapply(m$display$coarse, function(z) z[setdiff(names(z), "seconds")]),
      predict(m, d$x, type = "numeric", nthreads = nt)
    )
  }
  expect_identical(run(d$yb, "exponential", 1L), run(d$yb, "exponential", 2L))
  expect_identical(run(d$ym, "exponential", 1L), run(d$ym, "exponential", 2L))
  a <- run(d$ym, "logit", 1L, n_splits = 80)
  expect_identical(a, run(d$ym, "logit", 2L, n_splits = 80))
  # the display coarsened, so threads did not change it either (whether it coarsens depends on the platform's forest)
  if (on_macos_arm64()) {
    expect_gt(length(a[[8]]), 0)
  }
})

test_that("rpf(prune = TRUE) and rpf(tune = ) of a classification forest are the pruning of rpf(prune = FALSE)", {
  skip_on_cran()
  d <- cp_data(300)
  for (y in list(d$yb, d$ym)) {
    set.seed(8)
    a <- rpf(d$x, y, max_interaction = 3, n_families = 10, n_splits = NULL, prune = TRUE, max_panels = 2)
    set.seed(8)
    f <- rpf(d$x, y, max_interaction = 3, n_families = 10, n_splits = NULL, prune = FALSE)
    b <- prune_components(f, d$x, y, max_panels = 2)
    expect_s3_class(a, "rpf_pruned")
    expect_identical(
      a[c("rule", "tau", "beta", "intercept", "oob", "criterion")],
      b[c("rule", "tau", "beta", "intercept", "oob", "criterion")]
    )
    expect_identical(predict(a, d$x), predict(b, d$x))
    set.seed(8)
    tu <- rpf(d$x, y, max_interaction = 3, n_families = 10, tune = "quick")
    expect_s3_class(tu, "rpf_pruned")
    expect_true(all(tu$tune$search$oob_rmse > 0))
    expect_identical(tu$criterion, "logloss")
  }
})

test_that("the criterion: C++ (the adaptive bound) and R agree; the squared error is the regression formula", {
  ns <- cp_ns()
  d <- cp_data(200)
  set.seed(1)
  for (y in list(d$yb, d$ym)) {
    for (loss in c("L2", "logit", "exponential")) {
      fit <- cp_fit(d, y, loss, n_families = 4, n_splits = 20)
      crit <- ns$oob_criterion(fit$mode, loss, y)
      V <- crit$V
      n <- nrow(d$x)
      G <- matrix(rnorm(n * V, sd = 1.5), n, V)
      D <- matrix(rnorm(n * V, sd = 0.3), n, V)
      m0 <- sample(1:3, n, TRUE)
      m1 <- sample(1:3, n, TRUE)
      inv <- 1 / m0 + 1 / m1
      coef <- 1 / (m0 + m1) - 1 / 4
      cpp <- fit$fit$oob_row_errors(G, D, inv, coef)
      # R's crit_rows() with the halves g0 = g + d / 2, g1 = g - d / 2 (their difference is d), B = 4
      r <- ns$crit_rows(crit, seq_len(n), G, G + D / 2, G - D / 2, m0, m1, 4, impute = FALSE)
      if (ns$crit_is_scalar_mse(crit)) {
        e <- (crit$y - G[, 1])^2
        expect_equal(r$corr, e - coef * D[, 1]^2 / inv, tolerance = 1e-12)
        expect_equal(cpp, r$corr, tolerance = 1e-12)
      } else {
        expect_equal(cpp, r$corr, tolerance = 1e-10)
      }
      # a row without d: no correction
      D[1, 1] <- NA
      expect_equal(fit$fit$oob_row_errors(G, D, inv, coef)[1], r$e[1], tolerance = 1e-12)
    }
  }
  # regression: the squared error with the correction
  set.seed(2)
  rf <- rpf_050(d$x, d$y01 + rnorm(nrow(d$x)), n_families = 2, n_splits = 5)
  n <- nrow(d$x)
  G <- matrix(rnorm(n), n, 1)
  D <- matrix(rnorm(n), n, 1)
  y <- rf$fit$get_data()$Y[, 1]
  expect_equal(
    rf$fit$oob_row_errors(G, D, rep(0.5, n), rep(0.1, n)),
    (y - G[, 1])^2 - 0.1 * (D[, 1]^2 / 0.5),
    tolerance = 1e-14
  )
})

test_that("the log-loss Hessian is the numerical one; logit's centred shares do not depend on the reference level", {
  ns <- cp_ns()
  set.seed(2)
  cases <- list(
    list(link = "sigmoid", V = 1L, K = 2L),
    list(link = "softmax_ref", V = 2L, K = 3L),
    list(link = "softmax", V = 3L, K = 3L)
  )
  for (case in cases) {
    crit <- structure(
      list(type = "logloss", link = case$link, V = case$V, K = case$K, cls = 2L, centre = "none", n = 1L),
      class = "rpf_crit"
    )
    g <- rnorm(case$V)
    f <- function(gg) ns$crit_logloss(crit, matrix(gg, 1), 2L)$e
    h <- 1e-4
    H <- matrix(0, case$V, case$V)
    for (i in seq_len(case$V)) {
      for (j in seq_len(case$V)) {
        ei <- replace(numeric(case$V), i, h)
        ej <- replace(numeric(case$V), j, h)
        H[i, j] <- (f(g + ei + ej) - f(g + ei - ej) - f(g - ei + ej) + f(g - ei - ej)) / (4 * h^2)
      }
    }
    p <- ns$crit_logloss(crit, matrix(g, 1), 2L)$P[1, ]
    expect_equal(H, diag(p, case$V) - tcrossprod(p), tolerance = 1e-5)
  }
  # the scores against level 1 and the same scores against level 3: the centred variance is the same
  crit <- structure(list(centre = "ref"), class = "rpf_crit")
  M1 <- matrix(rnorm(300), 100, 3) # logits of levels 2..4 against level 1
  M3 <- cbind(-M1[, 2], M1[, 1] - M1[, 2], M1[, 3] - M1[, 2]) # levels 1, 2, 4 against level 3
  expect_equal(ns$crit_centred_var(crit, M1), ns$crit_centred_var(crit, M3), tolerance = 1e-12)
})

test_that("n_splits = NULL does not over-grow on a flat out-of-bag curve (three-class logit)", {
  skip_on_cran()
  # A family pauses at most one doubling past the first family's estimate and grows on when the stop comes later, so
  # on a long, flat out-of-bag curve each grows to about the stop (before that fix, the estimate from part of the
  # families kept finding new minima, and the families grew to about 2.8 times the stop). The curve here is long and
  # flat on macOS arm64 at the greediness 1 / 0.025: the stop on a new minimum alone stops at 200 splits and chooses
  # 100, the plateau stop (the default) at 141 and 100. (The logit loss of N18 stops after 15 splits at 10 / 0.4,
  # where every family grows past its stop by the pause's doubling.) Other platforms grow other forests, with other
  # curves: there the bounds are checked at their stops
  d <- cp_data(300)
  old <- options(randomPlantedForest.exact_draw = TRUE, randomPlantedForest.adaptive_tolerance = 0)
  on.exit(options(old))
  set.seed(8)
  f <- rpf(
    d$x,
    d$ym,
    loss = "logit",
    max_interaction = 3,
    n_families = 10,
    n_try_split_points = 1,
    share_candidates = 0.025,
    prune = FALSE
  )
  st <- f$fit$adaptive_stats()
  stop_f <- f$splits_oob$stopped_at
  expect_lte(f$params$splits, stop_f)
  expect_lte(st$splits_grown, 1.2 * 10 * stop_f)
  # the same forest as the fixed fit at the chosen splits
  set.seed(8)
  g <- rpf(
    d$x,
    d$ym,
    loss = "logit",
    max_interaction = 3,
    n_families = 10,
    n_try_split_points = 1,
    share_candidates = 0.025,
    n_splits = f$params$splits,
    prune = FALSE
  )
  expect_identical(f$fit$get_model(), g$fit$get_model())
  # the plateau stop (the default tolerance): the families grow only about as far as its stop
  options(randomPlantedForest.adaptive_tolerance = 1)
  set.seed(8)
  p <- rpf(
    d$x,
    d$ym,
    loss = "logit",
    max_interaction = 3,
    n_families = 10,
    n_try_split_points = 1,
    share_candidates = 0.025,
    prune = FALSE
  )
  expect_lte(p$splits_oob$n_splits, p$splits_oob$stopped_at)
  expect_lte(p$fit$adaptive_stats()$splits_grown, 1.2 * 10 * p$splits_oob$stopped_at)
  # the flat curve of the comment above
  if (on_macos_arm64()) {
    expect_identical(f$params$splits, 100L)
    expect_identical(stop_f, 200L)
    expect_identical(p$splits_oob[c("n_splits", "stopped_at")], list(n_splits = 100L, stopped_at = 141L))
  }
})

test_that("the OOB path of a classification forest: prefix property, n_splits = NULL = the fixed fit", {
  skip_on_cran()
  d <- cp_data(300)
  tot <- function(f) {
    o <- f$fit$purify_oob(0L, 1L)
    V <- ncol(o$intercept0)
    S <- o$intercept0 + o$intercept1
    for (j in seq_len(ncol(o$sum0) / V)) {
      for (v in seq_len(V)) {
        S[, v] <- S[, v] + o$sum0[, (j - 1) * V + v] + o$sum1[, (j - 1) * V + v]
      }
    }
    S
  }
  for (y in list(d$yb, d$ym)) {
    for (loss in c("L2", "logit", "exponential")) {
      for (s in c("leaf", "target_tree_then_leaf", "tree_then_leaf")) {
        set.seed(5)
        big <- rpf_050(
          d$x,
          y,
          loss = loss,
          max_interaction = 3,
          n_families = 6,
          n_splits = 60,
          oob_path = TRUE,
          split_structure = s
        )
        if (loss != "L2" || nlevels(y) > 2) {
          expect_identical(attr(big$oob_path, "criterion"), if (loss == "L2") "brier" else "logloss")
          expect_true(all(c("oob_loss", "oob_loss_raw", "se_loss", "se_loss_raw") %in% names(big$oob_path)))
        }
        p <- big$fit$oob_path()
        V <- ncol(p$sum0) / length(p$checkpoints)
        expect_equal(V, big$fit$get_shape()$value_size)
        c3 <- 6L
        set.seed(5)
        fk <- rpf_050(
          d$x,
          y,
          loss = loss,
          max_interaction = 3,
          n_families = 6,
          n_splits = p$checkpoints[c3],
          split_structure = s
        )
        cols <- (c3 - 1) * V + seq_len(V)
        expect_equal(
          unname(tot(fk)),
          unname(p$sum0[, cols, drop = FALSE] + p$sum1[, cols, drop = FALSE]),
          tolerance = 1e-12
        )
        # the forest is the fit without the path
        set.seed(5)
        plain <- rpf_050(d$x, y, loss = loss, max_interaction = 3, n_families = 6, n_splits = 60, split_structure = s)
        expect_identical(big$fit$get_model(), plain$fit$get_model())
      }
      set.seed(9)
      tu <- rpf_050(d$x, y, loss = loss, max_interaction = 3, n_families = 8, n_splits = NULL, nthreads = 2)
      expect_s4_class(tu$fit, "Rcpp_ClassificationRPF")
      expect_identical(tu$splits_oob$n_splits, tu$params$splits)
      expect_output(print(tu), "chosen out of bag")
      set.seed(9)
      fx <- rpf_050(d$x, y, loss = loss, max_interaction = 3, n_families = 8, n_splits = tu$params$splits)
      expect_identical(tu$fit$get_model(), fx$fit$get_model())
      expect_identical(predict(tu, d$x, type = "numeric"), predict(fx, d$x, type = "numeric"))
    }
  }
})

test_that("a pruned classification model predicts the same after rpf_marshal(), saveRDS() and readRDS()", {
  d <- cp_data(300)
  fit <- cp_fit(d, d$ym, "logit", n_families = 12)
  pd <- prune_components(fit, d$x, d$ym, max_panels = 2, seed = 1)
  p0 <- predict(pd, d$x)
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(rpf_marshal(pd), f)
  back <- rpf_unmarshal(readRDS(f))
  expect_identical(predict(back, d$x), p0)
  expect_identical(predict_components(back, d$x)$m, predict_components(pd, d$x)$m)
})

test_that("cross-validation and given rules work for classification", {
  d <- cp_data(250)
  for (y in list(d$yb, d$ym)) {
    fit <- cp_fit(d, y, "exponential", n_families = 6, n_splits = 40)
    cv <- prune_components(fit, d$x, y, select_by = "cv", folds = 3, max_panels = Inf, seed = 2)
    expect_true(all(c("cv_loss", "se_loss") %in% names(cv$cv)))
    expect_identical(attr(cv$cv, "criterion"), "logloss")
    expect_equal(unname(cp_recon(cv, d$x)), unname(as.matrix(predict(cv, d$x, type = "numeric"))), tolerance = 1e-12)
    dd <- prune_components(fit, d$x, y, rule = "D", tau = 0, max_panels = Inf)
    expect_equal(
      as.matrix(predict(dd, d$x, type = "numeric")),
      as.matrix(predict(fit, d$x, type = "numeric")),
      tolerance = 1e-12,
      ignore_attr = TRUE
    )
  }
  expect_error(prune_components(fit, d$x, factor(rep("a", 250), levels = "a")), "factor outcome")
  expect_error(prune_components(fit, d$x, as.numeric(d$ym)), "factor outcome")
})

test_that("rule R's likelihood fit: recovers beta, never above D's deviance, finite on separable data", {
  ns <- cp_ns()
  set.seed(4)
  n <- 3000
  # a multinomial logit (reference level 1) with shared slopes on two components of V = 2 columns
  F <- matrix(rnorm(n * 4), n, 4) # component-major: component 1 columns 1:2, component 2 columns 3:4
  eta <- cbind(0.3 + 1.5 * F[, 1] + 0.5 * F[, 3], -0.2 + 1.5 * F[, 2] + 0.5 * F[, 4])
  P <- cbind(1, exp(eta)) / (1 + rowSums(exp(eta)))
  cls <- apply(P, 1, function(p) sample.int(3, 1, prob = p))
  crit <- structure(
    list(type = "logloss", link = "softmax_ref", V = 2L, K = 3L, cls = cls, centre = "ref", n = n),
    class = "rpf_crit"
  )
  co <- ns$crit_fit_R(crit, F, 1:2, seq_len(n), start_alpha = c(0, 0))
  expect_equal(co$beta, c(1.5, 0.5), tolerance = 0.1)
  expect_equal(co$alpha, c(0.3, -0.2), tolerance = 0.15)
  dev <- function(a, b) sum(ns$crit_logloss(crit, ns$crit_predict_R(list(alpha = a, beta = b), F, 1:2, 2L), cls)$e)
  expect_lte(dev(co$alpha, co$beta), dev(c(0, 0), c(1, 1)))
  # separable binary data: finite coefficients
  x1 <- matrix(c(seq(-1, -0.1, length.out = 50), seq(0.1, 1, length.out = 50)), ncol = 1)
  cb <- structure(
    list(type = "logloss", link = "sigmoid", V = 1L, K = 2L, cls = rep(1:2, each = 50), centre = "none", n = 100L),
    class = "rpf_crit"
  )
  co2 <- ns$crit_fit_R(cb, x1, 1L, 1:100, start_alpha = 0)
  expect_true(all(is.finite(c(co2$alpha, co2$beta))))
})

test_that("rule R's likelihood fit with more coefficients than rows takes minimum-norm steps", {
  ns <- cp_ns()
  set.seed(11)
  # Every distinct row several times, in more than one class: the rows are not separable, so the likelihood has a
  # maximum, at the scores of the class shares; many coefficients reach them, and the fit is the start plus the
  # minimum-norm change that does
  m <- 8
  k <- 40
  F <- matrix(rnorm(m * k), m, k)[rep(seq_len(m), 3), ]
  n <- 3L * m
  cb <- structure(
    list(
      type = "logloss",
      link = "sigmoid",
      V = 1L,
      K = 2L,
      cls = rep(c(2L, 2L, 1L), each = m),
      centre = "none",
      n = n
    ),
    class = "rpf_crit"
  )
  co <- ns$crit_fit_R(cb, F, seq_len(k), seq_len(n), start_alpha = 0.2)
  Z <- cbind(1, F)
  start <- c(0.2, rep(1, k))
  expect_equal(as.numeric(Z %*% c(co$alpha, co$beta)), rep(log(2), n), tolerance = 1e-6)
  expect_equal(c(co$alpha, co$beta), start + ns$prune_lstsq(Z, log(2) - as.numeric(Z %*% start)), tolerance = 1e-6)
  # the fit on some of the rows, as the folds of the out-of-bag choice fit it
  idx <- c(1:6, m + 1:6, 2 * m + 1:6)
  cs <- ns$crit_fit_R(cb, F, seq_len(k), idx, start_alpha = 0.2)
  expect_equal(as.numeric(Z[idx, ] %*% c(cs$alpha, cs$beta)), rep(log(2), length(idx)), tolerance = 1e-6)

  # three classes (logit: two scores against the first level, one beta per component for both)
  m <- 6
  k <- 60
  F <- matrix(rnorm(m * 2 * k), m, 2 * k)[rep(seq_len(m), 4), ]
  n <- 4L * m
  cm <- structure(
    list(
      type = "logloss",
      link = "softmax_ref",
      V = 2L,
      K = 3L,
      cls = rep(c(1L, 2L, 2L, 3L), each = m),
      centre = "ref",
      n = n
    ),
    class = "rpf_crit"
  )
  co <- ns$crit_fit_R(cm, F, seq_len(k), seq_len(n), start_alpha = c(0.1, -0.1))
  Z <- rbind(cbind(1, 0, F[, 2 * seq_len(k) - 1]), cbind(0, 1, F[, 2 * seq_len(k)]))
  start <- c(0.1, -0.1, rep(1, k))
  target <- rep(c(log(2), 0), each = n)
  expect_equal(as.numeric(Z %*% c(co$alpha, co$beta)), target, tolerance = 1e-6)
  expect_equal(c(co$alpha, co$beta), start + ns$prune_lstsq(Z, target - as.numeric(Z %*% start)), tolerance = 1e-6)
})
