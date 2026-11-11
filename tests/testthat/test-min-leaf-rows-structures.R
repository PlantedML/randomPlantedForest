# min_leaf_rows with every split structure, in regression and classification (2 and 3 classes), random (both draws of
# "leaf") and deterministic. test-min-leaf-rows.R and test-classification-leaves.R hold the first tests of "leaf".

mlrs_structures <- c("leaf", "leaf_binned", "tree_then_leaf", "tree_then_all_leaves", "target_tree_then_leaf")

mlrs_sim <- function(outcome = c("regression", "binary", "3 classes"), n = 300, seed = 1) {
  outcome <- match.arg(outcome)
  set.seed(seed)
  x <- cbind(x1 = runif(n), x2 = runif(n), x3 = round(runif(n) * 6), x4 = rnorm(n), x5 = rbinom(n, 1, 0.3))
  eta <- 2 * (x[, 1] - 0.5) + sin(4 * x[, 2]) + 3 * (x[, 1] - 0.5) * (x[, 2] - 0.5) + 0.3 * x[, 3] + 1.5 * x[, 5] - 1.3
  y <- switch(
    outcome,
    "regression" = eta + rnorm(n, sd = 0.5),
    "binary" = factor(rbinom(n, 1, stats::plogis(eta))),
    "3 classes" = factor(findInterval(eta + stats::rlogis(n), c(-0.5, 0.5)))
  )
  list(x = x, y = y)
}

mlrs_fit <- function(d, structure, ..., exact = FALSE, seed = 3) {
  old <- options(randomPlantedForest.exact_draw = exact)
  on.exit(options(old), add = TRUE)
  set.seed(seed)
  rpf_050(d$x, d$y, split_structure = structure, ...)
}

# Per family: the rows of the bootstrap sample (with their multiplicity) in every leaf of every tree except the
# constant one (leaves are [lower, upper) boxes), and, for max_interaction = 1, the number of splits made: the sum
# over the trees of their leaves minus their rows / n (as in test-classification-leaves.R)
mlrs_leaf_stats <- function(fit, x) {
  inbag <- fit$fit$inbag()
  model <- fit$fit$get_model()
  lapply(seq_along(model), function(f) {
    family <- model[[f]]
    trees <- which(!vapply(family$variables, identical, logical(1), 0L))
    rows <- lapply(trees, function(t) {
      vapply(
        family$intervals[[t]],
        function(iv) sum(inbag[colSums(t(x) >= iv[1, ] & t(x) < iv[2, ]) == ncol(x), f]),
        1
      )
    })
    list(rows = unlist(rows), splits = sum(lengths(rows)) - sum(unlist(rows)) / nrow(x))
  })
}
mlrs_min_rows <- function(fit, x) min(unlist(lapply(mlrs_leaf_stats(fit, x), `[[`, "rows")))

test_that("min_leaf_rows = 1 is the default fit of every structure, bit for bit", {
  for (outcome in c("regression", "binary", "3 classes")) {
    d <- mlrs_sim(outcome)
    for (s in mlrs_structures) {
      label <- paste(outcome, s)
      a <- mlrs_fit(d, s, n_families = 2, n_splits = 30)
      b <- mlrs_fit(d, s, n_families = 2, n_splits = 30, min_leaf_rows = 1L)
      expect_identical(b$fit$get_model(), a$fit$get_model(), label = label)
      a <- mlrs_fit(d, s, n_families = 1, n_splits = 30, deterministic = TRUE)
      b <- mlrs_fit(d, s, n_families = 1, n_splits = 30, deterministic = TRUE, min_leaf_rows = 1L)
      expect_identical(b$fit$get_model(), a$fit$get_model(), label = label)
    }
  }
})

test_that("every leaf keeps min_leaf_rows rows of the bootstrap sample, in every structure and mode", {
  for (outcome in c("regression", "binary", "3 classes")) {
    d <- mlrs_sim(outcome)
    for (s in mlrs_structures) {
      # the limit binds: without it some leaf has fewer rows
      a <- mlrs_fit(d, s, n_families = 3, n_splits = 60)
      expect_lt(mlrs_min_rows(a, d$x), 15, label = paste(outcome, s, "without the limit"))
      for (m in c(5L, 15L)) {
        for (exact in if (s == "leaf") c(TRUE, FALSE) else FALSE) {
          label <- sprintf("%s, %s, m = %d, exact draw %s", outcome, s, m, exact)
          b <- mlrs_fit(d, s, n_families = 3, n_splits = 60, min_leaf_rows = m, exact = exact)
          expect_gte(mlrs_min_rows(b, d$x), m, label = label)
        }
        # deterministic: all rows once
        b <- mlrs_fit(d, s, n_families = 1, n_splits = 60, deterministic = TRUE, min_leaf_rows = m)
        expect_gte(mlrs_min_rows(b, d$x), m, label = paste(outcome, s, m, "deterministic"))
      }
    }
  }
})

test_that("leaves without a threshold that keeps min_leaf_rows rows are never drawn", {
  # Order 1, one candidate and one try per round: every round has a leaf with a valid threshold (the constant leaf,
  # with all rows, stays), so every round splits unless a leaf without one is drawn. delete_leaves (the default) makes
  # small leaves by overwriting the split leaf. Classification with L2, which splits whenever it can: logit and
  # exponential make no split in a round whose candidate gains nothing, and the draw is the same for every loss.
  for (outcome in c("regression", "binary", "3 classes")) {
    d <- mlrs_sim(outcome)
    for (s in mlrs_structures) {
      for (exact in if (s == "leaf") c(TRUE, FALSE) else FALSE) {
        label <- sprintf("%s, %s, exact draw %s", outcome, s, exact)
        f <- mlrs_fit(
          d,
          s,
          n_families = 3,
          n_splits = 80,
          max_interaction = 1,
          max_candidates = 1,
          n_try_split_points = 1,
          min_leaf_rows = 20,
          exact = exact,
          loss = if (outcome == "regression") NULL else "L2"
        )
        for (st in mlrs_leaf_stats(f, d$x)) {
          expect_gte(min(st$rows), 20, label = label)
          expect_equal(st$splits, 80, label = label)
        }
      }
    }
  }
})

test_that("a binary predictor is split when both sides keep min_leaf_rows rows, in every structure", {
  set.seed(2)
  xb <- cbind(b = rbinom(200, 1, 0.5))
  yb <- 2 * xb[, 1] + rnorm(200, sd = 0.1)
  for (s in mlrs_structures) {
    db <- list(x = xb, y = yb)
    fb <- mlrs_fit(db, s, n_families = 3, n_splits = 5, min_leaf_rows = 20)
    pb <- predict(fb, cbind(b = c(0, 1)))$.pred
    expect_gt(pb[2] - pb[1], 1.5, label = s)
    # classification: the class of b
    dc <- list(x = xb, y = factor(xb[, 1]))
    fc <- mlrs_fit(dc, s, n_families = 3, n_splits = 5, min_leaf_rows = 20)
    pc <- predict(fc, cbind(b = c(0, 1)), type = "prob")$.pred_1
    expect_gt(pc[2] - pc[1], 0.5, label = s)
    # with x5 (30% ones) next to the other predictors, deterministic
    for (outcome in c("regression", "binary")) {
      d <- mlrs_sim(outcome)
      f <- mlrs_fit(d, s, n_families = 1, n_splits = 60, deterministic = TRUE, min_leaf_rows = 10)
      expect_true(any(vapply(f$fit$get_model()[[1]]$variables, function(v) 5L %in% v, logical(1))), label = s)
    }
  }
})

test_that("min_leaf_rows fits do not depend on the number of threads", {
  for (outcome in c("regression", "binary", "3 classes")) {
    d <- mlrs_sim(outcome)
    for (s in mlrs_structures) {
      for (exact in if (s == "leaf") c(TRUE, FALSE) else FALSE) {
        label <- sprintf("%s, %s, exact draw %s", outcome, s, exact)
        a <- mlrs_fit(d, s, n_families = 4, n_splits = 40, min_leaf_rows = 10, exact = exact)
        b <- mlrs_fit(d, s, n_families = 4, n_splits = 40, min_leaf_rows = 10, exact = exact, nthreads = 2)
        expect_identical(b$fit$get_model(), a$fit$get_model(), label = label)
        expect_false(identical(
          mlrs_fit(d, s, n_families = 4, n_splits = 40, exact = exact)$fit$get_model(),
          a$fit$get_model()
        ))
      }
    }
  }
})

test_that("loss L2 on a 0/1 response is the regression fit with min_leaf_rows, in every structure", {
  d <- mlrs_sim("binary")
  x <- d$x
  y <- matrix(as.numeric(d$y == "1"), ncol = 1)
  # max_interaction, n_families, splits, n_try_split_points, share_candidates, purify, deterministic,
  # nthreads, cv, split_decay_rate, max_candidates, delete_leaves, split_structure (0 = target_tree_then_leaf,
  # 1 = tree_then_all_leaves, 2 = tree_then_leaf, 3 = leaf, 4 = leaf_binned), then [min_leaf_rows, exact_draw]
  base <- c(2, 3, 50, 5, 0.4, 0, 0, 1, 0, 0.1, 50, 1, 3)
  for (mode in 0:4) {
    for (det in c(0, 1)) {
      for (m in c(5, 20)) {
        p <- base
        p[7] <- det
        p[13] <- mode
        opt <- c(m, 1)
        set.seed(5)
        reg <- new(RandomPlantedForest, y, x, c(p, opt))
        set.seed(5)
        cls <- new(ClassificationRPF, y, x, "L2", c(p, opt))
        label <- sprintf("structure %d, min_leaf_rows %d, deterministic %d", mode, m, det)
        expect_identical(cls$get_model(), reg$get_model(), label = label)
        expect_identical(cls$predict_matrix(x, 0, 1L), reg$predict_matrix(x, 0, 1L), label = label)
      }
    }
  }
})

test_that("tune searches min_leaf_rows with every structure", {
  skip_on_cran()
  d <- mlrs_sim("regression", n = 200)
  set.seed(1)
  f <- rpf(
    d$x,
    d$y,
    split_structure = "tree_then_leaf",
    n_families = 10,
    max_interaction = 2,
    tune = tune_control(search = 5, min_leaf_rows = c(1, 20), neighbours = FALSE)
  )
  expect_setequal(f$tune$search$min_leaf_rows, c(1L, 20L))
})
