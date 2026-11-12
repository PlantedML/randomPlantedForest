# The logit and exponential losses: every region's step minimises its rows' loss plus that of a = 1 pseudo-row (a / K1
# per class) at the region's mean score, and a split is scored against the parent's own step (cpf.hpp, DriftFree).
# Their properties on one split, through the test hook rpf_classif_check(), and in small fits.

# classes 1..K1 coded as the losses code them: logit K1 - 1 columns 0 / 1 (class 1, the first level, is the reference);
# exponential one column -1 / 1 (binary) or K1 columns -1 / 1 (SAMME)
df_y <- function(cls, K1, loss) {
  if (loss == "logit") {
    outer(cls, 2:K1, "==") * 1
  } else if (K1 == 2) {
    matrix(ifelse(cls == 2, 1, -1))
  } else {
    outer(cls, 1:K1, "==") * 2 - 1
  }
}
df_cols <- function(K1, loss) if (loss == "logit" || K1 == 2) K1 - 1 else K1
df_check <- function(loss, Y, F, side, exact_step = TRUE) {
  rpf_classif_check(loss, as.matrix(Y), as.matrix(F), as.integer(side), exact_step, 1)
}
# the class probabilities of scores (rows): logit softmax(0, f), exponential binary plogis(f), SAMME softmax(f)
df_prob <- function(F, loss, K1) {
  F <- as.matrix(F)
  if (loss == "logit") {
    F <- cbind(0, F)
  } else if (K1 == 2) {
    F <- cbind(0, F)
  }
  e <- exp(F - apply(F, 1, max))
  e / rowSums(e)
}
df_cases <- list(c("logit", 2), c("logit", 3), c("exponential", 2), c("exponential", 3))

test_that("a pure region converges to a finite level and is then not moved again (drift-free)", {
  n <- 100
  for (case in df_cases) {
    loss <- case[1]
    K1 <- as.integer(case[2])
    # 100 rows of class 2 below the threshold, one row of class 1 above it
    cls <- c(rep(2L, n), 1L)
    Y <- df_y(cls, K1, loss)
    side <- c(rep(0L, n), 1L)
    # the level of the exact minimiser: probability (n + 1 / K1) / (n + 1) for the class
    target <- (n + 1 / K1) / (n + 1)
    for (exact in c(TRUE, FALSE)) {
      if (loss == "exponential" && !exact) {
        next
      }
      label <- sprintf("%s, %d classes, exact_step %s", loss, K1, exact)
      F <- matrix(0, n + 1, df_cols(K1, loss))
      steps <- numeric(0)
      for (r in 1:40) {
        ch <- df_check(loss, Y, F, side, exact)
        F[side == 0, ] <- sweep(F[side == 0, , drop = FALSE], 2, ch$step_s, "+")
        steps[r] <- max(abs(ch$step_s))
      }
      p <- df_prob(F[1, , drop = FALSE], loss, K1)[1, 2]
      expect_equal(p, target, tolerance = 1e-9, label = label)
      if (exact) {
        # the first step reaches the level, every later one is 0
        expect_gt(steps[1], 1, label = label)
        expect_lt(max(steps[-1]), 1e-9, label = label)
      } else {
        # Newton's steps shrink towards it and never overshoot
        expect_true(all(diff(steps) < 0 | steps[-1] < 1e-12), label = label)
        expect_lt(steps[40], 1e-9, label = label)
      }
    }
  }
})

test_that("splitting a pure leaf gains nothing at any level, nor does repeating a split", {
  set.seed(1)
  n <- 100
  for (case in df_cases) {
    loss <- case[1]
    K1 <- as.integer(case[2])
    for (exact in c(TRUE, FALSE)) {
      if (loss == "exponential" && !exact) {
        next
      }
      label <- sprintf("%s, %d classes, exact_step %s", loss, K1, exact)
      # a pure leaf split in half, at the scores 0, at its level and beyond it
      Y <- df_y(rep(2L, n), K1, loss)
      side <- rep(0:1, each = n / 2)
      level <- df_check(loss, Y, matrix(0, n, df_cols(K1, loss)), c(rep(0L, n - 1), 1L), exact)$step_s
      for (mult in c(0, 1, 3)) {
        F <- matrix(mult * level, n, length(level), byrow = TRUE)
        ch <- df_check(loss, Y, F, side, exact)
        expect_identical(ch$score, Inf, label = paste(label, "pure, level x", mult))
        expect_gt(ch$raw, 0, label = paste(label, "pure, level x", mult))
      }
      # a mixed region whose classes depend on x: the split at x < 0.5 gains; after its exact steps, the same split
      # again does not (after one Newton step it may: the repeat finishes a step that was too short)
      x <- runif(2 * n)
      cls <- 1L + (runif(2 * n) < ifelse(x < 0.5, 0.2, 0.8)) + if (K1 == 3) (runif(2 * n) < 0.3) else 0L
      Y <- df_y(cls, K1, loss)
      side <- as.integer(x >= 0.5)
      F <- matrix(0, 2 * n, df_cols(K1, loss))
      first <- df_check(loss, Y, F, side, exact)
      expect_lt(first$score, 0, label = paste(label, "first split"))
      F[side == 0, ] <- sweep(F[side == 0, , drop = FALSE], 2, first$step_s, "+")
      F[side == 1, ] <- sweep(F[side == 1, , drop = FALSE], 2, first$step_b, "+")
      again <- df_check(loss, Y, F, side, exact)
      if (exact) {
        expect_identical(again$score, Inf, label = paste(label, "repeat"))
        expect_gt(again$raw, 0, label = paste(label, "repeat"))
      } else {
        expect_lt(abs(again$raw), abs(first$raw) / 10, label = paste(label, "repeat"))
      }
    }
  }
})

test_that("binary logit: the K x K arithmetic is the binary closed form", {
  set.seed(2)
  n <- 60
  y <- rbinom(n, 1, 0.4)
  F <- matrix(rnorm(n), n)
  side <- rbinom(n, 1, 0.5)
  a <- 1
  region <- function(i) {
    p <- stats::plogis(F[i])
    m <- mean(F[i])
    G <- sum(y[i] - p) + a * (0.5 - stats::plogis(m))
    H <- sum(p * (1 - p)) + a / 4
    list(c = G / H, gain = G^2 / (2 * H), lps = a / 2 * (log1p(exp(-m)) + log1p(exp(m))))
  }
  s <- region(side == 0)
  b <- region(side == 1)
  p <- region(rep(TRUE, n))
  ch <- df_check("logit", matrix(y), F, side, exact_step = FALSE)
  expect_equal(ch$newton_s, s$c, tolerance = 1e-13)
  expect_equal(ch$newton_b, b$c, tolerance = 1e-13)
  expect_equal(ch$raw, -(s$gain + b$gain - p$gain) + s$lps + b$lps - p$lps, tolerance = 1e-12)
  # the exact step: the root of sum (y - plogis(F + c)) + (a / 2)(1 - 2 plogis(m + c)) = 0
  root <- function(i) {
    m <- mean(F[i])
    stats::uniroot(
      function(c) sum(y[i] - stats::plogis(F[i] + c)) + a / 2 * (1 - 2 * stats::plogis(m + c)),
      c(-20, 20),
      tol = 1e-14
    )$root
  }
  ch <- df_check("logit", matrix(y), F, side, exact_step = TRUE)
  expect_equal(ch$step_s, root(side == 0), tolerance = 1e-10)
  expect_equal(ch$step_b, root(side == 1), tolerance = 1e-10)
})

test_that("logit, 3 classes: the Newton step is the full K x K one; exact_step reaches the minimiser", {
  set.seed(3)
  n <- 80
  K1 <- 3
  cls <- sample.int(K1, n, TRUE, prob = c(0.5, 0.3, 0.2))
  Y <- df_y(cls, K1, "logit")
  F <- matrix(rnorm(2 * n), n)
  side <- rbinom(n, 1, 0.4)
  a <- 1
  i <- side == 0
  P <- df_prob(F[i, ], "logit", K1)[, -1]
  m <- colMeans(F[i, ])
  pm <- exp(c(0, m)) / sum(exp(c(0, m)))
  G <- colSums(Y[i, ] - P) + a * (1 / K1 - pm[-1])
  H <- matrix(0, 2, 2)
  for (r in which(i)) {
    p <- df_prob(F[r, , drop = FALSE], "logit", K1)[1, -1]
    H <- H + diag(p) - tcrossprod(p)
  }
  H <- H + a / 2 * (diag(2) - matrix(1 / K1, 2, 2))
  ch <- df_check("logit", Y, F, side, exact_step = FALSE)
  expect_equal(ch$newton_s, solve(H, G), tolerance = 1e-12)
  expect_equal(ch$step_s, ch$newton_s)
  # exact: the gradient of the penalised loss is 0 at the step
  ch <- df_check("logit", Y, F, side, exact_step = TRUE)
  c <- ch$step_s
  Pc <- df_prob(sweep(F[i, ], 2, c, "+"), "logit", K1)[, -1]
  pmc <- exp(c(0, m + c)) / sum(exp(c(0, m + c)))
  grad <- colSums(Pc - Y[i, ]) + a * (pmc[-1] - 1 / K1)
  expect_lt(max(abs(grad)), 1e-8)
  expect_gt(max(abs(c - ch$newton_s)), 1e-4)
})

test_that("one split from the scores 0: the leaves get the exact minimisers, or one Newton step", {
  # the leaves' steps, before the level (test-classif-level.R)
  rlang::local_options(randomPlantedForest.level = FALSE)
  set.seed(4)
  n <- 200
  x <- data.frame(x1 = runif(n))
  for (K1 in 2:3) {
    y <- factor(letters[pmin(K1, 1 + (x$x1 > 0.4) + (runif(n) < 0.3))], levels = letters[1:K1])
    for (loss in c("logit", "exponential")) {
      for (exact in c(TRUE, FALSE)) {
        set.seed(1)
        f <- rpf(
          x,
          y,
          loss = loss,
          exact_step = exact,
          max_interaction = 1,
          n_families = 1,
          n_splits = 1,
          n_try_split_points = 10,
          share_candidates = 1,
          deterministic = TRUE,
          prune = FALSE
        )
        p <- as.matrix(predict(f, x, type = "prob"))
        sides <- split(seq_len(n), round(p[, 1], 12))
        expect_length(sides, 2)
        for (s in sides) {
          cnt <- tabulate(as.integer(y[s]), K1)
          want <- if (exact || loss == "exponential") {
            (cnt + 1 / K1) / (length(s) + 1)
          } else {
            K <- K1 - 1
            H <- length(s) * (diag(K) / K1 - matrix(1 / K1^2, K, K)) + 0.5 * (diag(K) - matrix(1 / K1, K, K))
            step <- solve(H, cnt[-1] - length(s) / K1)
            exp(c(0, step)) / sum(exp(c(0, step)))
          }
          expect_equal(unname(p[s[1], ]), want, tolerance = 1e-12, label = paste(loss, K1, exact))
        }
      }
    }
  }
})

test_that("separable classes: one split each, then no candidate gains and the fit stays as it is", {
  # the leaves' steps, before the level (test-classif-level.R)
  rlang::local_options(randomPlantedForest.level = FALSE)
  n <- 120
  x <- data.frame(x1 = seq(0.005, 0.995, length.out = n))
  y <- factor(ifelse(x$x1 < 0.3 | x$x1 >= 0.7, "a", "b"))
  for (loss in c("logit", "exponential")) {
    # deterministic, every threshold tried
    fit <- function(splits) {
      rpf(
        x,
        y,
        loss = loss,
        max_interaction = 1,
        n_families = 1,
        n_splits = splits,
        n_try_split_points = n,
        share_candidates = 1,
        deterministic = TRUE,
        prune = FALSE
      )
    }
    a <- fit(20)
    b <- fit(200)
    expect_identical(a$fit$get_model(), b$fit$get_model(), label = loss)
    # every leaf boundary is a class boundary: no pure leaf was split
    m <- a$fit$get_model()[[1]]
    t1 <- which(vapply(m$variables, identical, logical(1), 1L))
    ends <- unlist(lapply(m$intervals[[t1]], function(iv) iv[, 1]))
    inner <- unique(ends[is.finite(ends) & ends > min(x$x1) & ends < max(x$x1)])
    expect_setequal(inner, x$x1[c(which(x$x1 >= 0.3)[1], which(x$x1 >= 0.7)[1])])
    # each region at its level, whatever the steps before: (n + 1 / 2) / (n + 1) for its class, n its rows
    p <- predict(a, x, type = "prob")$.pred_a
    for (r in list(x$x1 < 0.3, x$x1 >= 0.3 & x$x1 < 0.7, x$x1 >= 0.7)) {
      level <- (sum(r) + 0.5) / (sum(r) + 1)
      expect_equal(p[r], rep(if (y[r][1] == "a") level else 1 - level, sum(r)), tolerance = 1e-12, label = loss)
    }
  }
})

test_that("the logit forest does not depend on which class is the reference", {
  set.seed(5)
  n <- 300
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  eta <- cbind(0, 3 * x$x1 - 1.5, 2 * (x$x2 - 0.5) + x$x3)
  cls <- apply(eta, 1, function(e) sample.int(3, 1, prob = exp(e)))
  # deterministic (no bootstrap): a bootstrap sample's duplicated rows make candidates that tie exactly, and the codings
  # add their sums in other orders, so they can break such a tie differently
  for (exact in c(TRUE, FALSE)) {
    probs <- lapply(list(c("a", "b", "c"), c("b", "c", "a"), c("c", "a", "b")), function(lev) {
      y <- factor(c("a", "b", "c")[cls], levels = lev)
      f <- rpf(
        x,
        y,
        loss = "logit",
        exact_step = exact,
        max_interaction = 2,
        n_families = 1,
        n_splits = 40,
        deterministic = TRUE,
        prune = FALSE
      )
      as.matrix(predict(f, x, type = "prob"))[, c(".pred_a", ".pred_b", ".pred_c")]
    })
    expect_lt(max(abs(probs[[2]] - probs[[1]])), 1e-10, label = paste("exact_step", exact))
    expect_lt(max(abs(probs[[3]] - probs[[1]])), 1e-10, label = paste("exact_step", exact))
  }
})

test_that("the running sums choose the splits of the per-candidate pass, in every split structure", {
  # a refit through the C++ object (set_parameters()) gets no level: compared without it
  rlang::local_options(randomPlantedForest.level = FALSE)
  # The searches score a split point from the sums of the rows' terms along the leaf's sorted rows (or per bin, or by
  # one pass per try); set_parameters("scan", 0) (internal) scores every try by a pass over its leaf's rows instead.
  # The leaves' steps come from the winner's rows either way, so on data without ties (no bootstrap: deterministic;
  # continuous predictors; at least 5 rows per side) the fits are the same bit for bit.
  structures <- c("leaf", "leaf_binned", "tree_then_leaf", "tree_then_all_leaves", "target_tree_then_leaf")
  set.seed(7)
  n <- 250
  x <- data.frame(x1 = runif(n), x2 = rnorm(n), x3 = runif(n))
  eta <- cbind(0, 3 * x$x1 - 1.5, 2 * x$x2 * x$x3)
  cls <- apply(eta, 1, function(e) sample.int(3, 1, prob = exp(e)))
  for (K1 in 2:3) {
    y <- factor(letters[pmin(cls, K1)])
    for (loss in c("logit", "exponential")) {
      for (s in structures) {
        label <- paste(loss, K1, s)
        f <- rpf(
          x,
          y,
          loss = loss,
          split_structure = s,
          max_interaction = 2,
          n_families = 1,
          n_splits = 40,
          min_leaf_rows = 5,
          deterministic = TRUE,
          prune = FALSE
        )
        scan <- f$fit$get_model()
        expect_gt(sum(lengths(scan[[1]]$values)), 10, label = label)
        f$fit$set_parameters("scan", 0)
        expect_identical(f$fit$get_model(), scan, label = label)
      }
    }
  }
})

test_that("logit keeps its rows' terms or computes them: the same fit", {
  # a refit through the C++ object (set_parameters()) gets no level: compared without it
  rlang::local_options(randomPlantedForest.level = FALSE)
  set.seed(8)
  n <- 200
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  for (K1 in c(2, 3, 12)) {
    y <- factor(sample(letters[1:K1], n, TRUE))
    y[x$x1 < 0.3] <- "a"
    set.seed(9)
    f <- rpf(x, y, loss = "logit", max_interaction = 2, n_families = 2, n_splits = 30, prune = FALSE)
    default <- f$fit$get_model()
    for (keep in c(0, 1)) {
      set.seed(9)
      f$fit$set_parameters("logit_terms", keep)
      expect_identical(f$fit$get_model(), default, label = paste(K1, "classes, kept", keep))
    }
  }
})

test_that("logit and exponential fits do not depend on the number of threads", {
  set.seed(10)
  n <- 200
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  y3 <- factor(sample(c("a", "b", "c"), n, TRUE))
  y3[x$x1 < 0.4] <- "a"
  y2 <- factor(ifelse(x$x2 + rnorm(n, sd = 0.3) > 0.5, "u", "v"))
  for (y in list(y2, y3)) {
    for (loss in c("logit", "exponential")) {
      for (s in c("leaf", "tree_then_leaf")) {
        fit <- function(nt, exact = TRUE) {
          set.seed(11)
          rpf(
            x,
            y,
            loss = loss,
            split_structure = s,
            exact_step = exact,
            n_families = 4,
            n_splits = 30,
            nthreads = nt,
            prune = FALSE
          )$fit$get_model()
        }
        label <- paste(loss, nlevels(y), s)
        expect_identical(fit(2), fit(1), label = label)
        if (loss == "logit") {
          expect_identical(fit(2, FALSE), fit(1, FALSE), label = paste(label, "Newton"))
        }
      }
    }
  }
})

test_that("classification fits made before these losses load, print and predict", {
  # fixtures/classif_before_n18.rds: rpf_marshal() of three fits of the version before (N16, 11395bae), with delta
  # and epsilon in their params, and their predictions then (fixtures/make-classif-before-n18.R), made on macOS arm64:
  # elsewhere libm's exp and log may round the last digits of a prediction differently, so they are compared to 1e-12
  fx <- readRDS(test_path("fixtures", "classif_before_n18.rds"))
  expect_false(is.null(fx$exponential2$params$delta))
  for (nm in c("exponential2", "logit3", "pruned3")) {
    r <- rpf_unmarshal(fx[[nm]])
    expect_equal(predict(r, fx$new_data, type = "prob"), fx$pred[[nm]], tolerance = 1e-12, label = nm)
    out <- paste(utils::capture.output(print(r)), collapse = "\n")
    expect_no_match(out, "delta|epsilon")
  }
  r <- rpf_unmarshal(fx$logit3)
  # its link as saved then (K - 1 scores against the first level) is reported now with the first level's 0 joined and
  # centred across the classes
  old <- cbind(0, as.matrix(fx$pred$logit3_link))
  link <- predict(r, fx$new_data, type = "link")
  expect_named(link, paste0(".pred_", levels(r$blueprint$ptypes$outcomes[[1]])))
  expect_equal(unname(as.matrix(link)), unname(old - rowMeans(old)), tolerance = 1e-12)
  expect_match(paste(format(r), collapse = "\n"), "loss: logit")
  # with its data, it can be purified, and its components add up to its predictions
  purify(r)
  pc <- predict_components(r, fx$new_data)
  M <- as.matrix(pc$m)
  V <- length(pc$intercept)
  recon <- sapply(seq_len(V), function(v) pc$intercept[v] + rowSums(M[, seq(v, ncol(M), by = V), drop = FALSE]))
  expect_equal(unname(as.matrix(predict(r, fx$new_data, type = "link"))), unname(recon), tolerance = 1e-12)
})

test_that("tuning and pruning carry exact_step", {
  skip_on_cran()
  set.seed(12)
  n <- 200
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- factor(ifelse(x$x1 + 0.5 * x$x2 + rnorm(n, sd = 0.3) > 0.8, "a", "b"))
  set.seed(13)
  tuned <- rpf(x, y, loss = "logit", exact_step = FALSE, n_families = 10, tune = "quick")
  expect_s3_class(tuned, "rpf_pruned")
  expect_false(tuned$fit$params$exact_step)
  set.seed(13)
  pruned <- rpf(x, y, loss = "logit", exact_step = FALSE, n_families = 10)
  expect_false(pruned$fit$params$exact_step)
  r <- rpf_unmarshal(rpf_marshal(pruned))
  expect_false(r$fit$params$exact_step)
  expect_identical(predict(r, x, type = "prob"), predict(pruned, x, type = "prob"))
})
