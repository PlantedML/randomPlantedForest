# Classification with split_structure = "leaves": min_leaf_rows and the fast
# draw (options(randomPlantedForest.exact_draw = FALSE)), by regression's draw

with_exact_draw <- function(value, code) {
  old <- options(randomPlantedForest.exact_draw = value)
  on.exit(options(old), add = TRUE)
  force(code)
}

cl_sim <- function(n = 300, seed = 1, classes = 2) {
  set.seed(seed)
  x <- cbind(x1 = runif(n), x2 = runif(n), x3 = round(runif(n) * 6), x4 = rnorm(n), x5 = rbinom(n, 1, 0.3))
  eta <- 2 * (x[, 1] - 0.5) + sin(4 * x[, 2]) + 3 * (x[, 1] - 0.5) * (x[, 2] - 0.5) + 0.3 * x[, 3] + 1.5 * x[, 5] - 1.3
  y <- if (classes == 2) {
    rbinom(n, 1, stats::plogis(eta))
  } else {
    findInterval(eta + stats::rlogis(n), c(-0.5, 0.5))
  }
  list(x = x, y = factor(y))
}

# Per family, the rows of the bootstrap sample in every leaf of every tree
# except the constant one (leaves are [lower, upper) boxes), and the number of
# splits made: with max_interaction = 1, the sum over the trees of their
# leaves minus their rows / n (as in test-fast-draw.R)
cl_leaf_stats <- function(fit, x) {
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

test_that("loss L2 on a 0/1 response is the regression fit, bit for bit, with min_leaf_rows and the fast draw", {
  d <- cl_sim()
  x <- d$x
  y <- matrix(as.numeric(d$y == "1"), ncol = 1)
  # max_interaction, ntrees, splits, split_try, t_try, purify, deterministic,
  # nthreads, cv, split_decay_rate, max_candidates, delete_leaves,
  # split_structure ("leaves")
  base <- c(2, 3, 60, 5, 0.4, 0, 0, 1, 0, 0.1, 50, 1, 3)
  # [min_leaf_rows, exact_draw]
  for (opt in list(c(1, 1), c(5, 1), c(1, 0), c(5, 0), c(20, 0))) {
    for (det in c(0, 1)) {
      for (nt in c(1, 2)) {
        p <- base
        p[7] <- det
        p[8] <- nt
        set.seed(5)
        reg <- new(RandomPlantedForest, y, x, c(p, opt))
        set.seed(5)
        cls <- new(ClassificationRPF, y, x, "L2", c(p, 0, 0.1, opt))
        label <- sprintf("min_leaf_rows %d, exact %d, deterministic %d, %d threads", opt[1], opt[2], det, nt)
        expect_identical(cls$get_model(), reg$get_model(), label = label)
        expect_identical(cls$predict_matrix(x, 0, 1L), reg$predict_matrix(x, 0, 1L), label = label)
      }
    }
  }
})

test_that("the defaults are the fit of before; the options change random fits", {
  d <- cl_sim()
  fit <- function(...) {
    set.seed(9)
    rpf(d$x, d$y, ntrees = 3, splits = 40, ...)$fit$get_model()
  }
  default <- fit()
  expect_identical(fit(min_leaf_rows = 1L), default)
  expect_identical(with_exact_draw(TRUE, fit()), default)
  expect_false(identical(fit(min_leaf_rows = 10L), default))
  expect_false(identical(with_exact_draw(FALSE, fit()), default))
  # the fast draw is reproducible and the same on any number of threads
  expect_identical(with_exact_draw(FALSE, fit(nthreads = 2L)), with_exact_draw(FALSE, fit()))
})

test_that("classification leaves keep min_leaf_rows rows, and no infeasible entry is drawn, by both draws", {
  for (classes in c(2, 3)) {
    d <- cl_sim(classes = classes)
    for (loss in c("L2", "logit", "exponential")) {
      # One candidate per round, and always one with a valid threshold (the
      # constant leaf), so every round splits unless an entry without one is
      # drawn; overwritten leaves (delete_leaves) make such entries.
      fit <- function(exact, nthreads = 1L) {
        set.seed(3)
        with_exact_draw(
          exact,
          rpf(
            d$x,
            d$y,
            loss = loss,
            ntrees = 3,
            splits = 100,
            max_interaction = 1,
            max_candidates = 1,
            min_leaf_rows = 5,
            nthreads = nthreads
          )
        )
      }
      for (exact in c(TRUE, FALSE)) {
        label <- sprintf("%d classes, %s, exact %s", classes, loss, exact)
        for (st in cl_leaf_stats(fit(exact), d$x)) {
          expect_gte(min(st$rows), 5, label = label)
          expect_equal(st$splits, 100, label = label)
        }
      }
      expect_identical(fit(FALSE)$fit$get_model(), fit(FALSE, nthreads = 2L)$fit$get_model())
    }
  }

  # the default growth (order 2, more candidates): no leaf below 5 rows either
  d <- cl_sim()
  for (exact in c(TRUE, FALSE)) {
    set.seed(4)
    f <- with_exact_draw(exact, rpf(d$x, d$y, ntrees = 3, splits = 60, min_leaf_rows = 5))
    expect_gte(min(unlist(lapply(cl_leaf_stats(f, d$x), `[[`, "rows"))), 5)
  }
  # deterministic: all rows once
  f <- rpf(d$x, d$y, ntrees = 1, splits = 60, deterministic = TRUE, min_leaf_rows = 15)
  expect_gte(min(unlist(lapply(cl_leaf_stats(f, d$x), `[[`, "rows"))), 15)
})

test_that("classification min_leaf_rows is validated, stored, printed and restored", {
  d <- cl_sim()
  expect_error(rpf(d$x, d$y, ntrees = 1, min_leaf_rows = 5, split_structure = "hist"), "leaves")
  set.seed(2)
  f <- with_exact_draw(FALSE, rpf(d$x, d$y, ntrees = 2, splits = 20, min_leaf_rows = 10))
  expect_identical(f$params$min_leaf_rows, 10L)
  expect_false(f$params$exact_draw)
  out <- paste(format(f), collapse = "\n")
  expect_match(out, "min_leaf_rows: 10")
  expect_match(out, "fast draw")
  r <- rpf_unmarshal(rpf_marshal(f))
  expect_identical(predict(r, d$x, type = "prob"), predict(f, d$x, type = "prob"))
  # [..., split_mode, delta, epsilon, min_leaf_rows, exact_draw]; models saved
  # before have neither: 1 and the exact draw
  pars <- rpf_param_vector(f$params, "classification")
  expect_length(pars, 17)
  expect_equal(pars[16:17], c(10, 0))
  params <- f$params
  params$min_leaf_rows <- NULL
  params$exact_draw <- NULL
  expect_equal(rpf_param_vector(params, "classification")[16:17], c(1, 1))
  expect_length(new(ClassificationRPF, "L2", pars[1:15])$get_model(), 0)
})
