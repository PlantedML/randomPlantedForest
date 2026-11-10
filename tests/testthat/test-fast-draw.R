# options(randomPlantedForest.exact_draw): the exact and the fast draw of
# candidates for split_structure = "leaves"

with_exact_draw <- function(value, code) {
  old <- options(randomPlantedForest.exact_draw = value)
  on.exit(options(old), add = TRUE)
  force(code)
}

# How often each pool entry is drawn (internal test function); rows_ok: empty,
# or 0 / 1 per entry as min_leaf_rows sets it
draw_counts <- function(ages, rate, k, reps, exact, rows_ok = integer()) {
  rpf_leaves_draw_counts(ages, rate, as.integer(k), as.integer(reps), exact, as.integer(rows_ok))
}

# Inclusion probabilities of drawing k of the weights w without replacement,
# one at a time in proportion to the weights of those not yet drawn
# (successive sampling), by enumerating the ordered draws.
inclusion_probs <- function(w, k) {
  probs <- numeric(length(w))
  walk <- function(left, p, depth) {
    if (depth == k) {
      return()
    }
    for (i in which(left)) {
      q <- p * w[i] / sum(w[left])
      probs[i] <<- probs[i] + q
      left_i <- left
      left_i[i] <- FALSE
      walk(left_i, q, depth + 1)
    }
  }
  walk(rep(TRUE, length(w)), 1, 0)
  probs
}

test_that("both draws pick each entry with probability exp(-rate * age) / sum", {
  ages <- c(0, 1, 2, 5, 10, 0, 3)
  rate <- 0.3
  w <- exp(-rate * ages)
  reps <- 20000
  for (exact in c(TRUE, FALSE)) {
    set.seed(1)
    counts <- draw_counts(ages, rate, 1L, reps, exact)
    expect_equal(sum(counts), reps)
    # chi-squared test of the counts against the weights
    p_value <- suppressWarnings(stats::chisq.test(counts, p = w / sum(w))$p.value)
    expect_gt(p_value, 1e-3)
  }
})

test_that("both draws have the same inclusion probabilities for several candidates", {
  ages <- c(0, 1, 2, 4, 8, 0)
  rate <- 0.4
  k <- 3L
  reps <- 20000
  pi_i <- inclusion_probs(exp(-rate * ages), k)
  expect_equal(sum(pi_i), k)
  for (exact in c(TRUE, FALSE)) {
    set.seed(2)
    counts <- draw_counts(ages, rate, k, reps, exact)
    expect_equal(sum(counts), k * reps)
    # each count is binomial(reps, pi_i): within 4.5 standard deviations
    z <- (counts - reps * pi_i) / sqrt(reps * pi_i * (1 - pi_i))
    expect_true(all(abs(z) < 4.5), label = paste("exact =", exact, "z =", toString(round(z, 2))))
  }
})

test_that("the fast draw is without replacement and skips zero weights", {
  ages <- c(0, 3, 1, 2)
  set.seed(3)
  # the whole pool: every entry once per draw
  expect_equal(draw_counts(ages, 0.5, 4L, 100L, FALSE), rep(100, 4))
  # weights that underflow to 0 are not drawn while others are positive
  counts <- draw_counts(c(0, 1e4, 1, 1e4), 0.5, 2L, 100L, FALSE)
  expect_equal(counts, c(100, 0, 100, 0))
  # no positive weight: uniformly, as the exact draw
  counts <- draw_counts(c(1e4, 1e4, 1e4), 0.5, 1L, 300L, FALSE)
  expect_equal(sum(counts), 300)
  expect_true(all(counts > 0))
})

test_that("the default is the fast draw and is stored with the fit", {
  set.seed(4)
  fit_default <- rpf_050(mpg ~ ., data = mtcars, ntrees = 3, splits = 20)
  set.seed(4)
  fit_fast <- with_exact_draw(FALSE, rpf_050(mpg ~ ., data = mtcars, ntrees = 3, splits = 20))
  expect_false(fit_default$params$exact_draw)
  expect_identical(fit_default$fit$get_model(), fit_fast$fit$get_model())
})

test_that("the fast draw gives other, reproducible fits, the same on any number of threads", {
  fit <- function(seed, nthreads = 1L) {
    set.seed(seed)
    with_exact_draw(FALSE, rpf_050(mpg ~ ., data = mtcars, ntrees = 3, splits = 20, nthreads = nthreads))
  }
  fast <- fit(5)
  expect_false(fast$params$exact_draw)
  expect_identical(fast$fit$get_model(), fit(5)$fit$get_model())
  expect_identical(fast$fit$get_model(), fit(5, nthreads = 2L)$fit$get_model())

  set.seed(5)
  exact <- with_exact_draw(TRUE, rpf_050(mpg ~ ., data = mtcars, ntrees = 3, splits = 20))
  expect_false(identical(fast$fit$get_model(), exact$fit$get_model()))
  # an equally good fit: in-sample error of the same order
  mse <- function(f) mean((predict(f, mtcars)$.pred - mtcars$mpg)^2)
  expect_lt(mse(fast), 2 * mse(exact))
})

test_that("the option does not change other structures or deterministic fits", {
  same_fit <- function(...) {
    set.seed(6)
    a <- with_exact_draw(TRUE, rpf_050(...))
    set.seed(6)
    b <- with_exact_draw(FALSE, rpf_050(...))
    expect_identical(a$fit$get_model(), b$fit$get_model())
  }
  for (structure in c("hist", "cur_trees_1", "cur_trees_2", "res_trees")) {
    same_fit(mpg ~ ., data = mtcars, ntrees = 2, splits = 15, split_structure = structure)
  }
  same_fit(mpg ~ ., data = mtcars, ntrees = 1, splits = 15, deterministic = TRUE)
  cls <- data.frame(y = factor(mtcars$am), mtcars[, c("mpg", "wt", "hp")])
  same_fit(y ~ ., data = cls, ntrees = 2, splits = 15, split_structure = "hist")
  same_fit(y ~ ., data = cls, ntrees = 1, splits = 15, deterministic = TRUE)
})

test_that("the option is checked", {
  with_exact_draw("yes", expect_error(rpf_050(mpg ~ ., data = mtcars, ntrees = 1), "exact_draw"))
  with_exact_draw(NA, expect_error(rpf_050(mpg ~ ., data = mtcars, ntrees = 1), "exact_draw"))
})

test_that("models stored without min_leaf_rows and exact_draw rebuild with 1 and the exact draw", {
  params <- list(
    max_interaction = 2,
    ntrees = 1,
    splits = 5,
    split_try = 10,
    t_try = 0.4,
    purify = FALSE,
    deterministic = FALSE,
    nthreads = 1,
    split_decay_rate = 0.1,
    max_candidates = 50,
    delete_leaves = TRUE,
    split_structure = "leaves"
  )
  # regression: [..., split_mode, min_leaf_rows, exact_draw]
  pars <- rpf_param_vector(params, "regression")
  expect_length(pars, 15)
  expect_equal(pars[14:15], c(1, 1))
  expect_equal(rpf_param_vector(c(params, exact_draw = FALSE), "regression")[14:15], c(1, 0))
  expect_equal(rpf_param_vector(c(params, min_leaf_rows = 5L, exact_draw = FALSE), "regression")[14:15], c(5, 0))
  # classification: [..., split_mode, delta, epsilon, min_leaf_rows, exact_draw]
  pars <- rpf_param_vector(c(params, delta = 0, epsilon = 0.1), "classification")
  expect_length(pars, 17)
  expect_equal(pars[14:17], c(0, 0.1, 1, 1))
})

test_that("print() mentions the exact draw (not the default)", {
  set.seed(7)
  exact <- with_exact_draw(TRUE, rpf_050(mpg ~ ., data = mtcars, ntrees = 1, splits = 5))
  expect_match(paste(format(exact), collapse = "\n"), "exact draw")
  set.seed(7)
  fast <- rpf_050(mpg ~ ., data = mtcars, ntrees = 1, splits = 5)
  expect_no_match(paste(format(fast), collapse = "\n"), "draw")
})

test_that("both draws leave out the entries that min_leaf_rows rules out", {
  ages <- c(0, 1, 2, 5, 10, 0, 3)
  ok <- c(1L, 0L, 1L, 1L, 0L, 1L, 1L)
  rate <- 0.3
  w <- exp(-rate * ages[ok == 1])
  for (exact in c(TRUE, FALSE)) {
    set.seed(8)
    counts <- draw_counts(ages, rate, 1L, 20000L, exact, ok)
    expect_equal(counts[ok == 0], c(0, 0))
    p_value <- suppressWarnings(stats::chisq.test(counts[ok == 1], p = w / sum(w))$p.value)
    expect_gt(p_value, 1e-3)
    # as many candidates as entries left: each of them once per draw
    expect_equal(draw_counts(ages, rate, 5L, 100L, exact, ok), 100 * ok)
    # none left: uniformly over the whole pool, as without the check
    none <- draw_counts(ages, rate, 1L, 700L, exact, rep(0L, 7))
    expect_equal(sum(none), 700)
    expect_true(all(none > 0))
  }
})

# Per family, the number of splits made. With max_interaction = 1 a split
# either adds two leaves to tree {j} that cover the whole space (a split of the
# constant leaf) or splits a leaf of tree {j} in two (one leaf more, the same
# rows), so the splits are the sum over the trees of their leaves minus their
# rows of x (leaves are [lower, upper) boxes) / n.
fast_splits <- function(fit, x) {
  vapply(
    fit$fit$get_model(),
    function(family) {
      trees <- which(!vapply(family$variables, identical, logical(1), 0L))
      rows <- lapply(trees, function(t) {
        vapply(family$intervals[[t]], function(iv) sum(colSums(t(x) >= iv[1, ] & t(x) < iv[2, ]) == ncol(x)), 1)
      })
      sum(lengths(rows)) - sum(unlist(rows)) / nrow(x)
    },
    1
  )
}

test_that("both draws split every round under min_leaf_rows", {
  set.seed(1)
  n <- 300
  x <- cbind(x1 = runif(n), x2 = runif(n), x3 = round(runif(n) * 6), x4 = rnorm(n), x5 = rbinom(n, 1, 0.3))
  y <- 2 * x[, 1] + sin(4 * x[, 2]) + 3 * x[, 1] * x[, 2] + 0.3 * x[, 3] + 1.5 * x[, 5] + rnorm(n, sd = 0.5)
  # One candidate per round, and always one with a valid threshold (the
  # constant leaf), so every round splits unless an entry without one is drawn.
  # Leaves below 2 * 5 rows make such entries, among them the entries of a
  # split leaf overwritten by its upper half, whose weights the fast draw must
  # set again. (Without that, 118 to 127 of the 150 splits per family; without
  # leaving out entries at all, 99 to 114.)
  fit <- function(exact, nthreads = 1L) {
    set.seed(3)
    with_exact_draw(
      exact,
      rpf_050(
        x,
        y,
        ntrees = 3,
        splits = 150,
        max_interaction = 1,
        max_candidates = 1,
        min_leaf_rows = 5,
        nthreads = nthreads
      )
    )
  }
  for (exact in c(TRUE, FALSE)) {
    expect_equal(fast_splits(fit(exact), x), rep(150, 3))
  }
  expect_identical(fit(FALSE)$fit$get_model(), fit(FALSE, nthreads = 2L)$fit$get_model())

  # N2's test with the fast draw: with delete_leaves = FALSE every split adds
  # two leaves, so every round splits
  set.seed(5)
  f <- with_exact_draw(
    FALSE,
    rpf_050(x, y, ntrees = 3, splits = 15, min_leaf_rows = 40, max_candidates = 1, delete_leaves = FALSE)
  )
  expect_equal(vapply(f$fit$get_model(), function(family) sum(lengths(family$values)), 1), rep(1 + 2 * 15, 3))
})
