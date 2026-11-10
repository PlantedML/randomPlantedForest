mlr_sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- cbind(x1 = runif(n), x2 = runif(n), x3 = round(runif(n) * 6), x4 = rnorm(n), x5 = rbinom(n, 1, 0.3))
  y <- 2 * x[, 1] + sin(4 * x[, 2]) + 3 * x[, 1] * x[, 2] + 0.3 * x[, 3] + 1.5 * x[, 5] + rnorm(n, sd = 0.5)
  list(x = x, y = y)
}

mlr_fit <- function(d, ..., n_splits = 60, n_families = 5) {
  set.seed(5)
  rpf_050(d$x, d$y, n_families = n_families, n_splits = n_splits, ...)
}

# Rows of x in every leaf of every tree except the constant one (leaves are [lower, upper) boxes)
leaf_rows <- function(fit, x) {
  unlist(lapply(fit$fit$get_model(), function(family) {
    lapply(seq_along(family$variables), function(t) {
      if (identical(family$variables[[t]], 0L)) {
        return(NULL)
      }
      vapply(family$intervals[[t]], function(iv) sum(colSums(t(x) >= iv[1, ] & t(x) < iv[2, ]) == ncol(x)), 1)
    })
  }))
}

# Leaves per family; with delete_leaves = FALSE every split adds two
n_leaves <- function(fit) {
  vapply(fit$fit$get_model(), function(family) sum(lengths(family$values)), 1)
}

test_that("`min_leaf_rows = 1` is the default fit, bit for bit", {
  d <- mlr_sim()
  a <- mlr_fit(d)
  b <- mlr_fit(d, min_leaf_rows = 1L)
  expect_identical(b$fit$get_model(), a$fit$get_model())
  expect_identical(predict(b, d$x), predict(a, d$x))
  a <- mlr_fit(d, n_families = 1, deterministic = TRUE)
  b <- mlr_fit(d, n_families = 1, deterministic = TRUE, min_leaf_rows = 1)
  expect_identical(b$fit$get_model(), a$fit$get_model())
})

test_that("every leaf keeps at least `min_leaf_rows` rows", {
  d <- mlr_sim()
  # deterministic fits use all rows once, so the rows in a leaf are its rows of x
  a <- mlr_fit(d, n_families = 1, deterministic = TRUE)
  expect_lt(min(leaf_rows(a, d$x)), 15) # the limit binds
  for (m in c(5L, 15L, 40L)) {
    b <- mlr_fit(d, n_families = 1, deterministic = TRUE, min_leaf_rows = m)
    expect_gte(min(leaf_rows(b, d$x)), m)
  }
})

test_that("a binary feature is split when both sides keep `min_leaf_rows` rows", {
  d <- mlr_sim()
  b <- mlr_fit(d, n_families = 1, deterministic = TRUE, min_leaf_rows = 10)
  expect_true(any(vapply(b$fit$get_model()[[1]]$variables, function(v) 5L %in% v, logical(1))))
  # a design with one binary feature only: the forest learns it
  set.seed(2)
  xb <- cbind(b = rbinom(200, 1, 0.5))
  yb <- 2 * xb[, 1] + rnorm(200, sd = 0.1)
  set.seed(1)
  fb <- rpf_050(xb, yb, n_families = 3, n_splits = 5, min_leaf_rows = 20)
  pb <- predict(fb, cbind(b = c(0, 1)))$.pred
  expect_gt(pb[2] - pb[1], 1.5)
})

test_that("leaves without a valid threshold are not drawn; too large a value gives no split", {
  d <- mlr_sim()
  # one candidate per round, always one with a valid threshold: every round splits
  b <- mlr_fit(d, min_leaf_rows = 40, max_candidates = 1, delete_leaves = FALSE, n_splits = 15, n_families = 3)
  expect_equal(n_leaves(b), rep(1 + 2 * 15, 3))
  big <- mlr_fit(d, min_leaf_rows = 200, n_families = 2)
  expect_equal(n_leaves(big), c(1, 1))
  p <- predict(big, d$x)$.pred
  expect_lt(diff(range(p)), 1e-12)
})

test_that("`min_leaf_rows` is validated, stored and restored", {
  d <- mlr_sim()
  expect_error(mlr_fit(d, min_leaf_rows = 0), "min_leaf_rows")
  expect_error(mlr_fit(d, min_leaf_rows = 2.5), "min_leaf_rows")
  expect_error(mlr_fit(d, min_leaf_rows = 5, split_structure = "leaf_binned"), "leaf")
  # classification: test-classification-leaves.R

  b <- mlr_fit(d, min_leaf_rows = 10, n_families = 2)
  expect_identical(b$params$min_leaf_rows, 10L)
  expect_match(format(b), "min_leaf_rows: 10", all = FALSE)
  expect_false(any(grepl("min_leaf_rows", format(mlr_fit(d, n_families = 1)))))
  r <- rpf_unmarshal(rpf_marshal(b))
  expect_identical(predict(r, d$x), predict(b, d$x))
  # params saved before `min_leaf_rows` existed mean 1
  params <- b$params
  params$min_leaf_rows <- NULL
  expect_equal(rpf_param_vector(params, "regression")[[14]], 1)
})
