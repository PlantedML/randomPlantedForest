# rpf()'s defaults: a plain call is the explicit call it stands for (splits chosen out of bag, pruned with the display,
# 3 / 0.1 for "leaves", order 3, 200 families, the fast draw; the loss is master's, exponential for classification),
# bit for bit; with deterministic = TRUE the NULL defaults resolve to 30 splits and no pruning.

def_data <- function(n = 120) {
  set.seed(31)
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  f <- sin(2 * pi * x$x1) + 2 * (x$x2 - 0.5) * (x$x3 - 0.5)
  x$y <- f + rnorm(n, sd = 0.3)
  x$cl <- factor(ifelse(f + rnorm(n, sd = 0.5) > 0, "b", "a"))
  x$cl3 <- factor(c("a", "b", "c")[1 + (f > -0.3) + (f > 0.3)])
  x
}

def_parts <- function(m, d, type = NULL) {
  list(
    class = class(m),
    model = m$fit$fit$get_model(),
    params = m$fit$params,
    oob_path = m$fit$oob_path,
    rule = m$rule,
    tau = m$tau,
    beta = m$beta,
    intercept = m$intercept,
    pred = if (is.null(type)) predict(m, d) else predict(m, d, type = type)
  )
}

test_that("the default call is the explicit call it stands for", {
  skip_on_cran()
  d <- def_data(70)
  explicit <- function(f, ...) {
    old <- options(randomPlantedForest.exact_draw = FALSE)
    on.exit(options(old))
    set.seed(5)
    rpf(
      f,
      data = d,
      splits = NULL,
      prune = TRUE,
      max_panels = 10,
      split_try = 3,
      t_try = 0.1,
      max_interaction = 3,
      ntrees = 200,
      ...
    )
  }
  set.seed(5)
  reg <- rpf(y ~ x1 + x2 + x3, data = d)
  expect_s3_class(reg, "rpf_pruned")
  expect_false(reg$fit$params$exact_draw)
  expect_identical(def_parts(reg, d), def_parts(explicit(y ~ x1 + x2 + x3), d))
  for (y in c("cl", "cl3")) {
    f <- stats::as.formula(paste(y, "~ x1 + x2 + x3"))
    set.seed(5)
    m <- rpf(f, data = d)
    expect_identical(m$fit$params$loss, "exponential")
    expect_identical(def_parts(m, d, "prob"), def_parts(explicit(f, loss = "exponential"), d, "prob"))
  }
})

test_that("deterministic = TRUE resolves the defaults to 30 splits and no pruning", {
  d <- def_data(60)
  expect_warning(m <- rpf(y ~ x1 + x2 + x3, data = d, deterministic = TRUE), "200 tree families are identical")
  expect_s3_class(m, "rpf")
  expect_identical(m$params$splits, 30L)
  m1 <- rpf(y ~ x1 + x2 + x3, data = d, deterministic = TRUE, ntrees = 1)
  m30 <- rpf(y ~ x1 + x2 + x3, data = d, deterministic = TRUE, ntrees = 1, splits = 30, prune = FALSE)
  expect_identical(m1$fit$get_model(), m30$fit$get_model())
  expect_error(rpf(y ~ x1, data = d, deterministic = TRUE, ntrees = 1, prune = TRUE), "not available")
})

test_that("split_try and t_try default by split structure; NULL is not 'given' for tune", {
  d <- def_data(60)
  for (s in c("leaves", "hist", "cur_trees_1", "cur_trees_2", "res_trees")) {
    m <- rpf(y ~ x1 + x2, data = d, split_structure = s, ntrees = 1, splits = 5, prune = FALSE)
    expect_identical(c(m$params$split_try, m$params$t_try), if (s == "leaves") c(3, 0.1) else c(10, 0.4))
  }
  # giving a value is still possible
  m <- rpf(y ~ x1 + x2, data = d, ntrees = 1, splits = 5, split_try = 7, t_try = 0.3, prune = FALSE)
  expect_identical(c(m$params$split_try, m$params$t_try), c(7, 0.3))
  # NULL (the default) given explicitly with tune is not an error; a value is
  ns <- asNamespace("randomPlantedForest")
  given <- (function(splits = NULL, split_try = NULL, t_try = NULL, min_leaf_rows = 1, prune = NULL) {
    ns$tune_given(environment())
  })(split_try = NULL, prune = NULL)
  expect_false(any(given))
  expect_error(rpf(y ~ x1 + x2, data = d, ntrees = 10, split_try = 3, tune = "quick"), "split_try")
})
