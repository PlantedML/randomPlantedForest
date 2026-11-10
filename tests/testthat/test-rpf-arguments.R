test_that("invalid arguments are rejected", {
  fit_reg <- function(...) rpf_050(mpg ~ cyl + wt, data = mtcars, n_families = 2, ...)

  expect_error(fit_reg(max_interaction = -1))
  expect_error(fit_reg(n_families = 0))
  expect_error(fit_reg(n_splits = 0))
  expect_error(fit_reg(split_structure = "branches"))
  expect_error(fit_reg(n_try_split_points = 0))
  expect_error(fit_reg(share_candidates = 1.5))
  expect_error(fit_reg(max_candidates = 0))
  expect_error(fit_reg(split_decay_rate = -1))
  expect_error(fit_reg(delete_leaves = NA))
  expect_error(fit_reg(loss = "logit"))
  expect_error(fit_reg(purify = "yes"))
  expect_error(fit_reg(nthreads = 0))
  expect_error(fit_reg(export_forest = NA))
  expect_error(fit_reg(deterministic = 1))

  binary <- data.frame(x = mtcars$wt, y = factor(mtcars$am))
  fit_classif <- function(...) rpf_050(y ~ x, data = binary, n_families = 2, loss = "logit", ...)
  expect_error(fit_classif(delta = 2))
  expect_error(fit_classif(epsilon = -0.1))
})

test_that("split search parameters are passed through and stored", {
  set.seed(1)
  fit <- rpf_050(
    mpg ~ cyl + wt + hp,
    data = mtcars,
    n_families = 3,
    n_splits = 12,
    n_try_split_points = 3,
    share_candidates = 0.8,
    max_candidates = 5,
    split_decay_rate = 0,
    delete_leaves = FALSE
  )

  expect_identical(
    fit$params[c("splits", "split_try", "t_try", "max_candidates", "split_decay_rate", "delete_leaves")],
    list(splits = 12, split_try = 3, t_try = 0.8, max_candidates = 5, split_decay_rate = 0, delete_leaves = FALSE)
  )
  expect_true(all(is.finite(predict(fit, mtcars)$.pred)))
})

test_that("delete_leaves changes the fitted forest", {
  fit_with <- function(delete_leaves) {
    set.seed(1)
    fit <- rpf_050(mpg ~ cyl + wt + hp, data = mtcars, n_families = 3, n_splits = 30, delete_leaves = delete_leaves)
    predict(fit, mtcars)$.pred
  }

  expect_false(identical(fit_with(TRUE), fit_with(FALSE)))
})

test_that("delta and epsilon are used by the logit loss", {
  set.seed(1)
  dat <- data.frame(x1 = runif(100), x2 = runif(100))
  dat$y <- factor(ifelse(dat$x1 > 0.5, "a", "b"))
  fit_with <- function(...) {
    set.seed(2)
    fit <- rpf_050(y ~ ., data = dat, n_families = 3, loss = "logit", ...)
    predict(fit, dat, type = "link")$.pred
  }

  default <- fit_with()
  expect_false(identical(default, fit_with(delta = 0.1)))
  expect_false(identical(default, fit_with(epsilon = 0.3)))
})

test_that("max_interaction above the number of predictors is capped silently", {
  expect_no_message(fit <- rpf_050(mpg ~ cyl + wt, data = mtcars, n_families = 2, max_interaction = 5))
  expect_identical(fit$params$max_interaction, 2L)

  expect_no_message(fit_single <- rpf_050(mpg ~ wt, data = mtcars, n_families = 2))
  expect_identical(fit_single$params$max_interaction, 1L)
})

test_that("logical predictors are used like 0/1 integers", {
  set.seed(1)
  dat <- data.frame(x = runif(100), flag = runif(100) > 0.5)
  dat$y <- dat$x + 2 * dat$flag + rnorm(100, sd = 0.1)
  dat_int <- transform(dat, flag = as.integer(flag))

  set.seed(2)
  fit_lgl <- rpf_050(y ~ ., data = dat, n_families = 3)
  set.seed(2)
  fit_int <- rpf_050(y ~ ., data = dat_int, n_families = 3)

  expect_identical(predict(fit_lgl, dat), predict(fit_int, dat_int))
  expect_named(fit_lgl$blueprint$ptypes$predictors, c("x", "flag"))

  set.seed(2)
  fit_xy <- rpf_050(x = dat[c("x", "flag")], y = dat$y, n_families = 3)
  expect_identical(predict(fit_xy, dat), predict(fit_lgl, dat))
})

test_that("multiclass fits handle a single-level factor predictor", {
  set.seed(1)
  dat <- data.frame(x = runif(60), constant = factor("a"))
  dat$y <- factor(cut(dat$x, 3, labels = c("lo", "mid", "hi")))

  fit <- rpf_050(y ~ ., data = dat, n_families = 3)
  pred <- predict(fit, dat, type = "class")

  expect_gt(mean(pred$.pred_class == dat$y), 0.8)
})

test_that("purify() on a purified forest returns it unchanged", {
  fit <- rpf_050(mpg ~ cyl + wt, data = mtcars, max_interaction = 2, n_families = 3, purify = TRUE)
  before <- predict_components(fit, mtcars)

  expect_invisible(purify(fit))
  expect_true(is_purified(fit))
  expect_identical(predict_components(fit, mtcars), before)
})

test_that("the default loss depends on the outcome", {
  expect_identical(rpf_050(mpg ~ cyl + wt, data = mtcars, n_families = 2)$params$loss, "L2")
  expect_identical(rpf_050(Species ~ ., data = iris, n_families = 2)$params$loss, "exponential")
  expect_error(rpf_050(Species ~ ., data = iris, n_families = 2, loss = c("L1", "L2")))
})
