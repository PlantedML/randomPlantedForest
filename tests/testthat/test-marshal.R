test_that("params-only constructors create empty, unfitted objects", {
  # 13-element parameter vector: max_interaction, ntrees, splits, split_try,
  # t_try, purify, deterministic, nthreads, cv, split_decay_rate,
  # max_candidates, delete_leaves, split_mode
  pars <- c(1, 50, 30, 10, 0.4, 0, 0, 1, 0, 0.1, 50, 1, 3)
  fit <- new(RandomPlantedForest, pars)
  expect_length(fit$get_model(), 0)
  expect_false(fit$is_purified())

  cfit <- new(ClassificationRPF, "L2", c(pars, 0, 0.1))
  expect_length(cfit$get_model(), 0)

  expect_error(new(ClassificationRPF, "L2", pars), "15 parameters")
})

test_that("get_data returns training data from a fitted forest", {
  X <- as.matrix(mtcars[, c("wt", "cyl")])
  Y <- as.matrix(mtcars$mpg)
  fit <- new(RandomPlantedForest, Y, X, c(1, 10, 30, 10, 0.4, 0, 0, 1, 0, 0.1, 50, 1, 3))
  dat <- fit$get_data()
  expect_equal(unname(as.matrix(dat$X)), unname(X))
  expect_equal(unname(as.matrix(dat$Y)), unname(Y))
})

expect_roundtrip <- function(fit, new_data) {
  p_before <- predict(fit, new_data)
  blob <- rpf_marshal(fit)
  tmp <- tempfile(fileext = ".rds")
  saveRDS(blob, tmp)
  restored <- rpf_unmarshal(readRDS(tmp))
  expect_s3_class(restored, "rpf")
  expect_identical(predict(restored, new_data), p_before)
  restored
}

test_that("regression round-trip, no data", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 10)
  expect_roundtrip(fit, mtcars)
})

test_that("binary classification round-trip, no data", {
  # NOTE: `am` is pre-converted to a factor column rather than using
  # `as.factor(am)` inline in the formula - hardhat::mold() does not capture
  # an inline outcome transform in the blueprint ptype, which breaks
  # `predict()` independently of marshaling (pre-existing hardhat/formula
  # limitation, not specific to this task).
  dat <- transform(mtcars, am = factor(am))
  fit <- rpf(am ~ wt + hp, data = dat, ntrees = 10, loss = "logit")
  restored <- expect_roundtrip(fit, dat)
  expect_identical(
    predict(restored, dat, type = "prob"),
    predict(fit, dat, type = "prob")
  )
})

test_that("multiclass round-trip, no data", {
  fit <- rpf(Species ~ ., data = iris, ntrees = 10, loss = "L2")
  restored <- expect_roundtrip(fit, iris)
  expect_identical(
    predict(restored, iris, type = "prob"),
    predict(fit, iris, type = "prob")
  )
})

test_that("purified regression round-trip without data", {
  fit <- rpf(mpg ~ wt + cyl + hp, data = mtcars, ntrees = 10, max_interaction = 2)
  purify(fit)
  comp_before <- predict_components(fit, mtcars)
  restored <- expect_roundtrip(fit, mtcars)
  expect_true(is_purified(restored))
  expect_identical(predict_components(restored, mtcars), comp_before)
})

test_that("purified multiclass round-trip without data", {
  fit <- rpf(Species ~ ., data = iris, ntrees = 10, loss = "L2")
  purify(fit)
  expect_roundtrip(fit, iris)
})

test_that("marshaled blob contains no external pointers", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5)
  blob <- rpf_marshal(fit)
  has_extptr <- function(x) {
    if (typeof(x) == "externalptr") {
      return(TRUE)
    }
    if (is.list(x)) {
      return(any(vapply(x, has_extptr, logical(1))))
    }
    FALSE
  }
  expect_false(has_extptr(unclass(blob)))
})

test_that("include_data = TRUE allows purify after restore", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 10, max_interaction = 2)

  restored <- rpf_unmarshal(rpf_marshal(fit, include_data = TRUE))
  purify(restored)

  # restored shares the identical forest and training data with fit,
  # so purifying each independently must give identical predictions
  purify(fit)
  expect_identical(predict(restored, mtcars), predict(fit, mtcars))
})

test_that("purify on a data-free restored forest errors informatively", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5)
  restored <- rpf_unmarshal(rpf_marshal(fit))
  expect_error(purify(restored), "include_data")
})

test_that("restored-without-marshal rpf gives actionable error", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5)
  tmp <- tempfile(fileext = ".rds")
  saveRDS(fit, tmp)
  zombie <- readRDS(tmp)
  expect_false(rpf_is_valid(zombie))
  expect_error(predict(zombie, mtcars), "rpf_marshal")
  expect_error(purify(zombie), "rpf_marshal")
  expect_error(predict_components(zombie, mtcars), "rpf_marshal")
  expect_true(rpf_is_valid(fit))
})
