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
  comp_before <- predict_components(fit, iris)
  restored <- expect_roundtrip(fit, iris)
  expect_true(is_purified(restored))
  expect_identical(predict_components(restored, iris), comp_before)
})

test_that("marshaled blob contains no external pointers", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5)
  blob <- rpf_marshal(fit)
  has_extptr <- function(x) {
    if (typeof(x) == "externalptr") {
      return(TRUE)
    }
    # unclass so classed list-likes (e.g. package_version, whose `[[` returns
    # another package_version) recurse over their raw contents and terminate
    x <- unclass(x)
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
  expect_error(rpf_marshal(zombie), "readRDS")
  expect_true(rpf_is_valid(fit))
})

test_that("rpf_unmarshal errors on malformed blob missing fit_state", {
  expect_error(rpf_unmarshal(structure(list(), class = "rpf_marshaled")), "Unsupported")
})

test_that("corrupt blobs error instead of reading out of bounds", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5)

  # Interval matrix narrower than feature_size
  blob <- rpf_marshal(fit)
  blob$fit_state$model[[1]]$intervals[[1]][[1]] <- matrix(0, nrow = 2, ncol = 1)
  expect_error(rpf_unmarshal(blob), "interval matrix")

  # Purified grid: values matrix with wrong dimensions
  pfit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5)
  purify(pfit)
  # trees[[1]] is the null tree whose 1x1 values matrix is legitimately tiny;
  # corrupt a real component tree instead
  pblob <- rpf_marshal(pfit)
  pblob$fit_state$grid[[1]]$trees[[2]]$values <- matrix(0, nrow = 1, ncol = 1)
  expect_error(rpf_unmarshal(pblob), "values matrix")

  # Purified grid: non-positive grid dimension
  pblob2 <- rpf_marshal(pfit)
  pblob2$fit_state$grid[[1]]$trees[[2]]$dims <- -1L
  expect_error(rpf_unmarshal(pblob2), "non-positive")

  # Training data with wrong number of feature columns
  dblob <- rpf_marshal(fit, include_data = TRUE)
  dblob$fit_state$data$X <- dblob$fit_state$data$X[, 1, drop = FALSE]
  expect_error(rpf_unmarshal(dblob), "columns")
})

test_that("exported forest is rebuilt instead of duplicated in the blob", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5, export_forest = TRUE)
  blob <- rpf_marshal(fit)
  expect_null(blob$forest)
  restored <- rpf_unmarshal(blob)
  expect_s3_class(restored$forest, "rpf_forest")
  expect_identical(restored$forest, fit$forest)
})

test_that("blob records package version and warns when saved with a newer one", {
  fit <- rpf(mpg ~ wt + cyl, data = mtcars, ntrees = 5)
  blob <- rpf_marshal(fit)
  expect_identical(blob$fit_state$pkg_version, utils::packageVersion("randomPlantedForest"))

  blob$fit_state$pkg_version <- package_version("999.0.0")
  expect_warning(restored <- rpf_unmarshal(blob), "999.0.0")
  expect_identical(predict(restored, mtcars), predict(fit, mtcars))
})
