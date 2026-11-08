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

test_that("include_data = TRUE stores the in-bag counts", {
  set.seed(3)
  fit <- rpf(mpg ~ wt + cyl + hp, data = mtcars, ntrees = 6, max_interaction = 2)
  ib <- fit$fit$inbag()
  blob <- rpf_marshal(fit, include_data = TRUE)
  expect_type(blob$fit_state$inbag, "raw")
  restored <- rpf_unmarshal(blob)
  expect_identical(restored$fit$inbag(), ib)
  # the out-of-bag components after restoring are those of the original forest
  expect_identical(restored$fit$purify_oob(0L, 1L), fit$fit$purify_oob(0L, 1L))

  # without the data, no counts: a restored forest cannot be purified anyway
  plain <- rpf_marshal(fit)
  expect_null(plain$fit_state$inbag)
  expect_error(rpf_unmarshal(plain)$fit$inbag(), "include_data = TRUE")
  # a blob with data from a version without the field restores without counts
  old <- blob
  old$fit_state$inbag <- NULL
  expect_error(rpf_unmarshal(old)$fit$inbag(), "no in-bag counts")
  expect_error(rpf_unmarshal(old)$fit$purify_oob(0L, 1L), "no in-bag counts")
  # corrupt counts are rejected
  bad <- blob
  bad$fit_state$inbag <- bad$fit_state$inbag[-1]
  expect_error(rpf_unmarshal(bad), "in-bag counts do not match")
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
  purify(pfit, mode = 1L)
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

test_that("a purification capped by maxp_interaction survives the round trip", {
  set.seed(2)
  fit <- rpf(mpg ~ wt + cyl + hp + disp, data = mtcars, ntrees = 10, max_interaction = 3)
  p0 <- predict(fit, mtcars)
  purify(fit, maxp_interaction = 1)
  restored <- expect_roundtrip(fit, mtcars)
  expect_identical(predict(restored, mtcars), p0)
  expect_true(restored$fit$is_purify_capped())
  expect_identical(predict_components(restored, mtcars), predict_components(fit, mtcars))
})

test_that("purified round-trip keeps mode 1 grids and mode 2 rectangles apart", {
  fit <- rpf(mpg ~ wt + cyl + hp, data = mtcars, ntrees = 10, max_interaction = 3)
  purify(fit, mode = 1L)
  blob <- rpf_marshal(fit)
  expect_null(blob$fit_state$rects)
  expect_false(is.null(blob$fit_state$grid))
  comp1 <- predict_components(fit, mtcars)
  expect_identical(predict_components(rpf_unmarshal(blob), mtcars), comp1)

  fit$fit$purify_threads(0L, 0L, 2L) # purify() is a no-op on a purified fit
  blob <- rpf_marshal(fit)
  expect_null(blob$fit_state$grid)
  expect_false(is.null(blob$fit_state$rects))
  comp2 <- predict_components(fit, mtcars)
  expect_equal(as.matrix(comp2$m), as.matrix(comp1$m), tolerance = 1e-10)
  restored <- rpf_unmarshal(blob)
  expect_identical(predict_components(restored, mtcars), comp2)
  expect_identical(predict(restored, mtcars), predict(fit, mtcars))
  expect_identical(restored$fit$purify_stats(), fit$fit$purify_stats())

  # rpf_boxes() computes a family's grid from the training data
  expect_error(rpf_boxes(restored), "include_data = TRUE")
  with_data <- rpf_unmarshal(rpf_marshal(fit, include_data = TRUE))
  expect_identical(rpf_boxes(with_data, family = 2L), rpf_boxes(fit, family = 2L))

  # a corrupt store is rejected: point an index leaf of a two-feature group past its rectangles
  bad <- blob
  comps <- bad$fit_state$rects$components
  i <- which(lengths(lapply(comps, `[[`, "S")) == 2)[1]
  j <- which(vapply(comps[[i]]$groups, function(g) length(g$ids) > 0, logical(1)))[1]
  bad$fit_state$rects$components[[i]]$groups[[j]]$ids[1] <- 1e6L
  expect_error(rpf_unmarshal(bad), "Corrupt purified components")
})

test_that("a restored capped purification keeps its order, and is purified again only with data", {
  set.seed(2)
  fit <- rpf(mpg ~ wt + cyl + hp + disp, data = mtcars, ntrees = 10, max_interaction = 3)
  purify(fit, maxp_interaction = 2)

  restored <- expect_roundtrip(fit, mtcars)
  expect_identical(restored$fit$get_purify_cap(), 2L)
  # at the cap: nothing to do
  expect_invisible(purify(restored, maxp_interaction = 2))
  # above it, without training data: a clear error, and the forest is left as it was
  expect_error(purify(restored), "cannot be purified to higher orders")
  expect_identical(restored$fit$get_purify_cap(), 2L)
  expect_identical(predict_components(restored, mtcars), predict_components(fit, mtcars))

  # with the training data it is purified again
  with_data <- rpf_unmarshal(rpf_marshal(fit, include_data = TRUE))
  purify(with_data)
  purify(fit)
  expect_false(with_data$fit$is_purify_capped())
  expect_identical(predict_components(with_data, mtcars), predict_components(fit, mtcars))
})

test_that("a blob without the cap's order restores a capped purification", {
  set.seed(2)
  fit <- rpf(mpg ~ wt + cyl + hp + disp, data = mtcars, ntrees = 10, max_interaction = 3)
  purify(fit, maxp_interaction = 2)
  blob <- rpf_marshal(fit, include_data = TRUE)
  # as written before the field existed: the order is unknown and counts as 1
  blob$fit_state$purify_cap <- NULL
  restored <- rpf_unmarshal(blob)
  expect_true(restored$fit$is_purify_capped())
  expect_identical(restored$fit$get_purify_cap(), 1L)
  expect_identical(predict(restored, mtcars), predict(fit, mtcars))
  purify(restored, maxp_interaction = 2)
  expect_identical(restored$fit$get_purify_cap(), 2L)
  expect_identical(predict_components(restored, mtcars), predict_components(fit, mtcars))
})

test_that("a capped purification keeps its store and its cap through the round trip, in both modes", {
  for (mode in 1:2) {
    set.seed(6)
    fit <- rpf(mpg ~ wt + cyl + hp + disp, data = mtcars, ntrees = 10, max_interaction = 3)
    p0 <- predict(fit, mtcars)
    purify(fit, maxp_interaction = 2, mode = mode)
    blob <- rpf_marshal(fit, include_data = TRUE)
    # mode 2 stores rectangles, mode 1 grids; either way with the cap's order
    if (mode == 2) {
      expect_false(is.null(blob$fit_state$rects))
      expect_null(blob$fit_state$grid)
    } else {
      expect_null(blob$fit_state$rects)
      expect_false(is.null(blob$fit_state$grid))
    }
    expect_identical(blob$fit_state$purify_cap, 2L)

    restored <- rpf_unmarshal(blob)
    expect_identical(restored$fit$get_purify_cap(), 2L)
    expect_true(restored$fit$is_purify_capped())
    expect_identical(predict(restored, mtcars), p0)
    expect_identical(predict_components(restored, mtcars), predict_components(fit, mtcars))

    # purified again above the cap from the restored data, as the original forest is
    purify(restored, mode = mode)
    purify(fit, mode = mode)
    expect_identical(restored$fit$get_purify_cap(), 0L)
    expect_identical(predict_components(restored, mtcars), predict_components(fit, mtcars))
    expect_null(predict_components(restored, mtcars)$remainder)

    # and the full purification survives a second round trip uncapped
    again <- rpf_unmarshal(rpf_marshal(restored))
    expect_identical(again$fit$get_purify_cap(), 0L)
    expect_identical(predict_components(again, mtcars), predict_components(fit, mtcars))
  }
})

prune_marshal_fit <- function(n = 300) {
  set.seed(11)
  x <- matrix(runif(n * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
  y <- 1 + 2 * x[, 1] + 4 * x[, 2] * x[, 3] + rnorm(n, sd = 0.3)
  set.seed(5)
  fit <- rpf(x, y, max_interaction = 3, ntrees = 10, splits = 40)
  list(x = x, y = y, fit = fit)
}

test_that("a pruned forest survives the round trip and predicts the same", {
  d <- prune_marshal_fit()
  for (rule in c("D", "R")) {
    pr <- prune_components(d$fit, d$x, d$y, rule = rule, tau = c(0, 0.01), folds = 3, seed = 1)
    blob <- rpf_marshal(pr)
    expect_s3_class(blob, "rpf_pruned_marshaled")
    expect_s3_class(blob$fit, "rpf_marshaled")
    expect_null(blob$fit$fit_state$data)
    tmp <- tempfile(fileext = ".rds")
    saveRDS(blob, tmp)
    restored <- rpf_unmarshal(readRDS(tmp))
    expect_s3_class(restored, "rpf_pruned")
    expect_s3_class(restored$fit, "rpf")
    expect_true(rpf_is_valid(restored))
    expect_identical(predict(restored, d$x), predict(pr, d$x))
    expect_identical(predict(restored, d$x, nthreads = 2), predict(pr, d$x))
    expect_identical(predict_components(restored, d$x), predict_components(pr, d$x))
    expect_identical(
      predict_components(restored, d$x, max_interaction = 1),
      predict_components(pr, d$x, max_interaction = 1)
    )
    expect_identical(format(restored), format(pr))
    # everything but the forest is kept as it is
    keep <- setdiff(names(pr), "fit")
    expect_identical(unclass(restored)[keep], unclass(pr)[keep])
    expect_identical(predict(restored$fit, d$x), predict(d$fit, d$x))
  }
})

test_that("a pruned forest marshaled with its data can be pruned again", {
  d <- prune_marshal_fit()
  pr <- prune_components(d$fit, d$x, d$y, rule = "R", tau = 0.01)
  blob <- rpf_marshal(pr, include_data = TRUE)
  expect_false(is.null(blob$fit$fit_state$data))
  restored <- rpf_unmarshal(blob)
  again <- prune_components(restored$fit, d$x, d$y, rule = "R", tau = 0.01)
  expect_identical(again$beta, pr$beta)
  expect_identical(again$components, pr$components)
  # the in-bag counts come along: the out-of-bag choice and R's out-of-bag coefficients are the same
  oob_again <- prune_components(restored$fit, d$x, d$y, rule = "R", tau = c(0, 0.01), seed = 1)
  oob_ref <- prune_components(d$fit, d$x, d$y, rule = "R", tau = c(0, 0.01), seed = 1)
  expect_identical(oob_again$select_by, "oob")
  expect_identical(oob_again$oob, oob_ref$oob)
  expect_identical(oob_again$beta, oob_ref$beta)
  expect_identical(oob_again$intercept, oob_ref$intercept)
  # without its data, the restored forest still predicts, but cannot be pruned again
  plain <- rpf_unmarshal(rpf_marshal(pr))
  expect_identical(predict(plain, d$x), predict(pr, d$x))
  expect_error(prune_components(plain$fit, d$x, d$y, rule = "R", tau = 0.01), "no training data")
})

test_that("a pruned forest restored without marshal gives an actionable error", {
  d <- prune_marshal_fit(n = 100)
  pr <- prune_components(d$fit, d$x, d$y, rule = "D", tau = 0.01)
  tmp <- tempfile(fileext = ".rds")
  saveRDS(pr, tmp)
  zombie <- readRDS(tmp)
  expect_false(rpf_is_valid(zombie))
  expect_error(predict(zombie, d$x), "rpf_marshal")
  expect_error(predict_components(zombie, d$x), "rpf_marshal")
  expect_error(rpf_marshal(zombie), "readRDS")
  expect_true(rpf_is_valid(pr))
  expect_error(rpf_marshal(pr, include_data = NA))
})
