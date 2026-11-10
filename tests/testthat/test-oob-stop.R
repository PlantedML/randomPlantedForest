# rpf(splits = NULL): the adaptive stop along the out-of-bag path

stop_sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- cbind(x1 = runif(n), x2 = runif(n), x3 = round(runif(n) * 6), x4 = rnorm(n))
  y <- 2 * x[, 1] + sin(4 * x[, 2]) + 3 * x[, 1] * x[, 2] + 0.3 * x[, 3] + rnorm(n, sd = 0.5)
  list(x = x, y = y)
}

stop_fit <- function(d, splits, ...) {
  set.seed(5)
  rpf(d$x, d$y, max_interaction = 2, ntrees = 6, splits = splits, ...)
}

# The OOB path at one checkpoint, as the C++ fit passes it to the stopper: rows i with OOB families of both halves
# predicting v[i] (sums 1 * v[i]), so the corrected error is mean((y - v)^2) (no variance between the halves)
path_at <- function(v) {
  n <- length(v)
  list(checkpoints = 1L, sum0 = matrix(v), sum1 = matrix(v), count = matrix(1L, n, 2))
}

test_that("splits = NULL stops one doubling past the lowest corrected OOB error and returns that forest", {
  d <- stop_sim()
  ad <- stop_fit(d, NULL)
  path <- ad$oob_path
  j <- which.min(path$oob_rmse)
  # stopped once the lowest error was one doubling (4 checkpoints) back, and chose it
  expect_identical(nrow(path), j + 4L)
  expect_identical(ad$params$splits, path$splits[j])
  expect_identical(ad$splits_oob, list(splits = path$splits[j], stopped_at = max(path$splits), capped = FALSE))
  # the default checkpoints, four per doubling, up to where it stopped
  cp <- oob_path_checkpoints(25600)
  expect_identical(path$splits, cp[seq_len(nrow(path))])
  # the path is that of a fixed fit grown to the stop; the forest that of a fixed fit with the chosen splits
  expect_identical(path, stop_fit(d, max(path$splits), oob_path = path$splits)$oob_path)
  ref <- stop_fit(d, path$splits[j])
  expect_identical(ad$fit$get_model(), ref$fit$get_model())
  expect_identical(predict(ad, d$x), predict(ref, d$x))
  # the C++ forest has the chosen number of splits and no OOB path
  expect_error(ad$fit$oob_path(), "no OOB path")
  # `oob_path = TRUE` changes nothing
  expect_identical(
    stop_fit(d, NULL, oob_path = TRUE)[c("params", "oob_path", "splits_oob")],
    ad[c("params", "oob_path", "splits_oob")]
  )
})

test_that("the stopper: patience, floor, the cap and ties", {
  y <- rep(0, 3)
  run <- function(st, v) {
    out <- 0L
    for (c in seq_along(v)) {
      out <- st$step(c, path_at(rep(v[c], 3)))
      if (out > 0) break
    }
    list(sel = out, c = c, res = st$result())
  }
  # rmse 0.5 0.4 0.3 0.35 0.36 0.37 0.38 0.39 ...: the lowest at 3, stop at 3 + patience
  v <- c(0.5, 0.4, 0.3, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41)
  r <- run(randomPlantedForest:::oob_stopper(y, 2, 1:10, patience = 4L), v)
  expect_identical(r[c("sel", "c")], list(sel = 3L, c = 7L))
  expect_identical(r$res$splits, 3L)
  expect_identical(r$res$stopped_at, 7L)
  expect_false(r$res$capped)
  expect_equal(r$res$path$oob_rmse, v[1:7])
  # a floor of 9 splits: no stop before 9
  r <- run(randomPlantedForest:::oob_stopper(y, 2, 1:10, patience = 4L, floor = 9), v)
  expect_identical(r[c("sel", "c")], list(sel = 3L, c = 9L))
  # still falling at the cap: capped; the minimum far enough back at the cap (only the floor held): not capped
  r <- run(randomPlantedForest:::oob_stopper(y, 2, 1:6, patience = 4L), c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1))
  expect_identical(r[c("sel", "c")], list(sel = 6L, c = 6L))
  expect_true(r$res$capped)
  r <- run(randomPlantedForest:::oob_stopper(y, 2, 1:6, patience = 4L, floor = 100), c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
  expect_identical(r[c("sel", "c")], list(sel = 1L, c = 6L))
  expect_false(r$res$capped)
  # ties (relative 1e-12) go to fewer splits
  sel <- randomPlantedForest:::oob_select
  expect_identical(sel(c(5, 3, 2 * (1 + 1e-13), 2, 2.1)), 3L)
  expect_identical(sel(c(5, 3, 2.01, 2, 2.1)), 4L)
  expect_identical(sel(c(1, 2, 3)), 1L)
  # a corrected error at or below 0 (its correction overshot) is left out unless all are; then the first
  expect_identical(sel(c(0.5, -0.01, 0.3, 0.31)), 3L)
  expect_identical(sel(c(0.5, 0, 0.2, -0.3, 0.1)), 5L)
  expect_identical(sel(c(0, -1e-3, -0.3)), 1L)
  expect_identical(randomPlantedForest:::oob_lowest(c(0.2, 0, 0.1, 0.1 * (1 + 1e-13), -1)), c(3L, 4L))
})

test_that("growing stops at the cap with a warning while the error still falls", {
  d <- stop_sim()
  ref <- stop_fit(d, 30)
  data <- ref$fit$get_data()
  set.seed(5)
  expect_warning(
    tu <- randomPlantedForest:::tune_splits_oob(data$Y, data$X, ref$params, checkpoints = 1:9),
    "cap \\(9 splits\\)"
  )
  expect_true(tu$capped)
  expect_identical(tu$stopped_at, 9L)
  expect_identical(tu$path$splits, 1:9)
  # the engine with the default checkpoints is rpf(splits = NULL)
  set.seed(5)
  tu <- randomPlantedForest:::tune_splits_oob(data$Y, data$X, ref$params)
  ad <- stop_fit(d, NULL)
  expect_identical(tu$splits, ad$params$splits)
  expect_identical(tu$path, ad$oob_path)
  expect_identical(tu$fit$get_model(), ad$fit$get_model())
})

test_that("splits = NULL is 30 for classification and deterministic forests", {
  d <- stop_sim()
  yc <- factor(d$y > median(d$y))
  set.seed(5)
  a <- rpf(d$x, yc, ntrees = 3, splits = NULL)
  set.seed(5)
  b <- rpf(d$x, yc, ntrees = 3, splits = 30)
  expect_identical(a$params$splits, 30L)
  expect_identical(a$fit$get_model(), b$fit$get_model())
  expect_null(a$splits_oob)
  det <- rpf(d$x, d$y, ntrees = 1, splits = NULL, deterministic = TRUE)
  expect_identical(det$params$splits, 30L)
  expect_identical(det$fit$get_model(), rpf(d$x, d$y, ntrees = 1, deterministic = TRUE)$fit$get_model())
})

test_that("splits = NULL: arguments, options and errors", {
  d <- stop_sim()
  expect_error(stop_fit(d, NULL, oob_path = c(5, 10)), "must be `TRUE` or `FALSE`")
  expect_error(stop_fit(d, 0), "splits")
  old <- options(randomPlantedForest.adaptive_pause_gb = -1)
  expect_error(stop_fit(d, NULL), "adaptive_pause_gb")
  options(old)
  # an error in the stopper stops the fit cleanly: the threads are joined and the next fit runs normally
  ref <- stop_fit(d, 30)
  data <- ref$fit$get_data()
  boom <- list(
    checkpoints = 1:20,
    stopper = function(c, path) if (c == 3) stop("boom") else 0L,
    patience = 4L,
    floor = 0,
    pause_gb = 1
  )
  params <- ref$params
  params$splits <- 20L
  params$nthreads <- 2L
  expect_error(randomPlantedForest:::rpf_impl(data$Y, data$X, "regression", params, adaptive = boom), "boom")
  bad <- modifyList(boom, list(stopper = function(c, path) 99L))
  expect_error(randomPlantedForest:::rpf_impl(data$Y, data$X, "regression", params, adaptive = bad), "stopper")
  expect_identical(stop_fit(d, 30)$fit$get_model(), ref$fit$get_model())
})

test_that("splits = NULL with purify, export_forest, print and marshal", {
  d <- stop_sim()
  ad <- stop_fit(d, NULL, purify = TRUE, export_forest = TRUE)
  ref <- stop_fit(d, ad$params$splits, purify = TRUE, export_forest = TRUE)
  expect_true(is_purified(ad))
  expect_identical(ad$forest, ref$forest)
  expect_identical(predict_components(ad, d$x[1:20, ]), predict_components(ref, d$x[1:20, ]))
  out <- paste(capture.output(print(ad)), collapse = " ")
  expect_match(out, sprintf("%d splits each", ad$params$splits), fixed = TRUE)
  expect_match(out, sprintf("chosen out of bag .*stopped at %d splits", ad$splits_oob$stopped_at))
  expect_no_match(paste(capture.output(print(ref)), collapse = " "), "out of bag")
  m <- rpf_unmarshal(rpf_marshal(ad, include_data = TRUE))
  expect_identical(m$splits_oob, ad$splits_oob)
  expect_identical(m$oob_path, ad$oob_path)
  expect_identical(predict(m, d$x), predict(ad, d$x))
})
