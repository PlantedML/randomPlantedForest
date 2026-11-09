# rpf(n_splits = NULL): the adaptive stop along the out-of-bag path

stop_sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- cbind(x1 = runif(n), x2 = runif(n), x3 = round(runif(n) * 6), x4 = rnorm(n))
  y <- 2 * x[, 1] + sin(4 * x[, 2]) + 3 * x[, 1] * x[, 2] + 0.3 * x[, 3] + rnorm(n, sd = 0.5)
  list(x = x, y = y)
}

stop_fit <- function(d, splits, ...) {
  set.seed(5)
  rpf_050(d$x, d$y, max_interaction = 2, n_families = 6, n_splits = splits, ...)
}

# The OOB path at one checkpoint, as the C++ fit passes it to the stopper: rows i with OOB families of both halves
# predicting v[i] (sums 1 * v[i]), so the corrected error is mean((y - v)^2) (no variance between the halves)
path_at <- function(v) {
  n <- length(v)
  list(checkpoints = 1L, sum0 = matrix(v), sum1 = matrix(v), count = matrix(1L, n, 2))
}

test_that("n_splits = NULL stops one doubling past the lowest corrected OOB error and returns that forest", {
  # the patience alone (no plateau stop)
  old <- options(randomPlantedForest.adaptive_tolerance = 0)
  on.exit(options(old))
  d <- stop_sim()
  ad <- stop_fit(d, NULL)
  path <- ad$oob_path
  j <- which.min(path$oob_rmse)
  # stopped once the lowest error was one doubling (4 checkpoints) back, and chose it
  expect_identical(nrow(path), j + 4L)
  expect_identical(ad$params$splits, path$n_splits[j])
  expect_identical(ad$splits_oob, list(n_splits = path$n_splits[j], stopped_at = max(path$n_splits), capped = FALSE))
  # the default checkpoints, four per doubling, up to where it stopped
  cp <- oob_path_checkpoints(25600)
  expect_identical(path$n_splits, cp[seq_len(nrow(path))])
  # the path is that of a fixed fit grown to the stop; the forest that of a fixed fit with the chosen splits
  expect_identical(path, stop_fit(d, max(path$n_splits), oob_path = path$n_splits)$oob_path)
  ref <- stop_fit(d, path$n_splits[j])
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

test_that("the stopper's plateau: a gain within `tolerance` paired standard errors stops; tolerance 0 is the patience", {
  # 4 rows with y = 0 whose corrected errors at checkpoint c are 1 - g * c plus +-0.05 in a pattern that changes from
  # checkpoint to checkpoint and sums to 0: the mean error falls by g per checkpoint, so it never stops by the patience,
  # and the paired differences of the rows have a standard error of about 0.05
  y <- rep(0, 4)
  pat <- list(c(1, -1, 1, -1), c(1, 1, -1, -1), c(1, -1, -1, 1))
  rows_at <- function(c, g) (1 - g * c) + 0.05 * pat[[(c %% 3) + 1]]
  run <- function(st, g, C = 10) {
    out <- 0L
    for (c in seq_len(C)) {
      out <- st$step(c, path_at(sqrt(rows_at(c, g))))
      if (out > 0) break
    }
    list(sel = out, c = c, res = st$result())
  }
  stopper <- function(...) randomPlantedForest:::oob_stopper(y, 2, 1:10, patience = 4L, ...)
  # a gain of 4 x 0.001 per doubling (4 checkpoints) is within one paired SE: P stops at the first chance (c = 5) and
  # chooses the lowest error (c = 5); the patience alone grows to the cap
  p <- run(stopper(tolerance = 1), 0.001)
  expect_identical(p[c("sel", "c")], list(sel = 5L, c = 5L))
  expect_false(p$res$capped)
  se <- sd(rows_at(1, 0.001) - rows_at(5, 0.001)) / 2
  expect_lt(4 * 0.001, se)
  t0 <- run(stopper(tolerance = 0), 0.001)
  expect_identical(t0[c("sel", "c")], list(sel = 10L, c = 10L))
  expect_true(t0$res$capped)
  # tolerance 0 is the patience alone, on this path and on the paths of the test above
  expect_identical(t0, run(stopper(), 0.001))
  # a floor of 8 splits: the plateau stops at 8
  expect_identical(run(stopper(tolerance = 1, floor = 8), 0.001)[c("sel", "c")], list(sel = 8L, c = 8L))
  # a gain of 4 x 0.03 is not within one SE (about 0.05), nor within 2; it is within 3
  expect_identical(run(stopper(tolerance = 1), 0.03)[c("sel", "c")], list(sel = 10L, c = 10L))
  expect_identical(run(stopper(tolerance = 2), 0.03)[c("sel", "c")], list(sel = 10L, c = 10L))
  expect_identical(run(stopper(tolerance = 3), 0.03)[c("sel", "c")], list(sel = 5L, c = 5L))
  # where the patience stops, the plateau changes nothing (the lowest error 4 checkpoints back)
  v <- c(0.5, 0.4, 0.3, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41)
  run1 <- function(st) {
    for (c in seq_along(v)) {
      out <- st$step(c, path_at(rep(v[c], 3)))
      if (out > 0) break
    }
    list(sel = out, c = c, res = st$result())
  }
  st <- function(tol) randomPlantedForest:::oob_stopper(rep(0, 3), 2, 1:10, patience = 4L, tolerance = tol)
  expect_identical(run1(st(1)), run1(st(0)))
})

test_that("the plateau keeps the rows it compares when a corrected error is below 0", {
  # checkpoint 2's correction overshoots (OOB means 0.5 and -0.5 around y = 0, 100 families: corrected MSE -0.245), so
  # it is left out of the choice; the plateau compares the lowest of the others (checkpoint 3) with the best and must
  # keep its rows: the stop runs to the last checkpoint and chooses 7, the lowest error above 0
  neg <- list(checkpoints = 1L, sum0 = matrix(rep(0.5, 3)), sum1 = matrix(rep(-0.5, 3)), count = matrix(1L, 3, 2))
  st <- randomPlantedForest:::oob_stopper(rep(0, 3), 100, 1:10, patience = 4L, tolerance = 1)
  v <- c(0.5, NA, 0.4, 0.45, 0.46, 0.47, 0.39, 0.395, 0.396, 0.397)
  for (c in seq_along(v)) {
    out <- st$step(c, if (is.na(v[c])) neg else path_at(rep(v[c], 3)))
    if (out > 0) break
  }
  expect_identical(c(out = out, c = c), c(out = 7L, c = 10L))
  expect_lt(st$points[2, 1], 0)
})

test_that("splits = NULL with the plateau stop: no later than the patience alone, the same forest at the same choice", {
  d <- stop_sim()
  fit_tol <- function(tol, ...) {
    old <- options(randomPlantedForest.adaptive_tolerance = tol)
    on.exit(options(old))
    stop_fit(d, NULL, ...)
  }
  p <- fit_tol(1)
  expect_identical(p$splits_oob, stop_fit(d, NULL)$splits_oob) # the default is tolerance 1
  t0 <- fit_tol(0)
  expect_lte(p$splits_oob$stopped_at, t0$splits_oob$stopped_at)
  # the path is the first checkpoints of today's, the choice its lowest error
  path <- p$oob_path
  expect_identical(path$splits, t0$oob_path$splits[seq_len(nrow(path))])
  expect_identical(path$oob_rmse, t0$oob_path$oob_rmse[seq_len(nrow(path))])
  expect_identical(p$params$splits, path$splits[which.min(path$oob_rmse)])
  expect_identical(p$fit$get_model(), stop_fit(d, p$params$splits)$fit$get_model())
  if (p$params$splits == t0$params$splits) {
    expect_identical(p$fit$get_model(), t0$fit$get_model())
  }
  # the families grow no further than with the patience alone
  expect_lte(p$fit$adaptive_stats()$splits_grown, t0$fit$adaptive_stats()$splits_grown)
  old <- options(randomPlantedForest.adaptive_tolerance = -1)
  expect_error(stop_fit(d, NULL), "adaptive_tolerance")
  options(randomPlantedForest.adaptive_tolerance = NA_real_)
  expect_error(stop_fit(d, NULL), "adaptive_tolerance")
  options(old)
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
  expect_identical(tu$path$n_splits, 1:9)
  # the engine with the default checkpoints is rpf(n_splits = NULL)
  set.seed(5)
  tu <- randomPlantedForest:::tune_splits_oob(data$Y, data$X, ref$params)
  ad <- stop_fit(d, NULL)
  expect_identical(tu$splits, ad$params$splits)
  expect_identical(tu$path, ad$oob_path)
  expect_identical(tu$fit$get_model(), ad$fit$get_model())
})

test_that("n_splits = NULL is chosen out of bag for classification too, and is 30 for deterministic forests", {
  d <- stop_sim()
  yc <- factor(d$y > median(d$y))
  set.seed(5)
  a <- rpf_050(d$x, yc, n_families = 3, n_splits = NULL)
  set.seed(5)
  b <- rpf_050(d$x, yc, n_families = 3, n_splits = a$params$splits)
  expect_identical(a$splits_oob$n_splits, a$params$splits)
  expect_identical(attr(a$oob_path, "criterion"), "logloss")
  expect_identical(a$fit$get_model(), b$fit$get_model())
  det <- rpf_050(d$x, d$y, n_families = 1, n_splits = NULL, deterministic = TRUE)
  expect_identical(det$params$splits, 30L)
  expect_identical(det$fit$get_model(), rpf_050(d$x, d$y, n_families = 1, deterministic = TRUE)$fit$get_model())
})

test_that("n_splits = NULL: arguments, options and errors", {
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
    tolerance = 1,
    floor = 0,
    pause_gb = 1
  )
  params <- ref$params
  params$splits <- 20L
  params$nthreads <- 2L
  expect_error(randomPlantedForest:::rpf_impl(data$Y, data$X, "regression", params, adaptive = boom), "boom")
  bad <- modifyList(boom, list(stopper = function(c, path) 99L))
  expect_error(randomPlantedForest:::rpf_impl(data$Y, data$X, "regression", params, adaptive = bad), "stopper")
  bad <- modifyList(boom, list(stopper = function(c, path) 0L, tolerance = -1))
  expect_error(randomPlantedForest:::rpf_impl(data$Y, data$X, "regression", params, adaptive = bad), "tolerance >= 0")
  expect_identical(stop_fit(d, 30)$fit$get_model(), ref$fit$get_model())
})

test_that("n_splits = NULL with purify, export_forest, print and marshal", {
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
