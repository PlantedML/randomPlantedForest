# purify(mode = 2) stores the purified components as weighted rectangles (no dense grids). It must give the same
# components as the dense reference, mode 1, at every order, with a cap below the order, for multiclass
# classification, on new points (also outside the training range), and when predict_components() purifies itself.

rect_data <- function(n = 250) {
  set.seed(4701)
  x <- data.frame(
    x1 = runif(n),
    x2 = sample(0:6, n, replace = TRUE),
    x3 = rnorm(n),
    x4 = sample(0:2, n, replace = TRUE)
  )
  x$y <- x$x1 + 0.5 * x$x2 * x$x3 - 0.3 * x$x1 * x$x4 + 0.2 * x$x2 * x$x3 * x$x4 + rnorm(n, sd = 0.5)
  x
}

# new points: resampled training values, plus some outside the training range
rect_newdata <- function(dat, n = 300) {
  set.seed(4702)
  nd <- as.data.frame(lapply(dat[, c("x1", "x2", "x3", "x4")], function(v) sample(v, n, replace = TRUE)))
  nd[1:10, "x1"] <- c(-1, 2, -1, 2, -1, 2, -1, 2, -1, 2)
  nd[11:20, "x3"] <- c(-9, 9, -9, 9, -9, 9, -9, 9, -9, 9)
  nd
}

fit_rect <- function(dat, max_interaction, ...) {
  set.seed(47)
  rpf_050(y ~ x1 + x2 + x3 + x4, data = dat, max_interaction = max_interaction, ntrees = 12, splits = 40, ...)
}

expect_same_components <- function(a, b, tol = 1e-10) {
  expect_equal(colnames(a$m), colnames(b$m))
  expect_equal(as.matrix(a$m), as.matrix(b$m), tolerance = tol)
  expect_equal(a$intercept, b$intercept, tolerance = tol)
}

for (order in c(2, 3, 4)) {
  test_that(paste0("mode 2 equals mode 1 on training rows and new points, order ", order), {
    dat <- rect_data()
    nd <- rect_newdata(dat)
    f1 <- fit_rect(dat, order)
    f2 <- fit_rect(dat, order)
    f0 <- fit_rect(dat, order)
    pred_unpurified <- predict(f2, nd)

    purify(f1, mode = 1L)
    purify(f2, mode = 2L)
    for (d in list(dat, nd)) {
      expect_same_components(predict_components(f2, d), predict_components(f1, d))
    }
    # predict() of the purified forest equals that of the unpurified one, at any number of threads
    expect_equal(predict(f2, nd), pred_unpurified, tolerance = 1e-10)
    expect_identical(predict(f2, nd, nthreads = 2L), predict(f2, nd, nthreads = 1L))
    expect_equal(predict(f1, nd), pred_unpurified, tolerance = 1e-10)

    # predict_components() purifies an unpurified fit itself with mode 2
    expect_false(is_purified(f0))
    expect_same_components(predict_components(f0, nd), predict_components(f1, nd))
  })
}

test_that("mode 2 equals mode 1 with maxp_interaction below the order", {
  dat <- rect_data()
  nd <- rect_newdata(dat)
  for (cap in 1:2) {
    f1 <- fit_rect(dat, 3)
    f2 <- fit_rect(dat, 3)
    purify(f1, maxp_interaction = cap, mode = 1L)
    purify(f2, maxp_interaction = cap, mode = 2L)
    m1 <- predict_components(f1, nd)
    m2 <- predict_components(f2, nd)
    expect_same_components(m2, m1)
    # the components above the cap are zero
    high <- lengths(strsplit(colnames(m2$m), ":")) > cap
    expect_true(all(as.matrix(m2$m)[, high] == 0))
  }
})

test_that("mode 2 equals mode 1 for multiclass classification (value_size > 1)", {
  set.seed(4703)
  n <- 240
  dat <- data.frame(x1 = runif(n), x2 = sample(0:4, n, replace = TRUE), x3 = rnorm(n))
  eta <- cbind(dat$x1 * dat$x2, dat$x3 - dat$x1, 0)
  dat$y <- factor(apply(eta + matrix(rnorm(3 * n), n), 1, which.max), labels = c("a", "b", "c"))
  nd <- dat[sample(n, 100), c("x1", "x2", "x3")]
  nd$x3[1:5] <- c(-9, 9, -9, 9, -9)
  for (loss in c("L2", "logit")) {
    fits <- lapply(1:2, function(i) {
      set.seed(48)
      rpf_050(y ~ x1 + x2 + x3, data = dat, max_interaction = 3, ntrees = 8, splits = 30, loss = loss)
    })
    pred_unpurified <- predict(fits[[2]], nd, type = "numeric")
    purify(fits[[1]], mode = 1L)
    purify(fits[[2]], mode = 2L)
    m1 <- predict_components(fits[[1]], nd)
    m2 <- predict_components(fits[[2]], nd)
    expect_true(any(grepl("__class:c", colnames(m2$m))))
    expect_same_components(m2, m1)
    expect_equal(predict(fits[[2]], nd, type = "numeric"), pred_unpurified, tolerance = 1e-10)
  }
})

test_that("mode 2 does not depend on the number of threads and can be repeated or replaced by mode 1", {
  dat <- rect_data()
  nd <- rect_newdata(dat)
  f <- fit_rect(dat, 3)
  # purify_threads(maxp_interaction, nthreads, mode) purifies again also when the fit is purified already
  f$fit$purify_threads(0L, 1L, 2L)
  a <- predict_components(f, nd)
  f$fit$purify_threads(0L, 2L, 2L)
  b <- predict_components(f, nd)
  expect_identical(as.matrix(a$m), as.matrix(b$m))
  expect_identical(a$intercept, b$intercept)
  f$fit$purify_threads(0L, 1L, 1L)
  c1 <- predict_components(f, nd)
  expect_same_components(c1, a)
  f$fit$purify_threads(0L, 2L, 2L)
  expect_same_components(predict_components(f, nd), a)
})

test_that("mode 2 and its predictions do not depend on the number of threads, also at 3", {
  skip_on_cran() # more than 2 threads
  dat <- rect_data()
  nd <- rect_newdata(dat)
  f <- fit_rect(dat, 3)
  f$fit$purify_threads(0L, 1L, 2L)
  a <- predict_components(f, nd)
  p1 <- predict(f, nd, nthreads = 1L)
  f$fit$purify_threads(0L, 3L, 2L)
  b <- predict_components(f, nd)
  expect_identical(as.matrix(b$m), as.matrix(a$m))
  expect_identical(b$intercept, a$intercept)
  expect_identical(predict(f, nd, nthreads = 3L), p1)
})

test_that("rpf_boxes() of a family is the same after mode 1 and mode 2", {
  dat <- rect_data()
  f1 <- fit_rect(dat, 3)
  f2 <- fit_rect(dat, 3)
  purify(f1, mode = 1L)
  purify(f2, mode = 2L)
  for (family in c(1L, 7L)) {
    b1 <- rpf_boxes(f1, family = family)
    b2 <- rpf_boxes(f2, family = family)
    expect_identical(b2[names(b2) != "value"], b1[names(b1) != "value"])
    expect_equal(b2$value, b1$value, tolerance = 1e-10)
  }
  # one family's grid is that element of all families' grids
  expect_identical(f2$fit$get_grid_family(7L), f2$fit$get_grid_leaves()[[7]])
  expect_error(f2$fit$get_grid_family(13L), "No tree family")
})

test_that("purify_stats() describes the rectangle store of mode 2 only", {
  dat <- rect_data()
  f <- fit_rect(dat, 3)
  expect_identical(f$fit$purify_stats(), list())
  purify(f, mode = 2L)
  st <- f$fit$purify_stats()
  expect_named(st, c("n_components", "n_rects", "n_groups", "n_index_nodes", "n_table_cells", "bytes"))
  expect_equal(st$n_components, length(f$fit$get_purified_rects()$components))
  expect_true(all(unlist(st) > 0))
  f$fit$purify_threads(0L, 1L, 1L)
  expect_identical(f$fit$purify_stats(), list())
  expect_identical(f$fit$get_purified_rects(), list())
})
