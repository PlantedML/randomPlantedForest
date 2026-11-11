# keep_forest: a pruned model holds no forest by default (rpf(), prune_components()); keep_forest = TRUE keeps it.
# Everything a pruned model does works without the forest, the same bit for bit; what needs the forest is one error.

kf_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n), x4 = runif(n))
  eta <- 2 * d$x1 - 1 + 4 * d$x2 * d$x3 * d$x4
  d$y <- eta + rnorm(n, sd = 0.1)
  d$yb <- factor(ifelse(eta + rnorm(n, sd = 0.3) > 0.5, "hi", "lo"))
  d$y3 <- factor(c("a", "b", "c")[1 + (eta > 0.2) + (eta > 0.9)])
  d
}

# the model without `fit` and the caches (environments compare by identity), and with its times (`seconds`, of the
# display and the tuning) set to 0
kf_model_data <- function(m) {
  m <- unclass(m)
  m$fit <- NULL
  m$cache <- NULL
  if (!is.null(m$display)) {
    m$display$cache <- NULL
  }
  no_times <- function(x) {
    if (!is.list(x)) {
      return(x)
    }
    if (is.data.frame(x)) {
      if ("seconds" %in% names(x)) {
        x$seconds <- 0
      }
      return(x)
    }
    for (k in seq_along(x)) {
      if (identical(names(x)[k], "seconds")) {
        x[[k]] <- 0
      } else if (!is.null(x[[k]])) {
        x[[k]] <- no_times(x[[k]])
      }
    }
    x
  }
  no_times(m)
}

# the model's fit (the forest's description) without the forest
kf_fit_data <- function(m) {
  f <- unclass(m$fit)
  f$fit <- NULL
  f$forest <- NULL
  f
}

test_that("rpf() drops the forest by default and keeps it with keep_forest = TRUE: the same model", {
  d <- kf_data()
  cases <- list(
    regression = list(f = y ~ x1 + x2 + x3 + x4, types = "numeric"),
    binary = list(f = yb ~ x1 + x2 + x3 + x4, types = c("prob", "class", "link")),
    three = list(f = y3 ~ x1 + x2 + x3 + x4, types = c("prob", "class", "link"))
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    set.seed(2)
    a <- rpf(cs$f, data = d, n_families = 10, max_panels = 3)
    rng_a <- .Random.seed
    set.seed(2)
    b <- rpf(cs$f, data = d, n_families = 10, max_panels = 3, keep_forest = TRUE)
    # the same random numbers drawn
    expect_identical(.Random.seed, rng_a, label = nm)
    expect_s3_class(a, "rpf_pruned")
    expect_s3_class(a$fit, "rpf")
    expect_null(a$fit[["fit", exact = TRUE]])
    expect_false(is.null(b$fit$fit))
    expect_true(rpf_is_valid(a))
    expect_false(rpf_is_valid(a$fit))
    expect_true(rpf_is_valid(b))
    # the same model, the same description of the forest
    expect_identical(kf_model_data(a), kf_model_data(b), label = nm)
    expect_identical(kf_fit_data(a), kf_fit_data(b), label = nm)
    for (type in cs$types) {
      expect_identical(predict(a, d, type = type), predict(b, d, type = type), label = paste(nm, type))
    }
    expect_identical(predict(a, d, nthreads = 2), predict(b, d, nthreads = 2))
    expect_identical(predict_components(a, d), predict_components(b, d), label = nm)
    expect_identical(predict_components(a, d, max_interaction = 1), predict_components(b, d, max_interaction = 1))
    expect_identical(format(a), format(b), label = nm)
    expect_identical(format(a, n_show = 2), format(b, n_show = 2))
  }
  # a coarsened component is predicted from its stored form without the forest
  set.seed(2)
  a <- rpf(y ~ x1 + x2 + x3 + x4, data = d, n_families = 10, max_panels = 3)
  if (on_macos_arm64()) {
    expect_gt(length(a$display$coarse), 0)
  }
})

test_that("prune = FALSE returns the forest whatever keep_forest is", {
  d <- kf_data(n = 100)
  set.seed(3)
  f1 <- rpf(y ~ x1 + x2, data = d, n_families = 5, n_splits = 10, prune = FALSE)
  set.seed(3)
  f2 <- rpf(y ~ x1 + x2, data = d, n_families = 5, n_splits = 10, prune = FALSE, keep_forest = TRUE)
  expect_s3_class(f1, "rpf")
  expect_true(rpf_is_valid(f1))
  expect_identical(predict(f1, d), predict(f2, d))
  expect_error(rpf(y ~ x1, data = d, keep_forest = NA), "keep_forest")
})

test_that("export_forest goes with the forest", {
  d <- kf_data(n = 100)
  set.seed(4)
  a <- rpf(y ~ x1 + x2, data = d, n_families = 10, export_forest = TRUE)
  expect_null(a$fit[["forest", exact = TRUE]])
  set.seed(4)
  b <- rpf(y ~ x1 + x2, data = d, n_families = 10, export_forest = TRUE, keep_forest = TRUE)
  expect_s3_class(b$fit$forest, "rpf_forest")
  expect_identical(predict(a, d), predict(b, d))
})

test_that("tuning returns the chosen model without its forest, or with it", {
  skip_on_cran()
  d <- kf_data(n = 150)
  set.seed(5)
  a <- rpf(y ~ x1 + x2 + x3, data = d, n_families = 10, tune = "quick")
  set.seed(5)
  b <- rpf(y ~ x1 + x2 + x3, data = d, n_families = 10, tune = "quick", keep_forest = TRUE)
  expect_null(a$fit[["fit", exact = TRUE]])
  expect_false(is.null(b$fit$fit))
  expect_identical(kf_model_data(a), kf_model_data(b))
  expect_identical(predict(a, d), predict(b, d))
})

test_that("prune_components() drops the forest by default: the caller still has it", {
  d <- kf_data()
  set.seed(6)
  fit <- rpf(y ~ x1 + x2 + x3, data = d, n_families = 10, n_splits = 40, prune = FALSE)
  p1 <- prune_components(fit, d, seed = 1)
  p2 <- prune_components(fit, d, seed = 1, keep_forest = TRUE)
  expect_null(p1$fit[["fit", exact = TRUE]])
  expect_identical(p2$fit, fit)
  expect_identical(kf_model_data(p1), kf_model_data(p2))
  expect_identical(predict(p1, d), predict(p2, d))
  expect_identical(format(p1), format(p2))
  # the caller's forest is untouched and still predicts
  expect_true(rpf_is_valid(fit))
  expect_s3_class(predict(fit, d), "tbl_df")
  # a given rule and tau, and by CV
  p3 <- prune_components(fit, d, rule = "R", tau = 0.01)
  expect_identical(predict(p3, d), predict(prune_components(fit, d, rule = "R", tau = 0.01, keep_forest = TRUE), d))
  expect_error(prune_components(fit, d, keep_forest = "yes"), "keep_forest")
})

test_that("what needs the forest is an error naming keep_forest = TRUE", {
  d <- kf_data(n = 150)
  set.seed(7)
  m <- rpf(y ~ x1 + x2 + x3, data = d, n_families = 10)
  f <- m$fit
  msg <- "keep_forest = TRUE"
  expect_error(prune_components(f, d), msg, fixed = TRUE)
  expect_error(prune_components(m, d), msg, fixed = TRUE)
  expect_error(purify(f), msg, fixed = TRUE)
  expect_error(is_purified(f), msg, fixed = TRUE)
  expect_error(predict(f, d), msg, fixed = TRUE)
  expect_error(predict_components(f, d), msg, fixed = TRUE)
  expect_error(plot(f), msg, fixed = TRUE)
  expect_error(rpf_boxes(f), msg, fixed = TRUE)
  expect_error(rpf_marshal(f), msg, fixed = TRUE)
  # the forest's description prints, and says the forest was not kept
  out <- paste(format(f), collapse = "\n")
  expect_match(out, "not kept")
  expect_match(out, "tree families")
  # with the forest, pruning again works
  set.seed(7)
  k <- rpf(y ~ x1 + x2 + x3, data = d, n_families = 10, keep_forest = TRUE)
  expect_s3_class(prune_components(k$fit, d, rule = "D", tau = 0.01), "rpf_pruned")
  expect_true(is_purified(k$fit) %in% c(TRUE, FALSE))
})

test_that("a pruned model without its forest survives saveRDS() and marshals to its own size", {
  d <- kf_data()
  set.seed(8)
  m <- rpf(y3 ~ x1 + x2 + x3 + x4, data = d, n_families = 10, max_panels = 3)
  pred <- predict(m, d, type = "prob")
  pc <- predict_components(m, d)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(m, tmp)
  r <- readRDS(tmp)
  expect_true(rpf_is_valid(r))
  expect_identical(predict(r, d, type = "prob"), pred)
  expect_identical(predict(r, d, type = "link"), predict(m, d, type = "link"))
  expect_identical(predict_components(r, d), pc)
  expect_identical(format(r), format(m))
  # marshal: the model as it is, without the caches; no forest inside
  blob <- rpf_marshal(m)
  expect_s3_class(blob, "rpf_pruned_marshaled")
  expect_false(inherits(blob$fit, "rpf_marshaled"))
  expect_lt(as.numeric(utils::object.size(blob)), 1.1 * as.numeric(utils::object.size(m)))
  saveRDS(blob, tmp)
  u <- rpf_unmarshal(readRDS(tmp))
  expect_s3_class(u, "rpf_pruned")
  expect_identical(predict(u, d, type = "prob"), pred)
  expect_identical(predict_components(u, d), pc)
  expect_identical(kf_model_data(u), kf_model_data(m))
  # with the forest, the blob carries it, and the model still needs rpf_marshal()
  set.seed(8)
  k <- rpf(y3 ~ x1 + x2 + x3 + x4, data = d, n_families = 10, max_panels = 3, keep_forest = TRUE)
  bk <- rpf_marshal(k)
  expect_s3_class(bk$fit, "rpf_marshaled")
  expect_gt(as.numeric(utils::object.size(bk)), as.numeric(utils::object.size(blob)))
  expect_identical(predict(rpf_unmarshal(bk), d, type = "prob"), pred)
  saveRDS(k, tmp)
  expect_false(rpf_is_valid(readRDS(tmp)))
})

test_that("bundle round-trip of a pruned model without its forest", {
  skip_if_not_installed("bundle")
  d <- kf_data(n = 150)
  set.seed(9)
  m <- rpf(y ~ x1 + x2 + x3, data = d, n_families = 10)
  b <- bundle::bundle(m)
  expect_s3_class(b, "bundled_rpf_pruned")
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(b, tmp)
  restored <- bundle::unbundle(readRDS(tmp))
  expect_s3_class(restored, "rpf_pruned")
  expect_identical(predict(restored, d), predict(m, d))
})

test_that("pruned models saved with their forest by earlier versions load, print and predict", {
  # fixtures/pruned_with_forest.rds: rpf_marshal() of two pruned models of the version before keep_forest (D2,
  # c0c34ea8), which kept their forest, with their predictions, components and print then
  # (fixtures/make-pruned-with-forest.R); made on macOS arm64, so predictions are compared to 1e-12 (libm)
  fx <- readRDS(test_path("fixtures", "pruned_with_forest.rds"))
  for (nm in c("reg", "cls")) {
    r <- rpf_unmarshal(fx[[nm]])
    expect_s3_class(r, "rpf_pruned")
    expect_true(rpf_is_valid(r))
    expect_false(is.null(r$fit$fit))
    # (the first line is a rule, drawn in ASCII under testthat)
    expect_identical(format(r)[-1], fx$print[[nm]][-1], label = nm)
    expect_equal(
      as.matrix(predict_components(r, fx$new_data)$m),
      as.matrix(fx$components[[nm]]),
      tolerance = 1e-12,
      label = nm
    )
    # the restored forest predicts too
    expect_s3_class(predict(r$fit, fx$new_data), "tbl_df")
  }
  r <- rpf_unmarshal(fx$reg)
  expect_gt(length(r$display$coarse), 0)
  expect_equal(predict(r, fx$new_data), fx$pred$reg, tolerance = 1e-12)
  r <- rpf_unmarshal(fx$cls)
  expect_equal(predict(r, fx$new_data, type = "prob"), fx$pred$cls, tolerance = 1e-12)
  expect_equal(predict(r, fx$new_data, type = "link"), fx$pred$cls_link, tolerance = 1e-12)
  expect_identical(predict(r, fx$new_data, type = "class"), fx$pred$cls_class)
  # and the pruned model marshalled at N16 (fixtures/classif_before_n18.rds)
  n16 <- readRDS(test_path("fixtures", "classif_before_n18.rds"))
  r <- rpf_unmarshal(n16$pruned3)
  expect_false(is.null(r$fit$fit))
  expect_equal(predict(r, n16$new_data, type = "prob"), n16$pred$pruned3, tolerance = 1e-12)
  expect_match(paste(format(r), collapse = "\n"), "How it was chosen")
})
