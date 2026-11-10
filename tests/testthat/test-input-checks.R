input_data <- function(n = 60) {
  set.seed(42)
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  x$y <- x$x1 + x$x2 + rnorm(n, sd = 0.1)
  x
}

test_that("an infinite outcome is an error, infinite predictors are not", {
  dat <- input_data()
  dat$y[3] <- Inf
  expect_error(rpf_050(y ~ x1 + x2, data = dat, n_families = 2), "The outcome must be finite")
  expect_error(rpf_050(x = dat[c("x1", "x2")], y = -dat$y, n_families = 2), "1 infinite value")
  # a transformed outcome is checked as well
  dat <- input_data()
  dat$x1[1] <- 0
  expect_error(rpf_050(log(x1) ~ x2, data = dat, n_families = 2), "The outcome must be finite")

  dat <- input_data()
  dat$x1[c(2, 5)] <- c(Inf, -Inf)
  set.seed(1)
  fit <- rpf_050(y ~ x1 + x2, data = dat, n_families = 2)
  expect_true(all(is.finite(predict(fit, dat)$.pred)))
})

test_that("epsilon and delta must give a non-empty truncation range", {
  dat <- input_data()
  dat$cl <- factor(dat$y > median(dat$y))
  fit_cl <- function(...) rpf_050(cl ~ x1 + x2, data = dat, n_families = 2, n_splits = 5, ...)

  for (loss in c("logit", "exponential")) {
    for (eps in c(0, 0.5, 0.7, 1)) {
      expect_error(fit_cl(loss = loss, epsilon = eps), "epsilon. must be greater than 0 and less than 0.5")
    }
    expect_s3_class(fit_cl(loss = loss, epsilon = 0.49), "rpf")
    expect_error(fit_cl(loss = loss, delta = 0.5), "delta. must be")
    expect_s3_class(fit_cl(loss = loss, delta = 0.49), "rpf")
  }
  expect_error(fit_cl(loss = "logit", delta = 0), "delta. must be greater than 0")
  expect_s3_class(fit_cl(loss = "exponential", delta = 0), "rpf")

  # unused by the other losses, which keep accepting [0, 1]
  expect_s3_class(fit_cl(loss = "L2", epsilon = 0, delta = 0.7), "rpf")
  expect_s3_class(rpf_050(y ~ x1 + x2, data = dat, n_families = 2, epsilon = 0), "rpf")
})

test_that("share_candidates must be positive", {
  dat <- input_data()
  expect_error(
    rpf_050(y ~ x1 + x2, data = dat, n_families = 2, share_candidates = 0),
    "share_candidates. must be greater than 0"
  )
  expect_error(rpf_050(y ~ x1 + x2, data = dat, n_families = 2, share_candidates = -0.1))
  expect_s3_class(rpf_050(y ~ x1 + x2, data = dat, n_families = 2, share_candidates = 1e-9), "rpf")
})

test_that("unused levels of a factor outcome are dropped, also for predict()", {
  dat <- input_data()
  dat$cl <- factor(ifelse(dat$x1 > 0.5, "c", "a"), levels = c("a", "b", "c"))
  dropped <- transform(dat, cl = droplevels(cl))

  for (loss in c("exponential", "logit", "L2")) {
    set.seed(1)
    expect_warning(
      fit <- rpf_050(cl ~ x1 + x2, data = dat, n_families = 2, n_splits = 10, loss = loss),
      "Dropped 1 unused level of the outcome"
    )
    set.seed(1)
    ref <- rpf_050(cl ~ x1 + x2, data = dropped, n_families = 2, n_splits = 10, loss = loss)
    expect_identical(levels(fit$blueprint$ptypes$outcomes[[1]]), c("a", "c"))
    expect_identical(predict(fit, dat), predict(ref, dat))
    expect_identical(predict(fit, dat, type = "class"), predict(ref, dat, type = "class"))
  }

  # multiclass, and the x/y interface
  dat$cl3 <- factor(cut(dat$x1, c(0, 0.3, 0.6, 1), labels = c("a", "c", "d")), levels = c("a", "b", "c", "d"))
  set.seed(1)
  expect_warning(fit <- rpf_050(x = dat[c("x1", "x2")], y = dat$cl3, n_families = 2, n_splits = 10), "unused level")
  set.seed(1)
  ref <- rpf_050(x = dat[c("x1", "x2")], y = droplevels(dat$cl3), n_families = 2, n_splits = 10)
  expect_named(predict(fit, dat), c(".pred_a", ".pred_c", ".pred_d"))
  expect_identical(predict(fit, dat), predict(ref, dat))
})

test_that("the outcome's type and number of classes are checked with clear messages", {
  dat <- input_data()
  x <- dat[c("x1", "x2")]
  chr <- sample(c("a", "b"), nrow(dat), TRUE)
  lgl <- dat$x1 > 0.5
  msg <- "must be numeric \\(regression\\) or a factor \\(classification\\)"

  expect_error(rpf_050(x = x, y = chr, n_families = 2), msg)
  expect_error(rpf_050(x = as.matrix(x), y = lgl, n_families = 2), msg)
  expect_error(rpf_050(y ~ x1 + x2, data = transform(dat, y = chr), n_families = 2), msg)
  expect_error(rpf_050(y ~ x1 + x2, data = transform(dat, y = lgl), n_families = 2), msg)

  one <- factor(rep("a", nrow(dat)))
  expect_error(rpf_050(x = x, y = one, n_families = 2), "needs at least two classes")
  expect_error(rpf_050(y ~ x1 + x2, data = transform(dat, y = one), n_families = 2), "needs at least two classes")
  # one class in the data, the other one unused
  expect_error(rpf_050(x = x, y = factor(one, levels = c("a", "b")), n_families = 2), "Only \"a\" occurs")
})

test_that("predictor names must be present, unique and free of ':'", {
  dat <- input_data()
  x <- dat[c("x1", "x2")]

  unnamed <- as.matrix(x)
  colnames(unnamed) <- NULL
  expect_error(rpf_050(x = unnamed, y = dat$y, n_families = 2), "Every predictor must have a name")
  partly <- x
  names(partly) <- c("x1", "")
  expect_error(rpf_050(x = partly, y = dat$y, n_families = 2), "Column 2 has no name")
  duplicated_names <- x
  names(duplicated_names) <- c("a", "a")
  expect_error(rpf_050(x = duplicated_names, y = dat$y, n_families = 2), "Predictor names must be unique")
  with_colon <- x
  names(with_colon) <- c("x1", "a:b")
  expect_error(rpf_050(x = with_colon, y = dat$y, n_families = 2), "must not contain")
  colon_data <- dat
  colon_data[["a:b"]] <- runif(nrow(dat))
  expect_error(rpf_050(y ~ ., data = colon_data, n_families = 2), "must not contain")
})

test_that("formula interaction terms are an error that points to main terms", {
  dat <- input_data()
  dat$x3 <- runif(nrow(dat))
  for (formula in list(y ~ x1 * x2, y ~ x1:x2 + x3, y ~ (x1 + x2 + x3)^2, y ~ .^2)) {
    expect_error(rpf_050(formula, data = dat, n_families = 2), "must not contain interaction terms")
  }
  expect_error(rpf_050(y ~ x1 * x2, data = dat, n_families = 2), "y ~ x1 \\+ x2")
  expect_error(rpf_050(y ~ x1 * x2, data = dat, n_families = 2), "add it to the data as a column with its own name")
  # transformed terms of one predictor stay allowed, and so does a product added as a column
  expect_s3_class(rpf_050(y ~ log(x1) + I(x2^2), data = dat, n_families = 2), "rpf")
  # the error points at I(): a product wrapped in I() is one predictor, not an interaction term
  expect_s3_class(rpf_050(y ~ x1 + x2 + I(x1 * x2), data = dat, n_families = 2), "rpf")
  dat$x1_x2 <- dat$x1 * dat$x2
  expect_s3_class(rpf_050(y ~ x1 + x2 + x1_x2, data = dat, n_families = 2), "rpf")
})

test_that("zero-row new_data gives a zero-row prediction", {
  dat <- input_data()
  set.seed(1)
  fit <- rpf_050(y ~ x1 + x2, data = dat, n_families = 2)
  pred <- predict(fit, dat[0, ])
  expect_identical(pred, predict(fit, dat)[0, ])
  expect_named(pred, ".pred")

  dat$cl <- factor(dat$y > median(dat$y))
  dat$cl3 <- cut(dat$y, quantile(dat$y, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE, labels = c("a", "b", "c"))
  for (outcome in c("cl", "cl3")) {
    for (loss in c("exponential", "logit", "L1", "L2")) {
      set.seed(1)
      fit <- rpf_050(stats::reformulate(c("x1", "x2"), outcome), data = dat, n_families = 2, n_splits = 5, loss = loss)
      for (type in c("prob", "class", "numeric", "link")) {
        expect_identical(predict(fit, dat[0, ], type = type), predict(fit, dat, type = type)[0, ])
      }
    }
  }
})
