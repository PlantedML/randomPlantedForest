xdat <- data.frame(
  yint = sample(c(0L, 1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("hi", "mid", "lo"), 100, replace = TRUE)),
  ychar = sample(c("hi", "mid", "lo"), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = cut(runif(100), 3, labels = 1:3),
  x4 = cut(runif(100), 2, labels = 1:2)
)

# Basic model creation ----------------------------------------------------
test_that("Multiclass: All numeric", {
  classif_fit <- rpf(yfact ~ ., data = xdat)

  expect_s3_class(classif_fit, "rpf")
  expect_s4_class(classif_fit$fit, "Rcpp_ClassificationRPF")
})

# Classif task detection ---------------------------------------------------
test_that("Multiclass: Detection works", {
  # y 3-level factor
  y_fact <- rpf(yfact ~ ., xdat)
  expect_s4_class(y_fact$fit, "Rcpp_ClassificationRPF")

  # y is integer: should _not_ be treated as classif task
  y_int <- rpf(yint ~ ., xdat)
  expect_failure(expect_s4_class(y_int$fit, "Rcpp_ClassificationRPF"))

  # y 3-level character
  expect_error(rpf(ychar ~ x1 + x2, xdat), regexp = "must be numeric \\(regression\\) or a factor")
  expect_error(rpf(ychar ~ x3 + x4, xdat), regexp = "Ordering of factor columns only implemented")
})

# Multiclass logit learns signal (#40) -------------------------------------
test_that("Multiclass logit learns signal (#40)", {
  set.seed(42)
  n <- 400
  dat <- data.frame(x1 = runif(n), x2 = runif(n))
  dat$y <- factor(ifelse(dat$x1 > 0.5, "a", ifelse(dat$x2 > 0.5, "b", "c")))
  idx <- sample(n, 280)

  fit <- rpf(y ~ ., data = dat[idx, ], loss = "logit", ntrees = 10, max_interaction = 2)
  pred <- predict(fit, dat[-idx, ], type = "class")

  expect_gt(mean(pred$.pred_class == dat$y[-idx]), 0.85)
})

test_that("Remainder is calculcated correctly", {
  classif_fit <- rpf(yfact ~ ., data = xdat, max_interaction = 3)

  components <- predict_components(classif_fit, xdat, max_interaction = 2)

  expect_s3_class(components$m, "data.frame")
  expect_equal(nrow(components$m), nrow(components$remainder))
  expect_equal(ncol(components$remainder), length(components$target_levels))
  expect_named(components$remainder, components$target_levels)
})

test_that("Multiclass: components and class-specific intercept sum to prediction", {
  for (loss in c("L2", "logit", "exponential")) {
    fit <- rpf(yfact ~ ., data = xdat, max_interaction = 2, ntrees = 10, loss = loss)
    components <- predict_components(fit, xdat)
    pred <- predict(fit, xdat, type = "numeric")

    expect_named(components$intercept, components$target_levels)
    for (level in components$target_levels) {
      m_level <- components$m[, endsWith(names(components$m), paste0("__class:", level)), with = FALSE]
      expect_equal(
        rowSums(m_level) + components$intercept[[level]],
        pred[[paste0(".pred_", level)]],
        info = paste(loss, level)
      )
    }
  }
})
