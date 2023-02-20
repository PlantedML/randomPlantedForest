set.seed(124)
xdat <- data.frame(
  y01 = sample(c(0L, 1L), 100, replace = TRUE),
  y12 = sample(c(1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("pos", "neg"), 100, replace = TRUE)),
  ychar = sample(c("pos", "neg"), 100, replace = TRUE),
  ylogi = sample(c(TRUE, FALSE), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = cut(runif(100), 3, labels = 1:3),
  x4 = cut(runif(100), 2, labels = 1:2)
)

test_that("Default: L2 with 'prob'", {
  bin_fit <- rpf(yfact ~ ., data = xdat)
  bin_pred <- predict(bin_fit, new_data = xdat)

  expect_identical(bin_fit$params$loss, "L2")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

# Sanity ----
test_that("Predictions are not constant: L1", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "L1")
  bin_pred <- predict(bin_fit, new_data = xdat)

  # More than one unique prediction
  expect_gt(nrow(unique(bin_pred)), 1)

  # unique predictions should be different than all 0 or 1
  expect_failure(expect_equal(
    sort(unique(bin_pred[[1]])), c(0, 1)
  ))
})

test_that("Predictions are not constant: L2", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "L2")
  bin_pred <- predict(bin_fit, new_data = xdat)

  # More than one unique prediction
  expect_gt(nrow(unique(bin_pred)), 1)

  # unique predictions should be different than all 0 or 1
  expect_failure(expect_equal(
    sort(unique(bin_pred[[1]])), c(0, 1)
  ))
})

test_that("Predictions are not constant: logit", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")
  bin_pred <- predict(bin_fit, new_data = xdat)

  # More than one unique prediction
  expect_gt(nrow(unique(bin_pred)), 1)

  # unique predictions should be different than all 0 or 1
  expect_failure(expect_equal(
    sort(unique(bin_pred[[1]])), c(0, 1)
  ))
})

test_that("Predictions are not constant: exponential", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "exponential")
  bin_pred <- predict(bin_fit, new_data = xdat)

  # More than one unique prediction
  expect_gt(nrow(unique(bin_pred)), 1)

  # unique predictions should be different than all 0 or 1
  expect_failure(expect_equal(
    sort(unique(bin_pred[[1]])), c(0, 1)
  ))
})

# L2 loss -----------------------------------------------------------------
test_that("L2: Probability prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "L2")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")

  expect_identical(bin_fit$params$loss, "L2")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

test_that("L2: Class prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "L2")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("L2: Numeric prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "L2")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_gt(max(bin_pred$.pred), 1)
  expect_lt(min(bin_pred$.pred), 0)
})

# L1 loss -----------------------------------------------------------------
test_that("L1: Probability prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "L1")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")

  expect_identical(bin_fit$params$loss, "L1")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

test_that("L1: Class prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "L1")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("L1: Numeric prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "L1")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  # FIXME: Not sure how to test the possibility of <0 | >1 here
})

# logit loss --------------------------------------------------------------
test_that("logit: Probability prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")

  expect_identical(bin_fit$params$loss, "logit")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

test_that("logit: Class prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("logit: Numeric/link prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")

  bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")
  bin_pred_lnk <- predict(bin_fit, new_data = xdat, type = "link")

  expect_identical(bin_pred, bin_pred_lnk)
  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
})

# exponential loss --------------------------------------------------------
test_that("exponential: Probability prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "exponential")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")

  expect_identical(bin_fit$params$loss, "exponential")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

test_that("exponential: Class prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "exponential")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("exponential: Numeric/link prediction", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "exponential")

  bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")
  bin_pred_lnk <- predict(bin_fit, new_data = xdat, type = "link")

  expect_identical(bin_pred, bin_pred_lnk)
  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
})


# Classif and prob agree with each other ----------------------------------

test_that("prob and classif yield same result", {
  bin_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")
  pred_class <- predict(bin_fit, new_data = xdat, type = "class")
  pred_prob <- predict(bin_fit, new_data = xdat, type = "prob")

  pred_max <- apply(pred_prob, 1, which.max)

  # Fail early if only one unique value is predicted
  expect_equal(length(unique(pred_max)), 2)

  pred_prob$pred_prob_class <- factor(pred_max, labels = levels(xdat$yfact))

  pred_both <- cbind(pred_prob, pred_class)

  expect_equal(pred_both$pred_prob_class, pred_both$.pred_class)
})

