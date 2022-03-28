xdat <- data.frame(
  y01 = sample(c(0L, 1L), 100, replace = TRUE),
  y12 = sample(c(1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("pos", "neg"), 100, replace = TRUE)),
  ychar = sample(c("pos", "neg"), 100, replace = TRUE),
  ylogi = sample(c(TRUE, FALSE), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100)
)

test_that("Default: L2 with 'prob'", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat)
  bin_pred <- predict(bin_fit, new_data = xdat)

  expect_identical(bin_fit$loss, "L2")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

# L2 loss -----------------------------------------------------------------
test_that("L2: Probability prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L2")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")

  expect_identical(bin_fit$loss, "L2")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

test_that("L2: Class prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L2")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("L2: Numeric prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L2")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_gt(max(bin_pred$.pred), 1)
  expect_lt(min(bin_pred$.pred), 0)
})

# L1 loss -----------------------------------------------------------------
test_that("L1: Probability prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L1")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")

  expect_identical(bin_fit$loss, "L1")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

test_that("L1: Class prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L1")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("L1: Numeric prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L1")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  # FIXME: Not sure how to test the possibility of <0 | >1 here
  # expect_gt(max(bin_pred$.pred), 1)
  # expect_lt(min(bin_pred$.pred), 0)
})

# logit loss --------------------------------------------------------------
test_that("logit: Probability prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "logit")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")

  expect_identical(bin_fit$loss, "logit")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  # FIXME: Not sure how to test the possibility of <0 | >1 here
  # expect_gte(min(bin_pred), 0)
  # expect_lte(max(bin_pred), 1)
})

test_that("logit: Class prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "logit")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("logit: Numeric/link prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "logit")

  bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")
  bin_pred_lnk <- predict(bin_fit, new_data = xdat, type = "link")

  expect_identical(bin_pred, bin_pred_lnk)
  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
})

# exponential loss --------------------------------------------------------
test_that("exponential: Probability prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "exponential")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")

  expect_identical(bin_fit$loss, "exponential")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

test_that("exponential: Class prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "exponential")
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")

  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("exponential: Numeric/link prediction", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "exponential")

  bin_pred <- predict(bin_fit, new_data = xdat, type = "numeric")
  bin_pred_lnk <- predict(bin_fit, new_data = xdat, type = "link")

  expect_identical(bin_pred, bin_pred_lnk)
  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
})

