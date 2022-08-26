xdat <- data.frame(
  yint = sample(c(0L, 1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("hi", "mid", "lo"), 100, replace = TRUE)),
  ychar = sample(c("hi", "mid", "lo"), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = cut(runif(100), 3, labels = 1:3),
  x4 = cut(runif(100), 2, labels = 1:2)
)

test_that("Default: L2 with 'prob'", {
  classif_fit <- rpf(yfact ~ ., data = xdat)

  expect_identical(classif_fit$loss, "L2")

  classif_pred <- predict(classif_fit, xdat)
  classif_pred_prob <- predict(classif_fit, new_data = xdat, type = "prob")

  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(classif_pred), 0)
  expect_lte(min(classif_pred), 1)
  expect_identical(classif_pred, classif_pred_prob)
})

# L2 loss -----------------------------------------------------------------
test_that("L2: Probability prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "L2")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "prob")

  expect_identical(classif_fit$loss, "L2")
  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(classif_pred), 0)
  expect_lte(max(classif_pred), 1)
})

test_that("L2: Class prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "L2")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "class")

  expect_identical(classif_fit$loss, "L2")
  expect_equal(dim(classif_pred), c(nrow(xdat), 1))
  expect_true(all(classif_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("L2: Numeric prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "L2")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "numeric")

  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gt(max(classif_pred), 1)
  expect_lt(min(classif_pred), 0)
})

# L1 loss -----------------------------------------------------------------
test_that("L1: Probability prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "L1")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "prob")

  expect_identical(classif_fit$loss, "L1")
  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(classif_pred), 0)
  expect_lte(max(classif_pred), 1)
})

test_that("L1: Class prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "L1")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "class")

  expect_identical(classif_fit$loss, "L1")
  expect_equal(dim(classif_pred), c(nrow(xdat), 1))
  expect_true(all(classif_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("L1: Numeric prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "L1")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "numeric")

  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  # Hard to test for the possibility of pred > 1 | pred < 0 w/o fixed data
  # expect_lte(min(classif_pred), 0)
  # expect_gte(max(classif_pred), 1)
})

# logit loss --------------------------------------------------------------
test_that("logit: Probability prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "prob")

  # FIXME: Should this be stored as logit_2 or logit?
  expect_identical(classif_fit$loss, "logit_2")
  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(classif_pred), 0)
  expect_lte(max(classif_pred), 1)
})

test_that("logit: Class prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "class")

  expect_equal(dim(classif_pred), c(nrow(xdat), 1))
  expect_true(all(classif_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("logit: Numeric/link prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")

  classif_pred <- predict(classif_fit, new_data = xdat, type = "numeric")
  classif_pred_lnk <- predict(classif_fit, new_data = xdat, type = "link")

  expect_identical(classif_pred, classif_pred_lnk)
  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
})

# exponential loss --------------------------------------------------------
test_that("exponential: Probability prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "exponential")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "prob")

  expect_identical(classif_fit$loss, "exponential_2")
  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(classif_pred), 0)
  expect_lte(max(classif_pred), 1)
})

test_that("exponential: Class prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "exponential")
  classif_pred <- predict(classif_fit, new_data = xdat, type = "class")

  expect_equal(dim(classif_pred), c(nrow(xdat), 1))
  expect_true(all(classif_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("exponential: Numeric/link prediction", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "exponential")

  classif_pred <- predict(classif_fit, new_data = xdat, type = "numeric")
  classif_pred_lnk <- predict(classif_fit, new_data = xdat, type = "link")

  expect_identical(classif_pred, classif_pred_lnk)
  expect_equal(dim(classif_pred), c(nrow(xdat), nlevels(xdat$yfact)))
})

# Classif and prob agree with each other ----------------------------------

test_that("prob and classif yield same result", {
  classif_fit <- rpf(yfact ~ ., data = xdat, loss = "logit")
  pred_class <- predict(classif_fit, new_data = xdat, type = "class")
  pred_prob <- predict(classif_fit, new_data = xdat, type = "prob")

  pred_max <- apply(pred_prob, 1, which.max)
  pred_prob$pred_prob_class <- factor(pred_max, labels = levels(xdat$yfact))

  pred_both <- cbind(pred_prob, pred_class)

  expect_equal(pred_both$pred_prob_class, pred_both$.pred_class)
})
