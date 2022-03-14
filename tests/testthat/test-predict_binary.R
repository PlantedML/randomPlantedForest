xdat <- data.frame(
  y01 = sample(c(0L, 1L), 100, replace = TRUE),
  y12 = sample(c(1L, 2L), 100, replace = TRUE),
  yfact = factor(sample(c("pos", "neg"), 100, replace = TRUE)),
  ychar = sample(c("pos", "neg"), 100, replace = TRUE),
  ylogi = sample(c(TRUE, FALSE), 100, replace = TRUE),
  x1 = rnorm(100),
  x2 = rnorm(100)
)

test_that("Probability prediction: y = factor / L2 loss", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L2")
  
  bin_pred <- predict(bin_fit, new_data = xdat, type = "prob")
  expect_equal(dim(bin_pred), c(nrow(xdat), nlevels(xdat$yfact)))
  expect_gte(min(bin_pred), 0)
  expect_lte(max(bin_pred), 1)
})

test_that("Class prediction: y = factor / L2 loss", {
  bin_fit <- rpf(yfact ~ x1 + x2, data = xdat, loss = "L2")
  
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")
  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% levels(xdat$yfact)))
})

test_that("Class prediction: y = 0,1 / L2 loss", {
  bin_fit <- suppressWarnings(rpf(y01 ~ x1 + x2, data = xdat, loss = "L2"))
  
  bin_pred <- predict(bin_fit, new_data = xdat, type = "class")
  expect_equal(dim(bin_pred), c(nrow(xdat), 1))
  expect_true(all(bin_pred$.pred_class %in% c("0", "1")))
})
