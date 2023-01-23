set.seed(124)
xdat <- data.frame(
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = cut(runif(100), 3, labels = 1:3),
  x4 = cut(runif(100), 2, labels = 1:2),
  x5 = factor(sample(letters[1:4], 100, replace = TRUE))
)

xdat$yfact <- factor(rbinom(100, 1, prob = plogis(xdat$x1 + xdat$x2)), labels = c("neg", "pos"))

test_that("Purification works", {
  bin_fit <- rpf(yfact ~ x1 + x2 + x3 + x4 + x5, data = xdat, max_interaction = 3)

  expect_false(is_purified(bin_fit))
  expect_identical(is_purified(bin_fit), bin_fit$fit$is_purified())

  purify(bin_fit)

  expect_true(is_purified(bin_fit))
  expect_identical(is_purified(bin_fit), bin_fit$fit$is_purified())

  bin_pred <- predict(bin_fit, new_data = xdat)
  m <- extract_components(bin_fit, xdat)
  msum <- rowSums(m)

  # Purification seems to affect prediction results, but hard to gauge how much is considered acceptable.
  expect_lte(sqrt(mean((msum - bin_pred$.pred_pos)^2)), expected = 0.015)
  # expect_equal(msum, bin_pred$.pred_pos, tolerance = 1e-13)
})
