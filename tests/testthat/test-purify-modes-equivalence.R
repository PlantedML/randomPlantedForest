set.seed(2025)

test_that("single component predictions match across purify modes (non-capped)", {
  rp1 <- rpf(mpg ~ cyl + disp + hp + wt, data = mtcars, max_interaction = 3, ntrees = 30, deterministic = TRUE)
  rp2 <- rpf(mpg ~ cyl + disp + hp + wt, data = mtcars, max_interaction = 3, ntrees = 30, deterministic = TRUE)

  expect_false(is_purified(rp1))
  expect_false(is_purified(rp2))

  purify(rp1, mode = 1L)
  purify(rp2, mode = 2L)

  expect_true(is_purified(rp1))
  expect_true(is_purified(rp2))

  m1 <- predict_components(rp1, mtcars)
  m2 <- predict_components(rp2, mtcars)

  expect_equal(colnames(m1$m), colnames(m2$m))
  expect_equal(as.matrix(m1$m), as.matrix(m2$m), tolerance = 1e-8)
  expect_equal(m1$intercept, m2$intercept, tolerance = 1e-10)
})

test_that("single component predictions match across purify modes (capped)", {
  rp1 <- rpf(mpg ~ cyl + disp + hp + wt, data = mtcars, max_interaction = 3, ntrees = 30, deterministic = TRUE)
  rp2 <- rpf(mpg ~ cyl + disp + hp + wt, data = mtcars, max_interaction = 3, ntrees = 30, deterministic = TRUE)

  expect_false(is_purified(rp1))
  expect_false(is_purified(rp2))

  purify(rp1, maxp_interaction = 2L, mode = 1L)
  purify(rp2, maxp_interaction = 2L, mode = 2L)

  expect_true(is_purified(rp1))
  expect_true(is_purified(rp2))

  m1 <- predict_components(rp1, mtcars, max_interaction = 2L)
  m2 <- predict_components(rp2, mtcars, max_interaction = 2L)

  expect_equal(colnames(m1$m), colnames(m2$m))
  expect_equal(as.matrix(m1$m), as.matrix(m2$m), tolerance = 1e-8)
  expect_equal(m1$intercept, m2$intercept, tolerance = 1e-10)
})
