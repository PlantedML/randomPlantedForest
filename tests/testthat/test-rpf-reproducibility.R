test_that("Setting seed in R works", {
  set.seed(13)
  rpf_fit1 <- rpf(mpg ~ wt + cyl, data = mtcars)
  pred1 <- predict(rpf_fit1, mtcars[1:5, ])

  set.seed(13)
  rpf_fit2 <- rpf(mpg ~ wt + cyl, data = mtcars)
  pred2 <- predict(rpf_fit2, mtcars[1:5, ])

  # No seed set, so should be different
  rpf_fit3 <- rpf(mpg ~ wt + cyl, data = mtcars)
  pred3 <- predict(rpf_fit3, mtcars[1:5, ])

  expect_equal(pred1, pred2)
  expect_failure(expect_equal(pred1, pred3))
})

test_that("Rcpp RNG/R RNG interference", {
  set.seed(1)
  r11 <- runif(1)
  r12 <- runif(1)

  set.seed(1)
  r21 <- runif(1)
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)
  r22 <- runif(1)

  # If this fails R is broken
  expect_equal(r11, r21)
  # If this fails Rcpp does not properly affect R's RNG
  expect_failure(expect_equal(r12, r22))
})


test_that("Seeding from R works with threading", {
  set.seed(13)
  rpf_fit1 <- rpf(mpg ~ wt + cyl, nthreads = 2, data = mtcars)
  pred1 <- predict(rpf_fit1, mtcars[1:5, ])

  set.seed(13)
  rpf_fit2 <- rpf(mpg ~ wt + cyl, nthreads = 2, data = mtcars)
  pred2 <- predict(rpf_fit2, mtcars[1:5, ])

  expect_equal(pred1, pred2)
})
