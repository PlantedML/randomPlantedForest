test_that("params-only constructors create empty, unfitted objects", {
  # 13-element parameter vector: max_interaction, ntrees, splits, split_try,
  # t_try, purify, deterministic, nthreads, cv, split_decay_rate,
  # max_candidates, delete_leaves, split_mode
  pars <- c(1, 50, 30, 10, 0.4, 0, 0, 1, 0, 0.1, 50, 1, 3)
  fit <- new(RandomPlantedForest, pars)
  expect_length(fit$get_model(), 0)
  expect_false(fit$is_purified())

  cfit <- new(ClassificationRPF, "L2", c(pars, 0, 0.1))
  expect_length(cfit$get_model(), 0)
})

test_that("get_data returns training data from a fitted forest", {
  X <- as.matrix(mtcars[, c("wt", "cyl")])
  Y <- as.matrix(mtcars$mpg)
  fit <- new(RandomPlantedForest, Y, X, c(1, 10, 30, 10, 0.4, 0, 0, 1, 0, 0.1, 50, 1, 3))
  dat <- fit$get_data()
  expect_equal(unname(as.matrix(dat$X)), unname(X))
  expect_equal(unname(as.matrix(dat$Y)), unname(Y))
})
