# Fitting -----------------------------------------------------------------
test_that("Basic fitting", {
  # set parameters ------------------------
  n_splits <- 15
  max_inter <- 2
  n_trees <- 50
  split_try <- 10
  t_try <- 0.5
  deterministic_forest <- TRUE
  parallel <- TRUE
  purify_forest <- FALSE
  loss <- "logit"
  delta <- 0.1
  epsilon <- 0

  # train models ------------------------
  rpf_fit <- rpf(
    y_train, x_train,
    max_interaction = max_inter, t_try = t_try,
    ntrees = n_trees, splits = n_splits, split_try = split_try,
    deterministic = deterministic_forest, parallel = parallel
  )

  expect_s4_class(rpf_fit, "Rcpp_RandomPlantedForest")

  y_train_bin <- as.integer(y_train > 0)

  rpf_fit_binary <- rpf(
    y_train_bin, x_train,
    max_interaction = max_inter, t_try = t_try,
    ntrees = n_trees, splits = n_splits, split_try = split_try,
    deterministic = deterministic_forest, parallel = parallel
  )

  expect_s4_class(rpf_fit_binary, "Rcpp_RandomPlantedForest")
})


# Prediction --------------------------------------------------------------
test_that("Predicting works", {
  # set parameters -------------------------
  n_splits <- 15
  max_inter <- 2
  n_trees <- 50
  split_try <- 10
  t_try <- 0.5
  deterministic_forest <- TRUE
  parallel <- TRUE
  purify_forest <- FALSE
  loss <- "logit"
  delta <- 0.1
  epsilon <- 0

  # train models ----------------------------
  rpf_fit <- rpf(
    y_train, x_train,
    max_interaction = max_inter, t_try = t_try,
    ntrees = n_trees, splits = n_splits, split_try = split_try,
    deterministic = deterministic_forest, parallel = parallel
  )

  # predict ---------------------------------
  pred <- predict(rpf_fit, x_test[, 1:2], c(1,2))
  expect_equal(length(pred), nrow(x_test))

  # Vector input
  pred <- predict(rpf_fit, x_test[, 1], c(1))
  expect_equal(length(pred), length(x_test[, 1]))

  pred <- predict(rpf_fit, x_test[1, ], c(0))
  expect_length(pred, 1)

  pred <- predict(rpf_fit, x_test, c(0))
  expect_equal(length(pred), nrow(x_test))
})
