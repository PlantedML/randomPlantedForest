# Calling miscellaneous C++ functions
test_that("C-level functionality works", {
  train <- mtcars[1:20, ]
  test <- mtcars[21:32, ]

  ntrees <- 35
  set.seed(23)
  rpfit <- rpf(mpg ~ ., data = train, max_interaction = 3, ntrees = ntrees)
  pred <- predict(rpfit, test)

  mse <- rpfit$fit$MSE(as.matrix(pred$.pred), as.matrix(test$mpg))
  checkmate::expect_number(mse, lower = 5, upper = 30)
  # FIXME: Even with seed, results change between platforms. Yikes.
  # MSE here is ~14 on non-macOS platforms.
  # Difference probably compiler related, so can't reasonably expect fixed values I think.
  # if (Sys.info()[["sysname"]] == "Darwin") expect_equal(mse, 17.1905463)

  mod <- rpfit$fit$get_model()
  expect_true(inherits(mod, "list"))
  expect_equal(length(mod), ntrees) # length is number of trees

  # FIXME: Cosmetic issue but should be addressed at some point
  # This should be assignable instead of printing to stdout
  # An R-vector or list would probably be ideal
  params_captured <- capture.output(rpfit$fit$get_parameters())
  # params_return <- rpfit$fit$get_parameters()
  # expect_equal(params_captured, params_return)

  # could test if setting params work
  rpfit$fit$set_parameters("n_trees", 60)
  params_captured <- capture.output(rpfit$fit$get_parameters())
  expect_equal(
    strsplit(params_captured, "(, )|(: )")[[1]][[2]],
    "n_trees=60"
  )
})

test_that("set_data loads without fitting and fit() matches the fitting constructor", {
  X <- as.matrix(mtcars[, c("cyl", "wt", "hp")])
  Y <- as.matrix(mtcars$mpg)
  pars <- c(2, 10, 20, 10, 0.4, 0, 0, 1, 0, 0.1, 50, 1, 3)

  set.seed(1)
  direct <- new(RandomPlantedForest, Y, X, pars)

  manual <- new(RandomPlantedForest, pars)
  manual$set_data(Y, X)
  expect_length(manual$get_model(), 0)
  set.seed(1)
  manual$fit()

  expect_identical(manual$predict_matrix(X, 0, 1L), direct$predict_matrix(X, 0, 1L))
})

test_that("unknown classification loss errors", {
  pars <- c(1, 50, 30, 10, 0.4, 0, 0, 1, 0, 0.1, 50, 1, 3, 0, 0.1)
  expect_error(new(ClassificationRPF, "nope", pars), "Unknown loss")
})
