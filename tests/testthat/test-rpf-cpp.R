# Calling miscellaneous C++ functions
test_that("C-level functionality works", {
  train <-  mtcars[1:20, ]
  test <-  mtcars[21:32, ]

  set.seed(23)
  rpfit <- rpf(mpg ~., data = train, max_interaction = 3)
  pred <- predict(rpfit, test)

  mse <- rpfit$fit$MSE(as.matrix(pred$.pred), as.matrix(test$mpg))
  expect_true(inherits(mse, "numeric"))
  # With fixed seed, result should not change over time
  expect_equal(mse, 17.1905463)

  mod <- rpfit$fit$get_model()
  expect_true(inherits(mod, "list"))
  expect_equal(length(mod), 50)

  # FIXME: Cosmetic issue but should be adressed at some point
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

  # Not sure what to test for here
  rpfit$fit$new_purify()

  # large output, also not assignable
  print_captured <- capture.output(rpfit$fit$print())
  expect_length(print_captured, 6864)
})
