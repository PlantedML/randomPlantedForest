# Calling miscellaneous C++ functions
test_that("C-level functionality works", {
  train <-  mtcars[1:20, ]
  test <-  mtcars[21:32, ]
  
  rpfit <- rpf(mpg ~., data = train, max_interaction = 3)
  pred <- predict(rpfit, test)
  
  mse <- rpfit$fit$MSE(as.matrix(pred$.pred), as.matrix(test$mpg))
  expect_true(inherits(mse, "numeric"))
  expect_gt(mse, 0)
  
  mod <- rpfit$fit$get_model()
  expect_true(inherits(mod, "list"))
  expect_equal(length(mod), 50)
  
  # This should be capturable instead of printing to stdout
  params <- rpfit$fit$get_parameters()
  
  # could test if setting params work
  rpfit$fit$set_parameters("n_trees", 60)
  
  # Not sure what to test for here
  rpfit$fit$new_purify()
  
  # large output, also not capturable
  # rpfit$fit$print()
})
