test_that("Formula interface", {
  rpf_fit <- rpf(mpg ~ wt + cyl, data = mtcars)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})

test_that("XY data.frame interface", {
  pred_df <- mtcars[, c(2, 6)]
  rpf_fit <- rpf(x = pred_df, y = mtcars$mpg)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})

test_that("XY matrix interface", {
  pred_mat <- as.matrix(mtcars[, c(2, 6)])
  rpf_fit <- rpf(x = pred_mat, y = mtcars$mpg)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})

test_that("recipe interface", {
  skip_if_not_installed("recipes")

  rec <- recipes::recipe(mpg ~ ., data = mtcars)
  rpf_fit <- rpf(rec, data = mtcars)

  expect_s3_class(rpf_fit, "rpf")
  expect_s4_class(rpf_fit$fit, "Rcpp_RandomPlantedForest")
})

test_that("Unsupported interface", {
  pred_tab <- as.table(as.matrix(mtcars[, c(2, 6)]))
  expect_error(rpf(x = pred_tab, y = mtcars$mpg))
})


# Parameter sets/combinations ---------------------------------------------

test_that("Setting max_interaction = 0 works", {
  set.seed(100)
  rpf_fit0 <- rpf(mpg ~ ., data = mtcars[1:20, ], max_interaction = 0)
  pred0 <- predict(rpf_fit0, new_data = mtcars[21:32, ])

  set.seed(100)
  rpf_fit10 <- rpf(mpg ~ ., data = mtcars[1:20, ], max_interaction = 10)
  pred10 <- predict(rpf_fit10, new_data = mtcars[21:32, ])

  set.seed(100)
  rpf_fit_default <- rpf(mpg ~ ., data = mtcars[1:20, ])
  pred_default <- predict(rpf_fit_default, new_data = mtcars[21:32, ])

  # Sanity: pred of default (max_interaction = 1) should be different
  expect_failure(expect_equal(pred_default, pred0))

  # prediction result should be identical if max_interaction = 0 or
  # max_interaction = p (= 10 for mtcars, ncol(mtcars) - 1)
  expect_equal(pred0, pred10)
})
