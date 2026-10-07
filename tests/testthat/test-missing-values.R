test_that("missing values in training data are rejected with the affected columns", {
  with_na_predictor <- mtcars
  with_na_predictor$wt[3] <- NA
  expect_error(rpf(mpg ~ ., data = with_na_predictor, ntrees = 2), "Missing values in `wt`")
  expect_error(
    rpf(x = with_na_predictor[, -1], y = with_na_predictor$mpg, ntrees = 2),
    "must not contain missing values"
  )

  with_na_outcome <- mtcars
  with_na_outcome$mpg[1] <- NA
  expect_error(rpf(mpg ~ ., data = with_na_outcome, ntrees = 2), "The outcome must not contain missing values")
})

test_that("missing values in new data are rejected instead of predicted", {
  fit <- rpf(mpg ~ cyl + wt, data = mtcars, ntrees = 2)
  new_data <- mtcars
  new_data$wt[3] <- NA

  expect_error(predict(fit, new_data), "`new_data` must not contain missing values")
  expect_error(predict_components(fit, new_data), "`new_data` must not contain missing values")
})

test_that("a recipe can impute missing values before fitting and predicting", {
  skip_if_not_installed("recipes")
  with_na <- mtcars
  with_na$wt[c(3, 10)] <- NA
  rec <- recipes::step_impute_median(recipes::recipe(mpg ~ cyl + wt, data = with_na), wt)

  fit <- rpf(rec, data = with_na, ntrees = 2)
  pred <- predict(fit, with_na)

  expect_false(anyNA(pred$.pred))
})
