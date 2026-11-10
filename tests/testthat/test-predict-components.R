test_that("predict_components returns correct structure", {
  rp <- rpf(mpg ~ ., data = mtcars, max_interaction = 2)

  m <- predict_components(rp, mtcars)

  expect_named(
    m$m,
    c(
      "cyl",
      "disp",
      "hp",
      "drat",
      "wt",
      "qsec",
      "vs",
      "am",
      "gear",
      "carb",
      "cyl:disp",
      "cyl:hp",
      "cyl:drat",
      "cyl:wt",
      "cyl:qsec",
      "cyl:vs",
      "am:cyl",
      "cyl:gear",
      "carb:cyl",
      "disp:hp",
      "disp:drat",
      "disp:wt",
      "disp:qsec",
      "disp:vs",
      "am:disp",
      "disp:gear",
      "carb:disp",
      "drat:hp",
      "hp:wt",
      "hp:qsec",
      "hp:vs",
      "am:hp",
      "gear:hp",
      "carb:hp",
      "drat:wt",
      "drat:qsec",
      "drat:vs",
      "am:drat",
      "drat:gear",
      "carb:drat",
      "qsec:wt",
      "vs:wt",
      "am:wt",
      "gear:wt",
      "carb:wt",
      "qsec:vs",
      "am:qsec",
      "gear:qsec",
      "carb:qsec",
      "am:vs",
      "gear:vs",
      "carb:vs",
      "am:gear",
      "am:carb",
      "carb:gear"
    )
  )

  expect_s3_class(m$m, "data.frame")
  expect_equal(nrow(mtcars), nrow(m$m))
})

test_that("predict_components returns requested predictors, in order", {
  rp <- rpf(mpg ~ ., data = mtcars, max_interaction = 2)

  m <- predict_components(rp, mtcars, predictors = c("cyl", "am", "vs"))

  expect_named(m$m, c("cyl", "am", "vs", "am:cyl", "cyl:vs", "am:vs"))

  expect_s3_class(m$m, "data.frame")
  expect_length(m$m, 6)
  expect_equal(nrow(mtcars), nrow(m$m))
})

test_that("predict_components purifies if needed", {
  rp <- rpf(mpg ~ ., data = mtcars[, c(1:4)], max_interaction = 2)

  expect_false(is_purified(rp))

  predict_components(rp, mtcars)

  expect_true(is_purified(rp))
})


test_that("predict_components sums to prediction", {
  rp <- rpf(mpg ~ cyl + am + gear, data = mtcars, max_interaction = 3, purify = TRUE)

  m <- predict_components(rp, mtcars)

  expect_equal(
    predict(rp, mtcars)[[1]],
    rowSums(m$m) + m$intercept
  )
})

test_that("predict_components errors appropriately", {
  rp <- rpf(mpg ~ cyl + am + gear, data = mtcars, max_interaction = 3)

  # Unknown predictors requested
  expect_error(predict_components(rp, mtcars, predictors = c("x1", "x2")))

  # Predictors not in fit on
  expect_error(predict_components(rp, mtcars, predictors = c("cyl", "disp")))

  # new_data does not contain all original predictors
  expect_error(predict_components(rp, mtcars[, c("cyl", "am")]))

  # max_interaction higher than original fit
  expect_error(predict_components(rp, mtcars, max_interaction = 30))

  # max_interaction implausible
  expect_error(predict_components(rp, mtcars, max_interaction = -1))

  # predictors is wrong type
  expect_error(predict_components(rp, mtcars, predictors = list("cyl")))
})

test_that(".predict_single_component is consistent with predictor order", {
  rp <- rpf(mpg ~ cyl + am + gear, data = mtcars, max_interaction = 3, purify = TRUE)

  # Internal data preprocessing only done in predict_components to save time
  processed <- hardhat::forge(mtcars, rp$blueprint)
  new_data <- preprocess_predictors_predict(rp, processed$predictors)

  expect_equal(
    .predict_single_component(rp, new_data, c("cyl", "am")),
    .predict_single_component(rp, new_data, c("am", "cyl"))
  )

  expect_equal(
    Reduce(
      expect_equal,
      list(
        .predict_single_component(rp, new_data, c("cyl", "am", "gear")),
        .predict_single_component(rp, new_data, c("gear", "am", "cyl")),
        .predict_single_component(rp, new_data, c("am", "cyl", "gear")),
        .predict_single_component(rp, new_data, c("gear", "cyl", "am")),
        .predict_single_component(rp, new_data, c("am", "gear", "cyl")),
        .predict_single_component(rp, new_data, c("cyl", "gear", "am"))
      )
    ),
    .predict_single_component(rp, new_data, c("cyl", "gear", "am"))
  )
})

test_that("duplicated rows get the same components as their originals", {
  rp <- rpf(mpg ~ ., data = mtcars, max_interaction = 2, ntrees = 10)
  idx <- c(3, 1, 3, 2, 1, 1)

  once <- predict_components(rp, mtcars[1:3, ])
  repeated <- predict_components(rp, mtcars[idx, ])

  expect_equal(repeated$m, once$m[idx], ignore_attr = TRUE)
  expect_equal(repeated$intercept, once$intercept)
})

test_that("predict_components() takes a matrix and leaves a data.frame new_data unchanged", {
  set.seed(54)
  x <- as.matrix(mtcars[, 2:6])
  fit <- rpf(x, mtcars$mpg, max_interaction = 2, ntrees = 8, splits = 30)
  nd <- mtcars[1:12, 2:6]
  nd_before <- nd

  pc_df <- predict_components(fit, nd)
  # the caller's data.frame is neither converted to a data.table nor modified
  expect_identical(class(nd), "data.frame")
  expect_identical(nd, nd_before)

  pc_mat <- predict_components(fit, x[1:12, ])
  expect_identical(pc_mat$m, pc_df$m)
  expect_equal(pc_mat$x, pc_df$x, ignore_attr = TRUE)
  expect_equal(pc_mat$intercept, pc_df$intercept)

  # with a remainder (max_interaction below the fit's) and a subset of predictors
  pc_mat1 <- predict_components(fit, x[1:12, ], max_interaction = 1, predictors = c("hp", "cyl"))
  pc_df1 <- predict_components(fit, nd, max_interaction = 1, predictors = c("hp", "cyl"))
  expect_identical(pc_mat1$m, pc_df1$m)
  expect_identical(pc_mat1$remainder, pc_df1$remainder)
  expect_identical(names(pc_mat1$x), c("hp", "cyl"))
  expect_identical(class(nd), "data.frame")
})

# intercept + m + remainder minus predict(type = "numeric"), per class for multiclass
components_gap <- function(object, pc, new_data) {
  pred <- as.matrix(predict(object, new_data, type = "numeric"))
  m <- as.matrix(pc$m)
  rem <- if (is.null(pc$remainder)) 0 else as.matrix(pc$remainder)
  if (is.null(pc$target_levels)) {
    return(max(abs(pc$intercept + rowSums(m) + rem - pred[, 1])))
  }
  gaps <- vapply(
    seq_along(pc$target_levels),
    function(k) {
      cols <- endsWith(colnames(m), paste0("__class:", pc$target_levels[k]))
      rem_k <- if (is.null(pc$remainder)) 0 else rem[, k]
      max(abs(pc$intercept[[k]] + rowSums(m[, cols, drop = FALSE]) + rem_k - pred[, k]))
    },
    numeric(1)
  )
  max(gaps)
}

test_that("a subset of predictors returns the other components as remainder", {
  set.seed(7)
  n <- 150
  d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  lin <- d$x1 + 2 * d$x2 * d$x3
  d_reg <- cbind(d, y = lin + rnorm(n, sd = 0.2))
  d_bin <- cbind(d, y = factor(lin > stats::median(lin), labels = c("lo", "hi")))
  d_3 <- cbind(d, y = cut(lin, stats::quantile(lin, 0:3 / 3), include.lowest = TRUE, labels = c("a", "b", "c")))

  for (dd in list(d_reg, d_bin, d_3)) {
    set.seed(1)
    fit <- rpf(y ~ ., data = dd, max_interaction = 2, ntrees = 5, splits = 30)
    for (purified in c(FALSE, TRUE)) {
      # unpurified: predict_components() purifies in place; then on a forest purified with mode = 1
      if (purified) {
        purify(fit, mode = 1)
      }
      for (preds in list(c("x1", "x2"), "x3", c("x3", "x1"))) {
        # at the fit's order (the default), and below it
        for (k in list(NULL, 1L)) {
          pc <- predict_components(fit, dd[1:40, ], max_interaction = k, predictors = preds)
          expect_false(is.null(pc$remainder))
          expect_lt(components_gap(fit, pc, dd[1:40, ]), 1e-10)
        }
      }
      # all predictors at the fit's order: nothing left out, no remainder
      pc <- predict_components(fit, dd[1:40, ])
      expect_null(pc$remainder)
      expect_lt(components_gap(fit, pc, dd[1:40, ]), 1e-10)
    }
  }
})

test_that("predict_components() works with formula terms that are not raw columns", {
  set.seed(8)
  n <- 120
  d <- data.frame(x1 = runif(n, 0.1, 1), x2 = runif(n), x3 = runif(n))
  d$y <- log(d$x1) + d$x2^2 * d$x3 + rnorm(n, sd = 0.1)

  cases <- list(
    list(formula = y ~ log(x1) + x2 + x3, names = c("log(x1)", "x2", "x3")),
    list(formula = y ~ I(x2^2) + x1, names = c("I(x2^2)", "x1")),
    list(formula = y ~ . - x3, names = c("x1", "x2"))
  )
  for (case in cases) {
    set.seed(1)
    fit <- rpf(case$formula, data = d, max_interaction = 2, ntrees = 5, splits = 20)
    pc <- predict_components(fit, d[1:30, ], max_interaction = 1)
    expect_named(pc$m, case$names)
    expect_named(pc$x, case$names)
    expect_lt(components_gap(fit, pc, d[1:30, ]), 1e-10)

    pc <- predict_components(fit, d[1:30, ], predictors = case$names[1])
    expect_named(pc$m, case$names[1])
    expect_lt(components_gap(fit, pc, d[1:30, ]), 1e-10)
  }
  # x holds the term's values, as the forest sees them
  set.seed(1)
  fit <- rpf(y ~ log(x1) + x2, data = d, ntrees = 5, splits = 20)
  expect_equal(predict_components(fit, d[1:30, ])$x[["log(x1)"]], log(d$x1[1:30]))
})

test_that("components follow the fit's column order when factors come first", {
  set.seed(9)
  n <- 150
  d <- data.frame(f = factor(sample(c("u", "v", "w"), n, TRUE)), x1 = runif(n), x2 = runif(n))
  d$y <- as.integer(d$f) * d$x1 + d$x2 + rnorm(n, sd = 0.2)
  set.seed(1)
  fit <- rpf(y ~ ., data = d, max_interaction = 2, ntrees = 5, splits = 30)

  # hardhat::mold() puts the factor last
  pc <- predict_components(fit, d[1:40, ])
  expect_named(pc$m, c("x1", "x2", "f", "x1:x2", "f:x1", "f:x2"))
  expect_null(pc$remainder)
  expect_lt(components_gap(fit, pc, d[1:40, ]), 1e-10)
  expect_identical(pc$x$f, d$f[1:40])

  pc_f <- predict_components(fit, d[1:40, ], predictors = c("f", "x2"))
  expect_named(pc_f$m, c("f", "x2", "f:x2"))
  expect_equal(pc_f$m$f, pc$m$f)
  expect_equal(pc_f$m[["f:x2"]], pc$m[["f:x2"]])
  expect_lt(components_gap(fit, pc_f, d[1:40, ]), 1e-10)
})
