# Sum the values of all boxes containing each row, using the model's encoding
predict_from_boxes <- function(fit, boxes, data) {
  X <- preprocess_predictors_predict(fit, hardhat::forge(data, fit$blueprint)$predictors)
  inside <- matrix(TRUE, nrow(X), nrow(boxes))
  for (variable in colnames(X)) {
    lower <- boxes[[paste0(variable, "_lower")]]
    upper <- boxes[[paste0(variable, "_upper")]]
    used <- !is.na(lower)
    inside[, used] <- inside[, used] & outer(X[, variable], lower[used], `>=`) & outer(X[, variable], upper[used], `<`)
  }
  drop(inside %*% boxes$value)
}

test_that("boxes reproduce predictions before and after purification", {
  set.seed(1)
  dat <- data.frame(x1 = runif(100), x2 = runif(100), g = factor(sample(c("a", "b", "c"), 100, TRUE)))
  dat$y <- dat$x1 * dat$x2 + (dat$g == "b") + rnorm(100, sd = 0.1)
  fit <- rpf(y ~ ., data = dat, ntrees = 1, splits = 20)

  raw <- rpf_boxes(fit)
  expect_equal(predict_from_boxes(fit, raw, dat), predict(fit, dat)$.pred)

  purify(fit)
  purified <- rpf_boxes(fit)
  expect_equal(predict_from_boxes(fit, purified, dat), predict(fit, dat)$.pred)
  expect_gt(nrow(purified), nrow(raw))
})

test_that("boxes are named and bounded by the variables of their tree", {
  set.seed(1)
  fit <- rpf(mpg ~ wt + hp + cyl, data = mtcars, ntrees = 2, splits = 10)
  boxes <- rpf_boxes(fit, family = 2L)

  expect_named(
    boxes,
    c("tree", "order", "box", "value", paste0(rep(c("wt", "hp", "cyl"), each = 2), c("_lower", "_upper")))
  )
  expect_identical(boxes$tree[boxes$order == 0], "(Intercept)")
  pairs <- boxes[boxes$tree == "hp:wt", ]
  expect_true(all(!is.na(pairs$wt_lower) & !is.na(pairs$hp_upper) & is.na(pairs$cyl_lower)))
  expect_error(rpf_boxes(fit, family = 3L))
})

test_that("tree names follow the model's columns when factors are reordered", {
  set.seed(1)
  dat <- data.frame(y = mtcars$mpg, f = factor(mtcars$gear), x = mtcars$wt)
  fit <- rpf(y ~ f + x, data = dat, ntrees = 1, splits = 10)
  boxes <- rpf_boxes(fit)

  x_boxes <- boxes[boxes$tree == "x", ]
  expect_gt(nrow(x_boxes), 0)
  expect_true(all(x_boxes$x_lower >= min(mtcars$wt)))
  expect_true(all(is.na(x_boxes$f_lower)))
})

test_that("multiclass boxes have one value column per class", {
  set.seed(1)
  fit <- rpf(Species ~ ., data = iris, ntrees = 1, splits = 5)
  boxes <- rpf_boxes(fit)

  expect_contains(names(boxes), paste0("value_", levels(iris$Species)))
  expect_false("value" %in% names(boxes))
})

test_that("box plots draw one panel per tree up to max_interaction", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  set.seed(1)
  dat <- data.frame(x1 = runif(50), x2 = runif(50), x3 = runif(50), g = factor(sample(c("a", "b"), 50, TRUE)))
  dat$y <- dat$x1 * dat$x2 * dat$x3 + rnorm(50, sd = 0.1)
  fit <- rpf(y ~ ., data = dat, ntrees = 1, splits = 30, max_interaction = 3)
  boxes <- rpf_boxes(fit)
  n_panels <- length(unique(boxes$tree[boxes$order %in% 1:2]))

  p <- plot(fit, type = "boxes")
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, n_panels - 1)
  expect_no_error(ggplot2::ggplot_build(p[[1]]))
  expect_match(p$patches$annotation$caption, "Remainder: \\d+ trees? on 3 variables")
  expect_error(plot(fit, type = "boxes", max_interaction = 3))

  some <- unique(boxes$tree[boxes$order %in% 1:2])[1:2]
  p_some <- plot(fit, type = "boxes", trees = some)
  expect_length(p_some$patches$plots, 1)
  expect_match(p_some$patches$annotation$caption, "Remainder")
  expect_error(plot(fit, type = "boxes", trees = "x9"))

  p_main <- plot(fit, type = "boxes", max_interaction = 1)
  expect_length(p_main$patches$plots, length(unique(boxes$tree[boxes$order == 1])) - 1)

  purify(fit)
  expect_s3_class(plot(fit, type = "boxes", factor_levels = FALSE), "patchwork")
})

test_that("plot() selects the class for multiclass fits only", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  set.seed(1)
  fit_multi <- rpf(Species ~ Sepal.Length + Petal.Length, data = iris, ntrees = 1, splits = 5)
  fit_reg <- rpf(mpg ~ wt + hp, data = mtcars, ntrees = 1, splits = 5)

  expect_s3_class(plot(fit_multi, type = "boxes", class = "virginica"), "patchwork")
  expect_s3_class(plot(fit_multi, class = "virginica"), "ggplot")
  expect_error(plot(fit_multi, class = "rose"))
  expect_error(plot(fit_reg, class = "virginica"), "multiclass")
})

# Leaves of a reconstructed split tree, checking along the way that the
# children of every split exactly tile the split's region
expect_tiles <- function(node, variables) {
  if (node$kind == "leaf") {
    return(node$leaf$.id)
  }
  ids <- lapply(node$children, expect_tiles, variables = variables)
  if (node$kind == "split") {
    low <- node$children[[1]]$region[[node$variable]]
    high <- node$children[[2]]$region[[node$variable]]
    expect_identical(c(low[1], high[2]), node$region[[node$variable]])
    expect_identical(low[2], high[1])
  }
  unlist(ids)
}

test_that("split reconstruction uses every leaf once and tiles regions", {
  set.seed(1)
  dat <- data.frame(x1 = runif(200), x2 = runif(200), g = factor(sample(c("a", "b", "c"), 200, TRUE)))
  dat$y <- dat$x1 + 2 * dat$x1 * dat$x2 + (dat$g == "b") + rnorm(200, sd = 0.1)
  fit <- rpf(y ~ ., data = dat, ntrees = 1, splits = 30)
  boxes <- rpf_boxes(fit)
  boxes$.value <- boxes$value

  for (tree in unique(boxes$tree[boxes$order > 0])) {
    leaves <- boxes[boxes$tree == tree, ]
    leaves$.id <- seq_len(nrow(leaves))
    variables <- strsplit(tree, ":", fixed = TRUE)[[1]]
    region <- lapply(variables, \(v) range(leaves[[paste0(v, "_lower")]], leaves[[paste0(v, "_upper")]]))
    names(region) <- variables

    ids <- expect_tiles(split_leaves(leaves, variables, region), variables)
    expect_setequal(ids, leaves$.id)
    expect_length(ids, nrow(leaves))
  }
})

test_that("repeated splits of a region become separate tilings", {
  # two splits of the root on x: at 0.5, and at 0.3 refined at 0.7
  leaves <- data.frame(
    x_lower = c(0, 0.5, 0, 0.3, 0.7),
    x_upper = c(0.5, 1, 0.3, 0.7, 1),
    .id = 1:5
  )
  node <- split_leaves(leaves, "x", list(x = c(0, 1)))

  expect_identical(node$kind, "overlap")
  expect_length(node$children, 2)
  expect_setequal(lapply(node$children, tiling_ids), list(1:2, 3:5))
})

test_that("tree plots need an unpurified forest", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  fit <- rpf(mpg ~ wt + hp + cyl, data = mtcars, ntrees = 1, splits = 10)

  p <- plot(fit)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  purify(fit)
  expect_error(plot(fit), "purified")
  expect_error(plot(fit, type = "branches"))
})

test_that("remainder summarizes the trees above max_interaction", {
  boxes <- data.frame(tree = c("(Intercept)", "a", "a:b", "a:b:c", "a:b:c", "a:b:d"), order = c(0, 1, 2, 3, 3, 3))

  expect_null(summarize_remainder(boxes[boxes$order > 3, ]))
  expect_identical(
    summarize_remainder(boxes[boxes$order > 2, ]),
    "Remainder: 2 trees on 3 variables with 3 leaves, not drawn"
  )
  expect_identical(
    summarize_remainder(boxes[boxes$order > 1, ]),
    "Remainder: 3 trees on 2-3 variables with 4 leaves, not drawn"
  )
})
