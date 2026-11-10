no_split_data <- function(n = 100) {
  set.seed(42)
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  x$y <- 10 + x$x1 + rnorm(n, sd = 0.1)
  x$cl <- factor(ifelse(runif(n) < 0.7, "a", "b"))
  x$cl3 <- factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.5, 0.3, 0.2)))
  x
}

test_that("a family without a split predicts the mean of its rows", {
  dat <- no_split_data()
  constant <- dat
  constant[c("x1", "x2")] <- 1
  for (structure in c("leaves", "hist", "cur_trees_1", "cur_trees_2", "res_trees")) {
    fit <- rpf(y ~ x1 + x2, data = constant, ntrees = 1, splits = 5, split_structure = structure, deterministic = TRUE)
    expect_equal(predict(fit, dat[1:3, ])$.pred, rep(mean(dat$y), 3), label = structure)
  }

  # a single row: every bootstrap sample is that row
  set.seed(1)
  fit <- rpf(y ~ x1 + x2, data = dat[1, ], ntrees = 5, splits = 5)
  expect_equal(predict(fit, dat[1:3, ])$.pred, rep(dat$y[1], 3))
})

test_that("a classification family without a split predicts the class shares", {
  dat <- no_split_data()
  dat[c("x1", "x2")] <- 1
  shares <- function(y) as.numeric(prop.table(table(y)))
  for (loss in c("exponential", "logit", "L1", "L2")) {
    fit <- rpf(cl ~ x1 + x2, data = dat, ntrees = 1, splits = 5, deterministic = TRUE, loss = loss)
    expect_equal(unlist(predict(fit, dat[1, ]), use.names = FALSE), shares(dat$cl), label = paste(loss, "binary"))
    fit <- rpf(cl3 ~ x1 + x2, data = dat, ntrees = 1, splits = 5, deterministic = TRUE, loss = loss)
    expect_equal(unlist(predict(fit, dat[1, ]), use.names = FALSE), shares(dat$cl3), label = paste(loss, "3 classes"))
  }

  # logit and exponential truncate the shares to [epsilon, 1 - epsilon], as in every leaf
  dat$rare <- factor(ifelse(seq_len(nrow(dat)) <= 97, "a", "b"))
  for (loss in c("exponential", "logit")) {
    fit <- rpf(rare ~ x1 + x2, data = dat, ntrees = 1, splits = 5, deterministic = TRUE, loss = loss, epsilon = 0.1)
    expect_equal(unlist(predict(fit, dat[1, ]), use.names = FALSE), c(0.9, 0.1), label = loss)
  }
})

test_that("in a forest where only some families split, the others predict their rows' mean", {
  # a bootstrap sample of these three rows has one distinct x in 9 of 27 cases: then no split is possible
  x <- data.frame(x = c(0, 0, 1))
  y <- c(1, 10, 100)
  set.seed(3)
  fit <- rpf(x, y, ntrees = 40, splits = 3, export_forest = TRUE)
  null_value <- vapply(
    fit$forest,
    function(family) {
      null_tree <- which(vapply(family$variables, function(v) identical(as.integer(v), 0L), logical(1)))
      family$values[[null_tree]][[1]]
    },
    numeric(1)
  )
  split <- vapply(fit$forest, function(family) length(family$variables) > 1, logical(1))
  expect_true(any(split) && any(!split))
  # families with a split keep their null leaf at 0; the others hold the mean of three rows of one x value
  expect_true(all(null_value[split] == 0))
  expect_true(all(null_value[!split] %in% c(1, 4, 7, 10, 100)))

  # purification puts the null leaves' values into the intercept, in both modes
  pred <- predict(fit, x)$.pred
  for (mode in 1:2) {
    set.seed(3)
    fit <- rpf(x, y, ntrees = 40, splits = 3)
    purify(fit, mode = mode)
    components <- predict_components(fit, x)
    expect_equal(predict(fit, x)$.pred, pred)
    expect_equal(rowSums(components$m) + components$intercept, pred)
  }
})

test_that("min_leaf_rows above half the rows leaves every family without a split: they predict their rows' mean", {
  dat <- no_split_data()
  # 100 rows: no threshold keeps 51 rows on each side, with either draw of candidates
  for (exact in c(TRUE, FALSE)) {
    old <- options(randomPlantedForest.exact_draw = exact)
    set.seed(1)
    fit <- rpf(y ~ x1 + x2, data = dat, ntrees = 5, splits = 10, min_leaf_rows = 51, export_forest = TRUE)
    options(old)
    expect_true(all(vapply(fit$forest, function(family) length(family$variables) == 1, logical(1))))
    pred <- predict(fit, dat[1:3, ])$.pred
    expect_equal(pred, rep(pred[1], 3))
    # the mean of five bootstrap means
    expect_lt(abs(pred[1] - mean(dat$y)), 0.1)
  }
  # without the bootstrap, the mean itself
  fit <- rpf(y ~ x1 + x2, data = dat, ntrees = 1, splits = 10, min_leaf_rows = 51, deterministic = TRUE)
  expect_equal(predict(fit, dat[1:3, ])$.pred, rep(mean(dat$y), 3))
})
