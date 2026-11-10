# rpf(splits = NULL) grows the tree families one at a time and rebuilds the chosen forest from the records of its
# splits: it must be the forest of rpf() with the chosen splits and the same seed, and its path the first checkpoints
# of rpf(oob_path = ) grown to where it stopped, bitwise, whatever the threads, the pause budget or the split structure

alone_data <- function() {
  set.seed(157)
  n <- 300
  x <- matrix(stats::runif(n * 4), n, 4, dimnames = list(NULL, paste0("x", 1:4)))
  y <- sin(3 * x[, 1]) + x[, 2] * x[, 3] + (x[, 4] > 0.5) + stats::rnorm(n, sd = 0.3)
  list(x = x, y = y)
}

alone_fit <- function(d, splits = NULL, ..., pause_gb = NULL) {
  old <- options(randomPlantedForest.adaptive_pause_gb = pause_gb)
  on.exit(options(old))
  args <- utils::modifyList(
    list(max_interaction = 3, ntrees = 10, split_try = 5, t_try = 0.5, prune = FALSE),
    list(...)
  )
  set.seed(11)
  do.call(rpf, c(list(d$x, d$y, splits = splits), args))
}

# The adaptive fit against the fixed fit at its choice and the path of the fixed fit grown to where it stopped
expect_alone_equal <- function(d, ...) {
  a <- alone_fit(d, ...)
  k <- a$params$splits
  expect_identical(a$splits_oob$splits, k)
  expect_identical(k, a$oob_path$splits[which.min(a$oob_path$oob_rmse)])
  f <- alone_fit(d, splits = k, ...)
  expect_identical(a$fit$get_model(), f$fit$get_model())
  expect_identical(predict(a, d$x), predict(f, d$x))
  expect_identical(a$fit$inbag(), f$fit$inbag())
  p <- alone_fit(d, splits = a$splits_oob$stopped_at, oob_path = a$oob_path$splits, ...)
  expect_identical(a$oob_path, p$oob_path)
  invisible(a)
}

test_that("splits = NULL gives the fixed fit at the chosen splits and its path, every split structure", {
  d <- alone_data()
  for (s in c("leaves", "hist", "cur_trees_1", "cur_trees_2", "res_trees")) {
    expect_alone_equal(d, split_structure = s)
  }
  expect_alone_equal(d, delete_leaves = FALSE)
  expect_alone_equal(d, max_interaction = 1)
})

test_that("waiting families and families grown again give the same fit", {
  d <- alone_data()
  # a large budget: families wait at the estimated stop and grow on; none: they finish there and are grown again
  for (s in c("leaves", "cur_trees_1")) {
    big <- expect_alone_equal(d, split_structure = s, pause_gb = 10)
    expect_gt(big$fit$adaptive_stats()$pauses, 0)
    none <- expect_alone_equal(d, split_structure = s, pause_gb = 0)
    st <- none$fit$adaptive_stats()
    expect_identical(st$pauses, 0)
    expect_identical(none$fit$get_model(), big$fit$get_model())
  }
  # the first families stop early, so the later ones that waited or finished at that estimate grow further
  expect_gt(alone_fit(d, split_structure = "cur_trees_1", pause_gb = 0)$fit$adaptive_stats()$regrowths, 0)
})

test_that("the adaptive fit does not depend on the number of threads", {
  skip_on_cran()
  d <- alone_data()
  for (s in c("leaves", "cur_trees_2")) {
    a1 <- alone_fit(d, split_structure = s, nthreads = 1)
    a3 <- expect_alone_equal(d, split_structure = s, nthreads = 3)
    expect_identical(a3$fit$get_model(), a1$fit$get_model())
    expect_identical(a3$oob_path, a1$oob_path)
    expect_identical(a3$splits_oob, a1$splits_oob)
    expect_identical(alone_fit(d, split_structure = s, nthreads = 3, pause_gb = 0)$fit$get_model(), a1$fit$get_model())
  }
})
