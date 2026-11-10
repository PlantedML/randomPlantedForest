split_structures <- c("leaf", "leaf_binned", "tree_then_leaf", "tree_then_all_leaves", "target_tree_then_leaf")

test_that("every split structure learns a regression signal", {
  set.seed(1)
  n <- 200
  dat <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  dat$y <- 2 * dat$x1 + sin(3 * dat$x2) + rnorm(n, sd = 0.1)

  for (structure in split_structures) {
    set.seed(2)
    fit <- rpf_050(y ~ ., data = dat, max_interaction = 2, n_families = 5, n_splits = 20, split_structure = structure)
    pred <- predict(fit, dat)$.pred

    # the stored params keep the old code
    expect_identical(fit$params$split_structure, randomPlantedForest:::split_structure_code(structure))
    expect_gt(cor(pred, dat$y), 0.8, label = structure)
  }
})

test_that("every split structure learns a classification signal", {
  set.seed(1)
  n <- 200
  dat <- data.frame(x1 = runif(n), x2 = runif(n))
  dat$y <- factor(ifelse(dat$x1 + dat$x2 > 1, "a", "b"))

  for (structure in split_structures) {
    for (loss in c("L2", "exponential")) {
      set.seed(2)
      fit <- rpf_050(y ~ ., data = dat, n_families = 5, n_splits = 20, loss = loss, split_structure = structure)
      accuracy <- mean(predict(fit, dat, type = "class")$.pred_class == dat$y)

      expect_gt(accuracy, 0.8, label = paste(structure, loss))
    }
  }
})

test_that("every split structure is reproducible and independent of nthreads", {
  set.seed(1)
  dat <- data.frame(x1 = runif(100), x2 = runif(100))
  dat$y <- dat$x1 * dat$x2 + rnorm(100, sd = 0.1)

  for (structure in split_structures) {
    fit_with <- function(nthreads) {
      set.seed(3)
      fit <- rpf_050(
        y ~ .,
        data = dat,
        max_interaction = 2,
        n_families = 4,
        n_splits = 15,
        split_structure = structure,
        nthreads = nthreads
      )
      predict(fit, dat, nthreads = 1L)
    }

    expect_identical(fit_with(1L), fit_with(1L), label = structure)
    expect_identical(fit_with(1L), fit_with(2L), label = structure)
  }
})
