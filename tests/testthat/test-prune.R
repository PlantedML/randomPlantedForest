# prune_components(): effect-level pruning of a regression forest

prune_sim <- function(n = 400, seed = 11) {
  set.seed(seed)
  x <- matrix(runif(n * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
  # a main effect and a pure pair (no main effects of x2, x3); x4 is noise
  y <- 1 + 2 * x[, 1] + 4 * x[, 2] * x[, 3] + rnorm(n, sd = 0.3)
  list(x = x, y = y)
}

prune_fit <- function(d, ...) {
  set.seed(5)
  rpf(d$x, d$y, max_interaction = 2, ntrees = 20, splits = 40, ...)
}

test_that("CV pruning keeps the true components and drops the false ones", {
  d <- prune_sim()
  fit <- prune_fit(d)
  fold_id <- rep_len(1:5, nrow(d$x))
  pr <- prune_components(fit, d$x, d$y, select_by = "cv", fold_id = fold_id, seed = 3, select = "min")

  expect_s3_class(pr, "rpf_pruned")
  expect_setequal(names(pr$beta), c("x1", "x2:x3"))
  expect_true(all(pr$components$kept == pr$components$component %in% c("x1", "x2:x3")))
  # the unpruned forest has false components
  expect_gt(nrow(pr$components), 2)
  expect_identical(pr$seed, 3L)
  expect_identical(pr$fold_id, fold_id)
  expect_identical(pr$select, "min")
  expect_identical(lapply(pr$sets, function(s) colnames(d$x)[s]), strsplit(names(pr$beta), ":"))

  # the CV table: every candidate, one choice of each kind
  expect_equal(nrow(pr$cv), 16)
  expect_equal(sum(pr$cv$min), 1)
  expect_equal(sum(pr$cv$one_se), 1)
  chosen <- pr$cv[pr$cv$min, ]
  expect_equal(c(chosen$rule, chosen$tau), c(pr$rule, pr$tau))
  expect_equal(chosen$cv_rmse, min(pr$cv$cv_rmse))
  expect_gte(pr$cv$tau[pr$cv$one_se], pr$tau)

  # the 1-SE rule (the default) picks the marked candidate; R's coefficients under CV are fitted in sample
  pr1 <- prune_components(fit, d$x, d$y, select_by = "cv", fold_id = fold_id, seed = 3)
  expect_identical(pr1$select, "1se")
  expect_identical(pr1$cv, pr$cv)
  one <- pr$cv[pr$cv$one_se, ]
  expect_equal(c(pr1$rule, pr1$tau), c(one$rule, one$tau))
  direct <- prune_components(fit, d$x, d$y, rule = pr1$rule, tau = pr1$tau, select_by = "cv")
  expect_identical(pr1$beta, direct$beta)
  expect_identical(pr1$intercept, direct$intercept)
  expect_no_match(format(pr1), "fitted out of bag")

  # reproducible with the same folds and seed
  pr2 <- prune_components(fit, d$x, d$y, select_by = "cv", fold_id = fold_id, seed = 3, select = "min")
  expect_identical(pr2$cv, pr$cv)
  expect_identical(pr2$beta, pr$beta)
})

test_that("the result does not depend on nthreads", {
  d <- prune_sim()
  fold_id <- rep_len(1:3, nrow(d$x))
  pr1 <- prune_components(prune_fit(d), d$x, d$y, select_by = "cv", fold_id = fold_id, seed = 3, nthreads = 1)
  pr2 <- prune_components(prune_fit(d), d$x, d$y, select_by = "cv", fold_id = fold_id, seed = 3, nthreads = 2)
  expect_identical(pr2$cv, pr1$cv)
  expect_identical(pr2$beta, pr1$beta)
  expect_identical(pr2$components, pr1$components)
  expect_identical(predict(pr2, d$x), predict(pr1, d$x))
  expect_identical(predict(pr1, d$x, nthreads = 2), predict(pr1, d$x))
  expect_error(prune_components(prune_fit(d), d$x, d$y, rule = "D", tau = 0, nthreads = 0))
  expect_error(predict(pr1, d$x, nthreads = 0))
})

test_that("predict() is the intercept plus the scaled kept components", {
  d <- prune_sim()
  fit <- prune_fit(d)
  pr <- prune_components(fit, d$x, d$y, rule = "R", tau = 0.01)
  expect_null(pr$cv)
  expect_null(pr$oob)
  expect_null(pr$seed)
  expect_null(pr$fold_id)
  expect_identical(pr$rule, "R")
  # a given R: no choice, but its coefficients are fitted out of bag (test-prune-oob-refit.R)
  expect_identical(pr$select_by, "oob")
  expect_output(print(pr), "fitted out of bag")
  expect_output(print(pr), "Given, not chosen")

  new <- as.data.frame(d$x[1:50, ])
  comp <- predict_components(fit, new)$m
  manual <- pr$intercept + as.numeric(as.matrix(comp[, names(pr$beta), with = FALSE]) %*% pr$beta)
  pred <- predict(pr, new)
  expect_s3_class(pred, "tbl_df")
  expect_equal(pred$.pred, manual, tolerance = 1e-10)
  # a matrix works too
  expect_identical(predict(pr, d$x[1:50, ]), pred)

  # predict_components() of the pruned model: the kept components scaled by beta
  pc <- predict_components(pr, new)
  expect_s3_class(pc, "rpf_components")
  expect_identical(names(pc$m), names(pr$beta))
  expect_equal(pc$m[[names(pr$beta)[1]]], pr$beta[[1]] * comp[[names(pr$beta)[1]]], tolerance = 1e-10)
  expect_equal(rowSums(pc$m) + pc$intercept, pred$.pred, tolerance = 1e-10)
  expect_null(pc$remainder)
  expect_identical(names(pc$x), colnames(d$x))

  # a subset of the predictors or a lower order: the rest goes into the remainder
  pc1 <- predict_components(pr, new, predictors = "x1")
  expect_identical(names(pc1$m), intersect(names(pr$beta), "x1"))
  expect_identical(names(pc1$x), "x1")
  expect_equal(rowSums(pc1$m) + pc1$intercept + pc1$remainder, pred$.pred, tolerance = 1e-10)
  pc2 <- predict_components(pr, new, max_interaction = 1)
  expect_true(all(!grepl(":", names(pc2$m))))
  expect_equal(rowSums(pc2$m) + pc2$intercept + pc2$remainder, pred$.pred, tolerance = 1e-10)
  expect_error(predict_components(pr, new, max_interaction = 3))
  expect_error(predict_components(pr, new, predictors = "x9"))
  expect_error(predict_components(pr, new, max_interactions = 1), class = "rlib_error_dots_nonempty")
  expect_error(predict(pr, new, type = "numeric"), class = "rlib_error_dots_nonempty")

  # with select_by = "cv" the refit is least squares of y on the kept components (in sample)
  pin <- prune_components(fit, d$x, d$y, rule = "R", tau = 0.01, select_by = "cv")
  expect_null(pin$select_by)
  expect_null(pin$cv)
  expect_identical(names(pin$beta), names(pr$beta))
  expect_no_match(format(pin), "fitted out of bag")
  own <- predict_components(fit, as.data.frame(d$x))$m
  ls <- stats::lm.fit(cbind(1, as.matrix(own[, names(pin$beta), with = FALSE])), d$y)$coefficients
  expect_equal(unname(c(pin$intercept, pin$beta)), unname(ls), tolerance = 1e-8)
})

test_that("rule D with tau = 0 is the unpruned forest", {
  d <- prune_sim()
  fit <- prune_fit(d)
  pr <- prune_components(fit, d$x, d$y, rule = "D", tau = 0)
  expect_true(all(pr$components$kept))
  expect_true(all(pr$beta == 1))
  new <- as.data.frame(d$x[1:50, ])
  expect_equal(predict(pr, new)$.pred, predict(fit, new)$.pred, tolerance = 1e-10)
  expect_equal(pr$intercept, predict_components(fit, new)$intercept)
})

test_that("a threshold above every share keeps no component: the intercept or the mean", {
  d <- prune_sim()
  fit <- prune_fit(d)
  new <- d$x[1:5, ]
  pd <- prune_components(fit, d$x, d$y, rule = "D", tau = 0.99)
  expect_length(pd$beta, 0)
  expect_false(any(pd$components$kept))
  expect_equal(predict(pd, new)$.pred, rep(pd$intercept, 5))
  expect_equal(ncol(predict_components(pd, new)$m), 0)
  # R: the mean of y over the rows with out-of-bag families (fitted out of bag), or over all rows (in sample)
  pr <- prune_components(fit, d$x, d$y, rule = "R", tau = 0.99)
  use <- rowSums(fit$fit$purify_oob(0L, 1L)$count) > 0
  expect_equal(pr$intercept, mean(d$y[use]))
  expect_output(print(pr), "Kept: +0")
  pin <- prune_components(fit, d$x, d$y, rule = "R", tau = 0.99, select_by = "cv")
  expect_equal(pin$intercept, mean(d$y))
})

test_that("the formula interface takes the outcome from data, and print() works", {
  d <- prune_sim(n = 200)
  df <- data.frame(d$x, y = d$y)
  set.seed(5)
  fit <- rpf(y ~ ., data = df, max_interaction = 2, ntrees = 5, splits = 20)
  pr <- prune_components(fit, df, rule = "D", tau = c(0, 0.01), select_by = "cv", folds = 3, seed = 1)
  expect_equal(nrow(pr$cv), 2)
  expect_length(unique(pr$fold_id), 3)
  expect_identical(pr$seed, 1L)
  expect_equal(nrow(predict(pr, df[1:3, ])), 3)
  expect_output(print(pr), "Pruned Random Planted Forest")
  expect_output(print(pr), "Chosen by 3-fold CV")
  expect_output(print(pr), "Kept components")
  expect_output(print(pr, n_show = 1), "and [0-9]+ more|Kept components")
  output <- capture.output(result <- withVisible(print(pr)))
  expect_false(result$visible)
  expect_identical(result$value, pr)
  expect_identical(output, format(pr))

  # without a seed, one is drawn and stored, and it reproduces the CV
  pr1 <- prune_components(fit, df, rule = "R", tau = c(0, 0.01), select_by = "cv", folds = 3)
  pr2 <- prune_components(fit, df, rule = "R", tau = c(0, 0.01), select_by = "cv", folds = 3, seed = pr1$seed)
  expect_identical(pr2$fold_id, pr1$fold_id)
  expect_identical(pr2$cv, pr1$cv)

  # the outcome may also be given as y
  pr3 <- prune_components(
    fit,
    df[, colnames(d$x)],
    d$y,
    rule = "D",
    tau = c(0, 0.01),
    select_by = "cv",
    folds = 3,
    seed = 1
  )
  expect_identical(pr3$cv, pr$cv)
})

test_that("factor predictors are coded per fold, as rpf() codes them", {
  d <- prune_sim(n = 300)
  df <- data.frame(d$x, g = factor(sample(letters[1:4], 300, replace = TRUE)))
  df$y <- d$y + c(a = 0, b = 1, c = -1, d = 2)[as.character(df$g)]
  set.seed(5)
  fit <- rpf(y ~ ., data = df, max_interaction = 2, ntrees = 10, splits = 30)
  pr <- prune_components(fit, df, select_by = "cv", folds = 3, seed = 2)
  expect_true("g" %in% names(pr$beta))
  expect_equal(nrow(predict(pr, df)), 300)
  # a fold fit is rpf() on the fold's rows, with the same seed
  ns <- asNamespace("randomPlantedForest")
  rows <- pr$fold_id != 1
  tr <- ns$prune_training_data(fit, df, NULL)
  refit <- ns$prune_refit(fit, tr$predictors, tr$y, rows, seed = 7, nthreads = 1L)
  set.seed(7)
  direct <- rpf(y ~ ., data = df[rows, ], max_interaction = 2, ntrees = 10, splits = 30)
  expect_identical(refit$factor_levels, direct$factor_levels)
  expect_identical(predict(refit, df), predict(direct, df))
})

test_that("a fold fit keeps the level order of the factors the forest keeps it for", {
  d <- prune_sim(n = 300)
  df <- data.frame(d$x, g = factor(sample(letters[1:4], 300, replace = TRUE), levels = c("d", "b", "a", "c")))
  df$y <- d$y + c(a = 0, b = 1, c = -1, d = 2)[as.character(df$g)]
  set.seed(5)
  fit <- rpf(y ~ ., data = df, max_interaction = 2, ntrees = 10, splits = 30, keep_level_order = "g")
  ns <- asNamespace("randomPlantedForest")
  rows <- rep_len(1:3, 300) != 1
  tr <- ns$prune_training_data(fit, df, NULL)
  refit <- ns$prune_refit(fit, tr$predictors, tr$y, rows, seed = 7, nthreads = 1L)
  set.seed(7)
  direct <- rpf(y ~ ., data = df[rows, ], max_interaction = 2, ntrees = 10, splits = 30, keep_level_order = "g")
  expect_identical(refit$factor_levels, direct$factor_levels)
  expect_identical(refit$factor_levels$g, c("d", "b", "a", "c"))
  expect_identical(predict(refit, df), predict(direct, df))
})

test_that("predict_components() selects the right components with a factor before numeric columns", {
  # hardhat::mold() puts the factors last, so the fit's columns are x1, x2, g, not the data's g, x1, x2
  set.seed(3)
  df <- data.frame(g = factor(sample(letters[1:3], 300, replace = TRUE)), x1 = runif(300), x2 = runif(300))
  df$y <- 2 * (df$g == "b") + df$x1 + df$x1 * df$x2 + rnorm(300, sd = 0.2)
  set.seed(5)
  fit <- rpf(y ~ ., data = df, max_interaction = 2, ntrees = 10, splits = 30)
  pr <- prune_components(fit, df, rule = "D", tau = 0)
  all_m <- predict_components(pr, df)$m
  for (p in list("x1", "g", c("g", "x2"))) {
    pc <- predict_components(pr, df, predictors = p)
    comps <- names(pc$m)
    # exactly the components within p, by their names
    expect_identical(comps, names(all_m)[vapply(strsplit(names(all_m), ":"), function(f) all(f %in% p), TRUE)])
    expect_gt(length(comps), 0)
    expect_equal(rowSums(as.matrix(pc$m)) + pc$remainder, rowSums(as.matrix(all_m)))
  }
  pc1 <- predict_components(pr, df, max_interaction = 1)
  expect_setequal(names(pc1$m), c("x1", "x2", "g"))
})

test_that("the fold fits use all parameters of the forest", {
  d <- prune_sim(n = 300)
  args <- list(
    max_interaction = 3,
    ntrees = 4,
    splits = 25,
    split_structure = "res_trees",
    split_try = 4,
    t_try = 0.7,
    max_candidates = 20,
    split_decay_rate = 0.3,
    delete_leaves = FALSE
  )
  set.seed(5)
  fit <- do.call(rpf, c(list(d$x, d$y, nthreads = 2), args))
  ns <- asNamespace("randomPlantedForest")
  tr <- ns$prune_training_data(fit, d$x, d$y)
  # on all rows with the fit's seed, the refit is the forest itself (threads change no fit)
  refit <- ns$prune_refit(fit, tr$predictors, tr$y, rep(TRUE, 300), seed = 5, nthreads = 1L)
  expect_identical(refit$fit$get_model(), fit$fit$get_model())
  expect_identical(refit$params[names(args)], fit$params[names(args)])
  expect_identical(predict(refit, d$x), predict(fit, d$x))

  # with interaction constraints (an argument of rpf() once they exist): every fold fit keeps them
  if ("interaction_constraints" %in% names(formals(ns$rpf.matrix))) {
    set.seed(5)
    cfit <- rpf(
      d$x,
      d$y,
      max_interaction = 2,
      ntrees = 20,
      splits = 40,
      interaction_constraints = list("x1", c("x2", "x3"), "x4")
    )
    pr <- prune_components(cfit, d$x, d$y, select_by = "cv", fold_id = rep_len(1:3, nrow(d$x)), seed = 3)
    # at most x1, ..., x4 and x2:x3 in every fold
    expect_true(all(pr$cv$n_kept <= 5))
    expect_true(all(pr$components$component %in% c("x1", "x2", "x3", "x4", "x2:x3")))
  }
})

test_that("a deterministic forest is pruned by CV", {
  d <- prune_sim(n = 200)
  set.seed(5)
  fit <- rpf(d$x, d$y, max_interaction = 2, ntrees = 1, splits = 20, deterministic = TRUE)
  # it has no out-of-bag rows: the default scores by CV
  pr <- prune_components(fit, d$x, d$y, folds = 3, seed = 1)
  expect_identical(pr$select_by, "cv")
  expect_null(pr$oob)
  expect_equal(nrow(pr$cv), 16)
  expect_equal(nrow(predict(pr, d$x)), 200)
})

test_that("bad input gives clear errors", {
  d <- prune_sim(n = 100)
  fit <- prune_fit(d)
  expect_error(prune_components(fit), "training rows")
  expect_error(prune_components(fit, d$x), "Pass the outcome as")
  expect_error(prune_components(fit, d$x, d$y[-1]), "has 99 values but `data` has 100 rows")
  expect_error(prune_components(fit, d$x, as.character(d$y)), "numeric outcome")
  expect_error(prune_components(fit, as.list(as.data.frame(d$x)), d$y), "data.frame or matrix")
  expect_error(prune_components(fit, d$x, d$y, tau = 1), "below 1")
  expect_error(prune_components(fit, d$x, d$y, tau = -0.1))
  expect_error(prune_components(fit, d$x, d$y, tau = c(0, 0)))
  expect_error(prune_components(fit, d$x, d$y, rule = "X"))
  expect_error(prune_components(fit, d$x, d$y, select = "max"))
  expect_error(prune_components(fit, d$x, d$y, select_by = "bag"))
  expect_error(prune_components(fit, d$x, d$y, fold_id = rep(1, 100)), "at least two distinct folds")
  expect_error(prune_components(fit, d$x, d$y, fold_id = 1:3))
  expect_error(prune_components(fit, d$x, d$y, folds = 1))
  expect_error(prune_components(fit, d$x, d$y, seed = "a"))
  expect_error(prune_components(fit, d$x, d$y, "R"), class = "rlib_error_dots_nonempty")
  expect_error(prune_components(fit, d$x, d$y, rules = "R"), class = "rlib_error_dots_nonempty")
  expect_error(prune_components(list(), d$x, d$y))

  yc <- factor(d$y > 1)
  set.seed(1)
  cfit <- rpf(d$x, yc, ntrees = 5, splits = 10)
  expect_error(prune_components(cfit, d$x, yc), "regression forests only")
})
