# The arguments renamed after 0.5.0 keep their old names as deprecated aliases: the same fit, bit for bit, one
# deprecation warning per call, and an error when both names are given.

renamed_fit <- function(..., seed = 1) {
  set.seed(seed)
  rpf(mpg ~ ., data = mtcars, max_interaction = 2, prune = FALSE, ...)
}

# the whole fit: every family's trees, leaf values and intervals, and the predictions
renamed_numbers <- function(fit) {
  list(fit$fit$get_model(), predict(fit, mtcars), fit$params)
}

test_that("each old argument name gives the fit of its new name, with a deprecation warning", {
  new_args <- list(n_families = 3, n_splits = 20, n_try_split_points = 4, share_candidates = 0.3)
  old_names <- c(
    n_families = "ntrees",
    n_splits = "splits",
    n_try_split_points = "split_try",
    share_candidates = "t_try"
  )
  ref <- do.call(renamed_fit, new_args)
  for (new in names(new_args)) {
    args <- new_args
    names(args)[names(args) == new] <- old_names[[new]]
    expect_warning(
      fit <- do.call(renamed_fit, args),
      sprintf("Use `%s` instead of `%s`", new, old_names[[new]]),
      class = "rpf_deprecated"
    )
    expect_identical(renamed_numbers(fit), renamed_numbers(ref))
  }
})

test_that("all old names at once give one warning and the same fit", {
  ref <- renamed_fit(n_families = 3, n_splits = 20, n_try_split_points = 4, share_candidates = 0.3)
  warnings <- list()
  fit <- withCallingHandlers(
    renamed_fit(ntrees = 3, splits = 20, split_try = 4, t_try = 0.3, split_structure = "leaves"),
    warning = function(w) {
      warnings[[length(warnings) + 1]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1)
  expect_s3_class(warnings[[1]], "lifecycle_warning_deprecated")
  msg <- conditionMessage(warnings[[1]])
  for (old in c("ntrees", "splits", "split_try", "t_try", "\"leaves\"")) {
    expect_match(msg, old, fixed = TRUE)
  }
  expect_identical(renamed_numbers(fit), renamed_numbers(ref))
})

test_that("the new names give no warning", {
  expect_no_warning(renamed_fit(n_families = 3, n_splits = 20, n_try_split_points = 4, share_candidates = 0.3))
})

test_that("each old split structure name gives the fit of its new name", {
  old_new <- c(
    leaves = "leaf",
    hist = "leaf_binned",
    cur_trees_1 = "tree_then_leaf",
    cur_trees_2 = "tree_then_all_leaves",
    res_trees = "target_tree_then_leaf"
  )
  for (old in names(old_new)) {
    ref <- renamed_fit(n_families = 2, n_splits = 15, split_structure = old_new[[old]])
    expect_warning(
      fit <- renamed_fit(n_families = 2, n_splits = 15, split_structure = old),
      sprintf("split_structure = \"%s\"", old_new[[old]]),
      fixed = TRUE,
      class = "rpf_deprecated"
    )
    expect_identical(renamed_numbers(fit), renamed_numbers(ref))
    # the stored params keep the old code
    expect_identical(fit$params$split_structure, old)
  }
})

test_that("giving an argument under both names is an error", {
  both <- list(
    list(ntrees = 3, n_families = 3),
    list(splits = 20, n_splits = 20),
    list(split_try = 4, n_try_split_points = 4),
    list(t_try = 0.3, share_candidates = 0.3),
    # NULL counts as given, under either name
    list(splits = NULL, n_splits = 20),
    list(splits = 20, n_splits = NULL)
  )
  for (args in both) {
    expect_error(do.call(renamed_fit, args), "were given")
  }
})

test_that("an old name given as NULL is the new name's NULL", {
  set.seed(1)
  ref <- rpf(mpg ~ ., data = mtcars, n_families = 10, n_splits = NULL, prune = FALSE)
  set.seed(1)
  expect_warning(
    fit <- rpf(mpg ~ ., data = mtcars, n_families = 10, splits = NULL, prune = FALSE),
    class = "rpf_deprecated"
  )
  expect_identical(fit$params$splits, ref$params$splits)
  expect_identical(fit$oob_path, ref$oob_path)
  expect_identical(renamed_numbers(fit), renamed_numbers(ref))
})

test_that("the old names work in every rpf() method, and other unknown arguments are still an error", {
  x <- mtcars[, c("cyl", "wt", "hp")]
  set.seed(1)
  ref <- rpf(x, mtcars$mpg, n_families = 2, n_splits = 10, prune = FALSE)
  set.seed(1)
  expect_warning(fit <- rpf(x, mtcars$mpg, ntrees = 2, splits = 10, prune = FALSE), class = "rpf_deprecated")
  expect_identical(predict(fit, x), predict(ref, x))
  set.seed(1)
  expect_warning(fit <- rpf(as.matrix(x), mtcars$mpg, ntrees = 2, splits = 10, prune = FALSE), class = "rpf_deprecated")
  expect_identical(predict(fit, as.matrix(x)), predict(ref, x))
  skip_if_not_installed("recipes")
  rec <- recipes::recipe(mpg ~ cyl + wt + hp, data = mtcars)
  set.seed(1)
  expect_warning(fit <- rpf(rec, mtcars, ntrees = 2, splits = 10, prune = FALSE), class = "rpf_deprecated")
  expect_identical(predict(fit, mtcars), predict(ref, x))
  expect_error(rpf(x, mtcars$mpg, n_tree = 2), class = "rlib_error_dots_nonempty")
  expect_error(rpf(x, mtcars$mpg, ntrees = 2, foo = 1), class = "rlib_error_dots_nonempty")
})

test_that("tune: an old name counts as given", {
  expect_error(
    suppressWarnings(rpf(mpg ~ ., data = mtcars, n_families = 10, tune = TRUE, split_try = 3)),
    "n_try_split_points"
  )
})

test_that("oob_path_checkpoints() takes the old name", {
  expect_warning(
    cp <- oob_path_checkpoints(splits = 100),
    "Use `n_splits` instead of `splits`",
    class = "rpf_deprecated"
  )
  expect_identical(cp, oob_path_checkpoints(n_splits = 100))
  expect_error(suppressWarnings(oob_path_checkpoints(100, splits = 100)), "were given")
  expect_error(oob_path_checkpoints(100, ntrees = 5), "Unused argument")
})

test_that("tune_control() takes a ladder with the old column names", {
  old <- data.frame(split_try = c(3, 1), t_try = c(0.1, 0.025))
  new <- data.frame(n_try_split_points = c(3, 1), share_candidates = c(0.1, 0.025))
  expect_warning(tc <- tune_control(ladder = old, search = 1:2), class = "rpf_deprecated")
  expect_identical(tc, tune_control(ladder = new, search = 1:2))
  expect_error(tune_control(ladder = cbind(old, new), search = 1:2), "give one")
})

test_that("print() and the result tables use the new names", {
  set.seed(1)
  fit <- rpf(mpg ~ ., data = mtcars, n_families = 10, prune = FALSE)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "n_try_split_points: 3", fixed = TRUE)
  expect_match(out, "share_candidates: 0.1", fixed = TRUE)
  expect_match(out, "split_structure: leaf\n", fixed = TRUE)
  expect_match(out, "`n_splits` chosen out of bag", fixed = TRUE)
  expect_identical(names(fit$oob_path)[1], "n_splits")
  expect_identical(names(fit$splits_oob), c("n_splits", "stopped_at", "capped"))
  expect_named(greediness_ladder(), c("n_try_split_points", "share_candidates"))
})
